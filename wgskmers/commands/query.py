"""Functions for querying against reference sequences in database"""

import os
import logging
from csv import DictWriter
from collections import namedtuple

import click

from .util import with_db


logger = logging.getLogger()


query_metric_keys = ['hamming', 'jaccard', 'asym_jacc']


QueryMatch = namedtuple('QueryMatch', ['query', 'ref', 'metric', 'rank', 'score'])

def top_matches(queries, refs, metrics, scores, n):
	import numpy as np
	from wgskmers.query import query_metrics

	for i, query in enumerate(queries):
		for j, metric_name in enumerate(metrics):

			metric = query_metrics[metric_name]

			scores_slice = scores[j, :, i]

			# Reverse order for distances
			if metric.is_distance:
				best_idx = np.argsort(scores_slice)
			else:
				best_idx = np.argsort(-scores_slice)

			for rank, idx in enumerate(best_idx[:n]):
				yield QueryMatch(
					query=query,
					ref=refs[idx],
					metric=metric,
					rank=rank,
					score=scores_slice[idx],
				)


def matches_to_csv(fh, queries, refs, metrics, scores, topn):
	from wgskmers.genbank import get_record_url

	writer = DictWriter(fh, [
		'query_file',
		'metric',
		'rank',
		'score',
		'description',
		'organism',
		'genus',
		'species',
		'strain',
		'set',
		'accession',
		'database',
		'link',
	])
	writer.writeheader()

	for match in top_matches(queries, refs, metrics, scores, topn):

		genome = match.ref.genome

		if genome.gb_db is not None:
			if genome.gb_acc is not None:
				link = get_record_url(genome.gb_acc, genome.gb_db)
			elif genome.gb_id is not None:
				link = get_record_url(genome.gb_id, genome.gb_db)
			else:
				link = None
		else:
			link = None

		writer.writerow(dict(
			query_file=match.query,
			metric=match.metric.title,
			rank=match.rank + 1,
			score=match.score,
			description=genome.description,
			organism=genome.organism,
			genus=genome.tax_genus,
			species=genome.tax_species,
			strain=genome.tax_strain,
			set=genome.genome_sets[0].name if genome.genome_sets else None,
			accession=genome.gb_acc,
			database=genome.gb_db,
			link=link,
		))


def print_matches(queries, refs, metrics, scores, topn):
	from wgskmers.query import query_metrics

	for i, query in enumerate(queries):

		click.echo('\n\n>{}'.format(query))

		for j, metric_name in enumerate(metrics):

			metric = query_metrics[metric_name]

			click.echo('\nTop {} scores by {}:'.format(topn, metric.title))

			scores_slice = scores[j:j+1, :, i:i+1]

			for match in top_matches([query], refs, [metric_name], scores_slice, topn):
				click.echo('{} {}'.format(match.score, match.ref.genome.description))


@click.command(name='query')

@click.option('-q', '--q-threshold', type=int,
	help='Filter k-mers in query containing PHRED scores below this value')
@click.option('-c', '--c-threshold', type=int, default=1,
	help='Filter k-mers in query occuring less than this many times')
@click.option('-f', '--format', default=None,
	type=click.Choice(['fasta', 'fastq', 'fastq-sanger', 'fastq-solexa',
	                   'fastq-illumina']),
	help='File format, as argument to Bio.SeqIO.parse. If omitted, will infer '
	     'from file extensions.')
@click.option('-m', '--metric', type=click.Choice(query_metric_keys + ['all']),
              default='all', help='Query metric to use')
@click.option('-n', '--n-results', type=int, default=10,
              help='Number of results to return, sorted by best')
@click.option('--csv', type=click.File('w'), help='Write output to csv file')
@click.option('--no-check-ext', is_flag=True,
              help='Don\'t filter batch mode files by sequence file extension')
@click.option('--no-print', is_flag=True,
              help='Don\'t print results to stdout')
@click.argument('collection_id', type=int)
@click.argument('src', type=click.Path(exists=True))
@with_db(choose=True)
def query_command(ctx, db, collection_id, src, **kwargs):
	"""Query a sequence against the reference database"""

	from wgskmers.database import KmerSetCollection
	from wgskmers.kmers import KmerSpec
	from wgskmers.query import mp_query_coords
	from wgskmers.parse import find_seq_files, SeqFileInfo, parse_to_array

	src = os.path.abspath(src)

	# Get arguments
	q_threshold = kwargs.pop('q_threshold', None)
	c_threshold = kwargs.pop('c_threshold')
	file_format = kwargs.pop('format', None)
	metric_choice = kwargs.pop('metric', 'all')
	n_results = kwargs.pop('n_results')
	check_ext = not kwargs.pop('no_check_ext', False)
	no_print = kwargs.pop('no_print', False)
	csv_out = kwargs.pop('csv', None)

	# Get collection
	session = db.get_session()
	collection = session.query(KmerSetCollection).get(collection_id)
	if collection is None:
		raise click.ClickException(
			'No k-mer collection with id {}'
			.format(collection_id)
		)

	# Kmer spec
	spec = KmerSpec(k=collection.k, prefix=collection.prefix)

	# Get input files
	if os.path.isdir(src):
		files_info = find_seq_files(src, filter_ext=True, filter_contents=True,
		                            warn_contents=True)
		if not files_info:
			raise click.ClickException('No sequence files found in {}'.format(src))

	else:
		info = SeqFileInfo.get(src)
		files_info = [info]

		# Check format
		if file_format is not None:
			info.seq_format = file_format
		elif info.seq_format is None:
			raise RuntimeError('Couldn\'t infer format for {}'.format(src))

		# Check contents
		info.check_contents()
		if info.contents_ok is False:
			raise RuntimeError('Bad file contents in {}: {}'
			                   .format(src, info.contents_error))

	query_files = [info.path for info in files_info]

	# Get the query vectors
	pbar_args = dict(desc='Parsing query files', unit=' file(s)', leave='False')
	query_array = parse_to_array(files_info, spec, progress=pbar_args,
			                     q_threshold=q_threshold,
			                     c_threshold=c_threshold)

	# Reference sets
	ref_sets = collection.kmer_sets.all()

	# Metrics
	if metric_choice == 'all':
		metric_names = query_metric_keys
	else:
		metric_names = [metric_choice]

	# Make the query
	scores = mp_query_coords(query_array, db, collection, ref_sets, metric_names,
	                  progress=True)

	# Print the scores
	if not no_print:
		print_matches(query_files, ref_sets, metric_names, scores, n_results)

	# Scores to csv
	if csv_out is not None:
		matches_to_csv(csv_out, query_files, ref_sets, metric_names, scores,
		               n_results)
