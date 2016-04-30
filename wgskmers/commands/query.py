"""Functions for querying against reference sequences in database"""

import os
import logging

import click
from tqdm import tqdm
import numpy as np

from wgskmers.database import KmerSetCollection
from wgskmers.kmers import KmerSpec, KmerFinder, QualityKmerFinder
from wgskmers.query import query_metrics
from .util import with_db, ProgressSeqParser


logger = logging.getLogger()


def kmers_from_records(records, spec, quality_threshold=None):
	"""Generator yielding KmerFinders for a set of sequence records.

	Args:
		records: iterable of Bio.SeqRecord.SeqRecord, as output from
			Bio.SeqIO.parse. Records to find k-mers in.
		spec. KmerSpec. Spec defining k-mers to search for.
		quality_threshold: numeric|None. If not None, get quality scores from
			records and filter out k-mers containing score below this value.

	Yields:
		KmerFinder or QualityKmerFinder, depending if quality_threshold was
			None or not.
	"""

	# Parse file and iterate over sequences
	for record in records:

		# Upper case for search
		seq = record.seq.upper()

		# No quality
		if quality_threshold is None:
			yield spec.find(seq, revcomp=True)

		# With quality info
		else:
			phred_scores = record.letter_annotations['phred_quality']
			yield spec.find_quality(seq, revcomp=True, quality=phred_scores,
			                        threshold=quality_threshold)


def infer_format(path):
	"""Infers format for sequence file, returns format argument to
	Bio.SeqIO.parse

	All it does is check the extension.
	"""
	ext = os.path.splitext(path)[1]
	if ext in ['.fasta', '.fastq']:
		return ext[1:]
	else:
		return None


@click.command(name='query')

@click.option('-q', '--q-threshold', type=int, required=False,
	help='Filter k-mers in query containing PHRED scores below this value')
@click.option('-c', '--c-threshold', type=int, required=False, default=1,
	help='Filter k-mers in query occuring less than this many times')
@click.option('-f', '--format', default=None,
	type=click.Choice(['fasta', 'fastq', 'fastq-sanger', 'fastq-solexa',
	                   'fastq-illumina']),
	help='File format, as argument to Bio.SeqIO.parse. If omitted, will infer '
	     'from extension of first file encountered.')
@click.option('-b', '--batch', is_flag=True, help='Run in batch mode.')
@click.option('-n', '--n-results', type=int, default=25,
              help='Number of results to return, sorted by best')
@click.argument('collection_id', type=int)
@click.argument('src', type=click.Path(exists=True))
@with_db(choose=True)
def query_command(ctx, db, collection_id, src, **kwargs):
	"""Query a sequence against the reference database"""

	src = os.path.abspath(src)

	# Get arguments
	q_threshold = kwargs.pop('q_threshold', None)
	c_threshold = kwargs.pop('c_threshold')
	file_format = kwargs.pop('format', None)
	batch_mode = kwargs.pop('batch', False)
	n_results = kwargs.pop('n_results')

	# Get collection
	session = db.get_session()
	collection = session.query(KmerSetCollection).get(collection_id)
	if collection is None:
		raise click.ClickException(
			'No k-mer collection with id {}'
			.format(collection_id)
		)

	# Reference sets
	ref_sets = collection.kmer_sets.all()

	# Kmer spec
	spec = KmerSpec(k=collection.k, prefix=collection.prefix)

	# Get format
	if file_format is None:
		file_format = infer_format(src)
		if file_format is None:
			raise ValueError("Couldn't infer format for {}"
				.format(src))
		logger.debug('Inferring format "{}" from {}'
		             .format(file_format, src))

	# Allocate array for results
	metrics = query_metrics.values()
	scores = np.ndarray((len(ref_sets), len(metrics)))

	# Parse with progress bar
	records = ProgressSeqParser(src, fmt=file_format,
	                            desc='Parsing source file')

	# Find the k-mers with quality info
	if q_threshold is not None:
		finders = kmers_from_records(records, spec,
		                             quality_threshold=q_threshold)

	else:
		# Just the k-mers themselves
		finders = kmers_from_records(records, spec)

	# Get counts
	counts_vec = None
	for finder in finders:
		counts_vec = finder.counts_vec(out=counts_vec)

	# Count threshold
	query_vec = counts_vec >= c_threshold

	# Now loop through reference kmer sets
	iterator = db.load_kmer_sets_lazy(collection, ref_sets)
	iterator = tqdm(iterator, desc='Querying reference genomes',
	                total=len(ref_sets))
	for i, vec in enumerate(iterator):

		ref_vec = vec > 0

		for j, metric in enumerate(metrics):
			scores[i, j] = metric(query_vec, ref_vec)

	# Print the scores
	for i, metric in enumerate(metrics):

		best_idx =  scores[:, i].argsort()
		if not metric.is_distance:
			best_idx = best_idx[::-1]

		click.echo('\nTop {} scores by {}:'.format(n_results, metric.name))

		for idx in best_idx[:n_results]:
			click.echo('{} {}'.format(scores[idx, i],
				ref_sets[idx].genome.description))
