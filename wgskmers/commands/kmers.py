"""Commands for managed stored k-mer sets"""

from itertools import izip

import click
from tqdm import tqdm

from .util import choose_db, with_db


class RefCalculator(object):

	@classmethod
	def init(cls, db, spec):
		cls.db = db
		cls.spec = spec

	@classmethod
	def calc_ref(cls, genome):
		from Bio import SeqIO
		from wgskmers.parse import vec_from_records

		with cls.db.open_genome(genome) as fh:

			records = SeqIO.parse(fh, genome.file_format)

			# If assembled, get boolean vector. Otherwise get counts (why not)
			return vec_from_records(records, cls.spec, counts=not genome.is_assembled)


@click.group(name='refs',
             short_help='Calculate and manage reference k-mer collections')
@choose_db()
def kmers_group():
	"""Command group for managing pre-calculated reference k-mer collections"""
	pass


@kmers_group.command(short_help='List stored reference k-mer collections')
@with_db()
def listc(ctx, db):
	"""List stored reference k-mer collections

	The ID of the collection is listed at the beginning of the line.
	"""
	from wgskmers.database.models import KmerSetCollection

	session = db.get_session()
	for collection in session.query(KmerSetCollection).all():

		attrs = {a: getattr(collection, a) for a
		         in ['id', 'k', 'prefix', 'title']}
		count = collection.kmer_sets.count()
		click.echo(
			'{id}: [{k} - {prefix}] "{title}" ({0} calculated sets)'
			.format(count, **attrs)
		)


@kmers_group.command(short_help='Create new collection of reference k-mers')
@click.argument('k', type=int)
@click.argument('prefix', type=str)
@click.argument('title', type=str)
@with_db(confirm=True)
def makec(ctx, db, k, prefix, title):
	"""Create a new collection of reference k-mers with given parameters

	Args:
		[K]: Length of k-mers to find, INCLUDING prefix
		[PREFIX]: Nucleotide sequence that k-mers must start with
		[TITLE]: Unique title for k-mer collection
	"""
	from wgskmers.kmers import nucleotides
	from wgskmers.database.models import KmerSetCollection

	# Check k
	if k <= 0:
		raise click.ClickException('K must be positive')

	# Check prefix
	prefix = prefix.upper()
	if not prefix:
		raise click.ClickException('Prefix cannot be empty')
	if len(prefix) >= k:
		raise click.ClickException('Length of prefix must be <= K')
	if not set(prefix).issubset(nucleotides):
		raise click.ClickException('Prefix contains invalid characters')

	# Check title
	title = title.strip()
	if not title:
		raise click.ClickException('Title cannot be empty')

	session = db.get_session()
	existing = session.query(KmerSetCollection).filter_by(title=title)
	if existing.first() is not None:
		raise click.ClickException(
			'A k-mer collection already exists with this title')

	# Create it
	collection = db.create_kmer_collection(k=k, prefix=prefix, title=title,
	                                       format='coords')

	click.echo(
		'K-mer collection "{}" created with ID {}'
		.format(title, collection.id)
	)


@kmers_group.command(short_help='Calculate k-mer sets and add to collection')
@click.argument('collection_id', type=int)
@with_db(confirm=True)
def calc(ctx, db, collection_id):
	"""Caluclate k-mer sets and add to collection

	Currently this just calculates k-mer sets for all stored genomes and
	adds them to the collection.

	You can get the collection ID by running listc (it will be printed at the
	beginning of the line).
	"""
	import multiprocessing as mp

	from wgskmers.kmers import KmerSpec
	import wgskmers.multiprocess as kmp
	from wgskmers.database.models import Genome, KmerSet, KmerSetCollection

	# Get collection
	session = db.get_session()
	collection = session.query(KmerSetCollection).get(collection_id)
	if collection is None:
		raise click.ClickException(
			'No k-mer collection with id {}'
			.format(collection_id)
		)

	store_set = db.store_kmer_sets(collection)

	spec = KmerSpec(k=collection.k, prefix=collection.prefix)

	# Get genomes not already calculated in collection
	genome_query = session.query(Genome).filter(
		~Genome.kmer_sets.any(KmerSet.collection == collection)
	)
	genomes = genome_query.all()

	# Create pool
	init_args = (db, spec)
	pool = mp.Pool(initializer=RefCalculator.init, initargs=init_args,
	               maxtasksperchild=100)

	kmp.enable_method_pickling()

	try:

		# Start the workers
		results = pool.imap(RefCalculator.calc_ref, genomes)
		pool.close()

		# Iterate through results
		added, errors = 0, 0
		for vec, genome in tqdm(izip(results, genomes), total=len(genomes)):

			# Try adding the set
			try:
				store_set(vec, genome, has_counts=not genome.is_assembled)
				added += 1

			# Print exception and continue
			except Exception as e:
				click.secho(
					'Error finding k-mers for genome "{}": {}'
					.format(genome.description, e),
					err=True, fg='red'
				)
				errors += 1

		skipped = session.query(Genome).count() - added - errors
		click.echo(
			'Calculated {} sets, {} errors, {} already in collection'
			.format(added, errors, skipped)
		)

	finally:
		pool.terminate()
