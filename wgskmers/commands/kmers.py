""""""

import click
import sqlalchemy as sqla
from tqdm import tqdm
from Bio import SeqIO

from .util import choose_db, with_db
from wgskmers import models
from wgskmers.kmers import nucleotides, KmerFinder, KmerSpec


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

	session = db.get_session()
	for collection in session.query(models.KmerSetCollection).all():

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
	existing = session.query(models.KmerSetCollection).filter_by(title=title)
	if existing.first() is not None:
		raise click.ClickException(
			'A k-mer collection already exists with this title')

	# Create it
	collection = db.create_kmer_collection(k=k, prefix=prefix, title=title)

	click.echo(
		'K-mer collection "{}" created with ID {}'
		.format(title, collection.id)
	)


@kmers_group.command(short_help='Caluclate k-mer sets and add to collection')
@click.argument('collection_id', type=int)
@with_db(confirm=True)
def calc(ctx, db, collection_id):
	"""Caluclate k-mer sets and add to collection

	Currently this just calculates k-mer sets for all stored genomes and
	adds them to the collection.

	You can get the collection ID by running listc (it will be printed at the
	beginning of the line).
	"""

	# Get collection
	session = db.get_session()
	collection = session.query(models.KmerSetCollection).get(collection_id)
	if collection is None:
		raise click.ClickException(
			'No k-mer collection with id {}'
			.format(collection_id)
		)

	spec = KmerSpec(k=collection.k, prefix=collection.prefix)

	# Get genomes not already calculated in collection
	genome_query = session.query(models.Genome).filter(
		~models.Genome.kmer_sets.any(models.KmerSet.collection == collection)
	)

	# Iterate through genomes
	added, errors = 0, 0
	for genome in tqdm(genome_query.all()):

		# Parse it
		records = SeqIO.parse(db.open_genome(genome), genome.file_format)

		# Vector to store k-mers
		vec = None

		# Go through the records
		for record in records:
			finder = KmerFinder(spec, record.seq)

			# If assembled, get boolean vector
			if genome.assembled:
				vec = finder.bool_vec(out=vec)

			# Otherwise, get counts (why not)
			else:
				vec = finder.counts_vec(out=vec)

		# Try adding the set
		try:
			db.add_kmer_set(vec, collection, genome,
			                has_counts=not genome.assembled)
			added += 1

		# Print exception and continue
		except Exception as e:
			click.secho(
				'Error finding k-mers for genome "{}": {}'
				.format(genome.description, e),
				err=True, fg='red'
			)
			errors += 1

	skipped = session.query(models.Genome).count() - added - errors
	click.echo(
		'Calculated {} sets, {} errors, {} already in collection'
		.format(added, errors, skipped)
	)
