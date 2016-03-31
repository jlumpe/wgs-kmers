""""""

import os
import re
from glob import glob
from functools import wraps
from collections import OrderedDict
from textwrap import dedent
from csv import DictWriter, DictReader

import click
from tqdm import tqdm

from wgskmers import database, models


gb_header_re = re.compile(r'>gi\|(\d+)\|(?:gb|ref|emb)\|(.+)\|(.*)')


def parse_bool(string):
	string = string.lower().strip()

	if string in ['yes', 'true', 't']:
		return True
	elif string in ['no', 'false', 'f']:
		return False
	else:
		raise ValueError(string)


def guess_fasta_attrs(path):

	attrs = dict(file_format='fasta')

	with open(path) as fh:

		# Get first line, should be header
		# (otherwise, not much we can do...)
		first = fh.readline()
		if first.startswith('>'):

			# Match genbank header format
			match = gb_header_re.match(first.strip())
			if match is not None:
				gi_str, acc, desc = match.groups()
				attrs['description'] = desc.strip()
			else:
				gi_str, acc, desc = None, None, None

			# Check to see if other headers present
			# If so, attributes guessed from the first likely differ
			for line in fh:
				if line.startswith('>'):
					attrs['assembled'] = False
					break

			else:
				attrs['assembled'] = True

				if match is not None:
					attrs['ncbi_acc'] = acc.strip()

					try:
						attrs['ncbi_gi'] = int(gi_str)
					except ValueError:
						pass

	return attrs


genome_import_attrs = OrderedDict([
	('description', str),
	('ncbi_gi', int),
	('ncbi_acc', str),
	('assembled', parse_bool),
	('organism_name', str),
	('organism_taxonomy', str),
	('file_format', lambda value: str(value).lower()),
])

genome_show_attrs = ['id'] + genome_import_attrs.keys() + [
	'created_at',
	'updated_at',
]

genome_import_cols = genome_import_attrs.copy()
genome_import_cols['file'] = str


def parse_import_csv(fh, db):

	session = db.get_session()

	return_vals = []
	uq_vals = dict()

	# Read in template file

	reader = DictReader(fh)

	# Check missing columns
	missing_cols = set(genome_import_cols.keys()).difference(reader.fieldnames)
	if missing_cols:
		raise click.ClickException(
			'Missing column "{}" in import file'
			.format(missing_cols.pop())
		)

	for i, row in enumerate(reader):

		err_prefix = 'Error on row {}: '.format(i + 1)

		# Check file
		path = row['file']
		if not path.strip():
			raise click.ClickException(err_prefix + 'invalid file')
		elif not os.path.isfile(path):
			raise click.ClickException(
				err_prefix +
				'File {} does not exist'.format(path)
			)

		# Get attributes from row
		attrs = dict()
		for attrname, converter in genome_import_attrs.iteritems():
			val = row[attrname].strip()

			# Empty string is None
			if not val:
				attrs[attrname] = None

			# Try converting to correct format
			else:
				try:
					attrs[attrname] = converter(val)
				except ValueError:
					raise click.ClickException(
						err_prefix + 
						'Invalid value for column "{}"'.format(attrname)
					)

		# Check required attributes
		for attrname in ['description', 'assembled', 'file_format']:
			if attrs[attrname] is None:
				raise click.ClickException(
					err_prefix +
					'A value for column "{}" is required'.format(attrname)
				)

		# Check fasta format
		if attrs['file_format'] != 'fasta':
			raise click.ClickException(
				err_prefix + 
				'Only "fasta" file_format is currently supported'
			)

		# Check uniqueness of columns
		for uq_col in ['description', 'ncbi_gi', 'ncbi_acc']:
			val = attrs[uq_col]
			if val is None:
				continue

			# Check already in database
			found = (session.query(models.Genome)
			                .filter_by(**{uq_col: val})
			                .first())
			if found is not None:
				raise click.ClickException(
					err_prefix + 
					'Reference genome already exists in database with '
					'{} "{}"'
					.format(uq_col, val)
				)

			# Check already in file
			seen_vals = uq_vals.setdefault(uq_col, set())
			if val in seen_vals:
				raise click.ClickException(
					err_prefix + 
					'Duplicate value for {}'.format(uq_col, val)
				)

			seen_vals.add(val)

		# Should be ok
		return_vals.append((path, attrs))

	return return_vals


def with_db(confirm=False):
	"""
	Decorator for commands that use the database.

	Gets database from context and passes as second argument (after context).
	If the current database was automatically determined, print an explicit
	message or require confirmation (if confirm=True, used for write ops).
	"""
	def decorator(func):

		@wraps(func)
		def wrapper(ctx, *args, **kwargs):

			# Check if current database chosen automatically
			auto_method = ctx.obj['db_auto']
			if auto_method:

				# How was current database determined?
				if auto_method == 'cwd':
					method_str = 'current working directory'
				elif auto_method == 'config':
					method_str = 'default database in config file'
				elif auto_method == 'environ':
					method_str = ('from environment variable {}'
					              .format(database.DEFAULT_DB_PATH_VAR))
				elif auto_method == 'override':
					method_str = 'overridden in shell'
				else:
					method_str = 'huh?'

				# Require confirmation
				if confirm:
					click.confirm('Using database {} ({}), ok?'
					              .format(ctx.obj['db'].directory, method_str),
					              abort=True)

				# Print message
				else:
					click.secho('Using database {} ({})'
					            .format(ctx.obj['db'].directory, method_str),
					            fg='yellow', err=True)

			# Call wrapped function
			func(ctx, ctx.obj['db'], *args, **kwargs)

		# Return wrapper, with context
		return click.pass_context(wrapper)

	# Return decorator
	return decorator


@click.group(name='gen', short_help='Manage reference genomes')
@click.option('-N', '--db', type=str, help='Use registered database by name')
@click.option('-d', '--default', is_flag=True, help='Use default database')
@click.pass_context
def genomes_group(ctx, db=None, default=False):
	"""
	Commands to manage reference genomes in database.

	Specify database to use explicitly with "--db" or "--default" commands. If
	neither is given will use either database in the current working
	directory or the one set as default, with confirmation for write ops.
	"""
	# Use default database
	if default:
		if db is not None:
			raise click.ClickException('Cant use both --db and --default options')

		db_path = database.get_default_db()
		if db_path is None:
			raise click.ClickException('No default database set')

		ctx.obj['db_auto'] = False

	# Use registered database by name
	elif db is not None:
		db = db.lower()
		try:
			db_path = database.get_registered_dbs()[db]
		except KeyError:
			raise click.ClickException('Database "{}" not found'.format(db))

		ctx.obj['db_auto'] = False

	# Determine automatically
	else:
		db_path, method = database.get_current_db()
		if db_path is None:
			raise click.ClickException('No database currently active')

		ctx.obj['db_auto'] = method

	# Open database and add to context
	ctx.obj['db'] = database.Database.open(db_path)


@genomes_group.command()
@click.option('-c', '--csv', 'out_csv', is_flag=True)
@click.argument('dest', type=click.File('w'), default='-')
@with_db()
def list(ctx, db, dest, out_csv=False):
	"""List reference genomes"""

	session = db.get_session()
	genomes = session.query(models.Genome).all()

	if out_csv:

		writer = DictWriter(dest, genome_show_attrs)
		writer.writeheader()

		for g in genomes:
			writer.writerow({c: getattr(g, c) for c in genome_show_attrs})

	else:
		for g in genomes:
			click.echo(g.description)


@genomes_group.command(name='import', short_help='Import reference genomes')
@click.option('-c', '--csv', 'csv_out', type=click.Path(),
	help='Path to write import csv file to')
@click.option('-e', '--existing', type=click.File(),
	help='Existing import template .csv file to use')
@click.argument('directory', type=click.Path(exists=True))
@with_db(confirm=True)
def import_genomes(ctx, db, directory, csv_out=None, existing=None):

	directory = os.path.abspath(directory)
		
	# Create a new import template file
	if existing is None:

		if csv_out is None:
			csv_out = os.path.join(directory, 'genomes-import.csv')

		# Write import .csv template
		with open(csv_out, 'w') as fh:

			writer = DictWriter(fh, genome_import_cols)
			writer.writeheader()
			
			# Find fasta files in directory
			directory = os.path.abspath(directory)
			for path in glob(os.path.join(directory, '*.fasta')):

				# Guess attributes
				attrs = guess_fasta_attrs(path)
				attrs['file'] = path
				writer.writerow(attrs)

		# Open it if possible
		click.launch(csv_out)

		# Instructions
		click.echo(dedent('''
			Import file tempate has been written to {}

			You may edit fields in this file to your liking, or even add and
			remove rows (but keep the same columns).

			If this operation is canceled, you may re-use the same import
			file by running this command again with the --existing option.
			'''.format(csv_out)))

		# Wait for user to confirm
		click.confirm('Please confirm when you have finished editing',
		              abort=True)

		csv_fh = open(csv_out)

	else:
		csv_fh = existing

	# Read in template file and parse attributes
	added = 0
	errors = 0
	with csv_fh:
		for file_path, genome_attrs in tqdm(parse_import_csv(csv_fh, db)):

			# Try adding it
			try:
				db.store_genome(file_path, **genome_attrs)
				added += 1
			except Exception as e:
				click.echo(
					'Error adding file {}: {}'
					.format(file_path, e),
					err=True
				)

	click.echo(
		'Successfully imported {} genomes, with {} errors'
		.format(added, errors)
	)
