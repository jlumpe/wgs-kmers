"""Commands for managing stored reference genomes"""
from builtins import str

import os
import re
from collections import OrderedDict, namedtuple
from textwrap import dedent
from csv import DictWriter, DictReader

import click
from tqdm import tqdm

from .util import choose_db, with_db


gb_header_re = re.compile(r'>gi\|(\d+)\|(?:gb|ref|emb)\|(.+)\|(.*)')


def parse_bool(string):
	"""Parse a boolean value from string format"""
	string = string.lower().strip()

	if string in ['yes', 'y', 'true', 't']:
		return True
	elif string in ['no', 'n', 'false', 'f']:
		return False
	else:
		raise ValueError(string)


def import_str(obj, lower=False):
	"""Converts value to string and strips for import to string field"""
	if obj:
		s = str(obj).strip()
		if lower:
			s = s.lower()
		return s

	else:
		return None


def guess_fasta_attrs(info):
	from wgskmers.genbank import extract_acc

	assert info.seq_format == 'fasta'
	attrs = dict(file_format='fasta')

	with info.open() as fh:

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
					attrs['is_assembled'] = False
					break

			else:
				attrs['is_assembled'] = True

				if match is not None:
					attrs['gb_acc'] = acc.strip()

					try:
						attrs['gb_id'] = int(gi_str)
					except ValueError:
						pass

	# Try getting accession from file name
	fn_acc = extract_acc(info.basename)
	if attrs.get('gb_acc', None) is None and fn_acc is not None:
		attrs['gb_acc'] = fn_acc
		attrs['gb_id'] = None

	# Default description is file name
	if 'description' not in attrs:
		attrs['description'] = info.wo_ext

	return attrs


genome_import_attrs = OrderedDict([
	('description', import_str),
	('gb_db', import_str),
	('gb_id', int),
	('gb_acc', import_str),
	('is_assembled', parse_bool),
	('organism', import_str),
	('file_format', lambda value: import_str(value, lower=True)),
])

reqd_genome_import_attrs = ['description', 'is_assembled', 'file_format']
uq_genome_import_attrs = ['description', 'gb_id', 'gb_acc']

genome_show_attrs = ['id'] + list(genome_import_attrs.keys()) + [
	'created_at',
	'updated_at',
]

genome_import_cols = list(genome_import_attrs.keys())
genome_import_cols.append('file')
genome_import_cols.append('compression')


ImportFileItem = namedtuple('ImportFileItem', ['path', 'compression', 'attrs'])

def parse_import_csv(fh, db):
	from wgskmers.database.models import Genome

	session = db.get_session()

	return_vals = []
	uq_vals = dict()

	# Read in template file

	reader = DictReader(fh)

	# Check missing columns
	missing_cols = set(genome_import_cols).difference(reader.fieldnames)
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

		# Check compression
		compression = import_str(row['compression'], lower=True)
		if compression not in [None, 'gzip']:
			raise click.ClickException('Unknown compression type "{}"'
			                           .format(compression))

		# Get attributes from row
		attrs = dict()
		for attrname, converter in genome_import_attrs.items():
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
		for attrname in reqd_genome_import_attrs:
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
		attrs_ok = True
		for uq_col in uq_genome_import_attrs:
			val = attrs[uq_col]
			if val is None:
				continue

			# Check already in database
			found = (session.query(Genome)
			                .filter_by(**{uq_col: val})
			                .first())
			if found is not None:
				click.echo(
					err_prefix + 
					'Reference genome already exists in database with '
					'{} "{}" - skipping'
					.format(uq_col, val),
					err=True
				)
				attrs_ok = False
				continue

			# Check already in file
			seen_vals = uq_vals.setdefault(uq_col, set())
			if val in seen_vals:
				click.echo(
					err_prefix +  'Duplicate value {} for column {} - skipping'
					.format(uq_col, val),
					err=True
				)
				attrs_ok = False
				continue

			seen_vals.add(val)

		# Add if ok
		if attrs_ok:
			return_vals.append(ImportFileItem(path=path, attrs=attrs,
			                                  compression=compression))

	return return_vals


@click.group(name='gen', short_help='Manage reference genomes')
@choose_db()
def genomes_group():
	"""Commands to manage reference genomes in """
	pass


@genomes_group.command()
@click.option('-c', '--csv', 'out_csv', is_flag=True)
@click.argument('dest', type=click.File('w'), default='-')
@with_db()
def list(ctx, db, dest, out_csv=False):
	"""List reference genomes"""
	from wgskmers.database.models import Genome

	session = db.get_session()
	genomes = session.query(Genome).all()

	if out_csv:

		writer = DictWriter(dest, genome_show_attrs)
		writer.writeheader()

		for g in genomes:
			writer.writerow({c: getattr(g, c) for c in genome_show_attrs})

	else:
		for g in genomes:
			click.echo(g.description)


@genomes_group.command(short_help='List genome sets')
@click.option('-c', '--csv', 'out_csv', is_flag=True)
@click.argument('dest', type=click.File('w'), default='-')
@with_db()
def list_sets(ctx, db, dest, out_csv=False):
	from wgskmers.database.models import GenomeSet

	session = db.get_session()
	gsets = session.query(GenomeSet).all()

	if out_csv:

		writer = DictWriter(dest, ['id', 'name', 'genome_count'])
		writer.writeheader()

		for gset in gsets:
			writer.writerow(dict(id=gset.id, name=gset.name,
			                genome_count=gset.genomes.count()))

	else:
		for gset in gsets:
			click.echo('({0.id}) {0.name}'.format(gset))


@genomes_group.command(short_help='Create genome set')
@click.argument('name')
@click.argument('description', required=False)
@with_db(confirm=True)
def make_set(ctx, db, name, description=None):
	from wgskmers.database.models import GenomeSet

	session = db.get_session()

	if session.query(GenomeSet).filter_by(name=name).count() > 0:
		raise click.ClickException('Genome set already exists with name "{}"'
		                           .format(name))

	gset = GenomeSet(name=name, description=description)

	session.add(gset)
	session.commit()

	click.echo('Genome set "{0.name}" created with ID {0.id}'.format(gset))


@genomes_group.command(name='import', short_help='Import reference genomes')
@click.option('-c', '--csv', 'csv_out', type=click.Path(),
	help='Path to write import csv file to')
@click.option('-e', '--existing', type=click.File(),
	help='Existing import template .csv file to use')
@click.option('-s', '--gset', 'gset_id', type=int,
              help='ID of genome set to add to (check the list_sets command)')
@click.option('-r', '--recursive', is_flag=True,
              help='Search the directory for files recursively')
@click.option('--keep/--no-keep', default=True,
              help='Keep original files after they have been imported')
@click.option('--check-contents/--no-check-contents', default=True,
              help='Perform basic verification of file contents')
@click.argument('directory', type=click.Path(exists=True))
@with_db(confirm=True)
def import_genomes(ctx, db, directory, **kwargs):
	from wgskmers.database.models import GenomeSet
	from wgskmers.parse import find_seq_files

	csv_out = kwargs.pop('csv_out')
	existing = kwargs.pop('existing')
	gset_id = kwargs.pop('gset_id')
	recursive = kwargs.pop('recursive')
	keep = kwargs.pop('keep')
	check_contents = kwargs.pop('check_contents')
	assert not kwargs

	directory = os.path.abspath(directory)

	# Get genome set
	if gset_id is not None:
		session = db.get_session()
		gset = session.query(GenomeSet).get(gset_id)
		session.close()

		if gset is None:
			raise click.ClickException('No genome set with ID {}'.format(gset_id))

		gsets = [gset]

	else:
		gsets = []
		
	# Create a new import template file
	if existing is None:

		if csv_out is None:
			csv_out = os.path.join(directory, 'genomes-import.csv')

		# Find fasta files in directory
		directory = os.path.abspath(directory)
		files_info = find_seq_files(
			directory,
			filter_ext=True,
			filter_contents=True,
			check_contents=check_contents,
			warn_contents=check_contents,
			recursive=recursive,
			allow_compressed=True,
			tqdm=dict(desc='Finding sequence files')
		)

		# Write import .csv template
		with open(csv_out, 'w') as fh:

			writer = DictWriter(fh, genome_import_cols)
			writer.writeheader()

			# Guess attributes
			for info in tqdm(files_info, desc='Checking files'):

				# Make sure it's a supported format
				if info.seq_format != 'fasta':
					continue

				attrs = guess_fasta_attrs(info)
				attrs['file'] = info.path
				attrs['compression'] = info.compression
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
		for item in tqdm(parse_import_csv(csv_fh, db)):

			store_kwargs = dict(
				genome_sets=gsets,
				compression='gzip',
				src_compression=item.compression,
				keep_src=keep,
			)
			store_kwargs.update(item.attrs)

			# Try adding it
			try:
				db.store_genome(item.path, **store_kwargs)
				added += 1
			except Exception as e:
				click.echo(
					'Error adding file {}: {}: {}'
					.format(item.path, type(e), e),
					err=True
				)
				errors += 1

	click.echo(
		'Successfully imported {} genomes, with {} errors'
		.format(added, errors)
	)
