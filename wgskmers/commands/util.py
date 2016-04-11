"""Utilites for console commands"""

import os
import itertools
from functools import wraps
from textwrap import dedent

import click
from Bio import SeqIO
from tqdm import tqdm

from wgskmers import database


class ProgressSeqParser(object):
	"""Wraps generator from Bio.SeqIO.parse with tqdm progress bar."""

	def __init__(self, file_, fmt='fasta', **kwargs):
		self.file_ = file_
		self.fmt = fmt
		self.tqdm_args = kwargs

	def __iter__(self):

		# Open file if given as path
		if isinstance(self.file_, basestring):
			fh = open(self.file_)
			opened = True
		else:
			fh = self.file_
			opened = False

		try:
			# Starting position in file, probably zero but maybe not...
			last = fh.tell()

			# Size of file to parse
			total = os.fstat(fh.fileno()).st_size - last

			# Create progress bar
			pbar = tqdm(unit='B', unit_scale=True, total=total,
			            **self.tqdm_args)
			with pbar:

				# Parse and iterate over records
				for record in SeqIO.parse(fh, self.fmt):

					# Update progress bar
					current = fh.tell()
					pbar.update(current - last)
					last = current

					# Yield parsed record
					yield record

		# Close the file if we opened it
		finally:
			if opened:
				fh.close()


def iterator_empty(iterator):
	"""Checks if an iterator is empty, also returning a substitute iterator.

	Consumes first element of iterator, but stores it for the substitute
	iterator to yield first.
	"""

	try:
		first = [next(iterator)]
		empty = False
	except StopIteration:
		first = []
		empty = True

	substitute = itertools.chain(first, iterator)

	return substitute, empty


choose_db_doc = dedent("""

	Specify database to use explicitly with "--db" or "--default" commands.
	If neither is given will use either database in the current working
	directory or the one set as default, with confirmation for write ops.
""")


def choose_db(pass_opts=False, pass_context=False):
	"""Decorator for a command/group that chooses a database to work with

	Adds options to specify name of database to work with (or to use default).
	If none specified, will use the current default database. Will open
	database and add to context dict under key ['db']. Will also use the method
	for automatically picking the database under 'db_auto' if appropriate.
	"""
	def decorator(func):

		@wraps(func)
		def wrapper(ctx, *args, **kwargs):

			default_db = kwargs.get('default_db', None)
			db_name = kwargs.get('db_name', None)

			# Arguments for wrapped func
			if not pass_opts:
				kwargs.pop('default_db', None)
				kwargs.pop('db_name', None)

			if pass_context:
				args = tuple([ctx] + list(args))

			# Use default database
			if default_db:
				if db_name is not None:
					raise click.ClickException(
						'Cant use both --db and --default options')

				db_path = database.get_default_db()
				if db_path is None:
					raise click.ClickException('No default database set')

				ctx.obj['db_auto'] = False

			# Use registered database by name
			elif db_name is not None:
				db_name = db_name.lower()
				try:
					db_path = database.get_registered_dbs()[db_name]
				except KeyError:
					raise click.ClickException(
						'Database "{}" not found'
						.format(db_name)
					)

				ctx.obj['db_auto'] = False

			# Determine automatically
			else:
				db_path, method = database.get_current_db()
				if db_path is None:
					raise click.ClickException('No database currently active')

				ctx.obj['db_auto'] = method

			# Open database and add to context
			ctx.obj['db'] = database.Database.open(db_path)

			# Call wrapped func
			return func(*args, **kwargs)

		if wrapper.__doc__ is not None:
			wrapper.__doc__ += choose_db_doc

		# Add click decorators for options and context passing
		decorators = [
			click.option('-N', '--db', 'db_name', type=str, metavar='NAME',
			             help='Use registered database by name'),
			click.option('-d', '--default', 'default_db', is_flag=True,
			             help='Use default database'),
			click.pass_context,
		]
		for d in decorators:
			wrapper = d(wrapper)

		return wrapper

	return decorator


def with_db(confirm=False, choose=False):
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

		# Add context
		wrapped = click.pass_context(wrapper)

		# Choose if necessary
		if choose:
			wrapped = choose_db(pass_opts=False, pass_context=False)(wrapped)

		return wrapped

	# Return decorator
	return decorator
