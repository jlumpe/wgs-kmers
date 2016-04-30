"""Utilites for console commands"""

import os
from functools import wraps
from textwrap import dedent

import click

from wgskmers import database


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


def choose_db_path(pass_opts=False):
	"""Decorator for command/group that picks database path (but doesn't open)

	Used for upgrade commands, etc where database might be an older version.
	"""

	def decorator(func):

		@wraps(func)
		def wrapper(*args, **kwargs):

			path = kwargs.get('path', None)
			db_name = kwargs.get('db_name', None)
			default_db = kwargs.get('default_db', None)

			if not pass_opts:
				kwargs.pop('path', None)
				kwargs.pop('db_name', None)
				kwargs.pop('default_db', None)

			if path is not None:
				db_path = os.path.abspath(path)

			elif db_name is not None:
				try:
					db_path = database.get_registered_dbs()[db_name]
				except KeyError:
					raise click.ClickException(
						'No database registered with name {}'
						.format(name)
					)

			elif default_db:
				db_path = database.get_default_db()
				if db_path is None:
					raise click.ClickException('No default database set')

			else:
				db_path, method = database.get_current_db()
				if db_path is None:
					raise click.ClickException('No database currently active')
				else:
					click.confirm('Using database at {}, ok?'.format(db_path),
					              abort=True)

			func(db_path, *args, **kwargs)

		decorators = [
			click.option('-N', '--db', 'db_name', type=str, metavar='NAME',
			             help='Use registered database by name'),
			click.option('-d', '--default', 'default_db', is_flag=True,
			             help='Use default database'),
			click.option('-P', '--path', type=str, metavar='PATH',
			             help='Specify database by path'),
		]
		for d in decorators:
			wrapper = d(wrapper)

		return wrapper

	return decorator
