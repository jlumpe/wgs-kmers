"""Commands for managing registered databases"""

import os
import string

import click

from .util import choose_db_path


def check_valid_name(name):
	"""Raises a ClickException if the name is not a valid database name"""
	if not set(string.lowercase + string.digits + '_').issuperset(name):
		raise click.ClickException('Name characters must be alphanumeric or '
		                           'underscores')


@click.group(name='db')
def database_group():
	"""Manage databases"""
	pass


@database_group.command()
def list():
	"""List registered databases"""
	from wgskmers.config import get_config
	from wgskmers import database

	# Get config databases
	config = get_config()
	registered_dbs = config.get_registered_dbs()
	default_db = config.get_default_db()

	# Get current database
	current_path, method = database.get_current_db()

	# List default database first, remaining alphabetically
	databases = registered_dbs.items()
	if default_db is not None:
		databases.insert(0, (None, default_db))

	# None found
	if not databases:
		click.echo('No databases are currently registered', err=True)
		return

	# Print names and paths
	for name, path in databases:

		# Default with brackets
		print_name = '[default]' if name is None else name

		# Left-justify
		print_name = print_name.ljust(20)

		# Check if missing
		if not database.is_db_directory(path):
			name_str = click.style(print_name, fg='red', bold=True)
			path_str = click.style(path, fg='red')

		else:
			path_str = path

			# Check if current
			if current_path is not None and path == current_path:
				name_str = click.style(print_name, fg='cyan', bold=True)
			elif name is None:
				name_str = click.style(print_name, fg='green')
			else:
				name_str = print_name

		click.echo('{} {}'.format(name_str, path_str))


@database_group.command()
def which():
	"""Get the currently active database"""
	from wgskmers import database

	path, method = database.get_current_db()

	if path is not None:

		if database.is_db_directory(path):
			click.echo(path)
		else:
			click.secho(path, fg='red')

		if method == 'cwd':
			click.echo('(Current working directory)')
		elif method == 'config':
			click.echo('(Default database in config file)')
		elif method == 'environ':
			click.echo(
				'(From environment variable {})'
				.format(database.DEFAULT_DB_PATH_VAR)
			)
		elif method == 'override':
			click.echo('(Overridden in shell)')
		else:
			click.echo('(huh?)')

	else:
		click.echo('No database currently active (not in a database directory '
		           'and no default set in config file or environment variable).')


@database_group.command(short_help='Create a new database')
@click.option('-d', '--set-default', is_flag=True,
              help='Set as default database')
@click.option('-n', '--name', type=str, help='Register database by name')
@click.argument('directory', type=click.Path(), required=False)
def init(directory, set_default=False, name=None):
	"""
	Create a new database in the specified directory, or the current one if
	none specified.
	"""
	from wgskmers.config import get_config

	config = get_config()

	if directory is None:
		directory = os.getcwd()
	else:
		directory = os.path.abspath(directory)

	# Check if exists
	if os.path.exists(directory) and os.listdir(directory):
		raise click.ClickException('Directory exists and is not emptry, aborting')

	# Confirm if default exists
	if set_default:
		current_default = config.get_default_db()
		if current_default is not None:
			click.confirm(
				'Default database is currently set to {}, overwrite?'
				.format(current_default),
				abort=True
			)

	# Check name
	if name is not None:
		name = name.lower()
		check_valid_name(name)

		# Confirm if name exists
		current_named = config.get_registered_dbs().get(name, None)
		if current_named is not None:
			click.confirm(
				'Database in {} is already registered as {}, overwrite?'
				.format(current_named, name),
				abort=True
			)

	# Create database
	database.Database.create(directory)

	# Register
	if name is not None:
		config.register_db(directory, name, True)
	if set_default:
		config.set_default_db(directory, True)

	click.echo('Success!')


@database_group.command(short_help='Set the default database')
@click.option('-n', '--name', help='Name of currently registered database')
@click.argument('path', type=click.Path(exists=True), required=False)
def set_default(path=None, name=None):
	"""
	Sets the default database by either specifying a path or giving the name
	of another currently registered database. If no path given, will check if
	current working directory is within a database.
	"""
	from wgskmers.config import get_config
	from wgskmers import database

	config = get_config()
	registered_dbs = config.get_registered_dbs()

	# Set from explicit path
	if path is not None:
		path = os.path.abspath(path)

		if name is not None:
			raise click.ClickException('Can\'t give both name and path.')

		if not database.is_db_directory(path):
			raise click.ClickException(
				'{} does not contain a database.'
				.format(path)
			)

	# Set from another registered database
	elif name is not None:
		try:
			path = registered_dbs[name]
		except KeyError:
			raise click.ClickException('No database registered with name {}'
			                     .format(name))

	else:
		path = database.get_db_root(os.getcwd())
		if path is None:
			raise click.ClickException('Not currently in a database directory')


	# Check if database already exists
	current_default = config.get_default_db()
	if current_default is not None:

		# Check if it is already the default
		if current_default == path:
			click.echo('The default database already is {}'
			           .format(current_default))
			return

		else:
			click.confirm(
				'The default database is currently set to {}, are you sure '
				'you wish to overwrite it?'
				.format(current_default),
				abort=True
			)

	# Register it
	config.set_default_db(path, True)
	click.echo('Default database set to {}'.format(path))


@database_group.command(short_help='Register a database by name')
@click.argument('name', type=str)
@click.argument('path', type=click.Path(exists=True), required=False)
def register(name, path=None):
	"""
	Register an existing database in the global configuration file. If no
	path given, will check if current working directory is within a database.
	"""
	from wgskmers.config import get_config
	from wgskmers import database

	config = get_config()
	registered_dbs = config.get_registered_dbs()

	# Check name valid
	name = name.lower()
	check_valid_name(name)

	# Can't set default this way
	if name == 'default':
		click.echo('Use set_default to set the default database')
		return

	# Set from explicit path
	if path is not None:
		path = os.path.abspath(path)

		if not database.is_db_directory(path):
			raise click.ClickException(
				'{} does not contain a database.'
				.format(path)
			)

	else:
		path = database.get_db_root(os.getcwd())
		if path is None:
			raise click.ClickException('Not currently in a database directory')


	# Check if database already exists
	if name in registered_dbs:
		current_path = registered_dbs[name]

		# Check if it is already registered
		if current_path == path:
			click.echo('Database at {} is already registered as {}'
			           .format(current_path, name))
			return

		else:
			click.confirm(
				'Database at {} is currently registered as {}, are you sure '
				'you wish to overwrite it?'
				.format(name, current_path),
				abort=True
			)

	# Register it
	config.register_db(path, name, True)
	click.echo('Database {} registered as {}'.format(path, name))


@database_group.command(short_help='Remove a registered database')
@click.argument('name', type=str)
def unregister(name):
	"""
	Un-register a currently registered database. The database itself is left
	alone. Type "default" to un-register the default database.
	"""
	from wgskmers.config import get_config

	config = get_config()
	registered_dbs = config.get_registered_dbs()

	if name not in registered_dbs:
		raise click.ClickException('Database "{}" is not currently registered.'
		                           .format(name))

	# Confirm
	if name == 'default':
		click.confirm('Are you sure you wish to unset the default database?',
		              abort=True)
	else:
		click.confirm('Are you sure you wish to unregister database "{}"?'
		              .format(name),
		              abort=True)

	# Remove it
	config.unregister_db(name)
	if name =='default':
		click.echo('Default database removed from config file.')
	else:
		click.echo('Database "{}" removed from config file.'.format(name))


@database_group.command()
@choose_db_path()
def upgrade(db_path):
	"""Upgrade a database to the latest version"""
	from wgskmers.database.upgrade import upgrader

	upgrader.upgrade(db_path)
	click.echo('Success!')


@database_group.command()
@choose_db_path()
def version(db_path):
	"""Get version of existing database"""
	from wgskmers import database

	click.echo(database.get_db_version(db_path))
