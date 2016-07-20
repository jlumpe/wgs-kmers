"""Command line tools for development"""
from __future__ import print_function

import os
from functools import wraps

import click

from .util import with_db, choose_db_path


@click.group(name='dev')
def dev_group():
	"""Developer commands"""
	pass


@dev_group.command(short_help='Debug shell')
@click.option('--db', type=str)
@click.option('--no-db', is_flag=True, default=False)
def shell(db=None, no_db=False):
	"""Open up an IPython shell with modules imported"""

	import IPython
	from wgskmers.database import open_database
	from wgskmers.config import get_config

	import wgskmers
	from wgskmers import kmers, database, util
	from wgskmers.database import models

	config = get_config()

	ns = dict(
		wgskmers=wgskmers,
		kmers=kmers,
		models=models,
		database=database,
		config=config,
		util=util,
	)

	for name in models.__all__:
		ns[name] = getattr(models, name)

	if not no_db:

		if db is None:
			db_path = config.get_default_db()
		else:
			db_path = config.get_registered_dbs()[db]

		if db_path is not None:
			ns['db'] = open_database(db_path)
			ns['session'] = ns['db'].get_session()

	IPython.start_ipython(argv=[], user_ns=ns)


@dev_group.command(short_help='Remove .pyc files in project')
@click.option('-a', '--all', 'clean_all', is_flag=True,
              help='Remove all .pyc files instead of just those without a '
                   'corresponding .py file')
def clean(clean_all=False):
	"""
	Remove all .pyc files in project, by default only if they lack a
	corresponding source file.
	"""
	import wgskmers

	src_path = os.path.dirname(wgskmers.__file__)

	for dirpath, dirname, filenames in os.walk(src_path):
		for filename in filenames:
			if filename.endswith('.pyc'):
				if clean_all or (filename[-4:] + '.py') in filenames:
					file_path = os.path.join(dirpath, filename)
					click.echo('removed ' + file_path)
					os.unlink(file_path)


@dev_group.group(name='alembic', short_help='Alembic database migration tools')
def alembic_group():
	"""Alembic database migration tools"""
	pass


def with_alembic_config(use_db=False):
	"""Creates decorator for commands needing alembic config"""

	from wgskmers.database.upgrade import get_alembic_config, get_sqlite_path

	def decorator(func):

		if use_db:
			@wraps(func)
			def wrapper(db_path, *args, **kwargs):
				cfg = get_alembic_config(get_sqlite_path(db_path))
				return func(cfg, *args, **kwargs)

			return choose_db_path()(wrapper)

		else:
			@wraps(func)
			def wrapper(*args, **kwargs):
				cfg = get_alembic_config()
				return func(cfg, *args, **kwargs)

			return wrapper

	return decorator


@alembic_group.command()
@with_db(choose=True)
def diff(ctx, db):
	"""Show database diff for current SQLA models"""

	from pprint import pprint
	from alembic.migration import MigrationContext
	from alembic.autogenerate import compare_metadata
	from wgskmers.database.models import Base
	from wgskmers.database.upgrade import get_alembic_config, get_sqlite_path

	cfg = get_alembic_config(get_sqlite_path(db.directory))
	engine = db.engine

	mc = MigrationContext.configure(engine.connect())
	diff = compare_metadata(mc, Base.metadata)

	print(pprint(diff))


@alembic_group.command()
@click.option('--message', '-m', type=str,
              help='Message string to use with revision')
@click.option('--sql', is_flag=True,
              help='Don\'t emit SQL to database - dump to standard output/file '
                   'instead')
@with_alembic_config()
def revision(cfg, **kwargs):
	"""Create a new revision file"""
	from alembic import command as alembic_command
	alembic_command.revision(cfg, **kwargs)


@alembic_group.command()
@click.option('--message', '-m', type=str,
              help='Message string to use with revision')
@click.option('--sql', is_flag=True,
              help='Don\'t emit SQL to database - dump to standard output/file '
                   'instead')
@with_alembic_config(use_db=True)
def auto_revision(cfg, **kwargs):
	"""Create a new revision file"""
	from alembic import command as alembic_command
	alembic_command.revision(cfg, autogenerate=True, **kwargs)


@alembic_group.command()
@click.option('--verbose', '-v', is_flag=True, help='Use more verbose output')
@with_alembic_config(use_db=True)
def current(cfg, **kwargs):
	"""Display the current revision for a database."""
	from alembic import command as alembic_command
	alembic_command.current(cfg, **kwargs)


@alembic_group.command()
@click.option('--resolve_dependencies', is_flag=True,
	help='Treat dependency versions as down revisions')
@click.option('--verbose', '-v', is_flag=True, help='Use more verbose output')
@with_alembic_config()
def heads(cfg, **kwargs):
	"""Show current available heads in the script directory."""
	from alembic import command as alembic_command
	alembic_command.heads(cfg, **kwargs)


@alembic_group.command()
@click.option('--verbose', '-v', is_flag=True, help='Use more verbose output')
@with_alembic_config()
def branches(cfg, **kwargs):
	"""Show current branch points."""
	from alembic import command as alembic_command
	alembic_command.branches(cfg, **kwargs)


@alembic_group.command()
@click.argument('rev')
@with_alembic_config()
def show(cfg, rev):
	"""Show the revision(s) denoted by the given symbol."""
	from alembic import command as alembic_command
	alembic_command.show(cfg, rev)


@alembic_group.command()
@click.option('--sql', is_flag=True,
              help='Don\'t emit SQL to database - dump to standard output/file '
                   'instead')
@click.option('--tag', type=str,
              help='Arbitrary \'tag\' name - can be used by custom env.py '
                   'scripts')
@click.argument('revision')
@with_alembic_config(use_db=True)
def stamp(cfg, revision, **kwargs):
	"""'stamp' the revision table with the given revision, don't run any
	migrations.
	"""
	from alembic import command as alembic_command
	alembic_command.stamp(cfg, revision, **kwargs)


@alembic_group.command()
@click.option('--sql', is_flag=True,
              help='Don\'t emit SQL to database - dump to standard output/file '
                   'instead')
@click.option('--tag', type=str,
              help='Arbitrary \'tag\' name - can be used by custom env.py '
                   'scripts')
@click.argument('revision')
@with_alembic_config(use_db=True)
def upgrade(cfg, revision, **kwargs):
	"""Upgrade database to revision"""
	from alembic import command as alembic_command
	alembic_command.upgrade(cfg, revision, **kwargs)


@alembic_group.command()
@click.option('--sql', is_flag=True,
              help='Don\'t emit SQL to database - dump to standard output/file '
                   'instead')
@click.option('--tag', type=str,
              help='Arbitrary \'tag\' name - can be used by custom env.py '
                   'scripts')
@click.argument('revision')
@with_alembic_config(use_db=True)
def downgrade(cfg, revision, **kwargs):
	"""Upgrade database to revision"""
	from alembic import command as alembic_command
	alembic_command.downgrade(cfg, revision, **kwargs)

