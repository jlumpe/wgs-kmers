"""Command line tools for development"""

import os

import click
from pprint import pprint

from alembic import command as alembic_command
from alembic.migration import MigrationContext
from alembic.autogenerate import compare_metadata

import wgskmers
from wgskmers import database
from wgskmers.database.upgrade import get_alembic_config, get_sqlite_path
from .util import with_db, choose_db_path


@click.group(name='dev')
def dev_group():
	"""Developer commands"""
	pass


@dev_group.command(short_help='Debug shell')
@click.option('--db', type=str)
def shell(db=None):
	"""Open up an IPython shell with modules imported"""
	import IPython

	import wgskmers
	from wgskmers import kmers, database, config, util
	from wgskmers.database import models

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

	if db is not None:
		ns['db'] = database.Database(database.get_registered_dbs()[db])
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
	src_path = os.path.dirname(wgskmers.__file__)

	for dirpath, dirname, filenames in os.walk(src_path):
		for filename in filenames:
			if filename.endswith('.pyc'):
				if clean_all or (filename[-4:] + '.py') in filenames:
					file_path = os.path.join(dirpath, filename)
					click.echo('removed ' + file_path)
					os.unlink(file_path)


@dev_group.group(short_help='Alembic database migration tools')
def alembic():
	"""Alembic database migration tools"""
	pass


@alembic.command()
@choose_db_path()
def diff(ctx, db):
	"""Show database diff for current SQLA models"""
	cfg = get_alembic_config(get_sqlite_path(db_path))
	engine = db.engine

	mc = MigrationContext.configure(engine.connect())
	diff = compare_metadata(mc, database.Base.metadata)

	print pprint(diff)


@alembic.command()
@click.option('--message', '-m', type=str,
              help='Message string to use with revision')
@click.option('--sql', is_flag=True,
              help='Don\'t emit SQL to database - dump to standard output/file '
                   'instead')
def revision(**kwargs):
	"""Create a new revision file"""
	cfg = get_alembic_config()
	alembic_command.revision(cfg, **kwargs)


@alembic.command()
@click.option('--message', '-m', type=str,
              help='Message string to use with revision')
@click.option('--sql', is_flag=True,
              help='Don\'t emit SQL to database - dump to standard output/file '
                   'instead')
@choose_db_path()
def auto_revision(db_path, **kwargs):
	"""Create a new revision file"""
	cfg = get_alembic_config(get_sqlite_path(db_path))
	alembic_command.revision(cfg, autogenerate=True, **kwargs)


@alembic.command()
@click.option('--verbose', '-v', is_flag=True, help='Use more verbose output')
@choose_db_path()
def current(db_path, **kwargs):
	"""Display the current revision for a database."""
	cfg = get_alembic_config(get_sqlite_path(db_path))
	alembic_command.current(cfg, **kwargs)


@alembic.command()
@click.option('--resolve_dependencies', is_flag=True,
	help='Treat dependency versions as down revisions')
@click.option('--verbose', '-v', is_flag=True, help='Use more verbose output')
def heads(**kwargs):
	"""Show current available heads in the script directory."""
	cfg = get_alembic_config()
	alembic_command.heads(cfg, **kwargs)


@alembic.command()
@click.option('--verbose', '-v', is_flag=True, help='Use more verbose output')
def branches(**kwargs):
	"""Show current branch points."""
	cfg = get_alembic_config()
	alembic_command.branches(cfg, **kwargs)


@alembic.command()
@click.argument('rev')
def show(rev):
	"""Show the revision(s) denoted by the given symbol."""
	cfg = get_alembic_config()
	alembic_command.show(cfg, rev)


@alembic.command()
@click.option('--sql', is_flag=True,
              help='Don\'t emit SQL to database - dump to standard output/file '
                   'instead')
@click.option('--tag', type=str,
              help='Arbitrary \'tag\' name - can be used by custom env.py '
                   'scripts')
@click.argument('revision')
@choose_db_path()
def stamp(db_path, revision, **kwargs):
	"""'stamp' the revision table with the given revision, don't run any
	migrations.
	"""
	cfg = get_alembic_config(get_sqlite_path(db_path))
	alembic_command.stamp(cfg, revision, **kwargs)


@alembic.command()
@click.option('--sql', is_flag=True,
              help='Don\'t emit SQL to database - dump to standard output/file '
                   'instead')
@click.option('--tag', type=str,
              help='Arbitrary \'tag\' name - can be used by custom env.py '
                   'scripts')
@click.argument('revision')
@choose_db_path()
def upgrade(db_path, revision, **kwargs):
	"""Upgrade database to revision"""
	cfg = get_alembic_config(get_sqlite_path(db_path))
	alembic_command.upgrade(cfg, revision, **kwargs)


@alembic.command()
@click.option('--sql', is_flag=True,
              help='Don\'t emit SQL to database - dump to standard output/file '
                   'instead')
@click.option('--tag', type=str,
              help='Arbitrary \'tag\' name - can be used by custom env.py '
                   'scripts')
@click.argument('revision')
@choose_db_path()
def downgrade(db_path, revision, **kwargs):
	"""Upgrade database to revision"""
	cfg = get_alembic_config(get_sqlite_path(db_path))
	alembic_command.downgrade(cfg, revision, **kwargs)

