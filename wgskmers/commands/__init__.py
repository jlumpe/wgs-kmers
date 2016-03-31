"""This package defines console commands for interacting with the project."""

import logging

import click

from .find import find_command
from .config import config_group
from .database import database_group
from .genomes import genomes_group


# Configure logger for warnings and debug messeges
logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger()


# Top-level cli group
@click.group()
@click.option('--debug', is_flag=True, default=False,
              help='Print debug messages')
@click.pass_context
def cli(ctx, debug=False):
	# Context object as dictionary
	ctx.obj = dict()

	# Debug mode
	if debug:
		logger.setLevel(logging.DEBUG)


cli.add_command(find_command)
cli.add_command(config_group)
cli.add_command(database_group)
cli.add_command(genomes_group)


@cli.command(short_help='Debug shell')
@click.option('--db', type=str)
def shell(db=None):
	"""Open up an IPython shell with modules imported"""
	import IPython

	import wgskmers
	from wgskmers import kmers, models, database, config, util

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

	IPython.start_ipython([], user_ns=ns)
