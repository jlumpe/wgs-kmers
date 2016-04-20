"""This package defines console commands for interacting with the project."""

import logging

import click

from .find import find_command
from .config import config_group
from .database import database_group
from .genomes import genomes_group
from .kmers import kmers_group
from .query import query_command
from .dev import dev_group


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
cli.add_command(kmers_group)
cli.add_command(query_command)
cli.add_command(dev_group)
