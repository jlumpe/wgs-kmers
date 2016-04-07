"""This package defines console commands for interacting with the project."""

import logging

import click

from . import find


# Configure logger for warnings and debug messeges
logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger()


# Top-level cli group
@click.group()
@click.option('--debug', is_flag=True, default=False,
	help='Print debug messages')
def cli(debug=False):
	# Debug mode
	if debug:
		logger.setLevel(logging.DEBUG)


cli.add_command(find.find_command)
