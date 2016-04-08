""""""

import os

import click

from wgskmers import config



@click.group(name='config')
def config_group():
	"""View and edit global configuration"""
	pass


@config_group.command()
def locate():
	"""Print config directory path"""
	if config.config_exists():
		click.echo(config.config_file_path)
	else:
		click.secho('Config file does not exist', fg='red')


@config_group.command()
def reset():
	"""Reset the global config file"""
	if not config.config_exists():
		click.secho('Config file does not exist', fg='red')
		return

	click.confirm(
		click.style('WARNING', fg='red', bold=True) +
		' - this will clear all configuration, including registered datbases. '
		'The databases themselves will be left as-is. Continue?',
		abort=True
	)
	os.unlink(config.config_file_path)
	click.echo('Config file reset')


@config_group.command(name='print')
def print_config():
	"""Print contents of the config file"""
	if config.config_exists():
		with open(config.config_file_path) as fh:
			click.echo(fh.read())
	else:
		click.secho('Config file does not exist', fg='red')


@config_group.command(short_help='Open the config file for editing')
def manual_edit():
	"""Open the config file for editing. Use this at your own risk!"""
	click.edit(filename=config.config_file_path)
