"""Commands for editing configuration"""

import shutil

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
		click.echo(config.config_dir)
	else:
		click.secho('Config directory does not exist', err=True)


@config_group.command()
def reset():
	"""Reset the global config file(s)"""
	if not config.config_exists():
		click.secho('Config directory does not exist', err=True)
		return

	click.confirm(
		click.style('WARNING', fg='red', bold=True) +
		' - this will clear all configuration, including registered datbases. '
		'The databases themselves will be left as-is. Continue?',
		abort=True
	)
	shutil.rmtree(config.config_dir)
	click.echo('Configuration reset')


@config_group.command(name='open')
def open_config():
	"""Open the config directory"""
	if config.config_exists():
		click.launch(config.config_dir)
	else:
		click.secho('Config directory does not exist', err=True)
