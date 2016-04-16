"""Alembic migrations"""

from pkg_resources import resource_filename
from alembic.config import Config as AlembicConfig
import alembic.command as alembic_command


def get_alembic_config(sqlite_path=None):
	"""Create Alembic Config object for database (for migrations)"""

	cfg = AlembicConfig(resource_filename('wgskmers.database', 'alembic.ini'))
	cfg.set_main_option('script_location',
	                    resource_filename('wgskmers.database', 'alembic'))

	if sqlite_path is not None:
		cfg.set_main_option('sqlalchemy.url', 'sqlite:///' + sqlite_path)

	return cfg


def migrate(config, revision):
	"""Runs an alembic migration"""
	alembic_command.upgrade(config, revision)
