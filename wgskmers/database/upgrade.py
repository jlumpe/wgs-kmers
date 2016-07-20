"""Functions for performing database upgrades"""
from builtins import range
from builtins import object

import os

from .database import CURRENT_DB_VERSION, get_db_version, set_db_version
from .migrate import get_alembic_config, migrate


def get_sqlite_path(db_path):
	"""Gets sqlite path for database, including for previous versions"""

	# This will probably be more complicated in the future...
	return os.path.join(db_path, 'data.db')


class DatabaseUpgrader(object):
	"""Stores and applies upgrade scripts for databases"""

	def __init__(self):
		self.scripts = dict()

	def script(self, from_ver):
		def decorator(func):
			self.scripts[from_ver] = func
			return func
		return decorator

	def revision_script(self, from_ver, revision):
		def script(db_path, db_ver):
			cfg = get_alembic_config(os.path.join(db_path, 'data.db'))
			migrate(cfg, revision)
		self.scripts[from_ver] = script

	def upgrade(self, db_path):
		"""Upgrades database to the newest version"""
		db_ver = get_db_version(db_path)
		for v in range(db_ver, CURRENT_DB_VERSION):
			self.scripts[v](db_path, v)
			set_db_version(db_path, v + 1)


upgrader = DatabaseUpgrader()

# Just alembic migrations
upgrader.revision_script(1, '704356629cab')
upgrader.revision_script(2, '41c9af002856')
upgrader.revision_script(3, '40c711d276f0')
upgrader.revision_script(4, '0a1c81836e60')
