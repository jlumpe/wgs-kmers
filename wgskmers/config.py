from future import standard_library
standard_library.install_aliases()
from past.builtins import basestring
from builtins import object
import os
import json
from configparser import RawConfigParser

from appdirs import AppDirs


app_dirs = AppDirs('wgskmers', version='0.2')
config_dir = app_dirs.user_config_dir


# Don't load until requested
_config = None


def config_exists():
	return os.path.isdir(config_dir)


def get_config():
	if _config is None:
		reload_config()
	return _config


def reload_config():
	global _config
	if not os.path.isdir(config_dir):
		os.makedirs(config_dir)
	_config = SystemConfig(config_dir)
	return _config


class SystemConfig(object):
	"""Manages system-wide configuration for kmers project"""

	_paths = {
		'databases': 'databases.json',
	}

	def __init__(self, directory):
		self.directory = os.path.abspath(directory)

		self._db_conf = None

	def _path_for(self, item, *args):
		return os.path.join(self.directory, self._paths[item], *args)

	def _get_db_conf(self):
		if self._db_conf is None:
			try:
				with open(self._path_for('databases')) as fh:
					self._db_conf = json.load(fh)

			except (IOError, ValueError):
				self._db_conf = dict()

		return self._db_conf

	def _save_db_conf(self, conf):
		self. db_conf = conf
		with open(self._path_for('databases'), 'w') as fh:
			json.dump(conf, fh)

	def get_default_db(self):
		"""Add database to configuration file"""
		return self._get_db_conf().get('default', None)

	def set_default_db(self, path, overwrite=False):
		"""Remove database from configuration file"""
		if not isinstance(path, basestring) and path is not None:
			raise TypeError(path)

		db_conf = self._get_db_conf()

		if db_conf.get('default', None) is not None and not overwrite:
			raise RuntimeError('Default database exists, refusing to overwrite')

		db_conf['default'] = path

		self._save_db_conf(db_conf)

	def get_registered_dbs(self):
		"""Get list of databases in config file"""
		return dict(self._get_db_conf().get('registered', {}))

	def register_db(self, path, name, overwrite=False):
		"""Add database to configuration file"""
		if not isinstance(path, basestring):
			raise TypeError(path)
		if not isinstance(name, basestring):
			raise TypeError(name)

		db_conf = self._get_db_conf()
		databases = db_conf.setdefault('registered', {})

		if name in databases and not overwrite:
			raise RuntimeError('Database {} exists, refusing to overwrite')

		databases[name] = path
		self._save_db_conf(db_conf)

	def unregister_db(self, name):
		"""Remove database from configuration file"""
		db_conf = self._get_db_conf()
		databases = db_conf.setdefault('registered', {})

		if not name in databases:
			raise RuntimeError('Database {} not registered')

		del databases[name]
		self._save_db_conf(db_conf)
