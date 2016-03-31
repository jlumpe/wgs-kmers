import os
import shutil
import json
import re

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from .models import *
from .util import rmpath
from .sqla import ReadOnlySession
from .config import config, save_config


# Current database version number
DB_VERSION = 1


# This environment variable overrides all others to set current database
DB_PATH_OVERRIDE_VAR = 'KMERS_DB_PATH_OVERRIDE'


# Environment variable setting default database path
DEFAULT_DB_PATH_VAR = 'KMERS_DEFAULT_DB_PATH'


# The presence of a file with this name indicates a directory contains
# a k-mer database. It also has basic information about the datbase in JSON
# format.
INFO_FILE_NAME = '.kmer-db'


def is_db_directory(path):
	"""Checks if a directory contains a k-mer database"""
	return os.path.isfile(os.path.join(path, INFO_FILE_NAME))


def get_db_root(path):
	"""Get which of a directory or its parents contains a database"""
	path = os.path.abspath(path)

	while True:
		if is_db_directory(path):
			return path

		parent = os.path.dirname(path)
		if parent == path:
			return None

		path = parent


def get_current_db(cwd=None):
	"""Gets path to database that should currently be used"""
	if cwd is None:
		cwd = os.getcwd()

	# Check if currently overridden by environment variable
	if DB_PATH_OVERRIDE_VAR in os.environ:
		return os.environ[DB_PATH_OVERRIDE_VAR], 'override'

	# Otherwise, start by looking in working directory and parents
	parent_dir = get_db_root(cwd)
	if parent_dir is not None:
		return parent_dir, 'cwd'

	# Check config file
	if config.has_option('databases', 'default'):
		return config.get('databases', 'default'), 'config'

	# Check default environment variable
	if DEFAULT_DB_PATH_VAR in os.environ:
		return os.environ[DEFAULT_DB_PATH_VAR], 'environ'

	# Not found
	return None, None


def get_default_db():
	"""Get database registered as default"""

	# Check config file
	if config.has_option('databases', 'default'):
		return config.get('databases', 'default')

	# Check default environment variable
	if DEFAULT_DB_PATH_VAR in os.environ:
		return os.environ[DEFAULT_DB_PATH_VAR]

	# Not found
	return None


def get_registered_dbs():
	"""Get list of databases in config file"""
	return {name: config.get('databases', name) for name
	        in config.options('databases')}


def register_db(path, name, overwrite=False):
	"""Add database to configuration file"""
	if config.has_option('databases', name) and not overwrite:
		raise RuntimeError('Database {} exists, refusing to overwrite')

	config.set('databases', name, path)
	save_config()


def unregister_db(name):
	"""Remove database from configuration file"""
	if not config.has_option('databases', name):
		raise RuntimeError('Database {} not registered')

	config.remove_option('databases', name)
	save_config()


class Database(object):

	_rel_paths = {
		'sqlite': 'data.db',
		'kmer_sets': 'kmer_sets',
		'genomes': 'genomes',
	}

	def __init__(self, directory):
		self.directory = directory

		# SqlAlchemy engine
		self._engine = create_engine('sqlite:///' + self._get_path('sqlite'))

		# SqlAlchemy session classes
		self._Session = sessionmaker(bind=self._engine)
		self._ExpireSession = sessionmaker(bind=self._engine,
		                                   expire_on_commit=False)

	def get_session(self):
		"""Create a new READ-ONLY SQLAlchemy session"""
		return self._Session()

	def store_genome(self, *args, **kwargs):

		# Pass file path to store, create genome from kwargs
		if len(args) == 1:
			file_path, = args
			genome = Genome(**kwargs)

		# Pass file path and Genome object
		elif len(args) == 2:
			file_path, genome = args
			if kwargs:
				raise TypeError('No kwargs taken if Genome instance passed')

		# Bad positional arguments
		else:
			raise TypeError('Function takes one or two positional arguments')

		session = self._ExpireSession()

		# Determine the file name
		ext = os.path.splitext(file_path)[1]
		genome.filename = self._make_genome_file_name(genome, session, ext)

		# Copy to directory
		store_path = self._get_path('genomes', genome.filename)
		if os.path.exists(store_path):
			raise RuntimeError('{} already exists'.format(store_path))
		else:
			shutil.copyfile(file_path, store_path)

		# Try adding the genome
		try:
			session.add(genome)
			session.commit()
			session.close()

		# On error, remove the file
		except Exception as e:
			os.unlink(store_path)
			raise e

		return genome

	def remove_genome(self, genome):

		genome_path = self._get_path('genomes', genome.filename)

		if os.path.isfile(genome_path):
			os.unlink(genome_path)

		session = self._Session()
		session.delete(session.merge(genome))
		session.commit()
		session.close()

	def open_genome(self, genome):
		path = self._get_path('genomes', genome.filename)
		return open(path)

	def _get_path(self, which, *args):
		return os.path.join(self.directory, self._rel_paths[which], *args)

	@classmethod
	def open(cls, directory):
		directory = os.path.abspath(directory)

		info_path = os.path.join(directory, INFO_FILE_NAME)
		if not os.path.isfile(info_path):
			raise RuntimeError('Directory does not contain a database')
		with open(info_path) as fh:
			info = json.load(fh)

		if info['version'] != DB_VERSION:
			raise RuntimeError('Database is not of the current version')

		return cls(directory)

	@classmethod
	def create(cls, directory, overwrite=False):
		directory = os.path.abspath(directory)

		if not os.path.lexists(directory):
			os.makedirs(directory)

		elif overwrite:
			for relpath in cls._rel_paths.values():
				rmpath(os.path.join(directory, relpath))

		elif os.listdir(directory):
			raise RuntimeError('{} exists and is not empty'.format(directory))

		# Create info file
		info = dict(version=DB_VERSION)
		info_path = os.path.join(directory, INFO_FILE_NAME)
		with open(info_path, 'w') as fh:
			json.dump(info, fh)

		# Create sudirectories
		os.mkdir(os.path.join(directory, cls._rel_paths['kmer_sets']))
		os.mkdir(os.path.join(directory, cls._rel_paths['genomes']))

		db = cls(directory)
		
		# Create SQLite database tables
		Base.metadata.create_all(db._engine)

		return db

	def _make_genome_file_name(self, genome, session, ext):
		if genome.ncbi_acc is not None:
			val = genome.ncbi_acc
		else:
			val = genome.description

		base = re.sub(r'\W+', '_', val[:25])
		filename = base + ext
		i = 0

		session = self.get_session()
		while (session.query(Genome).filter_by(filename=filename).first()
				is not None):

			i += 1
			filename = '{}_{}'.format(base, i) + ext

		return filename
