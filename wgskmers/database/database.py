import os
import shutil
import json
import re

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from alembic import command as alembic_command
import numpy as np

from wgskmers.util import rmpath
from wgskmers.config import config, save_config
from .models import *
from .sqla import ReadOnlySession
from .migrate import get_alembic_config


# Current database version number
CURRENT_DB_VERSION = 2


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


def get_db_version(path):
	"""Gets version number of database located in path"""
	with open(os.path.join(path, INFO_FILE_NAME)) as fh:
		info = json.load(fh)
	return info['version']


def set_db_version(path, version):
	"""Set version number of database (to be used after upgrading)"""
	info_path = os.path.join(path, INFO_FILE_NAME)

	with open(info_path) as fh:
		info = json.load(fh)

	info['version'] = version

	with open(info_path, 'w') as fh:
		json.dump(info, fh)



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
		'kmer_collections': 'kmer_collections',
		'genomes': 'genomes',
	}

	def __init__(self, directory):
		self.directory = directory

		# SqlAlchemy engine
		self.engine = create_engine('sqlite:///' + self._get_path('sqlite'))

		# SqlAlchemy session classes
		self._Session = sessionmaker(bind=self.engine)
		self._ExpireSession = sessionmaker(bind=self.engine,
		                                   expire_on_commit=False)

	def get_session(self):
		"""Create a new SQLAlchemy session"""
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

	def create_kmer_collection(self, title, prefix, k, parameters=None):

		if parameters is None:
			parameters = dict()

		collection = KmerSetCollection(title=title, prefix=prefix, k=k,
		                               parameters=parameters)

		collection.directory = self._make_kmer_collection_dirname(collection)

		os.mkdir(self._get_path('kmer_collections', collection.directory))

		session = self._ExpireSession()
		session.add(collection)
		session.commit()

		return collection

	def add_kmer_set(self, vec, collection, genome, **kwargs):

		# Destination for file
		filename = 'gen-{}.npy'.format(genome.id)
		store_path = self._get_path('kmer_collections', collection.directory,
		                            filename)

		# Create k-mer set
		kmer_set = KmerSet(
			collection_id=collection.id,
			genome_id=genome.id,
			dtype_str=str(vec.dtype),
			count=(vec > 0).sum(),
			filename=filename,
			**kwargs
		)

		# Write to file
		with open(store_path, 'wb') as fh:
			np.save(fh, vec)

		# Try adding the set
		try:
			session = self._ExpireSession()
			session.add(kmer_set)
			session.commit()

		# On error, remove the file
		except Exception as e:
			os.unlink(store_path)
			raise e

		return kmer_set

	def get_kmer_set_vec(self, kmer_set):

		file_path = self._get_path(
			'kmer_collections',
			kmer_set.collection.directory,
			kmer_set.filename
		)

		with open(file_path, 'rb') as fh:
			return np.load(fh)

	def alembic_config(self):
		"""Creates Alembic configuration for sqlite database"""
		return get_alembic_config(self._get_path('sqlite'))

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

		if info['version'] != CURRENT_DB_VERSION:
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
		info = dict(version=CURRENT_DB_VERSION)
		info_path = os.path.join(directory, INFO_FILE_NAME)
		with open(info_path, 'w') as fh:
			json.dump(info, fh)

		# Create sudirectories
		os.mkdir(os.path.join(directory, cls._rel_paths['kmer_collections']))
		os.mkdir(os.path.join(directory, cls._rel_paths['genomes']))

		# Ready to oepn, create database object
		db = cls(directory)
		
		# Create SQLite database tables
		Base.metadata.create_all(db.engine)

		# Stamp with current alembic revision
		alembic_command.stamp(db.alembic_config(), 'head')

		return db

	def _make_genome_file_name(self, genome, session, ext):
		if genome.gb_acc is not None:
			val = genome.gb_acc
		else:
			val = genome.description

		base = re.sub(r'\W+', '_', val[:25])
		filename = base + ext
		i = 0

		session = self.get_session()
		q = session.query(Genome)
		while q.filter_by(filename=filename).first() is not None:
			i += 1
			filename = '{}_{}'.format(base, i) + ext

		return filename

	def _make_kmer_collection_dirname(self, kmer_collection):

		base = re.sub(r'\W+', '_', kmer_collection.title[:25])
		dirname = base
		i = 0

		session = self.get_session()
		q = session.query(KmerSetCollection)
		while q.filter_by(directory=dirname).first() is not None:
			i += 1
			dirname = '{}_{}'.format(base, i) + ext

		return dirname

