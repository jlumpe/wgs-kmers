import os
import shutil
import json
import re
import gzip

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from alembic import command as alembic_command
import numpy as np

from wgskmers.util import rmpath, kwargs_finished
from wgskmers.config import get_config
from wgskmers.kmers import KmerSpec, KmerCoordsCollection
from .models import *
from .sqla import ReadOnlySession
from .store import kmer_storage_formats
from .migrate import get_alembic_config


# Current database version number
CURRENT_DB_VERSION = 5


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
	config_default = get_config().get_default_db()
	if config_default is not None:
		return config_default, 'config'

	# Check default environment variable
	if DEFAULT_DB_PATH_VAR in os.environ:
		return os.environ[DEFAULT_DB_PATH_VAR], 'environ'

	# Not found
	return None, None


def get_default_db():
	"""Get database registered as default"""

	# Check config file
	config_default = get_config().get_default_db()
	if config_default is not None:
		return config_default

	# Check default environment variable
	if DEFAULT_DB_PATH_VAR in os.environ:
		return os.environ[DEFAULT_DB_PATH_VAR]

	# Not found
	return None


class Database(object):
	"""Stores genome files and pre-calculated k-mer sets, plus metadata.

	This class is pickleable, and is intended to be able to be used with
	multiprocessing.
	"""

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

	def __getstate__(self):
		"""For pickling"""
		return self.directory

	def __setstate__(directory):
		"""For unpickling"""
		self.__init__(directory)

	def get_session(self):
		"""Create a new SQLAlchemy session"""
		return self._Session()

	def store_genome(self, file_, **kwargs):

		session = self._ExpireSession()

		# Get function kwargs
		src_compression = kwargs.pop('src_compression', None)
		src_mode = kwargs.pop('src_mode', 'rb' if src_compression else None)
		keep_src = kwargs.pop('keep_src', True)

		# Create genome from remaining kwargs
		kwargs.setdefault('compression', 'gzip')
		genome = Genome(**kwargs)

		# Determine the file name
		genome.filename = self._make_genome_file_name(genome, session)

		# Path to store at
		store_path = self._get_path('genomes', genome.filename)
		if os.path.exists(store_path):
			raise RuntimeError('{} already exists'.format(store_path))

		# Open destination using correct mode and format
		if src_mode is not None:
			write_mode = 'wb' if 'b' in src_mode else 'w'
		elif not isinstance(file_, basestring) and hasattr(src_fh, 'mode') and 'b' in src_fh.mode:
			write_mode = 'wb'
		else:
			write_mode = 'w'

		# Check if source and destination are in the same format
		same_format = genome.compression == src_compression
		if same_format:

			# See if we are moving/copying a file in the same format, if so can
			# just use file system operations
			if isinstance(file_, basestring):

				# Copy it
				if keep_src:
					shutil.copyfile(file_, store_path)
					src_moved = False

				# Move it
				else:
					os.rename(file_, store_path)
					src_moved = True

			# Copying stream in same format
			else:
				src_moved = False

				with open(store_path, write_mode) as dest_fh:
					shutil.copyfileobj(file_, dest_fh)

			needs_delete = False

		# Different formats
		else:
			src_moved = False

			# Copying from file
			if isinstance(file_, basestring):

				if src_compression is None:
					src_fh = open(file_, 'rb')

				elif src_compression == 'gzip':
					src_fh = gzip.open(file_, src_mode)

				else:
					raise ValueError(
						'Unknown compression format "{}"'
						.format(src_compression)
					)

				needs_delete = not keep_src

			# Copying from file-like object
			else:
				if src_compression is None:
					src_fh = file_

				elif src_compression == 'gzip':
					src_fh = gzip.GzipFile(fileobj=file_, mode=src_mode)

				else:
					raise ValueError(
						'Unknown compression format "{}"'
						.format(src_compression)
					)

				needs_delete = False

			# Open destinaction file
			if genome.compression is None:
				dest_fh = open(store_path, write_mode)

			elif genome.compression == 'gzip':
				dest_fh = gzip.open(store_path, write_mode)

			else:
				raise ValueError(
						'Unknown compression format "{}"'
					.format(genome.compression)
				)

			# Copy data
			with dest_fh:
				shutil.copyfileobj(src_fh, dest_fh)

		# Try adding the genome
		try:
			session.add(genome)
			session.commit()

		# Try to reverse file operations on error
		except Exception:

			# Move stored file back or remove it
			if src_moved:
				os.rename(store_path, file_)
			else:
				os.unlink(store_path)

			# Propagate exception
			raise

		finally:
			# Always close session
			session.close()

		# Remove the original if needed
		if needs_delete:
			os.unlink(file_)

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

		if genome.compression is None:
			return open(path)

		elif genome.compression == 'gzip':
			return gzip.open(path, 'r')

		else:
			raise RuntimeError('Can\'t open genome with compression "{}"'
			                   .format(genome.compression))

	def create_kmer_collection(self, **kwargs):

		kwargs.setdefault('parameters', dict())
		kwargs.setdefault('format', 'coords')

		assert kwargs['format'] in kmer_storage_formats

		collection = KmerSetCollection(**kwargs)

		collection.directory = self._make_kmer_collection_dirname(collection)

		col_path = self._get_path('kmer_collections', collection.directory)
		os.mkdir(col_path)

		try:
			session = self._ExpireSession()
			session.add(collection)
			session.commit()

		except:
			os.rmdir(col_path)
			raise

		return collection

	def store_kmer_sets(self, collection):
		return KmerSetAdder(self, collection)

	def get_kmer_loader(self, collection):
		return KmerSetLoader(self, collection)

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

	def _make_genome_file_name(self, genome, session):
		if genome.gb_acc is not None:
			val = genome.gb_acc
		else:
			val = genome.description

		max_len = 25
		ext = '.' + genome.file_format
		if genome.compression == 'gzip':
			ext += '.gz'

		base = re.sub(r'\W+', '_', val[:max_len])
		filename = base + ext
		i = 0

		q = session.query(Genome)
		while q.filter_by(filename=filename).first() is not None:
			i += 1
			filename = '{}_{}'.format(base, i) + ext

		return filename

	def _make_kmer_collection_dirname(self, kmer_collection):

		base = re.sub(r'\W+', '_', kmer_collection.title[:25]).lower()
		dirname = base
		i = 0

		session = self.get_session()
		q = session.query(KmerSetCollection)
		while q.filter_by(directory=dirname).first() is not None:
			i += 1
			dirname = '{}_{}'.format(base, i) + ext

		return dirname


class KmerSetLoader(object):
	"""Loads k-mer set vectors from collection database

	This is created as a separate object so that some initialization can
	only be done once when loading a large number of sets.
	"""

	def __init__(self, db, collection):
		self.db = db
		self.collection = collection
		self.format = kmer_storage_formats[collection.format](collection)

	def __getstate__(self):
		"""For pickling"""
		return (self.db, self.collection)

	def __setstate__(state):
		"""For unpickling"""
		self.__init__(*state)

	def _path_for(self, kmer_set):
		return self.db._get_path(
			'kmer_collections',
			self.collection.directory,
			kmer_set.filename
		)

	def load(self, kmer_set):
		"""Load a single k-mer set vector from the collection"""
		assert kmer_set.collection_id == self.collection.id

		with open(self._path_for(kmer_set), 'rb') as fh:
			return self.format.load(fh, kmer_set)

	def load_array(self, kmer_sets, out=None, dtype=None):

		if dtype is None:
			dtype = np.dtype(bool)

		if out is None:
			spec = KmerSpec(self.collection.k, self.collection.prefix)
			out = np.ndarray((len(kmer_sets), spec.idx_len), dtype=dtype)

		for i, kmer_set in enumerate(kmer_sets):
			out[i, :] = self.load(kmer_set)

		return out

	def load_coords(self, kmer_set, counts=False):
		assert kmer_set.collection_id == self.collection.id

		with open(self._path_for(kmer_set), 'rb') as fh:
			array = self.format.load_coords(fh, kmer_set)

		if counts:
			if kmer_set.has_counts:
				return array
			else:
				return np.vstack(array, np.ones(len(array), dtype=array.dtype))

		else:
			if kmer_set.has_counts:
				return array[0, :]
			else:
				return array

	def load_coords_col(self, kmer_sets, cls=KmerCoordsCollection):
		assert all(ks.collection_id == self.collection.id for ks in kmer_sets)

		coords_col = cls.empty([ks.count for ks in kmer_sets])

		for i, ks in enumerate(kmer_sets):
			coords_col[i] = self.load_coords(ks, counts=False)

		return coords_col


class KmerSetAdder(object):

	def __init__(self, db, collection):
		self.db = db
		self.collection = collection
		self.format = kmer_storage_formats[collection.format](collection)

	def __call__(self, vec, genome, **kwargs):

		# Destination for file
		filename = 'gen-{}.npy'.format(genome.id)
		store_path = self.db._get_path(
			'kmer_collections',
			self.collection.directory,
			filename
		)

		# Create k-mer set
		kmer_set = KmerSet(
			collection_id=self.collection.id,
			genome_id=genome.id,
			dtype_str=str(vec.dtype),
			count=(vec > 0).sum(),
			filename=filename,
			**kwargs
		)

		# Write to file
		with open(store_path, 'wb') as fh:
			self.format.store(fh, vec, kmer_set)

		# Try adding the set
		try:
			session = self.db._ExpireSession()
			session.add(kmer_set)
			session.commit()
			session.close()

		# On error, remove the file
		except Exception as e:
			os.unlink(store_path)
			raise

		return kmer_set
