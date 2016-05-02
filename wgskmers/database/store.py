"""Storage formats for kmer sets"""

import numpy as np

from wgskmers.kmers import KmerSpec, vec_to_coords, coords_to_vec


class KmerSetStorageFormat(object):
	"""ABC defining a method of storing k-mer sets in database"""

	def __init__(self, collection):
		self.spec = KmerSpec(collection.k, collection.prefix)
		self.index_dtype = np.int64 if self.spec.idx_len >= 2**32 else np.int32

	def store(self, fh, vec, kmer_set):
		raise NotImplementedError()

	def load(self, fh, kmer_set):
		raise NotImplementedError()

	def store_coords(self, fh, coords, kmer_set):
		raise NotImplementedError()

	def load_coords(self, fh, kmer_set):
		raise NotImplementedError()


class RawFormat(KmerSetStorageFormat):
	"""Stores vector in raw numpy format"""

	def store(self, fh, vec, kmer_set):
		np.save(fh, vec)

	def load(self, fh, kmer_set):
		return np.load(fh)

	def store_coords(self, fh, coords, kmer_set):
		vec = coords_to_vec(coords, has_counts=kmer_set.has_counts,
		                    idx_len=self.spec.idx_len, dtype=kmer_set.dtype_str)
		self.store(fh, vec, kmer_set)

	def load_coords(self, fh, kmer_set):
		vec = np.load(fh)
		return vec_to_coords(vec, counts=kmer_set.has_counts,
		                     dtype=self.index_dtype)


class CoordsFormat(KmerSetStorageFormat):
	"""Stores vector by indices of nonzero elements"""

	def store(self, fh, vec, kmer_set):
		coords = vec_to_coords(vec, counts=kmer_set.has_counts,
		                       dtype=self.index_dtype)

		np.save(fh, coords)

	def load(self, fh, kmer_set):
		coords = np.load(fh)

		return coords_to_vec(coords, has_counts=kmer_set.has_counts,
		                     idx_len=self.spec.idx_len, dtype=kmer_set.dtype_str)

	def store_coords(self, fh, coords, kmer_set):
		np.save(fh, coords)

	def load_coords(self, fh, kmer_set):
		return np.load(fh)


kmer_storage_formats = {
	'raw': RawFormat,
	'coords': CoordsFormat,
}
