"""Storage formats for kmer sets"""

import numpy as np

from wgskmers.kmers import KmerSpec


class KmerSetStorageFormat(object):
	"""ABC defining a method of storing k-mer sets in database"""

	def __init__(self, collection):
		self.spec = KmerSpec(collection.k, collection.prefix)

	def store(self, fh, vec, kmer_set):
		raise NotImplementedError()

	def load(self, fh, kmer_set):
		raise NotImplementedError()


class RawFormat(KmerSetStorageFormat):
	"""Stores vector in raw numpy format"""

	def store(self, fh, vec, kmer_set):
		np.save(fh, vec)

	def load(self, fh, kmer_set):
		return np.load(fh)


class CoordsFormat(KmerSetStorageFormat):
	"""Stores vector by indices of nonzero elements"""

	def __init__(self, collection):
		super(CoordsFormat, self).__init__(collection)
		self.index_dtype = np.int64 if self.spec.idx_len >= 2**32 else np.int32

	def store(self, fh, vec, kmer_set):
		coords, = np.nonzero(vec)
		coords = coords.astype(self.index_dtype)

		if kmer_set.has_counts:
			counts = vec[coords]
			array = np.stack((coords, counts))
		else:
			array = coords

		np.save(fh, array)

	def load(self, fh, kmer_set):
		array = np.load(fh)

		vec = np.zeros(self.spec.idx_len, dtype=kmer_set.dtype_str)

		if kmer_set.has_counts:
			coords, counts = array
			vec[coords] = counts

		else:
			vec[array] = 1

		return vec


kmer_storage_formats = {
	'raw': RawFormat,
	'coords': CoordsFormat,
}
