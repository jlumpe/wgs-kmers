"""Functions for dealing with the k-mers themselves - finding them, converting
to and from their indices, etc.

This defines an indexing of all k-mers for a given k, determined by their
lexicographical ordering (nucleotides are ordered alphabetically).

Note that these functions are case-sensitive and expect sequence strings to
be upper case.
"""

import numpy as np

from Bio.Seq import Seq


# The four DNA nucleotides.
# Note the (alphabetical) order - this wil be used to define k-mer indices
nucleotides = sorted('ACGT')

# Dict mapping each nucleotide to its index in the above order
nucleotide_indices = dict((n, i) for i, n in enumerate(nucleotides))


def reverse_compliment(seq):
	"""A quick way of getting the reverse compliment of a sequence

	Args:
		seq: str|Bio.Seq.Seq. Sequence to get reverse compliment of.

	Returns:
		str|Bio.Seq.Seq. Reverse compliment of same type as input.
	"""
	if isinstance(seq, Seq):
		return seq.reverse_complement()
	else:
		return str(Seq(seq).reverse_complement())


def kmer_index(kmer):
	"""Gets the index of a k-mer

	Args:
		kmer: str|Bio.Seq.Seq. K-mer to find index of.

	Returns:
		int. Index of k-mer.
	"""
	index = 0
	for nuc in kmer:
		index <<= 2
		index += nucleotide_indices[nuc]

	return index


def kmer_at_index(index, k):
	"""Get the k-mer at the given index

	Args:
		index: int. Index of k-mer to get.
		k: int. Length of k-mer.

	Returns:
		str. K-mer at index.
	"""
	nucs_reversed = []
	for i in range(k):
		nucs_reversed.append(nucleotides[index % 4])
		index >>= 2
	return ''.join(reversed(nucs_reversed))


def locate_kmers(seq, k, prefix):
	"""Generator that finds locations of k-mers in sequence

	Args:
		seq: str|Bio.Seq.Seq. Sequence to search within.
		k: int. Length of k-mers to find, including prefix.
		prefix: str. Finds k-mers beginning with this subsequence.

	Yields:
		int. Start index of each match (beginning of prefix).
	"""
	start = 0
	end = len(seq) - k
	while True:
		p = seq.find(prefix, start, end)
		if p >= 0:
			yield p
			start = p + 1
		else:
			break


class KmerSpec(object):
	"""Specifications for a k-mer search operation.

	Contains methods for creating KmerFinder objects to find k-mers within
	sequences.

	Properties:
		k: Length of k-mers to find (including prefix).
		prefix: str. Find k-mers beginning with this subsequence. Will be
			converted to upper case.
		plen: int. Length of prefix.
		k_sfx: int. Length of suffix (non-constant part of k-mers).
		idx_len: int. Number of indices needed to index these k-mers - equal
			to 4 ** k_sfx.
	"""

	def __init__(self, k, prefix):
		"""
		Args:
			k: Length of k-mers to find (including prefix).
			prefix: str. Find k-mers beginning with this subsequence.
		"""
		self.k = k
		self.prefix = str(prefix).upper()
		self.plen = len(self.prefix)
		self.k_sfx = self.k - self.plen
		self.idx_len = 4 ** self.k_sfx

	def find(self, seq, **kwargs):
		"""Creates KmerFinder based on this spec that finds k-mers in sequence.

		Args:
			seq: str|Bio.Seq.Seq. Sequence to search within.

		**kwargs (forwarded to KmerFinder()):
			revcomp: bool. If true, search reverse compliment as well.
			circular: bool. If true, take sequence to be circular and wrap
				search around from the end to the beginning.
		"""
		return KmerFinder(self, seq, **kwargs)

	def find_quality(self, seq, quality, threshold, **kwargs):
		"""Creates QualityKmerFinder based on this spec that finds k-mers in
		sequence based on quality.

		Args:
			seq: str|Bio.Seq.Seq. Sequence to search within.
			quality. sequence of numeric. Quality scores, same length as
				sequence.
			threshold. numeric. K-mers found containing quality scores below
				this value with be discarded.

		**kwargs (forwarded to QualityKmerFinder()):
			revcomp: bool. If true, search reverse compliment as well.
			circular: bool. If true, take sequence to be circular and wrap
				search around from the end to the beginning.
		"""
		return QualityKmerFinder(self, seq, quality, threshold, **kwargs)


class KmerFinder(object):
	"""Finds and extracts k-mers from a specific sequence.

	This class should generally be used by calling KmerSpec.find() instead of
	instantiating directly.
	"""

	def __init__(self, spec, seq, revcomp=False, circular=False):
		"""
		Args:
			spec: KmerSpec. Spec defining how to search for k-mers.
			seq: str|Bio.Seq.Seq. Sequence to search within.
			revcomp: bool. If true, search reverse compliment as well.
			circular: bool. If true, take sequence to be circular and wrap
				search around from the end to the beginning.

		Note that the sequence is case-senstitive. Prefixes	are converted to
		upper case by KmerSpec, so sequences should be upper case as well.
		"""
		self.spec = spec
		self.seq = seq
		self.seqlen = len(seq)
		self.find_revcomp = revcomp
		self.seq_circular = circular

	def get_kmers(self):
		"""Generator that yields all k-mers found in the sequence.

		NOTE: this yields k-mers WITHOUT the prefix (it's redundant and
		usually not used anyways), so it's actually giving
		(k - len(prefix))-mers.

		Yields:
			str. Each k-mer found, excluding prefix.
		"""
		for kmer in self._get_kmers(revcomp=False):
			yield str(kmer)
		if self.find_revcomp:
			for kmer in self._get_kmers(revcomp=True):
				yield str(kmer)

	def get_indices(self):
		"""Generator that yields indices of all k-mers found in the sequence.

		Yields:
			int. Index of each found k-mer.
		"""
		for kmer in self.get_kmers():
			if set(kmer).issubset(nucleotides):
				yield kmer_index(kmer)

	def bool_vec(self, out=None, dtype=np.bool):
		"""Creates boolean vector indicating indices of k-mers found.

		Args:
			out: np.ndarray|None. Array to write values to. Indices
				corresponding to found k-mers will be set to True (1). If
				None, one will be created. Should be 1d of length idx_len
				of the KmerSpec.
			dtype: np.dtype. Dtype of output array, if created automatically.

		Returns:
			np.ndarray. Same as out argument if not None, otherwise array with
				dtype set by dtype argument.
		"""
		if out is None:
			out = np.zeros(self.spec.idx_len, dtype=dtype)

		for index in self.get_indices():
			out[index] = True

		return out

	def counts_vec(self, out=None, dtype=np.uint16):
		"""Creates vector giving counts of k-mers found at each index.

		Args:
			out: np.ndarray|None. Array to write values to. Indices
				corresponding to found k-mers will be incremented. If
				None, one will be created. Should be 1d of length idx_len
				of the KmerSpec.
			dtype: np.dtype. Dtype of output array, if created automatically.

		Returns:
			np.ndarray. Same as out argument if not None, otherwise array with
				dtype set by dtype argument.
		"""

		if out is None:
			out = np.zeros(self.spec.idx_len, dtype=dtype)

		try:
			old_err = np.seterr(over='raise')

			for index in self.get_indices():
				out[index] += 1

		finally:
			np.seterr(**old_err)

		return out

	def _get_kmers(self, revcomp=False):
		"""Internal generator method that extracts the k-mer sequences"""

		if revcomp:
			seq = reverse_compliment(self.seq)
		else:
			seq = self.seq

		# Extract in the forward direction as a linear sequence
		for loc in locate_kmers(seq, self.spec.k, self.spec.prefix):
			yield seq[loc + self.spec.plen : loc + self.spec.k]

		# Account for circular sequences
		if self.seq_circular:
			# Search from (k-1) from the end to (k-1) after the beginning
			# (the k-1 excludes matches we may have found before)
			wrap_seq = seq[-(self.spec.k-1):] + seq[:(self.spec.k-1)]
			for loc in locate_kmers(wrap_seq, self.spec.k, self.spec.prefix):
				yield wrap_seq[loc + self.spec.plen : loc + self.spec.k]


class QualityKmerFinder(KmerFinder):
	"""Finds and extracts k-mers from a sequence with quality scores.

	This class should generally be used by calling KmerSpec.find_quality()
	instead of instantiating directly.
	"""

	def __init__(self, spec, seq, quality, threshold, **kwargs):
		"""
		Args:
			spec: KmerSpec. Spec defining how to search for k-mers.
			seq: str|Bio.Seq.Seq. Sequence to search within.
			quality. sequence of numeric. Quality scores, same length as
				sequence.
			threshold. numeric. K-mers found containing quality scores below
				this value with be discarded.

		**kwargs:
			revcomp: bool. If true, search reverse compliment as well.
			circular: bool. If true, take sequence to be circular and wrap
				search around from the end to the beginning.

		Note that the sequence is case-senstitive. Prefixes	are converted to
		upper case by KmerSpec, so sequences should be upper case as well.
		"""
		super(QualityKmerFinder, self).__init__(spec, seq, **kwargs)

		self.quality = quality
		self.threshold = threshold

	def _get_kmers(self, revcomp=False):
		"""Internal generator method that extracts the k-mer sequences"""

		if revcomp:
			seq = reverse_compliment(self.seq)
			qual = self.quality[::-1]
		else:
			seq = self.seq
			qual = self.quality

		# Extract in the forward direction as a linear sequence
		for loc in locate_kmers(seq, self.spec.k, self.spec.prefix):
			s = slice(loc + self.spec.plen, loc + self.spec.k)
			if min(qual[s]) >= self.threshold:
				yield seq[s]

		# Account for circular sequences
		if self.seq_circular:
			# Search from (k-1) from the end to (k-1) after the beginning
			# (the k-1 excludes matches we may have found before)
			wrap_seq = seq[-(self.spec.k-1):] + seq[:(self.spec.k-1)]
			wrap_qual = qual[-(self.spec.k-1):] + qual[:(self.spec.k-1)]

			for loc in locate_kmers(wrap_seq, self.spec.k, self.spec.prefix):
				s = slice(loc + self.spec.plen, loc + self.spec.k)
				if min(wrap_qual[s]) >= self.threshold:
					yield wrap_seq[s]
