"""Functions for parsing sequence files into k-mer sets"""

import os

import numpy as np
from Bio import SeqIO
from tqdm import tqdm

from wgskmers.util import kwargs_finished


def infer_format(path):
	"""Infers format for sequence file, returns format argument to
	Bio.SeqIO.parse

	All it does is check the extension.
	"""
	ext = os.path.splitext(path)[1]

	if ext in ['.fasta', '.fna']:
		return 'fasta'

	elif ext in ['.fastq']:
		return 'fastq'

	else:
		return None


seq_file_exts = ['.fasta', '.fastq', '.fna']

def find_seq_files(directory, check_ext=False):
	"""Finds sequence files in a directory"""
	paths = []

	for fname in os.listdir(directory):

		if check_ext:
			if not any(fname.endswith(ext) for ext in seq_file_exts):
				continue

		if os.path.isfile(fname):
			paths.append(os.path.join(directory, fname))

	return paths


def vec_from_records(records, spec, counts=False, **kwargs):
	"""Create a k-mer vector from a set of sequence records.

	Args:
		records: iterable of Bio.SeqRecord.SeqRecord, as output from
			Bio.SeqIO.parse. Records to find k-mers in.
		spec. KmerSpec. Spec defining k-mers to search for.

	kwargs:
		quality_threshold: numeric|None. If not None, get quality scores from
			records and filter out k-mers containing score below this value.

	And more that I'm too lazy to write right now...

	Yields:
		np.ndarray ... TODO
	"""

	# Get kwargs
	q_threshold = kwargs.pop('q_threshold', None)
	c_threshold = kwargs.pop('c_threshold', None)
	out = kwargs.pop('out', None)
	kwargs_finished(kwargs)

	buf = out if counts or c_threshold is None else None

	# Parse file and iterate over sequences
	for record in records:

		# Upper case for search
		seq = record.seq.upper()

		# No quality
		if q_threshold is None:
			finder = spec.find(seq, revcomp=True)

		# With quality info
		else:
			phred_scores = record.letter_annotations['phred_quality']
			finder = spec.find_quality(seq, revcomp=True, quality=phred_scores,
			                           threshold=q_threshold)

		# Get kmer vectors
		if counts or c_threshold is not None:
			buf = finder.counts_vec(out=buf)

		else:
			buf = finder.bool_vec(out=buf)

	if c_threshold is not None:
		return np.greater_equal(buf, c_threshold, out=out)

	else:
		return buf


def parse_to_array(files, spec, out=None, **kwargs):
	"""Parse a set of files into a 2d stack of k-mer vectors"""

	file_format = kwargs.pop('file_format', 'fasta')
	progress_args = kwargs.pop('progress', None)

	if progress_args is True:
		progress_args = dict()

	if out is None:
		out = np.zeros((len(files), spec.idx_len), dtype=np.bool)

	if progress_args is not None:
		iterable = tqdm(files, **progress_args)
	else:
		iterable = files

	for i, file_ in enumerate(iterable):
		with open(file_) as fh:
			records = SeqIO.parse(fh, file_format)
			vec_from_records(records, spec, out=out[i, :], **kwargs)

	return out


class ProgressSeqParser(object):
	"""Wraps generator from Bio.SeqIO.parse with tqdm progress bar."""

	def __init__(self, file_, fmt='fasta', **kwargs):
		self.file_ = file_
		self.fmt = fmt
		self.tqdm_args = kwargs

	def __iter__(self):

		# Open file if given as path
		if isinstance(self.file_, basestring):
			fh = open(self.file_)
			opened = True
		else:
			fh = self.file_
			opened = False

		try:
			# Starting position in file, probably zero but maybe not...
			last = fh.tell()

			# Size of file to parse
			total = os.fstat(fh.fileno()).st_size - last

			# Create progress bar
			pbar = tqdm(unit='B', unit_scale=True, total=total,
			            **self.tqdm_args)
			with pbar:

				# Parse and iterate over records
				for record in SeqIO.parse(fh, self.fmt):

					# Update progress bar
					current = fh.tell()
					pbar.update(current - last)
					last = current

					# Yield parsed record
					yield record

		# Close the file if we opened it
		finally:
			if opened:
				fh.close()