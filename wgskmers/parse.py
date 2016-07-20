"""Functions for parsing sequence files into k-mer sets"""
from past.builtins import basestring
from builtins import object

import os
import gzip

import numpy as np
from Bio import SeqIO
from tqdm import tqdm

from wgskmers.util import kwargs_finished


# Mapping from file extension to sequnce format
seq_file_exts = [
	(['.fasta', '.fna', '.fas', '.ffn'], 'fasta'),
	(['.fastq'], 'fastq'),
]
seq_file_exts = {ext: fmt for exts, fmt in seq_file_exts for ext in exts}


# Named tuple to store info inferred from file
class SeqFileInfo(object):

	def __init__(self, path, **kwargs):
		self.path = path
		self.basename = os.path.basename(path)
		self.abspath = os.path.abspath(path)

		self.wo_ext = kwargs.pop('wo_ext', os.path.splitext(self.basename))[0]
		self.seq_ext = kwargs.pop('seq_ext', None)
		self.seq_format = kwargs.pop('seq_format', None)
		self.compression = kwargs.pop('compression', None)
		self.contents_ok = kwargs.pop('contents_ok', None)
		self.contents_error = kwargs.pop('contents_error', None)
		kwargs_finished(kwargs)

	def open(self, mode='r'):
		if self.compression is None:
			return open(self.abspath, mode)
		elif self.compression == 'gzip':
			return gzip.open(self.abspath, mode)
		else:
			raise RuntimeError(
				'Can\'t open file with compression "{}"'
				.format(self.compression)
			)

	def check_contents(self):

		# Check compression format
		if self.compression == 'gzip':

			# Check gzip magic number
			with open(self.abspath, 'rb') as gzfh:
				magic = gzfh.read(2)

			if magic != '\x1f\x8b':
				self.contents_ok = False
				self.contents_error = 'Does not appear to be valid gzip file'
				return

		else:
			assert self.compression is None

		# Read first character
		with self.open() as fh:
			first_char = fh.read(1)

		# Check sequence file contents
		if self.seq_format == 'fasta':

			# Should start with >
			if first_char == '>':
				self.contents_ok = True
			else:
				self.contents_ok = False
				self.contents_error = 'FASTA file should start with ">"'

		elif self.seq_format == 'fastq':

			# Should start with @
			if first_char == '@':
				self.contents_ok = True
			else:
				self.contents_ok = False
				self.contents_error = 'FASTQ file should start with "@"'

		else:
			# Unknown format
			self.contents_ok = False
			self.contents_error = 'Unknown format'

	@classmethod
	def get(cls, path, allow_compressed=False, check_contents=False):

		info = SeqFileInfo(path)

		# Get extension
		wo_ext, ext = os.path.splitext(info.basename)

		# Check compression
		if allow_compressed:
			if ext == '.gz':
				info.compression = 'gzip'
				wo_ext, ext = os.path.splitext(wo_ext)
			else:
				info.compression = None
		else:
			info.compression = None

		info.wo_ext = wo_ext

		# Check extension for file type
		info.seq_ext = ext or None
		info.seq_format = seq_file_exts.get(info.seq_ext, None)

		# Check file contents
		if check_contents:
			info.check_contents()

		return info


def find_seq_files(directory, **kwargs):
	"""Finds sequence files in a directory"""

	recursive = kwargs.pop('recursive', False)
	filter_ext = kwargs.pop('filter_ext', True)
	allow_compressed = kwargs.pop('allow_compressed', False)
	filter_contents = kwargs.pop('filter_contents', False)
	warn_contents = kwargs.pop('warn_contents', False)
	check_contents = kwargs.pop('check_contents', filter_contents or warn_contents)
	tqdm_args = kwargs.pop('tqdm', None)
	kwargs_finished(kwargs)

	if warn_contents:
		import click

	# Find file paths
	if recursive:
		paths = (os.path.join(dirpath, fn)
		         for dirpath, dirnames, filenames in os.walk(directory)
		         for fn in filenames)
	else:
		paths = (f for f in os.listdir(directory) if os.path.isfile(f))

	# Show progress
	if tqdm_args is True:
		paths = tqdm(paths)
	elif tqdm_args:
		paths = tqdm(paths, **tqdm_args)

	files_info = []
	for path in paths:

		# Get file info
		info = SeqFileInfo.get(path, allow_compressed=allow_compressed,
		                       check_contents=check_contents)

		# Filter bad extensions/unknown file type
		if filter_ext and info.seq_format is None:
			continue

		# Filter bad contents
		if filter_contents and info.contents_ok is False:
			if warn_contents:
				click.echo('Warning - bad file contents in {}: {}'
				           .format(path, info.contents_error),
				           err=True)
			continue

		files_info.append(info)

	return files_info


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

		if isinstance(file_, basestring):
			fh = open(file_)
			this_format = file_format
		elif isinstance(file_, SeqFileInfo):
			fh = file_.open()
			this_format = file_.seq_format
		else:
			fh = file_
			this_format = file_format

		with fh:
			records = SeqIO.parse(fh, this_format)
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
