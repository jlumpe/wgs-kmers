"""Utilites for console commands"""

import os

import itertools

from Bio import SeqIO


# Use tqdm for progress bars if installed, but don't fail if not
# Commands should check if it is None before using
try:
	from tqdm import tqdm
except ImportError:
	tqdm = None


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


def iterator_empty(iterator):
	"""Checks if an iterator is empty, also returning a substitute iterator.

	Consumes first element of iterator, but stores it for the substitute
	iterator to yield first.
	"""

	try:
		first = [next(iterator)]
		empty = False
	except StopIteration:
		first = []
		empty = True

	substitute = itertools.chain(first, iterator)

	return substitute, empty
