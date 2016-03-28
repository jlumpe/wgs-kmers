"""Script to find all k-mers in a FASTA file beginnig with a query sequence

By default, the script will output a sorted list of all k-mers found.
"""

import os
import sys
import argparse
import logging
from collections import Counter

from Bio import SeqIO as seqio


logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger()


nucs = set('ATGC')
nuc_idx = dict((n, i) for i, n in enumerate(sorted(nucs)))


default_k = 16
default_query = 'ATGAC'


def kmer_index(kmer):

	idx = 0
	for nuc in kmer:
		idx <<= 2
		idx += nuc_idx[nuc]

	return idx


def locate_kmers(query, k, seq):
	"""Generator that finds all k-mers in seq that begin with query

	Note that this is case-sensitive.

	Args:
		query: str|Bio.Seq.Seq. Sequence to search for.
		k: int. Length of k-mers to find (this much space at the end of the
			sequence won't be searched).
		seq: str|Bio.Seq.Seq. Sequence to search within.

	Yields:
		int. Start indices of all matches
	"""
	start = 0
	end = len(seq) - k
	while True:
		p = seq.find(query, start, end)
		if p >= 0:
			yield p
			start = p + 1
		else:
			break


def kmers_from_fasta(fh, query, k, fmt='fasta', with_qual=False):
	"""Generator yielding k-mers in FASTA/Q file beginning with query, plus
	quality information. 

	Args:
		fh: str|stream. File name or open file handle to read from (1st
			argument to Bio.SeqIO.parse).
		query: str|Bio.Seq.Seq. Sequence to search for.
		k: int. Length of k-mers to find (this much space at the end of the
			sequence won't be searched).
		fmt: str. File format (2nd argument to Bio.SeqIO.parse).

	Yields:
		str or (str, int). By default, yields matching k-mers as str. If
			with_qual is Ture, yields 2-tuples of k-mers plus maximum PHRED
			score within k-mer.
	"""
	has_seqs = False

	# Parse file and iterate over sequences
	for record in seqio.parse(fh, fmt):

		has_seqs = True

		# Forwards and backwards
		for revcomp in [False, True]:

			# Upper case for search
			seq = record.seq.upper()

			# Reverse compliment
			if revcomp:
				seq = seq.reverse_complement()

			# Find locations of kmers
			kmer_locs = locate_kmers(query, k, seq)

			# Get quality info
			if with_qual:
				phred_scores = record.letter_annotations['phred_quality']
				if revcomp:
					phred_scores = list(reversed(phred_scores))

			# Iterate kmer locations
			for loc in kmer_locs:

				# Get kmer as string
				kmer = str(seq[loc:loc + k])

				# Skip kmers containing bad characters
				if not set(kmer).issubset(nucs):
					continue

				# Yield kmer with quality info
				if with_qual:
					qual = max(phred_scores[loc:loc + k])
					yield kmer, qual

				# Kmer only
				else:
					yield kmer

	# None found
	if not has_seqs:
		logger.warn('No sequences found, bad file format?')


def write_kmer_list(stream, kmers):
	"""Write sorted list of unique k-mers to output stream

	Args:
		stream: writeable stream. Stream to write to.
		kmers: Iterable of str. K-mers to write.
	"""
	# Find unique k-mers
	kmer_set = set(kmers)
	logger.debug('Found {} matches'.format(len(kmer_set)))

	# Sort them
	kmers_sorted = sorted(kmer_set)

	# Write
	for kmer in kmers_sorted:
		stream.write(kmer + '\n')


def write_kmer_counts(stream, kmers):
	"""Write sorted list of unique k-mers plus counts to output stream

	Args:
		stream: writeable stream. Stream to write to.
		kmers: Iterable of str. K-mers to write.
	"""
	# Get k-mer counts
	kmer_counts = Counter(kmers)
	logger.debug('Found {} matches'.format(len(kmer_counts)))

	# Sort k-mers
	kmers_sorted = sorted(kmer_counts.keys())

	# Write
	for kmer in kmers_sorted:
		stream.write('{} {}\n'.format(kmer, kmer_counts[kmer]))


def write_kmer_hist(stream, kmers):
	"""Write histogram of k-mer counts to output stream

	Args:
		stream: writeable stream. Stream to write to.
		kmers: Iterable of str. K-mers to write.
	"""
	# Get k-mer counts
	kmer_counts = Counter(kmers)
	logger.debug('Found {} matches'.format(len(kmer_counts)))

	# Histogram of counts
	counts_hist = Counter(kmer_counts.values())

	# Write
	for n in sorted(counts_hist.keys()):
		stream.write('{} {}\n'.format(n, counts_hist[n]))


def write_kmer_vec(stream, kmers, k, plen):
	"""Write boolean vector of kmer occurrences to output stream

	All 4^k possible k-mers for a given k are assigned an index based on
	alphabetic order. 4^k bytes are written to the output stream. Each is one
	if the k-kmer with the associated index was present in the input,
	otherwise zero.

	Args:
		stream: writeable stream. Stream to write to.
		kmers: Iterable of str. K-mers to write.
		k: int. Size of k-mers.
		plen: int. Length of constant perfix to ignore in each k-mer.
	"""
	# Initialize byte array (yes, this is 8x bigger than it needs to be)
	vec = bytearray(4 ** (k - plen))

	# For each kmer, set the appropriate index to one
	for kmer in kmers:
		idx = kmer_index(kmer[plen:])
		vec[idx] = 1

	# Write to output
	stream.write(vec)


def make_dest_path(src_path, dest_dir, ext='.txt'):
	"""Creates file path with same name as file given in src_path, but in
	dest_dir and with a different extension
	"""
	src_name = os.path.splitext(os.path.basename(src_path))[0]
	return os.path.join(dest_dir, src_name + ext)


def infer_format(path):
	"""Infers format for sequence file, returns format argument to
	Bio.SeqIO.parse

	All it does is check the extension.
	"""
	ext = os.path.splitext(path)[1]
	if ext in ['.fasta', '.fastq']:
		return ext[1:]
	else:
		return None


# Argument parser
parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument('-k', '--k', type=int, default=default_k,
	help='Length of k-mers to find (default {})'.format(default_k))
parser.add_argument('-q', '--query', type=str, default=default_query,
	help='Target sequence to find (default {})'.format(default_query))
parser.add_argument('-t', '--threshold', type=int, required=False,
	help='Filter k-mers containing PHRED scores over this value')

output_group = parser.add_mutually_exclusive_group()
output_group.add_argument('-c', '--count', action='store_true',
	help='Write number of occurrences of each k-mer in ouput')
output_group.add_argument('-H', '--hist', action='store_true',
	help='Output histogram of counts')
output_group.add_argument('-B', '--bool', action='store_true',
	help='Output boolean vector, one 0/1 byte per possible k-mer ordered alphabetically')

parser.add_argument('-f', '--format', type=str, required=False,
	choices=['fasta', 'fastq', 'fastq-sanger', 'fastq-solexa',
	         'fastq-illumina'],
	help='File format, as argument to Bio.SeqIO.parse. If omitted, will infer '
	     'from extension of first file encountered.')

parser.add_argument('-b', '--batch', action='store_true',
	help='Run in batch mode. Will process each file in directory given by '
	     '"src" argument, writing output to directory given by "dest".')

parser.add_argument('--debug', action='store_true', help='Print debug messages')
parser.add_argument('-o', '--overwrite', action='store_true',
	help='Overwrite existing output files')

parser.add_argument('src', type=str,
	help='FASTA file to read from, or directory of files in batch mode')
parser.add_argument('dest', type=str, nargs='?',
	help='File or directory to write output to. If omitted will write to stdout')


def main(args=None):
	"""Runs the script"""

	# Parse args
	args = parser.parse_args(args)

	# Check query valid
	query = args.query.upper()
	if not set(query).issubset(nucs):
		raise ValueError('Query contains invalid characters')

	# Check k positive
	k = args.k
	if k <= 0:
		raise ValueError('K must be positive')

	# Absolute source and destination paths
	src_arg = os.path.realpath(args.src)
	dest_arg = os.path.realpath(args.dest) if args.dest is not None else None

	# Debug mode
	if args.debug:
		logger.setLevel(logging.DEBUG)

	# Output file extension
	if args.bool:
		out_ext = '.kmer_vec'
	elif args.hist:
		out_ext = '.hist.txt'
	elif args.count:
		out_ext = '.counts.txt'
	else:
		out_ext = '.kmers.txt'

	if args.threshold is not None:
		out_ext = '-t{}'.format(args.threshold) + out_ext

	# Batch mode
	if args.batch:

		# Check arguments compatible with batch mode
		if dest_arg is None:
			raise ValueError('Must give destination directory in batch mode')

		# Find files in source directory (raises OSError if doesn't exist)
		src_paths = [os.path.join(src_arg, n) for n in os.listdir(src_arg)]
		if not src_paths:
			raise RuntimeError('No files found in {}'.format(src_arg))

		# Create output directory if necessary
		if not os.path.isdir(dest_arg):
			os.mkdir(dest_arg)

		# Destination file paths
		dest_paths = [make_dest_path(src, dest_arg, ext=out_ext)
		              for src in src_paths]

	# Single-file mode
	else:

		src_paths = [src_arg]

		# Get destination path
		if dest_arg is not None:

			# If given a directory, output to file with same name
			if os.path.isdir(dest_arg):
				dest_paths = [make_dest_path(src_arg, dest_arg, ext=out_ext)]
			else:
				dest_paths = [dest_arg]

		else:
			dest_paths = [None]

	# Get format
	if args.format is None:
		fmt = infer_format(src_paths[0])
		if fmt is None:
			raise ValueError("Couldn't infer format for {}"
				.format(src_paths[0]))
		logger.debug('Inferring format "{}" from {}'.format(fmt, src_paths[0]))
	else:
		fmt = args.format

	# Should be ok, loop over the files
	for src_path, dest_path in zip(src_paths, dest_paths):
		
		logger.debug('Processing source file {}'.format(src_path))

		# Check if output file exists
		if dest_path is not None and os.path.exists(dest_path):
			if args.overwrite:
				logger.warn('Overwriting output file {}'.format(dest_path))
			else:
				logger.warn('Refusing to overwrite {}'.format(dest_path))
				continue

		# Find the k-mers with quality info
		if args.threshold is not None:
			kmers_q = kmers_from_fasta(src_path, query=query, k=k, fmt=fmt,
			                           with_qual=True)

			# Filter by quality threshold
			# (note generator)
			kmers = (kmer for kmer, qual in kmers_q if qual <= args.threshold)

		else:
			# Just the k-mers themselves
			kmers = kmers_from_fasta(src_path, query=query, k=k, fmt=fmt)

		# Write output in try block because we can't use a with statement here
		try:
			# Output stream - file (text/binary) or stdout
			if dest_path is None:
				out_stream = sys.stdout
			elif args.bool:
				out_stream = open(dest_path, 'wb')
			else:
				out_stream = open(dest_path, 'w')

			# Output based on arguments
			if args.bool:
				write_kmer_vec(out_stream, kmers, k=k, plen=len(query))
			if args.count:
				write_kmer_counts(out_stream, kmers)
			elif args.hist:
				write_kmer_hist(out_stream, kmers)
			else:
				write_kmer_list(out_stream, kmers)

		# Close file handle
		finally:
			if dest_path is not None:
				out_stream.close()


# If run as script
if __name__ == '__main__':
	main()
