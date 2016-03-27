"""Script to find all k-mers in a FASTA file beginnig with a query sequence"""

import os
import sys
import argparse
import logging

from Bio import SeqIO as seqio


logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger()


nucs = set('ATGC')


default_k = 16
default_query = 'ATGAC'


def locate_kmers(query, k, seq):
	"""Generator that finds all k-mers in seq that begin with query

	Note that this is case-sensitive.

	Args:
		query: str|Bio.Seq.Seq. Sequence to search for.
		k: int. Length of k-mers to find (this much space at the end of the
			sequence won't be searched).
		seq: str|Bio.Seq.Seq. Sequence to search within.

	Returns:
		Generator yielding start indices of all matches
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


def filter_kmers(locs, phred, k, thresh):
	"""Filters kmers based on PHRED scores

	Args:
		locs: sequence of int. Start indices of kmers.
		phred: sequence of int. PHRED scores.
		k: int. Length of kmers.
		thresh: int. kmers containing PHRES scores greater than this will be
			filtered out.

	Returns:
		Generator yielding filtered kmer start locations.
	"""
	for loc, score in zip(locs, phred):

		# Check if highest score in the range is under the threshold
		if max(phred[loc:loc + k]) <= thresh:
			yield loc


def process_file(fh, query, k, fmt='fasta', thresh=None):
	"""Creates sorted list of all unique k-mers in FASTA/Q file matching query

	Args:
		fh: str|stream. File name or open file handle to read from (1st
			argument to Bio.SeqIO.parse).
		query: str|Bio.Seq.Seq. Sequence to search for.
		k: int. Length of k-mers to find (this much space at the end of the
			sequence won't be searched).
		fmt: str. File format (2nd argument to Bio.SeqIO.parse).
		thresh: int|None. Filter any kmers containing PHRED scores higher
			than this value.

	Returns:
		list of str. Sorted list of kmers.
	"""
	kmers = set()

	# Parse file and iterate over sequences
	for record in seqio.parse(fh, fmt):

		# Forwards and backwards
		for revcomp in [False, True]:

			# Upper case for search
			seq = record.seq.upper()

			# Reverse compliment
			if revcomp:
				seq = seq.reverse_complement()

			# Find locations of kmers
			kmer_locs = locate_kmers(query, k, seq)

			# Filter by quality if needed
			if thresh is not None:
				phred = record.letter_annotations['phred_quality']
				kmer_locs = filter_kmers(kmer_locs, phred, k, thresh=thresh)

			# Add kmers to set of matches
			for loc in kmer_locs:
				kmer = str(seq[loc:loc + k])

				# Discard bad characters
				if set(kmer).issubset(nucs):
					kmers.add(kmer)

	# Return sorted list
	return sorted(kmers)


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
	help='Filter kmers containing PHRED scores over this value')

parser.add_argument('-f', '--format', type=str, required=False,
	choices=['fasta', 'fastq', 'fastq-sanger', 'fastq-solexa',
	         'fastq-illumina'],
	help='File format, as argument to Bio.SeqIO.parse. If omitted, will infer '
	     'from extension of first file encountered.')

parser.add_argument('-b', '--batch', action='store_true',
	help='Run in match mode. Will process each file in directory given by '
	     '"src" argument, writing output to directory given by "dest".')

parser.add_argument('--debug', action='store_true', help='Print debug messages')
parser.add_argument('-o', '--overwrite', action='store_true',
	help='Overwrite existing output files')

parser.add_argument('src', type=str,
	help='FASTA file to read from, or directory of files in batch mode')
parser.add_argument('dest', type=str, nargs='?',
	help='File or directory to write output to. If ommited will write to stdout')


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
		dest_paths = [make_dest_path(src, dest_arg) for src in src_paths]

	# Single-file mode
	else:

		src_paths = [src_arg]

		# Get destination path
		if dest_arg is not None:

			# If given a directory, output to file with same name
			if os.path.isdir(dest_arg):
				dest_paths = [make_dest_path(src_arg, dest_arg)]
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

		# Find the kmers
		kmers = process_file(src_path, query=query, k=k, fmt=fmt,
		                     thresh=args.threshold)

		# Warn if none found
		if not kmers:
			logger.warn('No kmers found in {}, bad file format?'
				.format(src_path))

		# Otherwise write to output
		else:
			logger.debug('Found {} matches'.format(len(kmers)))

			# Output to file
			if dest_path is not None:
				with open(dest_path, 'w') as fh:
					for kmer in kmers:
						fh.write(kmer + '\n')

			# Output to stdout
			else:
				for kmer in kmers:
					print kmer


# If run as script
if __name__ == '__main__':
	main()
