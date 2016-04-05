"""Stand-alone command to find k-mers in a file or set of files.

By default, the script will output a sorted list of all k-mers found.
"""

import os
import sys
import logging
from collections import Counter

from Bio import SeqIO

from .parser import subparsers
from .util import tqdm, ProgressSeqParser, iterator_empty
from wgskmers.kmers import KmerSpec, nucleotides


logger = logging.getLogger()


# Default arguments
default_k = 16
default_prefix = 'ATGAC'


#-----------------------------------------------------------------------------
# Argument parser
#-----------------------------------------------------------------------------

parser = subparsers.add_parser('find', description=__doc__,
	help='Find k-mers in a single file or set of files')

# Parameters
parser.add_argument('-k', '--k', type=int, default=default_k,
	help='Length of k-mers to find (default {})'.format(default_k))
parser.add_argument('-q', '--prefix', type=str, default=default_prefix,
	help='Target sequence to find (default {})'.format(default_prefix))
parser.add_argument('-t', '--threshold', type=int, required=False,
	help='Filter k-mers containing PHRED scores below this value')

# Output format
output_group = parser.add_mutually_exclusive_group()
output_group.add_argument('-c', '--count', action='store_true',
	help='Write number of occurrences of each k-mer in ouput')
output_group.add_argument('-H', '--hist', action='store_true',
	help='Output histogram of counts')
output_group.add_argument('-B', '--bool', action='store_true',
	help='Output boolean vector, one 0/1 byte per possible k-mer ordered alphabetically')

# Input/output options
parser.add_argument('-f', '--format', type=str, required=False,
	choices=['fasta', 'fastq', 'fastq-sanger', 'fastq-solexa',
	         'fastq-illumina'],
	help='File format, as argument to Bio.SeqIO.parse. If omitted, will infer '
	     'from extension of first file encountered.')
parser.add_argument('-b', '--batch', action='store_true',
	help='Run in batch mode. Will process each file in directory given by '
	     '"src" argument, writing output to directory given by "dest".')
parser.add_argument('-o', '--overwrite', action='store_true',
	help='Overwrite existing output files')

# Other
parser.add_argument('-p', '--progress', action='store_true',
	help='Display progress')

# Source and destination files
parser.add_argument('src', type=str,
	help='FASTA file to read from, or directory of files in batch mode')
parser.add_argument('dest', type=str, nargs='?',
	help='File or directory to write output to. If omitted will write to stdout')


#-----------------------------------------------------------------------------
# Code...
#-----------------------------------------------------------------------------

def kmers_from_records(records, spec, quality_threshold=None):
	"""Generator yielding KmerFinders for a set of sequence records.

	Args:
		records: iterable of Bio.SeqRecord.SeqRecord, as output from
			Bio.SeqIO.parse. Records to find k-mers in.
		spec. KmerSpec. Spec defining k-mers to search for.
		quality_threshold: numeric|None. If not None, get quality scores from
			records and filter out k-mers containing score below this value.

	Yields:
		KmerFinder or QualityKmerFinder, depending if quality_threshold was
			None or not.
	"""

	# Parse file and iterate over sequences
	for record in records:

		# Upper case for search
		seq = record.seq.upper()

		# No quality
		if quality_threshold is None:
			yield spec.find(seq, revcomp=True)

		# With quality info
		else:
			phred_scores = record.letter_annotations['phred_quality']
			yield spec.find_quality(seq, revcomp=True, quality=phred_scores,
			                        threshold=quality_threshold)


def write_kmer_list(stream, finders):
	"""Write sorted list of unique k-mers to output stream

	Args:
		stream: writeable stream. Stream to write to.
		finders: Iterable of KmerFinder. Finders yielding k-mers to write.

	Returns:
		int. Number of unique k-mers found.
	"""
	# Find unique k-mers
	kmer_set = {k for f in finders for k in f.get_kmers()}

	# Sort them
	kmers_sorted = sorted(kmer_set)

	# Write
	for kmer in kmers_sorted:
		stream.write(kmer + '\n')

	return len(kmer_set)


def write_kmer_counts(stream, finders):
	"""Write sorted list of unique k-mers plus counts to output stream

	Args:
		stream: writeable stream. Stream to write to.
		finders: Iterable of KmerFinder. Finders yielding k-mers to write.

	Returns:
		int. Number of unique k-mers found.
	"""
	# Get k-mer counts
	kmer_counts = Counter(k for f in finders for k in f.get_kmers())

	# Sort k-mers
	kmers_sorted = sorted(kmer_counts.keys())

	# Write
	for kmer in kmers_sorted:
		stream.write('{} {}\n'.format(kmer, kmer_counts[kmer]))

	return len(kmer_counts)


def write_kmer_hist(stream, finders):
	"""Write histogram of k-mer counts to output stream

	Args:
		stream: writeable stream. Stream to write to.
		finders: Iterable of KmerFinder. Finders yielding k-mers to write.

	Returns:
		int. Number of unique k-mers found.
	"""
	# Get k-mer counts
	kmer_counts = Counter(k for f in finders for k in f.get_kmers())

	# Histogram of counts
	counts_hist = Counter(kmer_counts.values())

	# Write
	for n in sorted(counts_hist.keys()):
		stream.write('{} {}\n'.format(n, counts_hist[n]))

	return len(kmer_counts)


def write_kmer_vec(stream, finders):
	"""Write boolean vector of kmer occurrences to output stream

	All 4^k possible k-mers for a given k are assigned an index based on
	alphabetic order. 4^k bytes are written to the output stream. Each is one
	if the k-kmer with the associated index was present in the input,
	otherwise zero.

	Args:
		stream: writeable stream. Stream to write to.
		finders: Iterable of KmerFinder. Finders yielding k-mers to write.

	Returns:
		int. Number of unique k-mers found.
	"""
	# Have the first KmerFinder.bool_vec create the output array for us,
	# afterwards re-use it so as to find the union.
	vec = None
	for finder in finders:
		vec = finder.bool_vec(out=vec)

	# Write to output
	vec.tofile(stream)

	return vec.sum()


def make_dest_path(src_path, dest_dir, ext='.txt'):
	"""Creates file path with same name as file given in src_path, but in
	dest_dir and with a different extension.
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


def find_batch_files(path):
	"""Finds source files in directory for batch processing"""
	return [f for f in os.listdir(path) if os.path.isfile(f)]


#-----------------------------------------------------------------------------
# Main function
#-----------------------------------------------------------------------------

def main(args):
	"""Runs the sub-command"""

	# Check prefix valid
	prefix = args.prefix.upper()
	if not set(prefix).issubset(nucleotides):
		raise ValueError('Prefix contains invalid characters')

	# Check k positive
	k = args.k
	if k <= 0:
		raise ValueError('K must be positive')

	# Kmer spec
	spec = KmerSpec(k, prefix)

	# Check tqdm present
	if args.progress:
		if tqdm is None:
			raise RuntimeError('tqdm package required to display progress')

	# Absolute source and destination paths
	src_arg = os.path.realpath(args.src)
	dest_arg = os.path.realpath(args.dest) if args.dest is not None else None

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
		src_paths = [os.path.join(src_arg, n) for n
		             in find_batch_files(src_arg)]
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
	paths_iter = zip(src_paths, dest_paths)
	if args.batch and args.progress:
		paths_iter = tqdm(paths_iter, unit='files')

	for src_path, dest_path in paths_iter:
		
		logger.debug('Processing source file {}'.format(src_path))

		# Check if output file exists
		if dest_path is not None and os.path.exists(dest_path):
			if args.overwrite:
				logger.warn('Overwriting output file {}'.format(dest_path))
			else:
				logger.warn('Refusing to overwrite {}'.format(dest_path))
				continue

		#
		if args.progress:
			records = ProgressSeqParser(src_path, fmt=fmt, nested=True)
		else:
			records = SeqIO.parse(src_path, fmt)

		# Find the k-mers with quality info
		if args.threshold is not None:
			finders = kmers_from_records(records, spec,
			                             quality_threshold=args.threshold)

		else:
			# Just the k-mers themselves
			finders = kmers_from_records(records, spec)

		# Check any sequences actually found
		finders, no_seqs = iterator_empty(finders)
		if no_seqs:
			logger.warn('No sequences found in {}, bad file format?'
			            .format(src_path))

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
				count = write_kmer_vec(out_stream, finders)
			if args.count:
				count = write_kmer_counts(out_stream, finders)
			elif args.hist:
				count = write_kmer_hist(out_stream, finders)
			else:
				count = write_kmer_list(out_stream, finders)

			logger.debug('Found {} unique k-mers'.format(count))

		# Close file handle
		finally:
			if dest_path is not None:
				out_stream.close()


# Run main function when sub-command is invoked
parser.set_defaults(func=main)
