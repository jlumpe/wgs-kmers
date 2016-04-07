"""Stand-alone command to find k-mers in a file or set of files.

By default, the script will output a sorted list of all k-mers found.
"""

import os
import sys
import logging
from collections import Counter

import click
from Bio import SeqIO

from .util import tqdm, ProgressSeqParser, iterator_empty
from wgskmers.kmers import KmerSpec, nucleotides


logger = logging.getLogger()


# Default arguments
default_k = 16
default_prefix = 'ATGAC'


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

	# Sort k-mers by count
	kmers_sorted = sorted(kmer_counts.items(), key=lambda i: i[1],
	                      reverse=True)

	# Write
	for kmer, count in kmers_sorted:
		stream.write('{} {}\n'.format(kmer, count))

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


def find_batch_files(directory):
	"""Finds source files in directory for batch processing"""
	paths = []

	for fname in os.listdir(directory):
		path = os.path.join(directory, fname)
		if os.path.isfile(path):
			paths.append(path)

	return paths


#-----------------------------------------------------------------------------
# Command
#-----------------------------------------------------------------------------

@click.command(name='find')

# Search parameters
@click.option('-k', '--k', default=default_k, show_default=True,
	help='Length of k-mers to find')
@click.option('-q', '--prefix', default=default_prefix, show_default=True,
	help='Target sequence to find', metavar='SEQUENCE')
@click.option('-t', '--threshold', type=int, required=False,
	help='Filter k-mers containing PHRED scores below this value')

# Output format
@click.option('-l', '--list', 'output_format', flag_value='list', default=True,
	help='Output sorted list of k-mers')
@click.option('-c', '--count', 'output_format', flag_value='counts',
	help='Output number of occurrences of each k-mer')
@click.option('-H', '--hist', 'output_format', flag_value='hist',
	help='Output histogram of counts')
@click.option('-B', '--bool', 'output_format', flag_value='bool',
	help='Output boolean vector')

# Input/output options
@click.option('-f', '--format', default=None,
	type=click.Choice(['fasta', 'fastq', 'fastq-sanger', 'fastq-solexa',
	                   'fastq-illumina']),
	help='File format, as argument to Bio.SeqIO.parse. If omitted, will infer '
	     'from extension of first file encountered.')
@click.option('-b', '--batch', is_flag=True, help='Run in batch mode.')
@click.option('-o', '--overwrite', is_flag=True,
	help='Overwrite existing output files')

# Other
@click.option('-p', '--progress', is_flag=True, help='Display progress')

# Source and destination files
@click.argument('src', type=click.Path(exists=True))
@click.argument('dest', type=click.Path(), required=False)

# The actual command...
def find_command(src, dest, prefix, k, **kwargs):
	"""
	Find k-mers in a file or set of files. In single-file mode, read sequence
	file [SRC] and write output to [DEST] (or stdout). In batch mode, process
	each file in directory given by [SRC], writing output to separate files
	in directory given by [DEST].
	"""

	show_progress = kwargs.pop('progress', False)
	output_format = kwargs.pop('output_format', 'list')
	threshold = kwargs.pop('threshold', None)
	batch_mode = kwargs.pop('batch', False)
	file_format = kwargs.pop('format', None)
	overwrite_output = kwargs.pop('overwrite', False)

	# Check prefix valid
	prefix = prefix.upper()
	if not set(prefix).issubset(nucleotides):
		raise ValueError('Prefix contains invalid characters')

	# Check k positive
	if k <= 0:
		raise ValueError('K must be positive')

	# Kmer spec
	spec = KmerSpec(k, prefix)

	# Check tqdm present
	if show_progress:
		if tqdm is None:
			raise RuntimeError('tqdm package required to display progress')

	# Absolute source and destination paths
	src = os.path.realpath(src)
	dest = os.path.realpath(dest) if dest is not None else None

	# Output file extension
	if output_format == 'list':
		out_ext = '.kmers.txt'
	elif output_format == 'bool':
		out_ext = '.kmer_vec'
	elif output_format == 'hist':
		out_ext = '.hist.txt'
	elif output_format == 'counts':
		out_ext = '.counts.txt'
	else:
		raise ValueError('Bad output format {}'.format(repr(output_format)))

	if threshold is not None:
		out_ext = '-t{}'.format(threshold) + out_ext

	# Batch mode
	if batch_mode:

		# Check arguments compatible with batch mode
		if dest is None:
			raise ValueError('Must give destination directory in batch mode')

		# Find files in source directory (raises OSError if doesn't exist)
		src_paths = find_batch_files(src)
		if not src_paths:
			raise RuntimeError('No files found in {}'.format(src))

		# Create output directory if necessary
		if not os.path.isdir(dest):
			os.mkdir(dest)

		# Destination file paths
		dest_paths = [make_dest_path(src, dest, ext=out_ext)
		              for src in src_paths]

	# Single-file mode
	else:

		src_paths = [src]

		# Get destination path
		if dest is not None:

			# If given a directory, output to file with same name
			if os.path.isdir(dest):
				dest_paths = [make_dest_path(src, dest, ext=out_ext)]
			else:
				dest_paths = [dest]

		else:
			dest_paths = [None]

	# Get format
	if file_format is None:
		file_format = infer_format(src_paths[0])
		if file_format is None:
			raise ValueError("Couldn't infer format for {}"
				.format(src_paths[0]))
		logger.debug('Inferring format "{}" from {}'
		             .format(file_format, src_paths[0]))

	# Should be ok, loop over the files
	paths_iter = zip(src_paths, dest_paths)
	if batch_mode and show_progress:
		paths_iter = tqdm(paths_iter, unit='files')

	for src_path, dest_path in paths_iter:
		
		logger.debug('Processing source file {}'.format(src_path))

		# Check if output file exists
		if dest_path is not None and os.path.exists(dest_path):
			if overwrite_output:
				logger.warn('Overwriting output file {}'.format(dest_path))
			else:
				logger.warn('Refusing to overwrite {}'.format(dest_path))
				continue

		# Wrap parse iterator in progress bar
		if show_progress:
			records = ProgressSeqParser(src_path, fmt=file_format, nested=True)
		else:
			records = SeqIO.parse(src_path, file_format)

		# Find the k-mers with quality info
		if threshold is not None:
			finders = kmers_from_records(records, spec,
			                             quality_threshold=threshold)

		else:
			# Just the k-mers themselves
			finders = kmers_from_records(records, spec)

		# Check any sequences actually found
		finders, no_seqs = iterator_empty(finders)
		if no_seqs:
			logger.warn('No sequences found in {}, bad file format?'
			            .format(src_path))

		# Output stream - file (text/binary) or StringIO
		if dest_path is None:
			from cStringIO import StringIO
			out_stream = StringIO()
		elif output_format == 'bool':
			out_stream = open(dest_path, 'wb')
		else:
			out_stream = open(dest_path, 'w')

		# Write output in try block because we can't use a with statement here
		try:

			# Output based on arguments
			if output_format == 'list':
				count = write_kmer_list(out_stream, finders)
			elif output_format == 'bool':
				count = write_kmer_vec(out_stream, finders)
			elif output_format == 'counts':
				count = write_kmer_counts(out_stream, finders)
			elif output_format == 'hist':
				count = write_kmer_hist(out_stream, finders)
			else:
				assert False, 'You shouldn\'t be here...'

			logger.debug('Found {} unique k-mers'.format(count))

			# Write to stdout
			if dest_path is None:
				if output_format == 'hist':
					click.echo(out_stream.getvalue())
				else:
					click.echo_via_pager(out_stream.getvalue())

		finally:
			# Close file handle
			out_stream.close()
