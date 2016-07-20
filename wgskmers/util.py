"""Misc utility functions for the project"""
from builtins import next

import os
import shutil
import itertools


def rmpath(path):
	"""Removes a file, link, or directory by path."""

	if os.path.islink(path) or os.path.isfile(path):
		os.unlink(path)
		return True

	elif os.path.isdir(path):
		shutil.rmtree(path)
		return True

	else:
		return False


def kwargs_finished(kwargs):
	"""Raise an error if all elements of kwargs dict haven't been popped"""
	if kwargs:
		raise TypeError('Unknown keyword argument {}'.format(repr(list(kwargs.keys())[0])))


def iterator_empty(iterator):
	"""Checks if an iterator is empty, also returning a substitute iterator.

	Consumes first element of iterator, but stores it for the substitute
	iterator to yield first.
	"""

	try:
		first = [next(iterator)]
		return itertools.chain(first, iterator), False

	except StopIteration:
		return iter([]), True
