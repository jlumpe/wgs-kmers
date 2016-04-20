"""Misc utility functions for the project"""

import os
import shutil


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
		raise TypeError('Unknown keyword argument {}'.format(repr(kwargs.keys()[0])))
