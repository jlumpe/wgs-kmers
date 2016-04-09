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
