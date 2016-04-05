


import logging


"""This package defines console commands for interacting with the project."""

from .parser import parser
from . import find


# Configure logger for warnings and debug messeges
logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger()



def main(args=None):
	"""Main entry point for the console commands.

	Deals with global arguments, then calls main function for correct
	sub-command.

	Args:
		args: list of str. Arguments to pass to main argument parser. If None,
			will use sys.argv.
	"""

	# Parse arguments
	args = parser.parse_args(args)

	# Debug mode
	if args.debug:
		logger.setLevel(logging.DEBUG)

	# Run command function
	args.func(args)
