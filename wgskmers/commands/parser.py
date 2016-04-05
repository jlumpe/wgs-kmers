"""Defines main argument parser"""

import argparse


# Create the parser
parser = argparse.ArgumentParser(description=__doc__)

# Global arguments
parser.add_argument('--debug', action='store_true', help='Print debug messages')

# Subparsers - imported by command modules to define sub-commands
subparsers = parser.add_subparsers(title='commands', dest='command')
