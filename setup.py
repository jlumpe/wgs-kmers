"""setuptools installation script for wgs-kmers project"""

from setuptools import setup, find_packages

setup(
	name='wgs-kmers',
	version='0.0.1',
	description='Kmer finding in whole genome sequencing data',
	author='Jared Lumpe',
	author_email='mjlumpe@gmail.com',
	packages=find_packages(),
	install_requires=[
		'numpy >= 1.10',
		'biopython >= 1.66',
		'click >= 6.6',
		'sqlalchemy >= 1.0.11',
		'tqdm >= 3.4.0',
		'appdirs >= 1.4.0',
		'alembic >=  0.8.5',
		'future >= 0.15',
	],
	entry_points={
		'console_scripts': [
			'kmers = wgskmers.commands:cli',
		]
	},
	include_package_data=True,
)
