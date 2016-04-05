"""setuptools setup script for wgs-kmers project"""

from setuptools import setup, find_packages

setup(
	name='wgs-kmers',
	version='0.0.1',
	description='Kmer finding in whole genome sequencing data',
	author='Jared Lumpe',
	author_email='jlumpe@primitybio.com',
	packages=find_packages(),
	install_requires=[
		'biopython >= 1.66',
	],
	extras_require={
		'tqdm': 'tqdm >= 3.4.0'
	},
	include_package_data=True,
	entry_points={
		'console_scripts': [
			'kmers = wgskmers.commands:main',
		]
	}
)
