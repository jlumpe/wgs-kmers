"""Functions for querying Genbank/Entrez and dealing with file formats"""

import os
from glob import glob
import re
import json
from xml.etree.cElementTree import parse as parse_xml
from cStringIO import StringIO

from Bio import Entrez


Entrez.email = 'mjlumpe@gmail.com' # Don't abuse this...


ncbi_url = 'http://www.ncbi.nlm.nih.gov'


# Regular expressions to match accession numbers
# http://www.ncbi.nlm.nih.gov/Sequin/acc.html
# Looks like there are between 5-10 digits in the numerical portion
acc_re = re.compile(r'[A-Z]+_?\d{5,10}(?:\.\d+)?')
acc_versioned_re = re.compile(r'[A-Z]+_?\d{5,10}\.\d+')
acc_unversioned_re = re.compile(r'[A-Z]+_?\d{5,10}')
acc_refseq_re = re.compile(r'[A-Z]{2}_\d{5,10}(?:\.\d+)?')

# Same as acc_re but meant to be used for searching within a larger string.
# Contains negative look(ahead|behind)s asserting that the number is not
# immediately preceded/followed by an alphanumeric character to prevent false
# positives.
acc_search_re = re.compile('(?<![A-Za-z0-9])' + acc_re.pattern +
                           '(?![A-Za-z0-9])')


def is_accession(string):
	"""Checks if a string looks like an accession number"""
	return acc_re.match(string) is not None


def is_refseq(accession):
	"""Checks if a given accession number is a refseq number"""
	return acc_refseq_re.match(accession) is not None


def extract_acc(string, one_only=False):
	"""Extracts a genbank accession number from a larger string

	Meant to be used to find accession numbers in file names or similar.
	"""
	matches = acc_search_re.findall(string)

	if len(matches) == 0:
		return None
	elif len(matches) == 1 or not one_only:
		return matches[0]
	else:
		return None


class StupidBiopythonUrlopenReplacer(object):
	"""It's stupid"""
	def __init__(self, timeout):
		self.timeout = timeout
		self.old_urlopen = None

	def __enter__(self):
		if self.timeout is not None:
			self.old_urlopen = Entrez._urlopen
			def new_urlopen(url, data=None, timeout=self.timeout):
				return self.old_urlopen(url, data=data, timeout=timeout)
			Entrez._urlopen = new_urlopen

	def __exit__(self, *excinfo):
		if self.timeout is not None:
			Entrez._urlopen = self.old_urlopen


def ezutils(util, *args, **kwargs):
	"""Entrez eutils with automatic parsing of response, other cool stuff"""

	parse = kwargs.pop('parse', True)
	timeout = kwargs.pop('timeout', None)

	if 'id' in kwargs:
		kwargs['id'] = str(kwargs['id'])

	with StupidBiopythonUrlopenReplacer(timeout):
		response = getattr(Entrez, util)(*args, **kwargs)
	
	if parse:
		mime = response.headers.type
		if mime == 'application/json':
			return json.load(response)
		elif mime == 'text/xml':
			return parse_xml(response)
		elif mime == 'text/plain':
			return response
		else:
			raise ValueError('Unknown mimetype "{}"'.format(mime))
	else:
		return response
		

def ids_for_accession(acc, db, field=None):
	"""Query nucleotide GI numbers by accession"""

	if field is None:
		if db == 'assembly':
			field = 'Assembly Accession'
		else:
			field = 'Accession'

	if field is not False:
		term = '"{}"[{}]'.format(acc, field)
	else:
		term = acc

	result = ezutils('esearch', db=db, retmode='json', term=term)
	return map(int, result['esearchresult']['idlist'])


def get_summary(id, db):
	"""Fetch document summary as JSON

	For whatever reason, JSON contains more information.
	"""
	sum_res = ezutils('esummary', db=db, id=id, retmode='json')
	return sum_res['result'][str(id)]


def genome_attrs_from_gi(id, db):
	"""Get best values for Genome model attributes given genbank GI"""

	attrs = dict(gb_id=int(id), gb_db=db)

	summary = get_summary(id=id, db=db)

	attrs['gb_summary'] = summary

	if 'accessionversion' in summary:
		attrs['gb_acc'] = summary['accessionversion']
	elif 'assemblyaccession' in summary:
		attrs['gb_acc'] = summary['assemblyaccession']

	try:
		attrs['organism'] = summary['organism']
	except KeyError:
		pass

	try:
		attrs['description'] = '[{}] {}'.format(attrs['gb_acc'], summary['title'])
	except KeyError:
		pass

	return attrs


def fetch_fasta(id, db, **kwargs):

	if db == 'nuccore':
		return ezutils('efetch', db='nuccore', id=id, rettype='fasta',
		               retmode='text')

	elif db == 'assembly':
		linkname = kwargs.pop('linkname', 'assembly_nuccore_refseq')
		link_resp = ezutils('elink', db='nuccore', dbfrom='assembly', id=id,
		                    linkname=linkname, retmode='json')
		linkset, = link_resp['linksets']
		linksetdb, = linkset['linksetdbs']

		contents = StringIO()

		for link_id in sorted(linksetdb['links']):
			fetch_resp = fetch_fasta(link_id, 'nuccore')
			contents.write(fetch_resp.read())

		contents.seek(0)
		return contents

	else:
		raise ValueError('DB not supported')


db_paths = dict(
	assembly='assembly',
	nuccore='nuccore',
)

def get_record_url(acc_or_id, db):
	"""Gets URL for a genbank record by database and either accession or id"""

	try:
		db_path = db_paths[db]
	except KeyError:
		raise ValueError('Don\'t know how to get url for db {}'.format(db))

	return '/'.join([ncbi_url, db, str(acc_or_id)])
