"""Functions for querying Genbank/Entrez and dealing with file formats"""

import os
from glob import glob
import re
import json
from xml.etree.cElementTree import parse as parse_xml
from cStringIO import StringIO

from Bio import Entrez


Entrez.email = 'mjlumpe@gmail.com' # Don't abuse this...


# Regular expressions to match accession numbers
# http://www.ncbi.nlm.nih.gov/Sequin/acc.html
accession_re = re.compile(r'[A-Z]+_?\d+(\.\d+)?')
accession_versioned_re = re.compile(r'[A-Z]+_?\d+\.\d+')
accession_unversioned_re = re.compile(r'[A-Z]+_?\d+')
accession_refseq_re = re.compile(r'[A-Z]{2}_\d+(\.\d+)?')


def is_accession(string):
	"""Checks if a string looks like an accession number"""
	return accession_re.match(string) is not None


def is_refseq(accession):
	"""Checks if a given accession number is a refseq number"""
	return accession_refseq_re.match(accession) is not None


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
		

def ids_for_accession(acc, db, field='Accession'):
	"""Query nucleotide GI numbers by accession"""

	if field is not None:
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

	attrs = dict(gb_id=int(id), gb_db=db, is_assembled=True)

	summary = get_summary(id=id, db=db)

	attrs['gb_summary'] = summary
	attrs['gb_acc'] = summary['accessionversion']
	attrs['organism'] = summary['organism']
	attrs['description'] = '[{accessionversion}] {title}'.format(**summary)

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
