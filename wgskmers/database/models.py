"""SQLAlchemy models"""

from sqlalchemy import Table, ForeignKey
from sqlalchemy import Column, Integer, String, Boolean, DateTime
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.declarative import declarative_base

from .sqla import TrackChangesMixin, JsonType, MutableJsonDict


__all__ = ['Base', 'KmerSetCollection', 'KmerSet', 'Genome', 'GenomeSet']


# SqlAlchemy declarative base
Base = declarative_base()


class Genome(Base, TrackChangesMixin):
	"""A reference genome"""

	__tablename__ = 'genomes'

	id = Column(Integer(), primary_key=True)
	description = Column(String(), nullable=False, unique=True)
	gb_db = Column(String())
	gb_id = Column(Integer(), unique=True)
	gb_acc = Column(Integer(), unique=True)
	gb_summary = Column(MutableJsonDict.as_mutable(JsonType))
	gb_taxid = Column(Integer())
	is_assembled = Column(Boolean(), nullable=False)
	organism = Column(String())
	filename = Column(String(), nullable=False, unique=True)
	file_format = Column(String(), nullable=False)
	extra = Column(MutableJsonDict.as_mutable(JsonType))


genome_set_assoc = Table(
	'genome_set_assoc',
	Base.metadata,
	Column('set_id', Integer(), ForeignKey('genome_sets.id'), primary_key=True),
	Column('genome_id', Integer(), ForeignKey('genomes.id'), primary_key=True),
)


class GenomeSet(Base, TrackChangesMixin):
	"""A set of genomes"""

	__tablename__ = 'genome_sets'

	id = Column(Integer(), primary_key=True)
	name = Column(String(), nullable=False, unique=True)
	description = Column(String())
	extra = Column(MutableJsonDict.as_mutable(JsonType))

	genomes = relationship('Genome', secondary=genome_set_assoc, lazy='dynamic',
	                       backref='genome_sets')


class KmerSetCollection(Base, TrackChangesMixin):
	"""A collection of k-mer counts/statistics for a set of genomes calculated
	with the same parameters.
	"""

	__tablename__ = 'kmer_collections'

	id = Column(Integer(), primary_key=True)
	title = Column(String(), nullable=False, unique=True)
	directory = Column(String(), nullable=False, unique=True)
	prefix = Column(String(), nullable=False)
	k = Column(Integer(), nullable=False)
	parameters = Column(MutableJsonDict.as_mutable(JsonType), nullable=False)
	format = Column(String(), nullable=False)
	extra = Column(MutableJsonDict.as_mutable(JsonType))


class KmerSet(Base):

	__tablename__ = 'kmer_sets'

	collection_id = Column(Integer(), ForeignKey('kmer_collections.id'),
	                       primary_key=True)
	genome_id = Column(Integer(), ForeignKey('genomes.id'),
	                   primary_key=True)
	dtype_str = Column(String(), nullable=False)
	has_counts = Column(Boolean(), nullable=False)
	count = Column(Integer(), nullable=False)
	filename = Column(String(), nullable=False)
	extra = Column(MutableJsonDict.as_mutable(JsonType))

	collection = relationship('KmerSetCollection',
	                          backref=backref('kmer_sets', lazy='dynamic'))
	genome = relationship('Genome', backref=backref('kmer_sets', lazy='dynamic'))
