"""SQLAlchemy custom types, utility stuff"""

import datetime
import json

from sqlalchemy import Column, String, DateTime
from sqlalchemy.event import listen
from sqlalchemy.types import TypeDecorator
from sqlalchemy.ext.mutable import Mutable
import sqlalchemy.orm.session


class ReadOnlySession(sqlalchemy.orm.session.Session):
	"""Session class that doesn't allow flushing/committing"""

	def _flush(self, *args, **kwargs):
		raise RuntimeError('This sessison is read-only')


class TrackChangesMixin(object):
	"""Mixin for SQLAlchemy models that tracks when updates are made"""

	created_at = Column(DateTime())
	updated_at = Column(DateTime())

	@classmethod
	def insert_time_callback(cls, mapper, connection, instance):
		now = datetime.datetime.utcnow()
		instance.created_at = now
		instance.updated_at = now

	@classmethod
	def update_time_callback(cls, mapper, connection, instance):
		now = datetime.datetime.utcnow()
		instance.updated_at = now

	@classmethod
	def __declare_last__(cls):
		"""Called after mapper configured, register listeners"""
		listen(cls, 'before_insert', cls.insert_time_callback)
		listen(cls, 'before_update', cls.update_time_callback)


class MethodFactoryMeta(type):
	"""Metaclass for classes than generate their methods through a function"""

	def __new__(cls, name, bases, dct):
		if 'make_methods' in dct:
			make_func = dct['make_methods']

			if isinstance(make_func, staticmethod):
				new_methods = make_func.__func__()
			else:
				raise TypeError('make_methods must be static method')

			dct.update(new_methods)
		return type.__new__(cls, name, bases, dct)


def wrap_mutation_method(method, name):
	def wrapper(self, *args, **kwargs):
		self.changed()
		return method(self, *args, **kwargs)

	wrapper.func_name = name
	wrapper.__doc__ = method.__doc__

	return wrapper


class MutableList(Mutable, list):
	__metaclass__ = MethodFactoryMeta

	@staticmethod
	def make_methods():

		list_methods = ['append', 'extend', 'index', 'insert',
		                'pop', 'remove', 'reverse', 'sort', '__setitem__',
		                '__setslice__', '__delitem__', '__delslice__',
		                '__iadd__', '__imul__']

		methods = dict()

		for method_name in list_methods:
			list_method = getattr(list, method_name)
			wrapped = wrap_mutation_method(list_method, method_name)
			methods[method_name] = wrapped

		return methods

	@classmethod
	def coerce(cls, key, value):
		if not isinstance(value, MutableList):
			if isinstance(value, list):
				return MutableList(value)
			else:
				return Mutable.coerce(key, value)
		else:
			return value


class StringList(TypeDecorator):
	"""
	SqlAlchemy type for encoding a list of strings and storing
	in a text field.
	"""
	impl = String

	def process_bind_param(self, value, dialect):
		if value is not None:
			assert all(isinstance(e, basestring) for e in value)
			return json.dumps(value, separators=(',', ':'))
		else:
			return None

	def process_result_value(self, value, dialect):
		if value is not None:
			return json.loads(value)
		else:
			return None

MutableList.associate_with(StringList)


class MutableDict(Mutable, dict):
	__metaclass__ = MethodFactoryMeta

	@staticmethod
	def make_methods():

		dict_methods = ['clear', 'pop', 'popitem', 'setdefault', 'update',
		                '__setitem__', '__delitem__']

		methods = dict()

		for method_name in dict_methods:
			dict_method = getattr(dict, method_name)
			wrapped = wrap_mutation_method(dict_method, method_name)
			methods[method_name] = wrapped

		return methods

	@classmethod
	def coerce(cls, key, value):
		if not isinstance(value, MutableDict):
			if isinstance(value, dict):
				return MutableDict(value)
			else:
				return Mutable.coerce(key, value)
		else:
			return value


class FlatDict(TypeDecorator):
	"""
	SqlAlchemy type for encoding a dict of strings, numbers, bools, and Nones
	and storing it in a text field.
	"""
	impl = String

	def process_bind_param(self, value, dialect):
		if value is not None:
			for k, v in value.iteritems():
				assert isinstance(k, basestring)
				assert isinstance(v, (basestring, int, long, float, None))
			return json.dumps(value, separators=(',', ':'))
		else:
			return None

	def process_result_value(self, value, dialect):
		if value is not None:
			return json.loads(value)
		else:
			return None

MutableDict.associate_with(FlatDict)
