"""SQLAlchemy custom types, utility stuff"""

import datetime
import json
import collections
import weakref

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


# Python types corresponding to non-collection types storable in JSON
jsonable_scalars = (int, long, float, basestring, bool, type(None))


class MutableJsonCollection(Mutable):
	"""ABC for SQLAlchemy JSON collection type, supporting mutation tracking

	When nested under another collection, keeps a weak reference to its
	parent so that mutation notifications can be propagated back up to the
	root collection.
	"""

	def __init__(self, parent=None):
		if parent is not None:
			self._parent = weakref.proxy(parent)
		else:
			self._parent = None

	def __getstate__(self):
		"""For pickling"""
		return self.as_builtin()

	def __setstate__(state):
		"""For unpickling"""
		self.__init__(state)

	def changed(self):
		super(MutableJsonCollection, self).changed()

		# Call changes() method on parent (if exists)
		if self._parent is not None:
			try:
				self._parent.changed()
			except weakref.ReferenceError:
				self._parent = None

	def as_builtin(self):
		"""Return version of collection as builtin Python type

		(for JSON serialization)
		"""
		raise NotImplementedError()

	def _transform_element(self, elem):
		"""Transforms python types into MutableJsonCollection where possible"""

		if isinstance(elem, MutableJsonCollection):
			return elem

		elif isinstance(elem, jsonable_scalars):
			return elem

		elif isinstance(elem, collections.Mapping):
			return MutableJsonDict(elem, parent=self)

		elif isinstance(elem, collections.Sequence):
			return MutableJsonList(elem, parent=self)

		else:
			raise TypeError('{} is not a JSONable type'.format(type(elem)))

	@classmethod
	def _element_as_builtin(cls, elem):
		"""Converts an element to builtin type"""
		if isinstance(elem, MutableJsonCollection):
			return elem.as_builtin()

		else:
			return elem



class MutableJsonList(MutableJsonCollection, collections.MutableSequence):
	"""List-like object corresponding to JSON array in SQLAlchemy.

	Mutations will be tracked in this object and all nested collections.
	"""

	def __init__(self, sequence, parent=None):
		# Recursively convert any nested collections to MutableJsonCollection
		self._list = map(self._transform_element, sequence)

		MutableJsonCollection.__init__(self, parent)

	def __getitem__(self, index):
		return self._list[index]

	def __setitem__(self, index, value):
		self._list[index] = self._transform_element(value)
		self.changed()

	def __delitem__(self, index):
		del self._list[index]
		self.changed()

	def __iter__(self):
		return iter(self._list)

	def __len__(self):
		return len(self._list)

	def __repr__(self):
		return repr(self._list)

	def insert(self, index, value):
		self._list.insert(index, self._transform_element(value))
		self.changed()

	def as_builtin(self):
		return map(self._element_as_builtin, self._list)

	@classmethod
	def coerce(cls, key, value):
		if not isinstance(value, MutableJsonList):
			if isinstance(value, collections.Sequence):
				return MutableJsonList(value)
			else:
				return Mutable.coerce(key, value)
		else:
			return value


class MutableJsonDict(MutableJsonCollection, collections.MutableMapping):
	"""Dict-like object corresponding to JSON object in SQLAlchemy.

	Mutations will be tracked in this object and all nested collections.
	"""

	def __init__(self, mapping, parent=None):
		# Recursively convert any nested collections to MutableJsonCollection
		self._dict = {k: self._transform_element(v) for k, v
		              in dict(mapping).iteritems()}

		MutableJsonCollection.__init__(self, parent)

	def __getitem__(self, key):
		return self._dict[key]

	def __setitem__(self, key, value):
		if not isinstance(key, basestring):
			raise TypeError('Key must be string')

		self._dict[key] = self._transform_element(value)
		self.changed()

	def __delitem__(self, key):
		del self._dict[key]
		self.changed()

	def __iter__(self):
		return iter(self._dict)

	def __len__(self):
		return len(self._dict)

	def __repr__(self):
		return repr(self._dict)

	def as_builtin(self):
		return {k: self._element_as_builtin(v) for k, v
		        in self._dict.iteritems()}

	@classmethod
	def coerce(cls, key, value):
		if not isinstance(value, MutableJsonDict):
			if isinstance(value, collections.Mapping):
				return MutableJsonDict(value)
			else:
				return Mutable.coerce(key, value)
		else:
			return value


class MutableJsonCollectionEncoder(json.JSONEncoder):
	"""JSON encoder capable of serializing MutableJsonCollection objects"""

	def default(self, obj):
		if isinstance(obj, MutableJsonCollection):
			return obj.as_builtin()
		else:
			return super(MutableJsonCollectionEncoder, self).default(obj)


class JsonType(TypeDecorator):
	"""SQLA column type for JSON data"""

	impl = String

	def process_bind_param(self, value, dialect):
		if value is not None:
			return json.dumps(value, separators=(',', ':'),
			                  cls=MutableJsonCollectionEncoder)
		else:
			return None

	def process_result_value(self, value, dialect):
		if value is not None:
			json_val = json.loads(value)

			if isinstance(json_val, dict):
				return MutableJsonDict(json_val)
			elif isinstance(json_val, list):
				return MutableJsonList(json_val)
			else:
				return json_val

		else:
			return None
