"""Utility functions and classes for multiprocessing"""

import sys
import traceback
import multiprocessing as mp
import copy_reg
import types
import ctypes

import numpy as np

from wgskmers.kmers import KmerSpec, KmerCoordsCollection


def _pickle_method(method):
	func_name = method.im_func.func_name
	obj = method.im_self
	cls = method.im_class

	if func_name.startswith('__') and not func_name.endswith('__'):
		cls_name = cls.__name__.lstrip('_')
		func_name = '_' + cls_name + func_name

	if cls is type:
		return _unpickle_method, (func_name, obj, obj)
	else:
		return _unpickle_method, (func_name, obj, cls)


def _unpickle_method(func_name, obj, cls):
	try:
		unbound_method = getattr(cls, func_name)
	except KeyError:
		raise ValueError('{} has no method {}'.format(cls, func_name))

	return unbound_method.__get__(obj, cls)


def enable_method_pickling():
	"""Needed to map instance/class methods with process pools"""
	copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)


class SharedNumpyArray(object):
	"""Pickleable numpy wrapper around shared-memory array"""

	def __init__(self, ctype, shape):
		self.ctype = ctype

		self.shared_array = mp.Array(ctype, np.prod(shape))
		self.np_array = self._make_np_array(self.shared_array.get_obj(), ctype,
		                                    shape)

	def __getattr__(self, attr):
		return getattr(self.np_array, attr)

	def __getitem__(self, index):
		return self.np_array[index]

	def __setitem__(self, index, value):
		self.np_array[index] = value

	def get_lock(self):
		return self.shared_array.get_lock()

	def __getstate__(self):
		return (self.shared_array, self.ctype, self.shape)

	def __setstate__(self, state):
		shared_array, ctype, shape = state

		self.shared_array = shared_array
		self.ctype = ctype

		self.np_array = self._make_np_array(shared_array, ctype, shape)

	@classmethod
	def _make_np_array(cls, shared_array, ctype, shape):
		np_array = np.frombuffer(shared_array, dtype=np.dtype(ctype))
		np_array.shape = shape
		return np_array

class SharedKmerCoordsCollection(KmerCoordsCollection):

	def __init__(self, shared_array, bounds):
		self.shared_array = shared_array

		super(SharedKmerCoordsCollection, self).__init__(shared_array.np_array,
		                                                 bounds)

	def __getstate__(self):
		return (self.shared_array, self.bounds)

	def __setstate__(self, state):
		self.shared_array, self.bounds = state

	@classmethod
	def empty(cls, lengths):
		bounds = cls._make_bounds(lengths)

		shared_array = SharedNumpyArray(ctypes.c_uint32, (int(bounds[-1]),))

		return cls(shared_array, bounds)
