from __future__ import division
from builtins import range
from past.utils import old_div
from builtins import object
import numpy as np
from tqdm import tqdm

import numba as nb
import numba.cuda as nb_cuda

from .kmers import vec_to_coords, KmerCoordsCollection
from .util import kwargs_finished




class QueryMetric(object):
	"""A distance or similarity metric for kmer sets

	May not actually obey the triangle inequality, see the true_metric
	property.
	"""

	def __init__(self, py_func, title, return_type, is_distance, **kwargs):
		self.title = title

		self._py_func = py_func
		self._nb_array_funcs = dict()
		self._nb_coords_single = None
		self._nb_coords_vectorized = None

		self.return_type = np.dtype(return_type)
		self.is_distance = is_distance

		self.key = kwargs.pop('key', py_func.__name__)
		self.true_metric = kwargs.pop('true_metric', False)
		kwargs_finished(kwargs)

	def __call__(self, query_vecs, ref_vecs, out=None, target=None):
		"""Evaluate the metric on array-based data

		The target argument specifies the numba compiled target to use
		(cpu, parallel, or cuda), or 'python' for the pure python function.
		If None it will default to 'cpu', if that is unavailable the python
		function wil be used.
		"""

		if target is None:
			func = self._nb_array_funcs.get('cpu', self.py_func)

		elif target == 'python':
			func = self._py_func

		else:
			try:
				func = self.nb_array_funcs[target]

			except KeyError:
				raise RuntimeError(
					'Target "{}" not supported on this system'
					.format(target)
				)

		return func(query_vecs, ref_vecs, out=out)

	def coords(self, query_coords, ref_coords):
		"""Evaluate the metric on two sets in coordinate format"""
		return self._nb_coords_single(query_coords, ref_coords)

	def coords_multi(self, query_sets, ref_sets, out=None):
		"""Vectorized version that uses collection of query and/or ref coords"""
		query_multi = isinstance(query_sets, KmerCoordsCollection)
		refs_multi = isinstance(ref_sets, KmerCoordsCollection)

		# Convert arguments to flat arrays of coords and lengths
		if query_multi:
			qc_flat = query_sets.coords_array
			q_bounds = query_sets.bounds
		else:
			qc_flat = query_sets
			q_bounds = np.asarray([0, len(query_sets)])

		if refs_multi:
			rc_flat = ref_sets.coords_array
			r_bounds = ref_sets.bounds
		else:
			rc_flat = ref_sets
			r_bounds = np.asarray([0, len(ref_sets)])

		# Allocate output array if needed
		if out is None:
			out_shape = []
			if query_multi:
				out_shape.append(len(q_bounds) - 1)
			if refs_multi:
				out_shape.append(len(r_bounds) - 1)
			out = np.empty(out_shape, dtype=self.return_type)

		# Reshape out into a 2d array for the vectorized function
		out_reshaped = out.reshape((len(q_bounds) - 1, len(r_bounds) - 1))

		# Call vectorized version of function
		self._nb_coords_vectorized(qc_flat, q_bounds, rc_flat, r_bounds,
		                           out_reshaped)

		# Return output array in original shape
		return out

	def nb_array_func(self, *args, **kwargs):
		"""Creates decorator to register the numba guvectorized array function"""

		nb_targets = kwargs.pop('targets', ['cpu', 'parallel', 'cuda'])
		kwargs_finished(kwargs)

		# Create signatures
		nb_return_type = nb.from_dtype(self.return_type)
		signatures = [(nb.boolean[:], nb.boolean[:], nb_return_type[:])]

		# Create decorator
		def decorator(func):
			for target in nb_targets:

				if target == 'cuda' and not nb_cuda.is_available():
					continue

				guv_dec = nb.guvectorize(signatures, '(n),(n)->()',
				                         target=target, nopython=True)
				self._nb_array_funcs[target] = guv_dec(func)

			return func

		# If function passed directory, call decorator on it and return
		if args:
			func, = args
			return decorator(func)

		# Otherwise return the decorator
		else:
			return decorator

	def nb_coords_func(self, single_func):
		"""Decorator to register the numba coordinate function"""

		# JIT single version
		single_jit = nb.jit(nopython=True)(single_func)
		self._nb_coords_single = nb.jit(nopython=True)(single_func)

		# JIT vectorized version
		@nb.jit(nopython=True)
		def vectorized_coords_func(qc_flat, q_bounds, rc_flat, r_bounds, out):
			"""
			Args:
				qc_flat: Coordinates of each query set flattened into a
					1d array.
				q_bounds: Slices of qc_flat yielding each set of coordinates.
					Last value should be len(qc_flat).
				rc_flat: Coordinates of each reference set flattened into a
					1d array.
				r_bounds: Slices of qc_flat yielding each set of coordinates.
					Last value should be len(rc_flat).
			"""

			for q_i in range(len(q_bounds) - 1):
				q_begin, q_end = q_bounds[q_i:q_i+2]
				query = qc_flat[q_begin:q_end]

				for r_i in range(len(r_bounds) - 1):
					r_begin, r_end = r_bounds[r_i:r_i+2]
					ref = rc_flat[r_begin:r_end]

					out[q_i, r_i] = single_jit(query, ref)

		self._nb_coords_vectorized = vectorized_coords_func

		return self._nb_coords_single


query_metrics = dict()

def metric(*args, **kwargs):
	"""Decorator to create and register metric by python function"""
	def decorator(py_func):
		qm = QueryMetric(py_func, *args, **kwargs)
		query_metrics[qm.key] = qm
		return qm
	return decorator


##### Hamming Distance ######

@metric('Hamming distance', return_type=np.uint32, is_distance=True)
def hamming(query, ref, out=None):
	return (query != ref).sum(axis=-1, out=out)

@hamming.nb_array_func
def nb_hamming(query, ref, out):

	dist = 0

	for i in range(query.shape[0]):
		dist += query[i] ^ ref[i]

	out[0] = dist

@hamming.nb_coords_func
def coords_hamming(query, ref):

	n = query.shape[0]
	m = ref.shape[0]

	dist = 0

	i = 0
	j = 0
	while i < n and j < m:
		q = query[i]
		r = ref[j]

		if q != r:
			dist += 1

		if q <= r:
			i += 1

		if r <= q:
			j += 1

	dist += n - i
	dist += m - j

	return dist


##### Jaccard Index #####

@metric('Jaccard Index', return_type=np.float32, is_distance=False)
def jaccard(query, ref, out=None):

	if out is None:
		out_shape = np.broadcast(query, ref).shape[:-1]
		out = np.ndarray(out_shape, dtype=np.float32)

	np.logical_and(query, ref).sum(axis=-1, out=out)
	out /= np.logical_or(query, ref).sum(axis=-1)

	if out_shape == ():
		return out.item()
	else:
		return out

@jaccard.nb_array_func
def nb_jaccard(query, ref, out):

	union = 0
	intersection = 0

	for i in range(query.shape[0]):
		union += query[i] | ref[i]
		intersection += query[i] & ref[i]

	out[0] = intersection
	out[0] /= union

@jaccard.nb_coords_func
def coords_jaccard(query, ref):

	N = query.shape[0]
	M = ref.shape[0]

	union = 0

	i = 0
	j = 0
	while i < N and j < M:
		q = query[i]
		r = ref[j]

		union += 1

		if q <= r:
			i += 1

		if r <= q:
			j += 1

	union += N - i
	union += M - j

	return nb.float32(N + M - union) / union


##### Asymmetrical Jaccard #####

@metric('Asymmetrical Jaccard', return_type=np.float32, is_distance=False)
def asym_jacc(query, ref, out=None):

	if out is None:
		out_shape = np.broadcast(query, ref).shape[:-1]
		out = np.ndarray(out_shape, dtype=np.float32)

	np.logical_and(query, ref).sum(axis=-1, out=out)
	out /= ref.sum(axis=-1)

	if out_shape == ():
		return out.item()
	else:
		return out

@asym_jacc.nb_array_func
def nb_asym_jacc(query, ref, out):

	ref_weight = 0
	intersection = 0

	for i in range(query.shape[0]):
		ref_weight += ref[i]
		intersection += query[i] & ref[i]

	out[0] = intersection
	out[0] /= ref_weight

@asym_jacc.nb_coords_func
def coords_asym_jacc(query, ref):

	N = query.shape[0]
	M = ref.shape[0]

	intersection = 0

	i = 0
	j = 0
	while i < N and j < M:
		q = query[i]
		r = ref[j]

		if q == r:
			intersection += 1

		if q <= r:
			i += 1

		if r <= q:
			j += 1

	return nb.float32(intersection) / M


##### Parallelized query functions #####

class QueryWorker(object):

	@classmethod
	def init(cls, query, dest, metric_names, loader):
		cls.query = query
		cls.dest = dest
		cls.metrics = [query_metrics[name] for name in metric_names]
		cls.loader = loader

		# Ignore floating point divide by zero
		np.seterr(divide='ignore', invalid='ignore')

	@classmethod
	def calc_scores(cls, args):
		index, ref_sets = args

		ref_vecs = cls.loader.load_array(ref_sets, dtype=bool)

		for i, metric in enumerate(cls.metrics):
			scores = metric(cls.query.np_array[:, None, :], ref_vecs[None, ...])
			cls.dest[i, index:index + len(ref_sets), :] = scores.T


def mp_query(query, db, collection, ref_sets, metrics, **kwargs):
	import multiprocessing as mp
	import ctypes
	import wgskmers.multiprocess as kmp

	kmp.enable_method_pickling()

	nworkers = kwargs.pop('nworkers', mp.cpu_count())
	chunk_size = kwargs.pop('chunk_size', 5)
	progress = kwargs.pop('progress', False)
	kwargs_finished(kwargs)

	loader = db.get_kmer_loader(collection)

	# Scores output as shared memory
	scores_shape = (len(metrics), len(ref_sets), query.shape[0])
	scores = kmp.SharedNumpyArray(ctypes.c_float, scores_shape)

	# Query array as shared memory
	assert query.dtype == np.bool
	query_arr = kmp.SharedNumpyArray(ctypes.c_bool, query.shape)
	query_arr[:] = query

	# Create pool
	init_args = (query_arr, scores, metrics, loader)
	pool = mp.Pool(processes=nworkers, initializer=QueryWorker.init,
	               initargs=init_args)

	# Start tasks
	chunks = [(i, ref_sets[i:i + chunk_size]) for i
	          in range(0, len(ref_sets), chunk_size)]
	results = pool.imap_unordered(QueryWorker.calc_scores, chunks)

	# Monitor progress
	if progress:
		pbar = tqdm(total=len(ref_sets), desc='Querying reference database')
		with pbar:
			for r in results:
				pbar.update(chunk_size)

	pool.close()
	pool.join()

	# Return array
	return scores.np_array


class CoordsQueryWorker(object):

	@classmethod
	def init(cls, query_coords, loader, metric_names, dest):
		cls.query_coords = query_coords
		cls.loader = loader
		cls.metrics = [query_metrics[name] for name in metric_names]
		cls.dest = dest

		# Ignore floating point divide by zero
		np.seterr(divide='ignore', invalid='ignore')

	@classmethod
	def calc_score(cls, args):
		ref_idx, ref_set = args

		ref_coords = cls.loader.load_coords(ref_set)

		for k, metric in enumerate(cls.metrics):
			for j, query_coords in enumerate(cls.query_coords):
				cls.dest[k, ref_idx, j] = metric.coords(query_coords, ref_coords)


def mp_query_coords(query, db, collection, ref_sets, metrics, **kwargs):
	import multiprocessing as mp
	import ctypes
	import wgskmers.multiprocess as kmp

	kmp.enable_method_pickling()

	nworkers = kwargs.pop('nworkers', mp.cpu_count())
	progress = kwargs.pop('progress', False)
	kwargs_finished(kwargs)

	loader = db.get_kmer_loader(collection)

	# Scores output as shared memory
	scores_shape = (len(metrics), len(ref_sets), query.shape[0])
	scores = kmp.SharedNumpyArray(ctypes.c_float, scores_shape)

	# Query coords in shared memory
	query_coords = kmp.SharedKmerCoordsCollection.empty(query.sum(axis=1))
	for i in range(query.shape[0]):
		query_coords[i] = vec_to_coords(query[i, :])

	# Create pool
	init_args = (query_coords, loader, metrics, scores)
	pool = mp.Pool(processes=nworkers, initializer=CoordsQueryWorker.init,
	               initargs=init_args)

	# Start workers
	tasks = list(enumerate(ref_sets))
	results = pool.imap_unordered(CoordsQueryWorker.calc_score, tasks)

	# Monitor progress
	if progress:
		pbar = tqdm(total=len(tasks), desc='Querying reference database')
		with pbar:
			for r in results:
				pbar.update(1)

	pool.close()
	pool.join()

	# Return array
	return scores.np_array
