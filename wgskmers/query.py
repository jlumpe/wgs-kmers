import numpy as np
from tqdm import tqdm

from wgskmers.util import kwargs_finished
import numba as nb
import numba.cuda as nb_cuda

from .kmers import vec_to_coords


# Numba function signatures
coords_nb_signature = [
	(nb.uint32[:], nb.uint32[:], nb.float32[:]),
	(nb.uint64[:], nb.uint64[:], nb.float32[:]),
	(nb.int32[:], nb.int32[:], nb.float32[:]),
	(nb.int64[:], nb.int64[:], nb.float32[:]),
]


class QueryMetric(object):
	def __init__(self, key, title, py_func, is_distance):
		self.key = key
		self.title = title
		self.py_func = py_func
		self.nb_funcs = dict()
		self.coords_funcs = dict()
		self.is_distance = is_distance

	def __call__(self, query_vec, ref_vec):
		return self.py_func(query_vec, ref_vec)

	def nb(self, query_vec, ref_vec, out=None):
		return self._call_nb('cpu', query_vec, ref_vec, out)

	def nb_parallel(self, query_vec, ref_vec, out=None):
		return self._call_nb('parallel', query_vec, ref_vec, out)

	def nb_cuda(self, query_vec, ref_vec, out=None):
		return self._call_nb('cuda', query_vec, ref_vec, out)

	def coords(self, query_vec, ref_vec, out=None):
		return self._call_coords('cpu', query_vec, ref_vec, out)

	def coords_parallel(self, query_vec, ref_vec, out=None):
		return self._call_coords('parallel', query_vec, ref_vec, out)

	def coords_cuda(self, query_vec, ref_vec, out=None):
		return self._call_coords('cuda', query_vec, ref_vec, out)

	def _call_nb(self, target, query_vec, ref_vec, out=None):
		func = self.nb_funcs[target]
		if func is None:
			raise RuntimeError(
				'Target "{}" not supported on this system'
				.format(target)
			)
		else:
			return func(query_vec, ref_vec, out)

	def _call_coords(self, target, query_vec, ref_vec, out=None):
		func = self.coords_funcs[target]
		if func is None:
			raise RuntimeError(
				'Target "{}" not supported on this system'
				.format(target)
			)
		else:
			return func(query_vec, ref_vec, out)

	def numba_func(self, signature, layout, **kwargs):
		nb_targets = kwargs.pop('nb_targets', ['cpu', 'parallel', 'cuda'])
		assert not kwargs

		def decorator(func):
			for target in nb_targets:

				if target == 'cuda' and not nb_cuda.is_available():
					self.nb_funcs[target] = None

				else:
					guv_dec = nb.guvectorize(signature, layout, target=target,
					                         nopython=True)
					self.nb_funcs[target] = guv_dec(func)

			return func

		return decorator

	def coords_func(self, signature, layout, **kwargs):
		nb_targets = kwargs.pop('nb_targets', ['cpu', 'parallel', 'cuda'])
		assert not kwargs

		def decorator(func):
			for target in nb_targets:

				if target == 'cuda' and not nb.cuda.is_available():
					self.coords_funcs[target] = None

				else:
					guv_dec = nb.guvectorize(signature, layout, target=target,
					                         nopython=True)
					self.coords_funcs[target] = guv_dec(func)

			return func

		return decorator


query_metrics = dict()

def metric(title, is_distance):
	def decorator(py_func):
		key = py_func.__name__
		qm = QueryMetric(key, title, py_func, is_distance)
		query_metrics[key] = qm
		return qm
	return decorator


@metric('Hamming distance', is_distance=True)
def hamming(query, ref, out=None):
	return (query != ref).sum(axis=-1, out=out)

@hamming.numba_func([(nb.bool_[:], nb.bool_[:], nb.uint32[:])],
                    '(n),(n)->()')
def nb_hamming(query, ref, out):

	dist = 0

	for i in range(query.shape[0]):
		dist += query[i] ^ ref[i]

	out[0] = dist

@hamming.coords_func(coords_nb_signature, '(n),(m)->()')
def coords_hamming(query, ref, out):

	n = query.shape[0]
	m = ref.shape[0]

	out[0] = 0

	i = 0
	j = 0
	while i < n and j < m:
		q = query[i]
		r = ref[j]

		if q != r:
			out[0] += 1

		if q <= r:
			i += 1

		if r <= q:
			j += 1

	out[0] += n - i
	out[0] += m - j


@metric('Jaccard Index', is_distance=False)
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

@jaccard.numba_func([(nb.bool_[:], nb.bool_[:], nb.float32[:])],
                    '(n),(n)->()')
def nb_jaccard(query, ref, out):

	union = 0
	intersection = 0

	for i in range(query.shape[0]):
		union += query[i] | ref[i]
		intersection += query[i] & ref[i]

	out[0] = intersection
	out[0] /= union

@jaccard.coords_func(coords_nb_signature, '(n),(m)->()')
def coords_jaccard(query, ref, out):

	n = query.shape[0]
	m = ref.shape[0]

	intersection = 0
	union = 0

	i = 0
	j = 0
	while i < n and j < m:
		q = query[i]
		r = ref[j]

		union += 1
		if q == r:
			intersection += 1

		if q <= r:
			i += 1

		if r <= q:
			j += 1

	union += n - i
	union += m - j

	out[0] = intersection
	out[0] /= union


@metric('Asymmetrical Jaccard', is_distance=False)
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

@asym_jacc.numba_func([(nb.bool_[:], nb.bool_[:], nb.float32[:])],
                      '(n),(n)->()')
def nb_asym_jacc(query, ref, out):

	ref_weight = 0
	intersection = 0

	for i in range(query.shape[0]):
		ref_weight += ref[i]
		intersection += query[i] & ref[i]

	out[0] = intersection
	out[0] /= ref_weight


@asym_jacc.coords_func(coords_nb_signature, '(n),(m)->()')
def coords_asym_jacc(query, ref, out):

	n = query.shape[0]
	m = ref.shape[0]

	intersection = 0

	i = 0
	j = 0
	while i < n and j < m:
		q = query[i]
		r = ref[j]

		if q == r:
			intersection += 1

		if q <= r:
			i += 1

		if r <= q:
			j += 1

	out[0] = intersection
	out[0] /= m


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
			scores = metric.nb(cls.query.np_array[:, None, :], ref_vecs[None, ...])
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
