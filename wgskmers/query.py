import numpy as np
from tqdm import tqdm

from wgskmers.util import kwargs_finished


class QueryMetric(object):
	def __init__(self, name, func, is_distance):
		self.name = name
		self.func = func
		self.is_distance = is_distance

	def __call__(self, query_vec, ref_vec):
		return self.func(query_vec, ref_vec)


query_metrics = dict()

def metric(name, is_distance):
	def decorator(func):
		km = QueryMetric(name, func, is_distance)
		query_metrics[func.__name__] = km
		return km
	return decorator


@metric('Hamming distance', is_distance=True)
def hamming(query, ref, out=None):
	return (query != ref).sum(axis=-1, out=out)


@metric('Jaccard Index', is_distance=False)
def jaccard(query, ref, out=None):

	if out is None:
		out_shape = np.broadcast(query, ref).shape[:-1]
		out = np.ndarray(out_shape, dtype=np.float32)

	np.logical_and(query, ref).sum(axis=-1, out=out)
	out /= np.logical_or(query, ref).sum(axis=-1)

	return out


@metric('Asymmetrical Jaccard', is_distance=False)
def asym_jacc(query, ref, out=None):

	if out is None:
		out_shape = np.broadcast(query, ref).shape[:-1]
		out = np.ndarray(out_shape, dtype=np.float32)

	np.logical_and(query, ref).sum(axis=-1, out=out)
	out /= ref.sum(axis=-1)

	return out


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
