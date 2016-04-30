import numpy as np


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
