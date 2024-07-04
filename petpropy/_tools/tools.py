import functools
from numpy import vectorize

def vectorize_decorator(func):
    vectorized_func = vectorize(func)
    return functools.wraps(func)(vectorized_func)
