""" Caching facility for SymPy """
from __future__ import print_function, division
from sympy.core.decorators import wraps
import types
from collections import namedtuple
# TODO: refactor CACHE & friends into class?

class _cache(list):
    """ List of cached functions """

    def print_cache(self):
        """print cache info"""

        for item in self:
            name = item.__name__
            myfunc = item
            while hasattr(myfunc,'__wrapped__'):
                if hasattr(myfunc,'cache_info'):
                    info = myfunc.cache_info()
                    break
                else:
                    myfunc = myfunc.__wrapped__
            else:
                info = None

            print(name,info)

    def clear_cache(self):
        """clear cache content"""
        for item in CACHE:
            myfunc = item
            while hasattr(myfunc,'__wrapped__'):
                if hasattr(myfunc,'cache_clear'):
                    myfunc.cache_clear()
                    break
                else:
                    myfunc = myfunc.__wrapped__

# global cache registry:
CACHE = _cache()
# make clear and print methods available
print_cache = CACHE.print_cache
clear_cache = CACHE.clear_cache


########################################

# lru_cache class compatible with python 2 and 3.  Original function based
# implementation is part of the functools standard library for Python 3
_CacheInfo = namedtuple("CacheInfo", ["hits", "misses", "fails","maxsize", "currsize"])

class _HashedSeq(list):
    """ This class guarantees that hash() will be called no more than once
        per element.  This is important because the lru_cache() will hash
        the key multiple times on a cache miss.

    """

    __slots__ = 'hashvalue'

    def __init__(self, tup, hash=hash):
        self[:] = tup
        self.hashvalue = hash(tup)

    def __hash__(self):
        return self.hashvalue

# from sympy.assumptions.assume import global_assumptions  # circular import
from sympy.core.evaluate import global_evaluate
_globals = global_evaluate

def _make_key(args, kwds, kwd_mark = (object(),),
              sorted=sorted, tuple=tuple, type=type, len=len):
    """Make a cache key from optionally typed positional and keyword arguments

    The key is constructed in a way that is flat as possible rather than
    as a nested structure that would take more memory.

    If there is only a single argument and its data type is known to cache
    its hash value, then that argument is returned without a wrapper.  This
    saves space and improves lookup speed.

    """
    key = args
    # type info
    key += tuple(type(v) for v in args)
    if kwds:
        sorted_items = sorted(kwds.items())
        key += kwd_mark
        for item in sorted_items:
            key += item
        # type info
        key += tuple(type(v) for k, v in sorted_items)
    if _globals:
        key += tuple(g for g in _globals)
    return _HashedSeq(key)

def __cacheit(maxsize):

    class _lru_cache(object):
        """caching decorator.

           important: the result of cached function must be *immutable*


           Examples
           ========

           >>> from sympy.core.cache import cacheit
           >>> @cacheit
           ... def f(a,b):
           ...    return a+b

           >>> @cacheit
           ... def f(a,b):
           ...    return [a,b] # <-- WRONG, returns mutable object

           to force cacheit to check returned results mutability and consistency,
           set environment variable SYMPY_USE_CACHE to 'debug'
        """

        def __init__(self,func):
            wraps(func)(self)
            if not hasattr(self,'__wrapped__'):
                setattr(self,'__wrapped__',func)
            # Users should only access the lru_cache through its public API:
            #       cache_info, cache_clear, and f.__wrapped__
            # The internals of the lru_cache are encapsulated for thread safety and
            # to allow the implementation to change (including a possible C version).

            # Constants shared by all lru cache instances:
            self._make_key = _make_key         # build a key from the function arguments
            self._PREV, self._NEXT, self._KEY, self._RESULT = 0, 1, 2, 3   # names for the link fields
            self._maxsize = maxsize
            self._cache = {}
            self._hits = self._misses = self._fails = 0
            self._full = False
            self._root = []              # root of the circular doubly linked list
            self._root[:] = [self._root, self._root, None, None]     # initialize by pointing to self
            # append to global cache list
            CACHE.append(self)

        def _cache_get(self,key):
            return self._cache.get(key)

        def __call__(self,*args,**kwargs):
            # Size limited caching that tracks accesses by recency
            try:
                key = self._make_key(args, kwargs)
                if key is None:
                    raise TypeError
            except TypeError:
                # args,kwargs not hashable. mark as fail and return result
                self._fails += 1
                return self.__wrapped__(*args, **kwargs)
            PREV,NEXT,KEY,RESULT = self._PREV,self._NEXT,self._KEY,self._RESULT

            link = self._cache_get(key)
            if link is not None:
                # Move the link to the front of the circular queue
                link_prev, link_next, _key, result = link
                link_prev[NEXT] = link_next
                link_next[PREV] = link_prev
                last = self._root[PREV]
                last[NEXT] = self._root[PREV] = link
                link[PREV] = last
                link[NEXT] = self._root
                self._hits += 1
                return result
            result = self.__wrapped__(*args, **kwargs)
            if self._full:
                # Use the old self._root to store the new key and result.
                oldroot = self._root
                oldroot[KEY] = key
                oldroot[RESULT] = result
                # Empty the oldest link and make it the new self._root.
                # Keep a reference to the old key and old result to
                # prevent their ref counts from going to zero during the
                # update. That will prevent potentially arbitrary object
                # clean-up code (i.e. __del__) from running while we're
                # still adjusting the links.
                self._root = oldroot[NEXT]
                oldkey = self._root[KEY]
                oldresult = self._root[RESULT]
                self._root[KEY] = self._root[RESULT] = None
                # Now update the cache dictionary.
                del self._cache[oldkey]
                # Save the potentially reentrant cache[key] assignment
                # for last, after the self._root and links have been put in
                # a consistent state.
                self._cache[key] = oldroot
            else:
                # Put result in a new link at the front of the queue.
                last = self._root[PREV]
                link = [last, self._root, key, result]
                last[NEXT] = self._root[PREV] = self._cache[key] = link
                self._full = (len(self._cache) >= self._maxsize)
            self._misses += 1
            return result


        def cache_info(self):
            """ Report cache statistics """
            return _CacheInfo(self._hits,self._misses,self._fails,
                              self._maxsize,len(self._cache))

        def cache_clear(self):
            """ Clear the cache and cache statistics """
            self._cache.clear()
            self._root[:] = [self._root, self._root, None, None]
            self._hits = self._misses = self._fails = 0
            self._full = False

        def __get__(self,instance,cls):
            """ Ensure proper bound method object creation """
            if instance is None:
                return self
            else:
                return types.MethodType(self,instance)
    return _lru_cache

def __cacheit_debug(maxsize):

    def cachit_debug(func):
        """cacheit + code to check cache consistency"""
        cfunc = __cacheit(maxsize=maxsize)(func)

        @wraps(func)
        def wrapper(*args, **kw_args):
            # always call function itself and compare it with cached version
            r1 = func(*args, **kw_args)
            r2 = cfunc(*args, **kw_args)

            # try to see if the result is immutable
            #
            # this works because:
            #
            # hash([1,2,3])         -> raise TypeError
            # hash({'a':1, 'b':2})  -> raise TypeError
            # hash((1,[2,3]))       -> raise TypeError
            #
            # hash((1,2,3))         -> just computes the hash
            hash(r1), hash(r2)

            # also see if returned values are the same
            if r1 != r2:
                raise RuntimeError("Returned values are not the same")
            return r1
        return wrapper
    return cachit_debug


def _getenv(key, default=None):
    from os import getenv
    return getenv(key, default)

# SYMPY_USE_CACHE=yes/no/debug
USE_CACHE = _getenv('SYMPY_USE_CACHE', 'yes').lower()
# SYMPY_CACHE_SIZE=<some integer>
try:
    SYMPY_CACHE_SIZE = int(_getenv('SYMPY_CACHE_SIZE',2000))
except ValueError:
    raise RuntimeError(
        'SYMPY_CACHE_SIZE must be a valid integer. Got: %s' % SYMPY_CACHE_SIZE)

if USE_CACHE == 'no':
    cacheit = __cacheit_nocache
elif USE_CACHE == 'yes':
    cacheit = __cacheit(SYMPY_CACHE_SIZE)
elif USE_CACHE == 'debug':
    cacheit = __cacheit_debug(SYMPY_CACHE_SIZE)   # a lot slower
else:
    raise RuntimeError(
        'unrecognized value for SYMPY_USE_CACHE: %s' % USE_CACHE)
