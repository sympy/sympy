""" Caching facility for SymPy """
from __future__ import print_function, division

from distutils.version import LooseVersion as V

class _cache(list):
    """ List of cached functions """

    def print_cache(self):
        """print cache info"""

        for item in self:
            name = item.__name__
            myfunc = item
            while hasattr(myfunc, '__wrapped__'):
                if hasattr(myfunc, 'cache_info'):
                    info = myfunc.cache_info()
                    break
                else:
                    myfunc = myfunc.__wrapped__
            else:
                info = None

            print(name, info)

    def clear_cache(self):
        """clear cache content"""
        for item in self:
            myfunc = item
            while hasattr(myfunc, '__wrapped__'):
                if hasattr(myfunc, 'cache_clear'):
                    myfunc.cache_clear()
                    break
                else:
                    myfunc = myfunc.__wrapped__

# global cache registry:
CACHE = _cache()
# make clear and print methods available
print_cache = CACHE.print_cache
clear_cache = CACHE.clear_cache

from sympy.core.compatibility import lru_cache
from functools import update_wrapper

try:
    import fastcache
    from warnings import warn
    # the version attribute __version__ is not present for all versions
    if not hasattr(fastcache, '__version__'):
        warn("fastcache version >= 0.4.0 required", UserWarning)
        raise ImportError
        # ensure minimum required version of fastcache is present
    if V(fastcache.__version__) < '0.4.0':
        warn("fastcache version >= 0.4.0 required, detected {}"\
             .format(fastcache.__version__), UserWarning)
        raise ImportError
    # Do not use fastcache if running under pypy
    import platform
    if platform.python_implementation() == 'PyPy':
        raise ImportError

except ImportError:

    def __cacheit(maxsize):
        """caching decorator.

           important: the result of cached function must be *immutable*


           Examples
           ========

           >>> from sympy.core.cache import cacheit
           >>> @cacheit
           ... def f(a, b):
           ...    return a+b

           >>> @cacheit
           ... def f(a, b):
           ...    return [a, b] # <-- WRONG, returns mutable object

           to force cacheit to check returned results mutability and consistency,
           set environment variable SYMPY_USE_CACHE to 'debug'
        """
        def func_wrapper(func):
            cfunc = lru_cache(maxsize, typed=True)(func)

            # wraps here does not propagate all the necessary info
            # for py2.7, use update_wrapper below
            def wrapper(*args, **kwargs):
                try:
                    retval = cfunc(*args, **kwargs)
                except TypeError:
                    retval = func(*args, **kwargs)
                return retval

            wrapper.__wrapped__ = cfunc.__wrapped__
            wrapper.cache_info = cfunc.cache_info
            wrapper.cache_clear = cfunc.cache_clear
            update_wrapper(wrapper, func)
            CACHE.append(wrapper)
            return wrapper

        return func_wrapper
else:

    def __cacheit(maxsize):
        """caching decorator.

           important: the result of cached function must be *immutable*


           Examples
           ========

           >>> from sympy.core.cache import cacheit
           >>> @cacheit
           ... def f(a, b):
           ...    return a+b

           >>> @cacheit
           ... def f(a, b):
           ...    return [a, b] # <-- WRONG, returns mutable object

           to force cacheit to check returned results mutability and consistency,
           set environment variable SYMPY_USE_CACHE to 'debug'
        """
        def func_wrapper(func):

            cfunc = fastcache.clru_cache(maxsize, typed=True, unhashable='ignore')(func)
            CACHE.append(cfunc)
            return cfunc

        return func_wrapper
########################################


def __cacheit_nocache(func):
    return func


def __cacheit_debug(maxsize):
    """cacheit + code to check cache consistency"""
    def func_wrapper(func):
        from .decorators import wraps

        cfunc = __cacheit(maxsize)(func)

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
    return func_wrapper


def _getenv(key, default=None):
    from os import getenv
    return getenv(key, default)

# SYMPY_USE_CACHE=yes/no/debug
USE_CACHE = _getenv('SYMPY_USE_CACHE', 'yes').lower()
# SYMPY_CACHE_SIZE=some_integer/None
# special cases :
#  SYMPY_CACHE_SIZE=0    -> No caching
#  SYMPY_CACHE_SIZE=None -> Unbounded caching
scs = _getenv('SYMPY_CACHE_SIZE', '1000')
if scs.lower() == 'none':
    SYMPY_CACHE_SIZE = None
else:
    try:
        SYMPY_CACHE_SIZE = int(scs)
    except ValueError:
        raise RuntimeError(
            'SYMPY_CACHE_SIZE must be a valid integer or None. ' + \
            'Got: %s' % SYMPY_CACHE_SIZE)

if USE_CACHE == 'no':
    cacheit = __cacheit_nocache
elif USE_CACHE == 'yes':
    cacheit = __cacheit(SYMPY_CACHE_SIZE)
elif USE_CACHE == 'debug':
    cacheit = __cacheit_debug(SYMPY_CACHE_SIZE)   # a lot slower
else:
    raise RuntimeError(
        'unrecognized value for SYMPY_USE_CACHE: %s' % USE_CACHE)


class _Desc(object):
    "Attribute descriptor for `@cached_property`."
    def __init__(self, func):
        self.func = func
    def __get__(self, instance, clazz=None):
        value = self.func(instance)
        instance.__dict__[self.func.__name__] = value
        return value

def cached_property(func):
    """A decorator that combines the effects of `@cacheit` (compute value and
    cache it so it does not need to be recomputed) and `@property` (implement
    a property through a function). However, this is far more efficient.

    It is recommended practice to use this decorator to initialize properties
    with code that needs to call into SymPy; such initialization tends to create
    import cycles.

    This decorator will not work with class or module properties: These do not
    offer the machinery needed to implement a property through a function.

    The usual `@cacheit` caveat applies: The function's result value must not
    depend on anything that can change over time.

    This can be efficient because it does not need to construct keys from
    function parameters and store them in a dict, as `@cacheit` needs to do
    so that it can associate results with previously-seen parameter
    combinations.

    The `@property` semantics is useful because parameterless, constant
    functions are indistinguishable from mere values."""
    return _Desc(func)
