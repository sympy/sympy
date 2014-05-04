""" Caching facility for SymPy """
from __future__ import print_function, division
from collections import deque
from weakref import WeakValueDictionary,WeakKeyDictionary
# TODO: refactor CACHE & friends into class?

# global cache registry:
CACHE = []  # [] of
            # (item, {} or tuple of {})

from sympy.core.decorators import wraps

def print_cache():
    """print cache content"""

    for item, cache in CACHE:
        item = str(item)
        head = '='*len(item)

        print(head)
        print(item)
        print(head)

        if not isinstance(cache, tuple):
            cache = (cache,)
            shown = False
        else:
            shown = True

        for i, kv in enumerate(cache):
            if shown:
                print('\n*** %i ***\n' % i)

            if hasattr(kv,'items'):
                for k, v in list(kv.items()):
                    print('  %s :\t%s' % (k, v))
            else:
                for k in kv:
                    print('  %s  ' % k )


def clear_cache():
    """clear cache content"""
    for item, cache in CACHE:
        if not isinstance(cache, tuple):
            cache = (cache,)

        for kv in cache:
            kv.clear()

########################################


def __cacheit_nocache(func):
    return func


# from sympy.assumptions.assume import global_assumptions  # circular import
from sympy.core.evaluate import global_evaluate
_globals = (global_evaluate,)


def __cacheit(func):
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
    # Many of the functions which we want to cache have arguments and return
    # values which are not weak-referenceable.  To get around this, the 
    # arguments of the function are used as the key in weakValueDict while 
    # the return values are stored as values in a weakKeyDict.  The dicts 
    # are linked via a dummy object which is weak-referenceable and hashable.
    # Messy...
    func._cache_it_arg2obj = func_cache_it_arg2obj = WeakValueDictionary()
    func._cache_it_obj2ret = func_cache_it_obj2ret = WeakKeyDictionary()
    # weak references will be deleted when values are removed from deque
    func._cache_it_cache = func_cache_it_cache = deque(maxlen=2000)

    CACHE.append((func, (func_cache_it_arg2obj,
                         func_cache_it_obj2ret,
                         func_cache_it_cache)))

    # dummy class
    class Object():
        pass

    @wraps(func)
    def wrapper(*args, **kw_args):
        """
        Assemble the args and kw_args to compute the hash.
        """
        k = [(x, type(x)) for x in args]
        if kw_args:
            keys = sorted(kw_args)
            k.extend([(x, kw_args[x], type(kw_args[x])) for x in keys])
        if _globals:
            k.extend([tuple(g) for g in _globals])
        k = tuple(k)
        try:
            hash(k)
        except TypeError: # k is unhashable
            return func(*args,**kw_args)
        
        # check if we are cached
        try:
            obj = func_cache_it_arg2obj[k]
            ret = func_cache_it_obj2ret[obj]
        except KeyError: # not cached
            ret = func(*args, **kw_args)
            obj = Object()
            func_cache_it_arg2obj[k] = obj
            func_cache_it_obj2ret[obj] = ret

        # move 'obj' to front of deque
        func_cache_it_cache.append(obj)
        return ret

    return wrapper


def __cacheit_debug(func):
    """cacheit + code to check cache consistency"""
    cfunc = __cacheit(func)

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


def _getenv(key, default=None):
    from os import getenv
    return getenv(key, default)

# SYMPY_USE_CACHE=yes/no/debug
USE_CACHE = _getenv('SYMPY_USE_CACHE', 'yes').lower()

if USE_CACHE == 'no':
    cacheit = __cacheit_nocache
elif USE_CACHE == 'yes':
    cacheit = __cacheit
elif USE_CACHE == 'debug':
    cacheit = __cacheit_debug   # a lot slower
else:
    raise RuntimeError(
        'unrecognized value for SYMPY_USE_CACHE: %s' % USE_CACHE)
