""" Caching facility for SymPy """
from __future__ import print_function, division

# TODO: refactor CACHE & friends into class?

# global cache registry:
CACHE = []  # [] of
            #    (item, {} or tuple of {})

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

            for k, v in list(kv.items()):
                print('  %s :\t%s' % (k, v))


def clear_cache():
    """clear cache content"""
    for item, cache in CACHE:
        if not isinstance(cache, tuple):
            cache = (cache,)

        for kv in cache:
            kv.clear()

########################################

def user_cacheit(func):
    """
    For user level functions, possibly allowing caching.

    If the user calls a wrapped function and caching is not on, then caching is
    turned on. If caching is turned on, then at the conclusion of the function
    caching is turned off and cleared.
    """
    @wraps(func)
    def wrapper(*args, **kw_args):
        if cacheit._use_cache == 'no':
            cache_flag = True
            cacheit._use_cache = 'yes'
        else:
            cache_flag = False

        calc_me = func(*args, **kw_args)

        if cache_flag is True:
            clear_cache()
            cacheit._use_cache = 'no'
        return calc_me
    return wrapper


def cacheit(func):
    """
    Dispatches to one of 3 cache functions.

    Goes to either off, normal, or debug. Basically uses previous cache
    functions, just wrapped into one now.
    """

    func._cache_it_cache = func_cache_it_cache = {}
    CACHE.append((func, func_cache_it_cache))
    @wraps(func)
    def wrapper(*args, **kw_args):
        if cacheit._use_cache == 'yes':
            k = [(x, type(x)) for x in args]
            if kw_args:
                keys = sorted(kw_args)
                k.extend([(x, kw_args[x], type(kw_args[x])) for x in keys])
            k = tuple(k)

            try:
                return func_cache_it_cache[k]
            except (KeyError, TypeError):
                pass
            r = func(*args, **kw_args)
            try:
                func_cache_it_cache[k] = r
            except TypeError: # k is unhashable
                # Note, collections.Hashable is not smart enough to be used here.
                pass
        elif cacheit._use_cache == 'no' or cacheit._use_cache == 'debug':
            r = func(*args, **kw_args)
        else:
            raise ValueError('Bad cacheing setting')
        return r
    return wrapper
"""

    if cacheit._use_cache == 'yes':
        func._cache_it_cache = func_cache_it_cache = {}
        CACHE.append((func, func_cache_it_cache))

        @wraps(func)
        def wrapper(*args, **kw_args):
            ""
            Assemble the args and kw_args to compute the hash.
            ""
            import sys
            sys.stdout.write(str(cacheit._use_cache))
            k = [(x, type(x)) for x in args]
            if kw_args:
                keys = sorted(kw_args)
                k.extend([(x, kw_args[x], type(kw_args[x])) for x in keys])
            k = tuple(k)

            try:
                return func_cache_it_cache[k]
            except (KeyError, TypeError):
                pass
            r = func(*args, **kw_args)
            try:
                func_cache_it_cache[k] = r
            except TypeError: # k is unhashable
                # Note, collections.Hashable is not smart enough to be used here.
                pass
            return r
        return wrapper
    elif cacheit._use_cache == 'debug':
        #cacheit + code to check cache consistency
        cfunc = __cacheit(func)

        @wraps(func)
        def wrapper(*args, **kw_args):
            # always call function itself and compare it with cached version
            r1 = func(*args, **kw_args)
            r2 = cfunc(*args, **kw_args)

            hash(r1), hash(r2)
            assert r1 == r2

            return r1
        return wrapper
    elif cacheit._use_cache == 'no':
        return func
        import sys
        sys.stdout.write(str(cacheit._use_cache))
    else:
        raise RuntimeError(
            'unrecognized value for SYMPY_USE_CACHE: %s' % use_cache)
"""
cacheit._use_cache = 'no'
