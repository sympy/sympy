all_caches = {}

def memoizer_immutable_args(name):
    def make_memoized(func):
        func._cache_it_cache = func_cache_it_cache = {}
        def wrapper(*args):
            try:
                return func_cache_it_cache[args]
            except KeyError:
                pass
            except TypeError, msg:
                if 'dict objects are unhashable'==str(msg):
                    return func(*args)
                raise
            func_cache_it_cache[args] = r = func(*args)
            return r
        all_caches[name] = func_cache_it_cache
        return wrapper
    return make_memoized

def memoizer_Symbol_new(func):
    func._cache_it_cache = func_cache_it_cache = {}
    def wrapper(cls, name, dummy=False, **options):
        if dummy:
            return func(cls, name, dummy=dummy, **options)
        try:
            return func_cache_it_cache[name]
        except KeyError:
            pass
        func_cache_it_cache[name] = r = func(cls, name, dummy=dummy, **options)
        return r
    all_caches['Symbol.__new__'] = func_cache_it_cache
    return wrapper

def memoizer_Interval_new(func):
    func._cache_it_cache = func_cache_it_cache = {}
    def wrapper(cls, a, b=None, **options):
        if b is None:
            # to ensure that Interval(a) is Interval(a,a)
            args = (a,a)
        else:
            args = (a,b)
        try:
            return func_cache_it_cache[args]
        except KeyError:
            pass
        func_cache_it_cache[args] = r = func(cls, a, b, **options)
        return r
    all_caches['Interval.__new__'] = func_cache_it_cache
    return wrapper

def clear_cache():
    """Clear all cached objects."""
    for cache in all_caches.values():
        cache.clear()
