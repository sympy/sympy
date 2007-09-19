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

def clear_cache():
    """Clear all cached objects."""
    for cache in all_caches.values():
        cache.clear()
