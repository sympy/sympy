
def memoizer_immutable_args(func):
    func._cache_it_cache = func_cache_it_cache = {}
    def wrapper(*args):
        try:
            return func_cache_it_cache[args]
        except KeyError:
            pass
        func_cache_it_cache[args] = r = func(*args)
        return r
    return wrapper
