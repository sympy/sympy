
def recurrence_memo(initial):
    """
    Memo decorator for sequences defined by recurrence
    
    See usage examples e.g. in the specfun/combinatorial module
    """
    cache = initial
    def decorator(f):
        def g(n):
            L = len(cache)
            if n <= L - 1:
                return cache[n]
            for i in xrange(L, n+1):
                cache.append(f(i, cache))
            return cache[-1]
        return g
    return decorator
