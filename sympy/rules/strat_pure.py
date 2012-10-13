# Generic strategies. No dependence on SymPy

def exhaust(rule):
    """ Apply a rule repeatedly until it has no effect """
    def exhaustive_rl(expr):
        new, old = rule(expr), expr
        while(new != old):
            new, old = rule(new), new
        return new
    return exhaustive_rl

def memoize(rule):
    """ Memoized version of a rule """
    cache = {}
    def memoized_rl(expr):
        if expr in cache:
            return cache[expr]
        else:
            result = rule(expr)
            cache[expr] = result
            return result
    return memoized_rl

def condition(cond, rule):
    """ Only apply rule if condition is true """
    def conditioned_rl(expr):
        if cond(expr):
            return rule(expr)
        else:
            return      expr
    return conditioned_rl

def chain(*rules):
    """
    Compose a sequence of rules so that they apply to the expr sequentially
    """
    def chain_rl(expr):
        for rule in rules:
            expr = rule(expr)
        return expr
    return chain_rl

def debug(rule):
    """ Print out before and after expressions each time rule is used """
    def debug_rl(expr):
        result = rule(expr)
        if result != expr:
            print "Rule: %s"%rule.func_name
            print "In: %s\nOut: %s\n"%(expr, result)
        return result
    return debug_rl

def null_safe(rule):
    """ Return original expr if rule returns None """
    def null_safe_rl(expr):
        result = rule(expr)
        if result is None:
            return expr
        else:
            return result
    return null_safe_rl

def try_safe(rule):
    """ Return original expr if rule raises exception """
    def try_rl(expr):
        try:
            return rule(expr)
        except:
            return expr
    return try_rl

def do_one(*rules):
    """ Try each of the rules until one works. Then stop. """
    def do_one_rl(expr):
        for rl in rules:
            result = rl(expr)
            if result != expr:
                return result
        return expr
    return do_one_rl
