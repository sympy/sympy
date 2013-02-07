""" Generic SymPy-Independent Strategies """
from functools import partial

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

def debug(rule, file=None):
    """ Print out before and after expressions each time rule is used """
    if file is None:
        from sys import stdout
        file = stdout
    def debug_rl(expr):
        result = rule(expr)
        if result != expr:
            file.write("Rule: %s\n"%rule.func_name)
            file.write("In:   %s\nOut:  %s\n\n"%(expr, result))
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

def tryit(rule):
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

def switch(key, ruledict):
    """ Select a rule based on the result of key called on the function """
    def switch_rl(expr):
        rl = ruledict.get(key(expr), identity)
        return rl(expr)
    return switch_rl

identity = lambda x: x

def minimize(*rules, **kwargs):
    """ Select result of rules that minimizes objective

    >>> from sympy.rules import minimize
    >>> inc = lambda x: x + 1
    >>> dec = lambda x: x - 1
    >>> rl = minimize(inc, dec)
    >>> rl(4)
    3

    >>> rl = minimize(inc, dec, objective=lambda x: -x)  # maximize
    >>> rl(4)
    5
    """

    objective = kwargs.get('objective', identity)
    def minrule(expr):
        return min([rule(expr) for rule in rules], key=objective)
    return minrule

join = {list: chain, tuple: minimize}
def treesearch(tree, join=join):
    """ Transform call-tree into function

    Each node in the tree can be a

    function - returned
    list     - all trees within list are chained together
    tuple    - optimal value from within tuple is returned

    Textual example
    ---------------

    Text:
    Try either ``expand`` then ``simplify`` or try ``factor`` then ``foosimp``.

    Code:
    tree = ([expand, simplify], [factor, foosimp])

    Example

    >>> from sympy.rules import treesearch
    >>> inc    = lambda x: x + 1
    >>> dec    = lambda x: x - 1
    >>> double = lambda x: 2*x

    >>> tree = (inc, [dec, double]) # either inc or dec-then-double
    >>> fn = treesearch(tree)
    >>> fn(4)  # lowest value comes from the inc
    5
    >>> fn(1)  # lowest value comes from dec then double
    0

    By default this function uses the strategies ``chain`` and ``minimize``.
    These can be changed with the join keyword

    >>> from functools import partial
    >>> from sympy.rules import chain, minimize
    >>> maximize = partial(minimize, objective=lambda x: -x)
    >>> d = {list: chain, tuple: maximize}
    >>> fn = treesearch(tree, join=d)
    >>> fn(4)  # highest value comes from the dec then double
    6
    >>> fn(1)  # highest value comes from the inc
    2
    """

    if callable(tree):
        return tree
    for typ, strat in join.iteritems():
        if type(tree) == typ:
            return strat(*map(partial(treesearch, join=join), tree))
    raise TypeError("Type of tree %s not found in known types:\n\t%s"%(
                str(type(tree)), "{" + ", ".join(map(str, join.keys()))))
