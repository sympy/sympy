""" Generic Unification algorithm for expression trees with lists of children

The implementation is a direct translation of

Artificial Intelligence: A Modern Approach
by
Stuart Russel and Peter Norvig

Second edition, section 9.2, page 276

It is modified in the following ways:

1.  We allow associative and commutative Compound expressions. This results in
    combinatorial blowup.
2.  We explore the tree lazily
3.  We provide generic interfaces to symbolic algebra libraries in Python.

A more traditional version can be found here
http://aima.cs.berkeley.edu/python/logic.html
"""

from collections import namedtuple
from itertools import combinations
Compound = namedtuple('Compound', 'op args')
Variable = namedtuple('Variable', 'arg')
CondVariable = namedtuple('Variable', 'arg valid')
from sys import stdout

def _unify(x, y, s, **fns):
    if x == y:
        yield s
    elif isinstance(x, (Variable, CondVariable)):
        for x in _unify_var(x, y, s, **fns):
            yield x
    elif isinstance(y, (Variable, CondVariable)):
        for x in _unify_var(y, x, s, **fns):
            yield x
    elif isinstance(x, Compound) and isinstance(y, Compound):
        is_commutative = fns.get('is_commutative', lambda x: False)
        is_associative = fns.get('is_associative', lambda x: False)
        for sop in _unify(x.op, y.op, s, **fns):
            if is_associative(x) and is_associative(y):
                a, b = (x, y) if len(x.args) < len(y.args) else (y, x)
                if is_commutative(x) and is_commutative(y):
                    combinations = allcombinations(a.args, b.args, False)
                else:
                    combinations = allcombinations(a.args, b.args, True)
                for aaargs, bbargs in combinations:
                    aa = [unpack(Compound(a.op, arg)) for arg in aaargs]
                    bb = [unpack(Compound(b.op, arg)) for arg in bbargs]
                    for x in _unify(aa, bb, sop, **fns):
                        yield x
            elif len(x.args) == len(y.args):
                for x in _unify(x.args, y.args, sop, **fns):
                    yield x

    elif is_args(x) and is_args(y) and len(x) == len(y):
        if len(x) == 0:
            yield s
        else:
            for shead in _unify(x[0], y[0], s, **fns):
                for x in _unify(x[1:], y[1:], shead, **fns):
                    yield x

def _unify_var(var, x, s, **fns):
    if var in s:
        for x in _unify(s[var], x, s, **fns):
            yield x
    elif occur_check(var, x):
        pass
    elif isinstance(var, CondVariable) and var.valid(x):
        yield assoc(s, var, x)
    elif isinstance(var, Variable):
        yield assoc(s, var, x)

def occur_check(var, x):
    """ var occurs in subtree owned by x? """
    if var == x:
        return True
    elif isinstance(x, Compound):
        return occur_check(var, x.args)
    elif is_args(x):
        if any(occur_check(var, xi) for xi in x): return True
    return False

def assoc(d, key, val):
    """ Return copy of d with key associated to val """
    d = d.copy()
    d[key] = val
    return d

def is_args(x):
    """ Is x a traditional iterable? """
    return type(x) in (tuple, list, set)

def unpack(x):
    if isinstance(x, Compound) and len(x.args) == 1:
        return x.args[0]
    else:
        return x

def allcombinations(A, B, ordered):
    """
    Restructure A and B to have the same number of elements

    Assuming either
    associativity - ordered == True
    commutativity - ordered == None

    A and B can be rearranged so that the larger of the two lists is
    reorganized into smaller sublists.

    >>> for x in allcombinations((1, 2, 3), (5, 6), True): print x
    (((1,), (2, 3)), (5, 6))
    (((1, 2), (3,)), (5, 6))

    >>> for x in allcombinations((1, 2, 3), (5, 6), None): print x
    (((1,), (2, 3)), (5, 6))
    (((2,), (3, 1)), (5, 6))
    (((3,), (1, 2)), (5, 6))
    (((1, 2), (3,)), (5, 6))
    (((2, 3), (1,)), (5, 6))
    (((3, 1), (2,)), (5, 6))
    """
    sm, bg = (A, B) if len(A) < len(B) else (B, A)
    for part in kbin(range(len(bg)), len(sm), ordered=ordered):
        if bg == B:
            yield tuple((a,) for a in A), partition(B, part)
        else:
            yield partition(A, part), tuple((b,) for b in B)

def partition(it, part):
    """ Partition a tuple/list into pieces defined by indices

    >>> partition((10, 20, 30, 40), [[0, 1, 2], [3]])
    ((10, 20, 30), (40,))
    """

    t = type(it)
    return t([index(it, ind) for ind in part])


def index(it, ind):
    """ Fancy indexing into an indexable iterable (tuple, list)

    >>> index([10, 20, 30], (1, 2, 0))
    [20, 30, 10]
    """
    return type(it)([it[i] for i in ind])

def kbin(l, k, ordered=True):
    """
    Return sequence ``l`` partitioned into ``k`` bins.
    If ordered is True then the order of the items in the
    flattened partition will be the same as the order of the
    items in ``l``; if False, all permutations of the items will
    be given; if None, only unique permutations for a given
    partition will be given.

    Examples
    ========

    >>> from sympy.utilities.iterables import kbin
    >>> for p in kbin(range(3), 2):
    ...     print p
    ...
    [[0], [1, 2]]
    [[0, 1], [2]]
    >>> for p in kbin(range(3), 2, ordered=False):
    ...     print p
    ...
    [(0,), (1, 2)]
    [(0,), (2, 1)]
    [(1,), (0, 2)]
    [(1,), (2, 0)]
    [(2,), (0, 1)]
    [(2,), (1, 0)]
    [(0, 1), (2,)]
    [(0, 2), (1,)]
    [(1, 0), (2,)]
    [(1, 2), (0,)]
    [(2, 0), (1,)]
    [(2, 1), (0,)]
    >>> for p in kbin(range(3), 2, ordered=None):
    ...     print p
    ...
    [[0], [1, 2]]
    [[1], [2, 0]]
    [[2], [0, 1]]
    [[0, 1], [2]]
    [[1, 2], [0]]
    [[2, 0], [1]]

    """
    from sympy.utilities.iterables import partitions
    from itertools import permutations
    def rotations(seq):
        for i in range(len(seq)):
            yield seq
            seq.append(seq.pop(0))
    if ordered is None:
        func = rotations
    else:
        func = permutations
    for p in partitions(len(l), k):
        if sum(p.values()) != k:
            continue
        for pe in permutations(p.keys()):
            rv = []
            i = 0
            for part in pe:
                for do in range(p[part]):
                    j = i + part
                    rv.append(l[i: j])
                    i = j
            if ordered:
                yield rv
            else:
                template = [len(i) for i in rv]
                for pp in func(l):
                    rvp = []
                    ii = 0
                    for t in template:
                        jj = ii + t
                        rvp.append(pp[ii: jj])
                        ii = jj
                    yield rvp
