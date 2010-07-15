from sympy.core import C

#XXX: When we drop Python 2.4 support, replace minkey, iff, all, and any
# with their builtin equivalents.
def minkey(sequence, key):
    """
    Implements the min() function with the key argument.

    This is because the key argument isn't supported in Python 2.4.  The
    argument works the same as the key argument to sorted. `key` should be a
    function on the elements of `sequence` that returns objects that are
    comparable with each other, such as integers.  If more than one element
    of the sequence is the smallest with respect to key, then the first one
    will be returned.  This only supports the sequence version of min().

    == Example ==

    >>> from sympy.utilities.iterables import minkey
    >>> minkey(['ab', 'a', 'abc', 'x'], key=len)
    'a'

    """
    smallest = sequence[0]
    smallestkey = key(sequence[0])
    for i in xrange(1, len(sequence)):
        keyi = key(sequence[i])
        if keyi < smallestkey:
            smallest = sequence[i]
            smallestkey = keyi
    return smallest

def iff(condition, result1, result2):
    """
    Return result1 if condition else result2

    This is a replacement for the conditional if statement that is part of
    python 2.5+. If the condition must should not be called unless the
    condition is met, then wrap the result in a lambda; it will be called
    to return the result:

    >>> from sympy import iff
    >>> x = 0.5
    >>> iff(x == 0, x, lambda: 1/x)
    2.0
    >>> x = 0
    >>> iff(x == 0, x, lambda: 1/x)
    0
    """

    if condition:
        rv = result1
    else:
        rv = result2
    # XXX this is fragile; is there a better way to tell if it's a lambda?
    if '<lambda>' in str(rv):
        return rv()
    else:
        return rv


from sympy.core.basic import Basic


def all(iterable):
    """
    Return True if all elements are set to True. This
    function does not support predicates explicitly,
    but this behavior can be simulated easily using
    list comprehension.

    >>> from sympy import all
    >>> all( [True, True, True] )
    True
    >>> all( [True, False, True] )
    False
    >>> all( [ x % 2 == 0 for x in [2, 6, 8] ] )
    True
    >>> all( [ x % 2 == 0 for x in [2, 6, 7] ] )
    False

    NOTE: Starting from Python 2.5 this a built-in.
    """
    for item in iterable:
        if not item:
            return False
    return True

def any(iterable):
    """
    Return True if at least one element is set to True.
    This function does not support predicates explicitly,
    but this behavior can be simulated easily using
    list comprehension.

    >>> from sympy import any
    >>> any( [False, False, False] )
    False
    >>> any( [False, True, False] )
    True
    >>> any( [ x % 2 == 1 for x in [2, 6, 8] ] )
    False
    >>> any( [ x % 2 == 1 for x in [2, 6, 7] ] )
    True

    NOTE: Starting from Python 2.5 this a built-in.
    """
    for item in iterable:
        if item:
            return True
    return False

def make_list(expr, kind):
    """
    Returns a list of elements taken from specified expression
    when it is of sequence type (Add or Mul) or singleton list
    otherwise (Rational, Pow etc.).

    >>> from sympy import Symbol, make_list, Mul, Add
    >>> x, y = map(Symbol, 'xy')

    >>> make_list(x*y, Mul)
    [x, y]
    >>> make_list(x*y, Add)
    [x*y]
    >>> set(make_list(x*y + y, Add)) == set([y, x*y])
    True

    """
    if isinstance(expr, kind):
        return list(expr.args)
    else:
        return [expr]

def flatten(iterable, levels=None, cls=None):
    """
    Recursively denest iterable containers.

    >>> from sympy.utilities.iterables import flatten

    >>> flatten([1, 2, 3])
    [1, 2, 3]
    >>> flatten([1, 2, [3]])
    [1, 2, 3]
    >>> flatten([1, [2, 3], [4, 5]])
    [1, 2, 3, 4, 5]
    >>> flatten([1.0, 2, (1, None)])
    [1.0, 2, 1, None]

    If you want to denest only a specified number of levels of
    nested containers, then set ``levels`` flag to the desired
    number of levels::

    >>> ls = [[(-2, -1), (1, 2)], [(0, 0)]]

    >>> flatten(ls, levels=1)
    [(-2, -1), (1, 2), (0, 0)]

    If cls argument is specified, it will only flatten instances of that
    class, for example:

    >>> from sympy.core import Basic
    >>> class MyOp(Basic):
    ...     pass
    ...
    >>> flatten([MyOp(1, MyOp(2, 3))], cls=MyOp)
    [1, 2, 3]

    adapted from http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks
    """
    if levels is not None:
        if not levels:
            return iterable
        elif levels > 0:
            levels -= 1
        else:
            raise ValueError("expected non-negative number of levels, got %s" % levels)

    if cls is None:
        reducible = lambda x: hasattr(x, "__iter__") and not isinstance(x, basestring)
    else:
        reducible = lambda x: isinstance(x, cls)

    result = []

    for el in iterable:
        if reducible(el):
            if hasattr(el, 'args'):
                el = el.args
            result.extend(flatten(el, levels=levels, cls=cls))
        else:
            result.append(el)

    return result

def group(container, multiple=True):
    """
    Splits a container into a list of lists of equal, adjacent elements.

    >>> from sympy.utilities.iterables import group

    >>> group([1, 1, 1, 2, 2, 3])
    [[1, 1, 1], [2, 2], [3]]

    >>> group([1, 1, 1, 2, 2, 3], multiple=False)
    [(1, 3), (2, 2), (3, 1)]

    """
    if not container:
        return []

    current, groups = [container[0]], []

    for elem in container[1:]:
        if elem == current[-1]:
            current.append(elem)
        else:
            groups.append(current)
            current = [elem]

    groups.append(current)

    if multiple:
        return groups

    for i, current in enumerate(groups):
        groups[i] = (current[0], len(current))

    return groups

def split(container, key, reverse=False):
    """
    Splits a container into a list of lists with elements equivalent wrt ``key``.

    >>> from sympy.utilities.iterables import split

    >>> split([16, 8, 3, 1, 2, 5, 7], key=lambda a: a % 3)
    [[3], [16, 1, 7], [8, 2, 5]]

    """
    spliter, result = {}, []

    for elem in container:
        _key = key(elem)

        if _key in spliter:
            spliter[_key].append(elem)
        else:
            spliter[_key] = [elem]

    keys = sorted(spliter.keys(), reverse=reverse)

    for _key in keys:
        result.append(spliter[_key])

    return result

def postorder_traversal(node):
    """
    Do a postorder traversal of a tree.

    This generator recursively yields nodes that it has visited in a postorder
    fashion. That is, it descends through the tree depth-first to yield all of
    a node's children's postorder traversal before yielding the node itself.

    Parameters
    ----------
    node : sympy expression
        The expression to traverse.

    Yields
    ------
    subtree : sympy expression
        All of the subtrees in the tree.

    Examples
    --------
    >>> from sympy import symbols
    >>> from sympy.utilities.iterables import postorder_traversal
    >>> from sympy.abc import x, y, z
    >>> set(postorder_traversal((x+y)*z)) == set([z, y, x, x + y, z*(x + y)])
    True

    """
    if isinstance(node, Basic):
        for arg in node.args:
            for subtree in postorder_traversal(arg):
                yield subtree
    elif hasattr(node, "__iter__"):
        for item in node:
            for subtree in postorder_traversal(item):
                yield subtree
    yield node

def preorder_traversal(node):
    """
    Do a pre-order traversal of a tree.

    This generator recursively yields nodes that it has visited in a pre-order
    fashion. That is, it yields the current node then descends through the tree
    breadth-first to yield all of a node's children's pre-order traversal.

    Parameters
    ----------
    node : sympy expression
        The expression to traverse.

    Yields
    ------
    subtree : sympy expression
        All of the subtrees in the tree.

    Examples
    --------
    >>> from sympy import symbols
    >>> from sympy.utilities.iterables import preorder_traversal
    >>> from sympy.abc import x, y, z
    >>> set(preorder_traversal((x+y)*z)) == set([z, x + y, z*(x + y), x, y])
    True

    """
    yield node
    if isinstance(node, Basic):
        for arg in node.args:
            for subtree in preorder_traversal(arg):
                yield subtree
    elif hasattr(node, "__iter__"):
        for item in node:
            for subtree in preorder_traversal(item):
                yield subtree

def subsets(M, k):
    """
    Generates all k-subsets of n-element set.

    A k-subset of n-element set is any subset of length exactly k. The
    number of k-subsets on n elements is given by binom(n, k), whereas
    there are 2**n subsets all together.

    >>> from sympy.utilities.iterables import subsets
    >>> list(subsets([1, 2, 3], 2))
    [[1, 2], [1, 3], [2, 3]]

    """
    def recursion(result, M, k):
        if k == 0:
            yield result
        else:
            for i, item in enumerate(M[:len(M) + 1 - k]):
                for elem in recursion(result + [item], M[i + 1:], k - 1):
                    yield elem

    M = list(M)

    for i, item in enumerate(M[:len(M) + 1 - k]):
        for elem in recursion([item], M[i + 1:], k - 1):
            yield elem


def cartes(seq0, seq1, modus='pair'):
    """
    Return the Cartesian product of two sequences

    >>> from sympy.utilities.iterables import cartes
    >>> cartes([1,2], [3,4])
    [[1, 3], [1, 4], [2, 3], [2, 4]]

    """
    if  modus == 'pair':
        return [[item0, item1] for item0 in seq0 for item1 in seq1]
    elif modus == 'triple':
        return [item0 + [item1] for item0 in seq0 for item1 in seq1]

def variations(seq, n, repetition=False):
    """
    Returns all the variations of the list of size n.

    variations(seq, n, True) will return all the variations of the list of
    size n with repetitions

    variations(seq, n, False) will return all the variations of the list of
    size n without repetitions

    >>> from sympy.utilities.iterables import variations
    >>> variations([1,2,3], 2)
    [[1, 2], [1, 3], [2, 1], [2, 3], [3, 1], [3, 2]]
    >>> variations([1,2,3], 2, repetition=True)
    [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]

    """
    def setrep(seq):  # remove sets with duplicates (repetition is relevant)
        def delrep(seq):  # remove duplicates while maintaining order
            result = []
            for item in seq:
                if item not in result:
                    result.append(item)
            return result
        return [item for item in seq if item == delrep(item)]

    if n == 1:
        return [[item] for item in seq]
    result = range(len(seq))
    cartesmodus = 'pair'
    for i in range(n-1):
        result = cartes(result, range(len(seq)), cartesmodus)
        if not repetition:
            result = setrep(result)
        cartesmodus = 'triple'
    return [[seq[index] for index in indices] for indices in result]

def numbered_symbols(prefix='x', cls=None, start=0, *args, **assumptions):
    """
    Generate an infinite stream of Symbols consisting of a prefix and
    increasing subscripts.

    Parameters
    ----------
    prefix : str, optional
        The prefix to use. By default, this function will generate symbols of
        the form "x0", "x1", etc.

    cls : class, optional
        The class to use. By default, it uses Symbol, but you can also use Wild.

    start : int, optional
        The start number.  By default, it is 0.

    Yields
    ------
    sym : Symbol
        The subscripted symbols.
    """
    if cls is None:
        cls = C.Symbol

    while True:
        name = '%s%s' % (prefix, start)
        yield cls(name, *args, **assumptions)
        start += 1

def take(iter, n):
    """Return ``n`` items from ``iter`` iterator. """
    return [ iter.next() for i in xrange(n) ]

def dict_merge(*dicts):
    """Merge dictionaries into a single dictionary. """
    merged = {}

    for dict in dicts:
        merged.update(dict)

    return merged

