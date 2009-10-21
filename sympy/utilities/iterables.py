from sympy.core.symbol import Symbol

def iff(condition, result1, result2):
    """return result1 if condition else result2

    This is a replacement for the conditional if statement that is part of
    python 2.5+. If the condition must should not be called unless the
    condition is met, then wrap the result in a lambda; it will be called
    to return the result:

    iff(x == 0, x, lambda: 1/x).
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

def all(iterable):
    """Return True if all elements are set to True. This
       function does not support predicates explicitely,
       but this behaviour can be simulated easily using
       list comprehension.

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
    """Return True if at least one element is set to True.
       This function does not support predicates explicitely,
       but this behaviour can be simulated easily using
       list comprehension.

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
    """Returns a list of elements taken from specified expresion
       when it is of sequence type (Add or Mul) or singleton list
       otherwise (Rational, Pow etc.).

       >>> from sympy import *
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


def flatten(iterable, cls=None):
    """Recursively denest iterable containers.

       >>> flatten([1, 2, 3])
       [1, 2, 3]
       >>> flatten([1, 2, [3]])
       [1, 2, 3]
       >>> flatten([1, [2, 3], [4, 5]])
       [1, 2, 3, 4, 5]
       >>> flatten( (1,2, (1, None)) )
       [1, 2, 1, None]

       If cls argument is specif, it will only flatten instances of that
       class, for example:

       >>> from sympy.core import Basic
       >>> class MyOp(Basic):
       ...     pass
       ...
       >>> flatten([MyOp(1, MyOp(2, 3))], cls=MyOp)
       [1, 2, 3]



    adapted from http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks
    """
    if cls is None:
        reducible = lambda x: hasattr(x, "__iter__") and not isinstance(x, basestring)
    else:
        reducible = lambda x: isinstance(x, cls)
    result = []
    for el in iterable:
        if reducible(el):
            if hasattr(el, 'args'):
                el = el.args
            result.extend(flatten(el, cls=cls))
        else:
            result.append(el)
    return result

def postorder_traversal(node):
    """ Do a postorder traversal of a tree.

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
    >>> x,y,z = symbols('xyz')
    >>> set(postorder_traversal((x+y)*z)) == set([z, y, x, x + y, z*(x + y)])
    True

    """
    for arg in node.args:
        for subtree in postorder_traversal(arg):
            yield subtree
    yield node

def preorder_traversal(node):
    """ Do a preorder traversal of a tree.

    This generator recursively yields nodes that it has visited in a preorder
    fashion. That is, it yields the current node then descends through the tree
    breadth-first to yield all of a node's children's preorder traversal.

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
    >>> x,y,z = symbols('xyz')
    >>> set(preorder_traversal((x+y)*z)) == set([z, x + y, z*(x + y), x, y])
    True

    """
    yield node
    for arg in node.args:
        for subtree in preorder_traversal(arg):
            yield subtree

def subsets(M, k):
    """Generates all k-subsets of n-element set.

       A k-subset of n-element set is any subset of length exactly k. The
       number of k-subsets on n elements is given by binom(n, k), whereas
       there are 2**n subsets all together.

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
    """Return the cartesian product of two sequences

    >>> cartes([1,2], [3,4])
    [[1, 3], [1, 4], [2, 3], [2, 4]]

    """
    if  modus == 'pair':
        return [[item0, item1] for item0 in seq0 for item1 in seq1]
    elif modus == 'triple':
        return [item0 + [item1] for item0 in seq0 for item1 in seq1]

def variations(seq, n, repetition=False):
    """Returns all the variations of the list of size n.

    variations(seq, n, True) will return all the variations of the list of
        size n with repetitions

    variations(seq, n, False) will return all the variations of the list of
        size n without repetitions

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

def numbered_symbols(prefix='x', function=Symbol, start=0, *args, **assumptions):
    """ Generate an infinite stream of Symbols consisting of a prefix and
    increasing subscripts.

    Parameters
    ----------
    prefix : str, optional
        The prefix to use. By default, this function will generate symbols of
        the form "x0", "x1", etc.

    function : function, optional
        The function to use. By default, it uses Symbol, but you can also use
        Wild.

    start : int, optional
        The start number.  By default, it is 0.

    Yields
    ------
    sym : Symbol
        The subscripted symbols.
    """

    while True:
        name = '%s%s' % (prefix, start)
        yield function(name, *args, **assumptions)
        start += 1
