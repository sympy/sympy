from sympy.core import Basic, C
from sympy.core.compatibility import minkey, iff, all, any #for backwards compatibility

import random

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

def interactive_traversal(expr):
    """Traverse a tree asking a user which branch to choose. """
    from sympy.printing import pprint

    RED, BRED = '\033[0;31m', '\033[1;31m'
    GREEN, BGREEN = '\033[0;32m', '\033[1;32m'
    YELLOW, BYELLOW = '\033[0;33m', '\033[1;33m'
    BLUE, BBLUE = '\033[0;34m', '\033[1;34m'
    MAGENTA, BMAGENTA = '\033[0;35m', '\033[1;35m'
    CYAN, BCYAN = '\033[0;36m', '\033[1;36m'
    END = '\033[0m'

    def cprint(*args):
        print "".join(map(str, args))

    def _interactive_traversal(expr, stage):
        if stage > 0:
            print

        cprint("Current expression (stage ", BYELLOW, stage, END, "):")
        print BCYAN
        pprint(expr)
        print END

        if isinstance(expr, Basic):
            args = expr.args
        elif hasattr(expr, "__iter__"):
            args = list(expr)
        else:
            return expr

        n_args = len(args)

        if not n_args:
            return expr

        for i, arg in enumerate(args):
            cprint(GREEN, "[", BGREEN, i, GREEN, "] ", BLUE, type(arg), END)
            pprint(arg)
            print

        if n_args == 1:
            choices = '0'
        else:
            choices = '0-%d' % (n_args-1)

        try:
            choice = raw_input("Your choice [%s,f,l,r,d]: " % choices)
        except EOFError:
            result = expr
            print
        else:
            if choice in ['d', '']:
                result = expr
            elif choice == 'f':
                result = _interactive_traversal(args[0], stage+1)
            elif choice == 'l':
                result = _interactive_traversal(args[-1], stage+1)
            elif choice == 'r':
                result = _interactive_traversal(random.choice(args), stage+1)
            else:
                try:
                    choice = int(choice)
                except ValueError:
                    result = _interactive_traversal(expr, stage)
                else:
                    if choice < 0 or choice >= n_args:
                        result = _interactive_traversal(expr, stage)
                    else:
                        result = _interactive_traversal(args[choice], stage+1)

        return result

    return _interactive_traversal(expr, 0)

def cartes(*seqs):
    """Return Cartesian product (combinations) of items from iterable
    sequences, seqs, as a generator.

    Examples::
    >>> from sympy import Add, Mul
    >>> from sympy.abc import x, y
    >>> from sympy.utilities.iterables import cartes
    >>> do=list(cartes([Mul, Add], [x, y], [2]))
    >>> for di in do:
    ...     print di[0](*di[1:])
    ...
    2*x
    2*y
    2 + x
    2 + y
    >>>

    >>> list(cartes([1, 2], [3, 4, 5]))
    [[1, 3], [1, 4], [1, 5], [2, 3], [2, 4], [2, 5]]
    """

    if not seqs:
        yield []
    else:
        for item in seqs[0]:
            for subitem in cartes(*seqs[1:]):
                yield [item] + subitem

def variations(seq, n, repetition=False):
    """Returns a generator of the variations (size n) of the list `seq` (size N).
    `repetition` controls whether items in seq can appear more than once;

    Examples:

    variations(seq, n) will return N! / (N - n)! permutations without
    repetition of seq's elements:
        >>> from sympy.utilities.iterables import variations
        >>> list(variations([1, 2], 2))
        [[1, 2], [2, 1]]

    variations(seq, n, True) will return the N**n permutations obtained
    by allowing repetition of elements:
        >>> list(variations([1, 2], 2, repetition=True))
        [[1, 1], [1, 2], [2, 1], [2, 2]]

    If you ask for more items than are in the set you get the empty set unless
    you allow repetitions:
        >>> list(variations([0, 1], 3, repetition=False))
        []
        >>> list(variations([0, 1], 3, repetition=True))[:4]
        [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1]]


    Reference:
        http://code.activestate.com/recipes/190465/
    """

    if n == 0:
        yield []
    else:
        if not repetition:
            for i in xrange(len(seq)):
                for cc in variations(seq[:i] + seq[i + 1:], n - 1, False):
                    yield [seq[i]] + cc
        else:
            for i in xrange(len(seq)):
                for cc in variations(seq, n - 1, True):
                    yield [seq[i]] + cc

def subsets(seq, k=None, repetition=False):
    """Generates all k-subsets (combinations) from an n-element set, seq.

       A k-subset of an n-element set is any subset of length exactly k. The
       number of k-subsets of an n-element set is given by binomial(n, k),
       whereas there are 2**n subsets all together. If k is None then all
       2**n subsets will be returned from shortest to longest.

       Examples:
           >>> from sympy.utilities.iterables import subsets

       subsets(seq, k) will return the n!/k!/(n - k)! k-subsets (combinations)
       without repetition, i.e. once an item has been removed, it can no
       longer be "taken":
           >>> list(subsets([1, 2], 2))
           [[1, 2]]
           >>> list(subsets([1, 2]))
           [[], [1], [2], [1, 2]]
           >>> list(subsets([1, 2, 3], 2))
           [[1, 2], [1, 3], [2, 3]]


       subsets(seq, k, repetition=True) will return the (n - 1 + k)!/k!/(n - 1)!
       combinations *with* repetition:
           >>> list(subsets([1, 2], 2, repetition=True))
           [[1, 1], [1, 2], [2, 2]]

       If you ask for more items than are in the set you get the empty set unless
       you allow repetitions:
           >>> list(subsets([0, 1], 3, repetition=False))
           []
           >>> list(subsets([0, 1], 3, repetition=True))
           [[0, 0, 0], [0, 0, 1], [0, 1, 1], [1, 1, 1]]
       """

    if type(seq) is not list:
        seq = list(seq)
    if k == 0:
        yield []
    elif k is None:
        yield []
        for k in range(1, len(seq) + 1):
            for s in subsets(seq, k, repetition=repetition):
                yield list(s)
    else:
        if not repetition:
            for i in xrange(len(seq)):
                for cc in subsets(seq[i + 1:], k - 1, False):
                    yield [seq[i]] + cc
        else:
            nmax = len(seq) - 1
            indices = [0] * k
            yield seq[:1] * k
            while 1:
                indices[-1] += 1
                if indices[-1] > nmax:
                    #find first digit that can be incremented
                    for j in range(-2, -k - 1, -1):
                        if indices[j] < nmax:
                            indices[j:] = [indices[j] + 1] * -j
                            break # increment and copy to the right
                    else:
                        break # we didn't for-break so we are done
                yield [seq[li] for li in indices]

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
        if 'dummy' in assumptions and assumptions.pop('dummy'):
            import warnings
            warnings.warn("\nuse cls=Dummy to create dummy symbols",
                          DeprecationWarning)
            cls = C.Dummy
        else:
            cls = C.Symbol

    while True:
        name = '%s%s' % (prefix, start)
        yield cls(name, *args, **assumptions)
        start += 1

def capture(func):
    """Return the printed output of func().

    `func` should be an argumentless function that produces output with
    print statements.

    >>> from sympy.utilities.iterables import capture
    >>> def foo():
    ...     print 'hello world!'
    ...
    >>> 'hello' in capture(foo) # foo, not foo()
    True
    """
    import StringIO
    import sys

    stdout = sys.stdout
    sys.stdout = file = StringIO.StringIO()
    func()
    sys.stdout = stdout
    return file.getvalue()

def take(iter, n):
    """Return ``n`` items from ``iter`` iterator. """
    return [ value for _, value in zip(xrange(n), iter) ]

def dict_merge(*dicts):
    """Merge dictionaries into a single dictionary. """
    merged = {}

    for dict in dicts:
        merged.update(dict)

    return merged

def prefixes(seq):
    """
    Generate all prefixes of a sequence.

    Example
    =======

    >>> from sympy.utilities.iterables import prefixes

    >>> list(prefixes([1,2,3,4]))
    [[1], [1, 2], [1, 2, 3], [1, 2, 3, 4]]

    """
    n = len(seq)

    for i in xrange(n):
        yield seq[:i+1]

def postfixes(seq):
    """
    Generate all postfixes of a sequence.

    Example
    =======

    >>> from sympy.utilities.iterables import postfixes

    >>> list(postfixes([1,2,3,4]))
    [[4], [3, 4], [2, 3, 4], [1, 2, 3, 4]]

    """
    n = len(seq)

    for i in xrange(n):
        yield seq[n-i-1:]
