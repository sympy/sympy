from sympy.core import Basic, C
from sympy.core.compatibility import is_sequence, iterable #logically, they belong here

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
    elif iterable(node):
        for item in node:
            for subtree in postorder_traversal(item):
                yield subtree
    yield node

class preorder_traversal(object):
    """
    Do a pre-order traversal of a tree.

    This iterator recursively yields nodes that it has visited in a pre-order
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
    def __init__(self, node):
        self._skip_flag = False
        self._pt = self._preorder_traversal(node)

    def _preorder_traversal(self, node):
        yield node
        if self._skip_flag:
            self._skip_flag = False
            return
        if isinstance(node, Basic):
            for arg in node.args:
                for subtree in self._preorder_traversal(arg):
                    yield subtree
        elif iterable(node):
            for item in node:
                for subtree in self._preorder_traversal(item):
                    yield subtree

    def skip(self):
        """
        Skip yielding current node's (last yielded node's) subtrees.

        Examples
        --------
        >>> from sympy import symbols
        >>> from sympy.utilities.iterables import preorder_traversal
        >>> from sympy.abc import x, y, z
        >>> pt = preorder_traversal((x+y*z)*z)
        >>> for i in pt:
        ...     print i
        ...     if i == x+y*z:
        ...             pt.skip()
        z*(x + y*z)
        z
        x + y*z
        """
        self._skip_flag = True

    def next(self):
        return self._pt.next()

    def __iter__(self):
        return self

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
        print "".join(map(str, args)) + END

    def _interactive_traversal(expr, stage):
        if stage > 0:
            print

        cprint("Current expression (stage ", BYELLOW, stage, END, "):")
        print BCYAN
        pprint(expr)
        print END

        if isinstance(expr, Basic):
            if expr.is_Add:
                args = expr.as_ordered_terms()
            elif expr.is_Mul:
                args = expr.as_ordered_factors()
            else:
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
            choice = raw_input("Your choice [%s,f,l,r,d,?]: " % choices)
        except EOFError:
            result = expr
            print
        else:
            if choice == '?':
                cprint(RED, "%s - select subexpression with the given index" % choices)
                cprint(RED, "f - select the first subexpression")
                cprint(RED, "l - select the last subexpression")
                cprint(RED, "r - select a random subexpression")
                cprint(RED, "d - done\n")

                result = _interactive_traversal(expr, stage)
            elif choice in ['d', '']:
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
                    cprint(BRED, "Choice must be a number in %s range\n" % choices)
                    result = _interactive_traversal(expr, stage)
                else:
                    if choice < 0 or choice >= n_args:
                        cprint(BRED, "Choice must be in %s range\n" % choices)
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
    x + 2
    y + 2
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

def sift(expr, keyfunc):
    """Sift the arguments of expr into a dictionary according to keyfunc.

    INPUT: expr may be an expression or iterable; if it is an expr then
    it is converted to a list of expr's args or [expr] if there are no args.

    OUTPUT: each element in expr is stored in a list keyed to the value
    of keyfunc for the element.

    EXAMPLES:

    >>> from sympy.utilities import sift
    >>> from sympy.abc import x, y
    >>> from sympy import sqrt, exp

    >>> sift(range(5), lambda x: x%2)
    {0: [0, 2, 4], 1: [1, 3]}

    It is possible that some keys are not present, in which case you should
    used dict's .get() method:

    >>> sift(x+y, lambda x: x.is_commutative)
    {True: [y, x]}
    >>> _.get(False, [])
    []

    Sometimes you won't know how many keys you will get:
    >>> sift(sqrt(x) + x**2 + exp(x) + (y**x)**2,
    ... lambda x: x.as_base_exp()[0])
    {E: [exp(x)], x: [x**(1/2), x**2], y: [y**(2*x)]}
    >>> _.keys()
    [E, x, y]

    """
    d = {}
    if hasattr(expr, 'args'):
        expr = expr.args or [expr]
    for e in expr:
        d.setdefault(keyfunc(e), []).append(e)
    return d

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

def topological_sort(graph, key=None):
    r"""
    Topological sort of graph's vertices.

    **Parameters**

    ``graph`` : ``tuple[list, list[tuple[T, T]]``
        A tuple consisting of a list of vertices and a list of edges of
        a graph to be sorted topologically.

    ``key`` : ``callable[T]`` (optional)
        Ordering key for vertices on the same level. By default the natural
        (e.g. lexicographic) ordering is used (in this case the base type
        must implement ordering relations).

    **Examples**

    Consider a graph::

        +---+     +---+     +---+
        | 7 |\    | 5 |     | 3 |
        +---+ \   +---+     +---+
          |   _\___/ ____   _/ |
          |  /  \___/    \ /   |
          V  V           V V   |
         +----+         +---+  |
         | 11 |         | 8 |  |
         +----+         +---+  |
          | | \____   ___/ _   |
          | \      \ /    / \  |
          V  \     V V   /  V  V
        +---+ \   +---+ |  +----+
        | 2 |  |  | 9 | |  | 10 |
        +---+  |  +---+ |  +----+
               \________/

    where vertices are integers. This graph can be encoded using
    elementary Python's data structures as follows::

        >>> V = [2, 3, 5, 7, 8, 9, 10, 11]
        >>> E = [(7, 11), (7, 8), (5, 11), (3, 8), (3, 10),
        ...      (11, 2), (11, 9), (11, 10), (8, 9)]

    To compute a topological sort for graph ``(V, E)`` issue::

        >>> from sympy.utilities.iterables import topological_sort

        >>> topological_sort((V, E))
        [3, 5, 7, 8, 11, 2, 9, 10]

    If specific tie breaking approach is needed, use ``key`` parameter::

        >>> topological_sort((V, E), key=lambda v: -v)
        [7, 5, 11, 3, 10, 8, 9, 2]

    Only acyclic graphs can be sorted. If the input graph has a cycle,
    then :py:exc:`ValueError` will be raised::

        >>> topological_sort((V, E + [(10, 7)]))
        Traceback (most recent call last):
        ...
        ValueError: cycle detected

    .. seealso:: http://en.wikipedia.org/wiki/Topological_sorting

    """
    V, E = graph

    L = []
    S = set(V)
    E = list(E)

    for v, u in E:
        S.discard(u)

    if key is None:
        key = lambda value: value

    S = sorted(S, key=key, reverse=True)

    while S:
        node = S.pop()
        L.append(node)

        for u, v in list(E):
            if u == node:
                E.remove((u, v))

                for _u, _v in E:
                    if v == _v:
                        break
                else:
                    kv = key(v)

                    for i, s in enumerate(S):
                        ks = key(s)

                        if kv > ks:
                            S.insert(i, v)
                            break
                    else:
                        S.append(v)

    if E:
        raise ValueError("cycle detected")
    else:
        return L

def rotate_left(x, y):
    """
    Left rotates a list x by the number of steps specified
    in y.

    Examples:
    >>> from sympy.utilities.iterables import rotate_left
    >>> a = [0, 1, 2]
    >>> rotate_left(a, 1)
    [1, 2, 0]
    """
    if len(x) == 0:
        return x
    y = y % len(x)
    return x[y:] + x[:y]

def rotate_right(x, y):
    """
    Left rotates a list x by the number of steps specified
    in y.

    Examples:
    >>> from sympy.utilities.iterables import rotate_right
    >>> a = [0, 1, 2]
    >>> rotate_right(a, 1)
    [2, 0, 1]
    """
    if len(x) == 0:
        return x
    y = len(x) - y % len(x)
    return x[y:] + x[:y]

def multiset_partitions(multiset, m):
    """
    This is the algorithm for generating multiset partitions
    as described by Knuth in TAOCP Vol 4.

    Given a multiset, this algorithm visits all of its
    m-partitions, that is, all partitions having exactly size m
    using auxiliary arrays as described in the book.

    Examples:
    >>> from sympy.utilities.iterables import multiset_partitions
    >>> list(multiset_partitions([1,2,3,4], 2))
    [[[1, 2, 3], [4]], [[1, 3], [2, 4]], [[1], [2, 3, 4]], [[1, 2], \
    [3, 4]], [[1, 2, 4], [3]], [[1, 4], [2, 3]], [[1, 3, 4], [2]]]
    >>> list(multiset_partitions([1,2,3,4], 1))
    [[[1, 2, 3, 4]]]
    >>> list(multiset_partitions([1,2,3,4], 4))
    [[[1], [2], [3], [4]]]
    """
    cache = {}

    def visit(n, a):
        ps = [[] for i in xrange(m)]
        for j in xrange(n):
            ps[a[j + 1]].append(multiset[j])
        canonical = tuple(tuple(j) for j in ps)
        if not canonical in cache:
            cache[canonical] = 1
            return ps

    def f(m_arr, n_arr, sigma, n, a):
        if m_arr <= 2:
            v = visit(n, a)
            if not v is None:
                yield v
        else:
            for v in f(m_arr - 1, n_arr - 1, (m_arr + sigma) % 2, n, a):
                yield v
        if n_arr == m_arr + 1:
            a[m_arr] = m_arr - 1
            v = visit(n, a)
            if not v is None:
                yield v
            while a[n_arr] > 0:
                a[n_arr] = a[n_arr] - 1
                v = visit(n, a)
                if not v is None:
                    yield v
        elif n_arr > m_arr + 1:
            if (m_arr + sigma) % 2 == 1:
                a[n_arr - 1] = m_arr - 1
            else:
                a[m_arr] = m_arr - 1
            func = [f, b][(a[n_arr] + sigma) % 2]
            for v in func(m_arr, n_arr - 1, 0, n, a):
                if v is not None:
                    yield v
            while a[n_arr] > 0:
                a[n_arr] = a[n_arr] - 1
                func = [f, b][(a[n_arr] + sigma) % 2]
                for v in func(m_arr, n_arr - 1, 0, n, a):
                    if v is not None:
                        yield v

    def b(m_arr, n_arr, sigma, n, a):
        if n_arr == m_arr + 1:
            v = visit(n, a)
            if not v is None:
                yield v
            while a[n_arr] < m_arr - 1:
                a[n_arr] = a[n_arr] + 1
                v = visit(n, a)
                if not v is None:
                    yield v
            a[m_arr] = 0
            v = visit(n, a)
            if not v is None:
                yield v
        elif n_arr > m_arr + 1:
            func = [f, b][(a[n_arr] + sigma) % 2]
            for v in func(m_arr, n_arr - 1, 0, n, a):
                if v is not None:
                    yield v
            while a[n_arr] < m_arr - 1:
                a[n_arr] = a[n_arr] + 1
                func = [f, b][(a[n_arr] + sigma) % 2]
                for v in func(m_arr, n_arr - 1, 0, n, a):
                    if v is not None:
                        yield v
            if (m_arr + sigma) % 2 == 1:
                a[n_arr - 1] = 0
            else:
                a[m_arr] = 0
        if m_arr <= 2:
            v = visit(n, a)
            if not v is None:
                yield v
        else:
            for v in b(m_arr - 1, n_arr - 1, (m_arr + sigma) % 2, n, a):
                if v is not None:
                    yield v

    n = len(multiset)
    a = [0] * (n + 1)
    for j in xrange(1, m + 1):
        a[n - m + j] = j - 1
    return f(m, n, 0, n, a)

def partitions(n, m=None, k=None):
    """Generate all partitions of integer n (>= 0).

    'm' limits the number of parts in the partition, e.g. if m=2 then
        partitions will contain no more than 2 numbers, while
    'k' limits the numbers which may appear in the partition, e.g. k=2 will
        return partitions with no element greater than 2.

    Each partition is represented as a dictionary, mapping an integer
    to the number of copies of that integer in the partition.  For example,
    the first partition of 4 returned is {4: 1}: a single 4.

    >>> from sympy.utilities.iterables import partitions

    Maximum key (number in partition) limited with k (in this case, 2):

    >>> for p in partitions(6, k=2):
    ...     print p
    {2: 3}
    {1: 2, 2: 2}
    {1: 4, 2: 1}
    {1: 6}

    Maximum number of parts in partion limited with m (in this case, 2):

    >>> for p in partitions(6, m=2):
    ...     print p
    ...
    {6: 1}
    {1: 1, 5: 1}
    {2: 1, 4: 1}
    {3: 2}

    Note that the _same_ dictionary object is returned each time.
    This is for speed:  generating each partition goes quickly,
    taking constant time independent of n.

    >>> [p for p in partitions(6, k=2)]
    [{1: 6}, {1: 6}, {1: 6}, {1: 6}]

    If you want to build a list of the returned dictionaries then
    make a copy of them:

    >>> [p.copy() for p in partitions(6, k=2)]
    [{2: 3}, {1: 2, 2: 2}, {1: 4, 2: 1}, {1: 6}]

    Reference:
        modified from Tim Peter's version to allow for k and m values:
        code.activestate.com/recipes/218332-generator-for-integer-partitions/
    """

    if n < 0:
        raise ValueError("n must be >= 0")
    m = min(m or n, n)
    if m < 1:
        raise ValueError("maximum numbers in partition, m, must be > 0")
    k = min(k or n, n)
    if k < 1:
        raise ValueError("maximum value in partition, k, must be > 0")

    if m*k < n:
        return

    q, r = divmod(n, k)
    ms = {k: q}
    keys = [k]  # ms.keys(), from largest to smallest
    if r:
        ms[r] = 1
        keys.append(r)
    room = m - q - bool(r)
    yield ms

    while keys != [1]:
        # Reuse any 1's.
        if keys[-1] == 1:
            del keys[-1]
            reuse = ms.pop(1)
            room += reuse
        else:
            reuse = 0

        while 1:
            # Let i be the smallest key larger than 1.  Reuse one
            # instance of i.
            i = keys[-1]
            newcount = ms[i] = ms[i] - 1
            reuse += i
            if newcount == 0:
                del keys[-1], ms[i]
            room += 1


            # Break the remainder into pieces of size i-1.
            i -= 1
            q, r = divmod(reuse, i)
            need = q + bool(r)
            if need > room:
                if not keys:
                    return
                continue

            ms[i] = q
            keys.append(i)
            if r:
                ms[r] = 1
                keys.append(r)
            break
        room -= need
        yield ms

def binary_partitions(n):
    """
    Generates the binary partition of n.

    A binary partition consists only of numbers that are
    powers of two. Each step reduces a 2**(k+1) to 2**k and
    2**k. Thus 16 is converted to 8 and 8.

    Reference: TAOCP 4, section 7.2.1.5, problem 64

    Examples:
    >>> from sympy.utilities.iterables import binary_partitions
    >>> for i in binary_partitions(5):
    ...     print i
    ...
    [4, 1]
    [2, 2, 1]
    [2, 1, 1, 1]
    [1, 1, 1, 1, 1]
    """
    from math import ceil, log
    pow = int(2**(ceil(log(n, 2))))
    sum = 0
    partition = []
    while pow:
        if sum + pow <= n:
            partition.append(pow)
            sum += pow
        pow >>= 1

    last_num = len(partition) - 1 - (n & 1)
    while last_num >= 0:
        yield partition
        if partition[last_num] == 2:
            partition[last_num] = 1
            partition.append(1)
            last_num -= 1
            continue
        partition.append(1)
        partition[last_num] >>= 1
        x = partition[last_num + 1] = partition[last_num]
        last_num += 1
        while x > 1:
            if x <= len(partition) - last_num - 1:
                del partition[-x + 1:]
                last_num += 1
                partition[last_num] = x
            else:
                x >>= 1
    yield [1]*n

def uniq(seq):
    '''
    Remove repeated elements from an iterable, preserving order of first
    appearance.

    Returns a sequence of the same type of the input, or a list if the input
    was not a sequence.

    Examples:
    --------
    >>> from sympy.utilities.iterables import uniq
    >>> uniq([1,4,1,5,4,2,1,2])
    [1, 4, 5, 2]
    >>> uniq((1,4,1,5,4,2,1,2))
    (1, 4, 5, 2)
    >>> uniq(x for x in (1,4,1,5,4,2,1,2))
    [1, 4, 5, 2]

    '''
    from sympy.core.function import Tuple
    seen = set()
    result = (s for s in seq if not (s in seen or seen.add(s)))
    if not hasattr(seq, '__getitem__'):
        return list(result)
    if isinstance(seq, Tuple):
        return Tuple(*tuple(result))
    return type(seq)(result)

def generate_bell(n):
    """
    Generates the bell permutations.

    In a Bell permutation, each cycle is a decreasing
    sequence of integers.

    Reference:
    [1] Generating involutions, derangements, and relatives by ECO
        Vincent Vajnovszki, DMTCS vol 1 issue 12, 2010

    Examples:
    >>> from sympy.utilities.iterables import generate_bell
    >>> list(generate_bell(3))
    [(0, 1, 2), (0, 2, 1), (1, 0, 2), (2, 0, 1), (2, 1, 0)]
    """
    P = [i for i in xrange(n)]
    T = [0]
    cache = set()
    def gen(P, T, t):
        if t == (n - 1):
            cache.add(tuple(P))
        else:
            for i in T:
                P[i], P[t+1] = P[t+1], P[i]
                if tuple(P) not in cache:
                    cache.add(tuple(P))
                    gen(P, T, t + 1)
                P[i], P[t+1] = P[t+1], P[i]
            T.append(t + 1)
            cache.add(tuple(P))
            gen(P, T, t + 1)
            T.remove(t + 1)
    gen(P, T, 0)
    return sorted(cache)

def generate_involutions(n):
    """
    Generates involutions.

    An involution is a permutation that when multiplied
    by itself equals the identity permutation. In this
    implementation the involutions are generated using
    Fixed Points.

    Alternatively, an involution can be considered as
    a permutation that does not contain any cycles with
    a length that is greater than two.

    Reference:
    http://mathworld.wolfram.com/PermutationInvolution.html

    Examples:
    >>> from sympy.utilities.iterables import \
    generate_involutions
    >>> generate_involutions(3)
    [(0, 1, 2), (0, 2, 1), (1, 0, 2), (2, 1, 0)]
    >>> len(generate_involutions(4))
    10
    """
    P = range(n) # the items of the permutation
    F = [0] # the fixed points {is this right??}
    cache = set()
    def gen(P, F, t):
        if t == n:
            cache.add(tuple(P))
        else:
            for j in xrange(len(F)):
                P[j], P[t] = P[t], P[j]
                if tuple(P) not in cache:
                    cache.add(tuple(P))
                    Fj = F.pop(j)
                    gen(P, F, t + 1)
                    F.insert(j, Fj)
                P[j], P[t] = P[t], P[j]
            t += 1
            F.append(t)
            cache.add(tuple(P))
            gen(P, F, t)
            F.pop()
    gen(P, F, 1)
    return sorted(cache)

def generate_derangements(perm):
    """
    Routine to generate derangements.

    TODO: This will be rewritten to use the
    ECO operator approach once the permutations
    branch is in master.

    Examples:
    >>> from sympy.utilities.iterables import generate_derangements
    >>> list(generate_derangements([0,1,2]))
    [[1, 2, 0], [2, 0, 1]]
    >>> list(generate_derangements([0,1,2,3]))
    [[1, 0, 3, 2], [1, 2, 3, 0], [1, 3, 0, 2], [2, 0, 3, 1], \
    [2, 3, 0, 1], [2, 3, 1, 0], [3, 0, 1, 2], [3, 2, 0, 1], \
    [3, 2, 1, 0]]
    >>> list(generate_derangements([0,1,1]))
    []
    """
    indices = range(len(perm))
    p = variations(indices, len(indices))
    for rv in \
            uniq(tuple(perm[i] for i in idx) \
                 for idx in p if all(perm[k] != \
                                     perm[idx[k]] for k in xrange(len(perm)))):
        yield list(rv)

def unrestricted_necklace(n, k):
    """
    A routine to generate unrestriced necklaces.

    Here n is the length of the necklace and k - 1
    is the maximum permissible element in the
    generated necklaces.

    Reference:
    http://mathworld.wolfram.com/Necklace.html

    Examples:
    >>> from sympy.utilities.iterables import unrestricted_necklace
    >>> [i[:] for i in unrestricted_necklace(3, 2)]
    [[0, 0, 0], [0, 1, 1]]
    >>> [i[:] for i in unrestricted_necklace(4, 4)]
    [[0, 0, 0, 0], [0, 0, 1, 0], [0, 0, 2, 0], [0, 0, 3, 0], \
    [0, 1, 1, 1], [0, 1, 2, 1], [0, 1, 3, 1], [0, 2, 2, 2], \
    [0, 2, 3, 2], [0, 3, 3, 3]]
    """
    a = [0] * n
    def gen(t, p):
        if (t > n - 1):
            if (n % p == 0):
                yield a
        else:
            a[t] = a[t - p]
            for necklace in gen(t + 1, p):
                yield necklace
            for j in xrange(a[t - p] + 1, k):
                a[t] = j
                for necklace in gen(t + 1, t):
                    yield necklace
    return gen(1, 1)

def generate_oriented_forest(n):
    """
    This algorithm generates oriented forests.

    An oriented graph is a directed graph having no symmetric pair of directed
    edges. A forest is an acyclic graph, i.e., it has no cycles. A forest can
    also be described as a disjoint union of trees, which are graphs in which
    any two vertices are connected by exactly one simple path.

    Reference:
    [1] T. Beyer and S.M. Hedetniemi: constant time generation of \
        rooted trees, SIAM J. Computing Vol. 9, No. 4, November 1980
    [2] http://stackoverflow.com/questions/1633833/
        oriented-forest-taocp-algorithm-in-python

    Examples:
    >>> from sympy.utilities.iterables import generate_oriented_forest
    >>> list(generate_oriented_forest(4))
    [[0, 1, 2, 3], [0, 1, 2, 2], [0, 1, 2, 1], [0, 1, 2, 0], \
    [0, 1, 1, 1], [0, 1, 1, 0], [0, 1, 0, 1], [0, 1, 0, 0], [0, 0, 0, 0]]
    """
    P = range(-1, n)
    while True:
        yield P[1:]
        if P[n] > 0:
            P[n] = P[P[n]]
        else:
            for p in xrange(n - 1, 0, -1):
                if P[p] != 0:
                    target = P[p] - 1
                    for q in xrange(p - 1, 0, -1):
                        if P[q] == target:
                            break
                    offset = p - q
                    for i in xrange(p, n + 1):
                        P[i] = P[i - offset]
                    break
            else:
                break
