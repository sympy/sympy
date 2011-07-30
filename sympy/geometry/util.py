"""Utility functions for geometrical entities.

Contains
--------
intersection
convex_hull
are_similar

"""
from sympy import Symbol, Function, solve

def idiff(eq, y, x, dep=None):
    """Return dy/dx assuming that y and any other variables given in dep
    depend on x.

    >>> from sympy.abc import x, y, a
    >>> from sympy.geometry.util import idiff

    >>> idiff(x**2 + y**2 - 4, y, x)
    -x/y
    >>> idiff(x + a + y, y, x)
    -1
    >>> idiff(x + a + y, y, x, [a])
    -Derivative(a, x) - 1

    """
    if not dep:
        dep = []
    dep = set(dep)
    dep.add(y)

    f = dict([(s, Function(s.name)(x)) for s in eq.atoms(Symbol) if s != x and s in dep])
    dydx = Function(y.name)(x).diff(x)
    return solve(eq.subs(f).diff(x), dydx)[0].subs(
        [(b, a) for a, b in f.iteritems()])

def _symbol(s, matching_symbol=None):
    """Return s if s is a Symbol, else return either a new Symbol (real=True)
    with the same name s or the matching_symbol if s is a string and it matches
    the name of the matching_symbol.

    >>> from sympy import Symbol
    >>> from sympy.geometry.util import _symbol
    >>> x = Symbol('x')
    >>> _symbol('y')
    y
    >>> _.is_real
    True
    >>> _symbol(x)
    x
    >>> _.is_real is None
    True
    >>> arb = Symbol('foo')
    >>> _symbol('arb', arb) # arb's name is foo so foo will not be returned
    arb
    >>> _symbol('foo', arb) # now it will
    foo

    NB: the symbol here may not be the same as a symbol with the same
    name defined elsewhere as a result of different assumptions.

    """
    if isinstance(s, basestring):
        if matching_symbol and matching_symbol.name == s:
            return matching_symbol
        return Symbol(s, real=True)
    elif isinstance(s, Symbol):
        return s
    else:
        raise ValueError('symbol must be string for symbol name or Symbol')

def intersection(*entities):
    """The intersection of a collection of GeometryEntity instances.

    Parameters
    ----------
    entities : sequence of GeometryEntity

    Returns
    -------
    intersection : list of GeometryEntity

    Raises
    ------
    NotImplementedError
        When unable to calculate intersection.

    Notes
    -----
    The intersection of any geometrical entity with itself should return
    a list with one item: the entity in question.
    An intersection requires two or more entities. If only a single
    entity is given then the function will return an empty list.
    It is possible for `intersection` to miss intersections that one
    knows exists because the required quantities were not fully
    simplified internally.
    Reals should be converted to Rationals, e.g. Rational(str(real_num))
    or else failures due to floating point issues may result.

    Examples
    --------
    >>> from sympy.geometry import Point, Line, Circle, intersection
    >>> p1, p2, p3 = Point(0, 0), Point(1, 1), Point(-1, 5)
    >>> l1, l2 = Line(p1, p2), Line(p3, p2)
    >>> c = Circle(p2, 1)
    >>> intersection(l1, p2)
    [Point(1, 1)]
    >>> intersection(l1, l2)
    [Point(1, 1)]
    >>> intersection(c, p2)
    []
    >>> intersection(c, Point(1, 0))
    [Point(1, 0)]
    >>> intersection(c, l2)
    [Point(-5**(1/2)/5 + 1, 2*5**(1/2)/5 + 1), Point(5**(1/2)/5 + 1, -2*5**(1/2)/5 + 1)]

    """
    from entity import GeometryEntity

    if len(entities) <= 1:
        return []

    for i, e in enumerate(entities):
        if not isinstance(e, GeometryEntity):
            try:
                entities[i] = Point(e)
            except NotImplementedError:
                raise ValueError('%s is not a GeometryEntity and cannot be made into Point' % str(e))

    res = entities[0].intersection(entities[1])
    for entity in entities[2:]:
        newres = []
        for x in res:
            newres.extend(x.intersection(entity))
        res = newres
    return res


def convex_hull(*args):
    """The convex hull of a collection of 2-dimensional points.

    Parameters
    ----------
    args : a collection of Points

    Returns
    -------
    convex_hull : Polygon

    Notes
    -----
    This can only be performed on a set of non-symbolic points.

    See Also
    --------
    Point

    References
    ----------
    [1] http://en.wikipedia.org/wiki/Graham_scan

    [2] Andrew's Monotone Chain Algorithm
    ( A.M. Andrew, "Another Efficient Algorithm for Convex Hulls in Two Dimensions", 1979)
    http://softsurfer.com/Archive/algorithm_0109/algorithm_0109.htm

    Examples
    --------
    >>> from sympy.geometry import Point, convex_hull
    >>> points = [(1,1), (1,2), (3,1), (-5,2), (15,4)]
    >>> convex_hull(*points)
    Polygon(Point(-5, 2), Point(1, 1), Point(3, 1), Point(15, 4))

    """
    from entity import GeometryEntity
    from point import Point
    from line import Segment
    from polygon import Polygon

    p = set()
    for e in args:
        if not isinstance(e, GeometryEntity):
            try:
                e = Point(e)
            except NotImplementedError:
                raise ValueError('%s is not a GeometryEntity and cannot be made into Point' % str(e))
        if isinstance(e, Point):
            p.add(e)
        elif isinstance(e, Segment):
            p.update(e.points)
        elif isinstance(e, Polygon):
            p.update(e.vertices)
        else:
            raise NotImplementedError('Convex hull for %s not implemented.' % type(e))

    p = list(p)
    if len(p) == 1:
        return p[0]
    elif len(p) == 2:
        return Segment(p[0], p[1])

    def orientation(p, q, r):
        '''Return positive if p-q-r are clockwise, neg if ccw, zero if
        collinear.'''
        return (q[1] - p[1])*(r[0] - p[0]) - (q[0] - p[0])*(r[1] - p[1])

    # scan to find upper and lower convex hulls of a set of 2d points.
    U = []
    L = []
    p.sort()
    for p_i in p:
        while len(U) > 1 and orientation(U[-2], U[-1], p_i) <= 0:
            U.pop()
        while len(L) > 1 and orientation(L[-2], L[-1], p_i) >= 0:
            L.pop()
        U.append(p_i)
        L.append(p_i)
    U.reverse()
    convexHull = tuple(L + U[1:-1])

    if len(convexHull) == 2:
        return Segment(convexHull[0], convexHull[1])
    return Polygon(*convexHull)


def are_similar(e1, e2):
    """Are two geometrical entities similar.

    Can one geometrical entity be uniformly scaled to the other?

    Parameters
    ----------
    e1 : GeometryEntity
    e2 : GeometryEntity

    Returns
    -------
    are_similar : boolean

    Raises
    ------
    GeometryError
        When `e1` and `e2` cannot be compared.

    Notes
    -----
    If the two objects are equal then they are similar.

    Examples
    --------
    >>> from sympy import Point, Circle, Triangle, are_similar
    >>> c1, c2 = Circle(Point(0, 0), 4), Circle(Point(1, 4), 3)
    >>> t1 = Triangle(Point(0, 0), Point(1, 0), Point(0, 1))
    >>> t2 = Triangle(Point(0, 0), Point(2, 0), Point(0, 2))
    >>> t3 = Triangle(Point(0, 0), Point(3, 0), Point(0, 1))
    >>> are_similar(t1, t2)
    True
    >>> are_similar(t1, t3)
    False

    """
    if e1 == e2:
        return True
    try:
        return e1.is_similar(e2)
    except AttributeError:
        try:
            return e2.is_similar(e1)
        except AttributeError:
            n1 = e1.__class__.__name__
            n2 = e2.__class__.__name__
            raise GeometryError("Cannot test similarity between %s and %s" % (n1, n2))
