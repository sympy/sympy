"""Utility functions for geometrical entities.

Contains
========
intersection
convex_hull
closest_points
farthest_points
are_coplanar
are_similar

"""
from __future__ import division, print_function

from sympy import Function, Symbol, solve, Dummy, Expr, Mul, S
from sympy.core.evalf import pure_complex
from sympy.solvers.solvers import unrad
from sympy.polys.polytools import real_roots
from sympy.solvers.solveset import solveset_real
from sympy.core.compatibility import (
    is_sequence, range, string_types, ordered)


def _ordered_points(p):
    """Return the tuple of points sorted numerically according to args"""
    return tuple(sorted(p, key=lambda x: x.args))


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

    See Also
    ========

    sympy.core.symbol.Symbol

    """
    if isinstance(s, string_types):
        if matching_symbol and matching_symbol.name == s:
            return matching_symbol
        return Symbol(s, real=True)
    elif isinstance(s, Symbol):
        return s
    else:
        raise ValueError('symbol must be string for symbol name or Symbol')


def _uniquely_named_symbol(xname, *exprs):
    """Return a symbol which, when printed, will have a name unique
    from any other already in the expressions given. The name is made
    unique by prepending underscores.
    """
    prefix = '%s'
    x = prefix % xname
    syms = set().union(*[e.free_symbols for e in exprs])
    while any(x == str(s) for s in syms):
        prefix = '_' + prefix
        x = prefix % xname
    return _symbol(x)


def are_coplanar(*e):
    """ Returns True if the given entities are coplanar otherwise False

    Parameters
    ==========

    e: entities to be checked for being coplanar

    Returns
    =======

    Boolean

    Examples
    ========

    >>> from sympy import Point3D, Line3D
    >>> from sympy.geometry.util import are_coplanar
    >>> a = Line3D(Point3D(5, 0, 0), Point3D(1, -1, 1))
    >>> b = Line3D(Point3D(0, -2, 0), Point3D(3, 1, 1))
    >>> c = Line3D(Point3D(0, -1, 0), Point3D(5, -1, 9))
    >>> are_coplanar(a, b, c)
    False

    """
    from .line import LinearEntity3D
    from .point import Point3D
    from .plane import Plane
    from .entity import GeometryEntity
    # XXX update tests for coverage

    e = set(e)
    # first work with a Plane if present
    for i in list(e):
        if isinstance(i, Plane):
            e.remove(i)
            return all(p.is_coplanar(i) for p in e)

    if all(isinstance(i, Point3D) for i in e):
        if len(e) < 3:
            return False

        # remove pts that are collinear with 2 pts
        a, b = e.pop(), e.pop()
        for i in list(e):
            if Point3D.are_collinear(a, b, i):
                e.remove(i)

        if not e:
            return False
        else:
            # define a plane
            p = Plane(a, b, e.pop())
            for i in e:
                if i not in p:
                    return False
            return True
    else:
        pt3d = []
        for i in e:
            if isinstance(i, Point3D):
                pt3d.append(i)
            elif isinstance(i, LinearEntity3D):
                pt3d.extend(i.args)
            elif isinstance(i, GeometryEntity):  # XXX we should have a GeometryEntity3D class so we can tell the difference between 2D and 3D -- here we just want to deal with 2D objects; if new 3D objects are encountered that we didn't hanlde above, an error should be raised
                # all 2D objects have some Point that defines them; so convert those points to 3D pts by making z=0
                for p in i.args:
                    if isinstance(p, Point):
                        pt3d.append(Point3D(*(p.args + (0,))))
        return are_coplanar(*pt3d)


def are_similar(e1, e2):
    """Are two geometrical entities similar.

    Can one geometrical entity be uniformly scaled to the other?

    Parameters
    ==========

    e1 : GeometryEntity
    e2 : GeometryEntity

    Returns
    =======

    are_similar : boolean

    Raises
    ======

    GeometryError
        When `e1` and `e2` cannot be compared.

    Notes
    =====

    If the two objects are equal then they are similar.

    See Also
    ========

    sympy.geometry.entity.GeometryEntity.is_similar

    Examples
    ========

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
    from .exceptions import GeometryError

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
            raise GeometryError(
                "Cannot test similarity between %s and %s" % (n1, n2))


def centroid(*args):
    """Find the centroid (center of mass) of the collection containing only Points,
    Segments or Polygons. The centroid is the weighted average of the individual centroid
    where the weights are the lengths (of segments) or areas (of polygons).
    Overlapping regions will add to the weight of that region.

    If there are no objects (or a mixture of objects) then None is returned.

    See Also
    ========

    sympy.geometry.point.Point, sympy.geometry.line.Segment,
    sympy.geometry.polygon.Polygon

    Examples
    ========

    >>> from sympy import Point, Segment, Polygon
    >>> from sympy.geometry.util import centroid
    >>> p = Polygon((0, 0), (10, 0), (10, 10))
    >>> q = p.translate(0, 20)
    >>> p.centroid, q.centroid
    (Point2D(20/3, 10/3), Point2D(20/3, 70/3))
    >>> centroid(p, q)
    Point2D(20/3, 40/3)
    >>> p, q = Segment((0, 0), (2, 0)), Segment((0, 0), (2, 2))
    >>> centroid(p, q)
    Point2D(1, -sqrt(2) + 2)
    >>> centroid(Point(0, 0), Point(2, 0))
    Point2D(1, 0)

    Stacking 3 polygons on top of each other effectively triples the
    weight of that polygon:

        >>> p = Polygon((0, 0), (1, 0), (1, 1), (0, 1))
        >>> q = Polygon((1, 0), (3, 0), (3, 1), (1, 1))
        >>> centroid(p, q)
        Point2D(3/2, 1/2)
        >>> centroid(p, p, p, q) # centroid x-coord shifts left
        Point2D(11/10, 1/2)

    Stacking the squares vertically above and below p has the same
    effect:

        >>> centroid(p, p.translate(0, 1), p.translate(0, -1), q)
        Point2D(11/10, 1/2)

    """

    from sympy.geometry import Polygon, Segment, Point
    if args:
        if all(isinstance(g, Point) for g in args):
            c = Point(0, 0)
            for g in args:
                c += g
            den = len(args)
        elif all(isinstance(g, Segment) for g in args):
            c = Point(0, 0)
            L = 0
            for g in args:
                l = g.length
                c += g.midpoint*l
                L += l
            den = L
        elif all(isinstance(g, Polygon) for g in args):
            c = Point(0, 0)
            A = 0
            for g in args:
                a = g.area
                c += g.centroid*a
                A += a
            den = A
        c /= den
        return c.func(*[i.simplify() for i in c.args])


def closest_points(*args):
    """Return the subset of points from a set of points that were
    the closest to each other in the 2D plane.

    Parameters
    ==========

    args : a collection of Points on 2D plane.

    Notes
    =====

    This can only be performed on a set of points whose coordinates can
    be ordered on the number line. If there are no ties then a single
    pair of Points will be in the set.

    References
    ==========

    [1] http://www.cs.mcgill.ca/~cs251/ClosestPair/ClosestPairPS.html

    [2] Sweep line algorithm
    https://en.wikipedia.org/wiki/Sweep_line_algorithm

    Examples
    ========

    >>> from sympy.geometry import closest_points, Triangle
    >>> Triangle(sss=(3, 4, 5)).args
    (Point2D(0, 0), Point2D(3, 0), Point2D(3, 4))
    >>> closest_points(*_)
    {(Point2D(0, 0), Point2D(3, 0))}

    """
    from collections import deque
    from math import hypot, sqrt as _sqrt
    from sympy.functions.elementary.miscellaneous import sqrt
    from sympy.geometry.point import Point, Point2D

    p = [Point2D(i) for i in set(args)]
    if len(p) < 2:
        raise ValueError('At least 2 distinct points must be given.')

    try:
        p.sort(key=lambda x: x.args)
    except TypeError:
        raise ValueError("The points could not be sorted.")

    if any(not i.is_Rational for j in p for i in j.args):
        def hypot(x, y):
            arg = x*x + y*y
            if arg.is_Rational:
                return _sqrt(arg)
            return sqrt(arg)

    rv = [(0, 1)]
    best_dist = hypot(p[1].x - p[0].x, p[1].y - p[0].y)
    i = 2
    left = 0
    box = deque([0, 1])
    while i < len(p):
        while left < i and p[i][0] - p[left][0] > best_dist:
            box.popleft()
            left += 1

        for j in box:
            d = hypot(p[i].x - p[j].x, p[i].y - p[j].y)
            if d < best_dist:
                rv = [(j, i)]
            elif d == best_dist:
                rv.append((j, i))
            else:
                continue
            best_dist = d
        box.append(i)
        i += 1

    return {tuple([p[i] for i in pair]) for pair in rv}


def convex_hull(*args, **kwargs):
    """The convex hull surrounding the Points contained in the list of entities.

    Parameters
    ==========

    args : a collection of Points, Segments and/or Polygons

    Returns
    =======

    convex_hull : Polygon if ``polygon`` is True else as a tuple `(U, L)` where ``L`` and ``U`` are the lower and upper hulls, respectively.

    Notes
    =====

    This can only be performed on a set of points whose coordinates can
    be ordered on the number line.

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Graham_scan

    [2] Andrew's Monotone Chain Algorithm
    (A.M. Andrew,
    "Another Efficient Algorithm for Convex Hulls in Two Dimensions", 1979)
    http://geomalgorithms.com/a10-_hull-1.html

    See Also
    ========

    sympy.geometry.point.Point, sympy.geometry.polygon.Polygon

    Examples
    ========

    >>> from sympy.geometry import Point, convex_hull
    >>> points = [(1, 1), (1, 2), (3, 1), (-5, 2), (15, 4)]
    >>> convex_hull(*points)
    Polygon(Point2D(-5, 2), Point2D(1, 1), Point2D(3, 1), Point2D(15, 4))
    >>> convex_hull(*points, **dict(polygon=False))
    ([Point2D(-5, 2), Point2D(15, 4)],
     [Point2D(-5, 2), Point2D(1, 1), Point2D(3, 1), Point2D(15, 4)])

    """
    from .entity import GeometryEntity
    from .point import Point
    from .line import Segment
    from .polygon import Polygon

    polygon = kwargs.get('polygon', True)
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
            raise NotImplementedError(
                'Convex hull for %s not implemented.' % type(e))

    # make sure all our points are of the same dimension
    if any(len(x) != 2 for x in p):
        raise ValueError('Can only compute the convex hull in two dimensions')

    p = list(p)
    if len(p) == 1:
        return p[0] if polygon else (p[0], None)
    elif len(p) == 2:
        s = Segment(p[0], p[1])
        return s if polygon else (s, None)

    def _orientation(p, q, r):
        '''Return positive if p-q-r are clockwise, neg if ccw, zero if
        collinear.'''
        return (q.y - p.y)*(r.x - p.x) - (q.x - p.x)*(r.y - p.y)

    # scan to find upper and lower convex hulls of a set of 2d points.
    U = []
    L = []
    try:
        p.sort(key=lambda x: x.args)
    except TypeError:
        raise ValueError("The points could not be sorted.")
    for p_i in p:
        while len(U) > 1 and _orientation(U[-2], U[-1], p_i) <= 0:
            U.pop()
        while len(L) > 1 and _orientation(L[-2], L[-1], p_i) >= 0:
            L.pop()
        U.append(p_i)
        L.append(p_i)
    U.reverse()
    convexHull = tuple(L + U[1:-1])

    if len(convexHull) == 2:
        s = Segment(convexHull[0], convexHull[1])
        return s if polygon else (s, None)
    if polygon:
        return Polygon(*convexHull)
    else:
        U.reverse()
        return (U, L)

def farthest_points(*args):
    """Return the subset of points from a set of points that were
    the furthest apart from each other in the 2D plane.

    Parameters
    ==========

    args : a collection of Points on 2D plane.

    Notes
    =====

    This can only be performed on a set of points whose coordinates can
    be ordered on the number line. If there are no ties then a single
    pair of Points will be in the set.

    References
    ==========

    [1] http://code.activestate.com/recipes/117225-convex-hull-and-diameter-of-2d-point-sets/

    [2] Rotating Callipers Technique
    https://en.wikipedia.org/wiki/Rotating_calipers

    Examples
    ========

    >>> from sympy.geometry import farthest_points, Triangle
    >>> Triangle(sss=(3, 4, 5)).args
    (Point2D(0, 0), Point2D(3, 0), Point2D(3, 4))
    >>> farthest_points(*_)
    {(Point2D(0, 0), Point2D(3, 4))}

    """
    from math import hypot, sqrt as _sqrt
    from sympy.geometry.point import Point, Point2D
    def rotatingCalipers(Points):
        U, L = convex_hull(*Points, **dict(polygon=False))

        if L is None:
            if isinstance(U, Point):
                raise ValueError('At least two distinct points must be given.')
            yield U.args
        else:
            i = 0
            j = len(L) - 1
            while i < len(U) - 1 or j > 0:
                yield U[i], L[j]
                # if all the way through one side of hull, advance the other side
                if i == len(U) - 1:
                    j -= 1
                elif j == 0:
                    i += 1
                # still points left on both lists, compare slopes of next hull edges
                # being careful to avoid divide-by-zero in slope calculation
                elif (U[i+1].y - U[i].y) * (L[j].x - L[j-1].x) > \
                        (L[j].y - L[j-1].y) * (U[i+1].x - U[i].x):
                    i += 1
                else:
                    j -= 1

    p = [Point2D(i) for i in set(args)]

    if any(not i.is_Rational for j in p for i in j.args):
        def hypot(x, y):
            arg = x*x + y*y
            if arg.is_Rational:
                return _sqrt(arg)
            return sqrt(arg)

    rv = []
    diam = 0
    for pair in rotatingCalipers(args):
        h, q = _ordered_points(pair)
        d = hypot(h.x - q.x, h.y - q.y)
        if d > diam:
            rv = [(h, q)]
        elif d == diam:
            rv.append((h, q))
        else:
            continue
        diam = d

    return set(rv)


def idiff(eq, y, x, n=1):
    """Return ``dy/dx`` assuming that ``eq == 0``.

    Parameters
    ==========

    y : the dependent variable or a list of dependent variables (with y first)
    x : the variable that the derivative is being taken with respect to
    n : the order of the derivative (default is 1)

    Examples
    ========

    >>> from sympy.abc import x, y, a
    >>> from sympy.geometry.util import idiff

    >>> circ = x**2 + y**2 - 4
    >>> idiff(circ, y, x)
    -x/y
    >>> idiff(circ, y, x, 2).simplify()
    -(x**2 + y**2)/y**3

    Here, ``a`` is assumed to be independent of ``x``:

    >>> idiff(x + a + y, y, x)
    -1

    Now the x-dependence of ``a`` is made explicit by listing ``a`` after
    ``y`` in a list.

    >>> idiff(x + a + y, [y, a], x)
    -Derivative(a, x) - 1

    See Also
    ========

    sympy.core.function.Derivative: represents unevaluated derivatives
    sympy.core.function.diff: explicitly differentiates wrt symbols

    """
    if is_sequence(y):
        dep = set(y)
        y = y[0]
    elif isinstance(y, Symbol):
        dep = {y}
    else:
        raise ValueError("expecting x-dependent symbol(s) but got: %s" % y)

    f = dict([(s, Function(
        s.name)(x)) for s in eq.free_symbols if s != x and s in dep])
    dydx = Function(y.name)(x).diff(x)
    eq = eq.subs(f)
    derivs = {}
    for i in range(n):
        yp = solve(eq.diff(x), dydx)[0].subs(derivs)
        if i == n - 1:
            return yp.subs([(v, k) for k, v in f.items()])
        derivs[dydx] = yp
        eq = dydx - yp
        dydx = dydx.diff(x)


def intersection(*entities):
    """The intersection of a collection of GeometryEntity instances.

    Parameters
    ==========

    entities : sequence of GeometryEntity

    Returns
    =======

    intersection : list of GeometryEntity

    Raises
    ======

    NotImplementedError
        When unable to calculate intersection.

    Notes
    =====

    The intersection of any geometrical entity with itself should return
    a list with one item: the entity in question.
    An intersection requires two or more entities. If only a single
    entity is given then the function will return an empty list.
    It is possible for `intersection` to miss intersections that one
    knows exists because the required quantities were not fully
    simplified internally.
    Reals should be converted to Rationals, e.g. Rational(str(real_num))
    or else failures due to floating point issues may result.

    See Also
    ========

    sympy.geometry.entity.GeometryEntity.intersection

    Examples
    ========

    >>> from sympy.geometry import Point, Line, Circle, intersection
    >>> p1, p2, p3 = Point(0, 0), Point(1, 1), Point(-1, 5)
    >>> l1, l2 = Line(p1, p2), Line(p3, p2)
    >>> c = Circle(p2, 1)
    >>> intersection(l1, p2)
    [Point2D(1, 1)]
    >>> intersection(l1, l2)
    [Point2D(1, 1)]
    >>> intersection(c, p2)
    []
    >>> intersection(c, Point(1, 0))
    [Point2D(1, 0)]
    >>> intersection(c, l2)
    [Point2D(-sqrt(5)/5 + 1, 2*sqrt(5)/5 + 1),
     Point2D(sqrt(5)/5 + 1, -2*sqrt(5)/5 + 1)]

    """
    from .entity import GeometryEntity
    from .point import Point

    if len(entities) <= 1:
        return []

    # entities may be an immutable tuple
    entities = list(entities)
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


def gsolve(ge1, ge2, x=None, y=None, check=True, simplify=True):
    """Return the set of solution to 2 equations in the form of
    A*x**2 + B*y**2 + C*x*y + D*x + E*y + F == 0 where A - F are
    numerical constants.
    """
    from sympy.geometry.entity import GeometryEntity
    from sympy.core.relational import Eq
    # ------- helpers ----------
    verbose = False
    def what(ge1):
        # return given expression or GeometryEntity's equation
        # which must be in terms of x and y
        if not isinstance(ge1, GeometryEntity):
            return S(ge1)
        if isinstance(ge1, Expr):
            return ge1
        assert isinstance(ge1, GeometryEntity), "expecting Expr or GeometryEntity"
        try:
            return ge1
        except AttributeError:
            raise AttributeError('Only objects with `equation` method are supported.')
    def verify(i):
        # verify that `i` is in the form
        # a*x**2 + b*y**2 + c*x*y + d*x + e*y + f
        con = []
        for p in (x**2, y**2, x*y, x, y):
            i, d = i.as_independent(p, as_Add=True)
            if not d:
                con.append(S.Zero)
                continue
            if d == p:
                con.append(S.One)
                continue
            assert d.is_Mul, "expected Mul not %s" % d
            from sympy.utilities.iterables import sift
            sifted = sift(d.args, lambda x: x.is_number)
            assert Mul(*sifted[False]) == p, "unexpected expression: %s" % _p
            con.append(Mul(*sifted[True]))
        assert i.is_number, "expected a number for f in a*x**2+b*y**2+c*x*y+d*x+e*y+f but got %s" % i
        return con + [i]
    # ------- end of helpers ----------
    e1 = what(ge1)
    e2 = what(ge2)
    if all(isinstance(g, GeometryEntity) for g in (e1, e2)):
        x, y = Dummy('x'), Dummy('y')
        e1 = e1.equation(x, y)
        e2 = e2.equation(x, y)
    elif not any(isinstance(g, GeometryEntity) for g in (e1, e2)):
        free = e1.free_symbols | e2.free_symbols
        if len(free) == 2:
            x, y = list(ordered(free))
    elif not isinstance(e1, GeometryEntity):
        free = e1.free_symbols
        if len(free) == 2:
            x, y = list(ordered(free))
            e2 = e2.equation(x, y)
    elif not isinstance(e2, GeometryEntity):
        free = e2.free_symbols
        if len(free) == 2:
            x, y = list(ordered(free))
            e1 = e1.equation(x, y)
    assert x and y, "When providing expressions of more than 2 variables, also give x and y."
    e1 = e1.expand()
    e2 = e2.expand()
    c1 = verify(e1)
    c2 = verify(e2)
    free1 = e1.free_symbols
    free2 = e2.free_symbols
    if len(free2) == 1 and not (len(free1) == 1 and y in free1):
        e1, e2 = e2, e1
        free1, free2 = free2, free1
        a, b, c, d, e, f = c2
    else:
        a, b, c, d, e, f = c1
    if x not in free1:
        y1 = solveset_real(e1, y)
        if x not in free2:
            y2 = solveset_real(e2, y)
            return set(
                [(x, i) for i in y1.intersection(y2)])  # f(y), g(y)
        elif y not in free2:
            x2 = solveset_real(e2, x)
            return set(
                [(i, j) for i in x2 for j in y1])       # f(y), g(x)
        rv = set()
        for i in y1:
            for j in solveset_real(e2.subs(y, i), x):
                rv.add((j, i))
        return rv                                       # f(y), g(x, y)
    elif y not in free1:
        x1 = solveset_real(e1, x)
        if y not in free2:
            x2 = solveset_real(e2, x)
            return set(
                [(i, y) for i in x1.intersection(x2)])  # f(x), g(x)
        rv = set()
        for i in x1:
            for j in solveset_real(e2.subs(x, i), y):
                rv.add((i, j))
        return rv                                       # f(x), g(x, y)
    del free1, free2
    # f(x, y), g(x, y)
    x1 = solve(e1, x)
    # e2y is referenced in the unit tests
    e2y = e2.subs(x, x1[0])  # it doesn't matter which one you pick
    uu = unrad(e2y)
    exclude = set()
    if uu:
        assert not uu[1], "not expecting a change of variables"
        u = uu[0]
    else:
        u, den = e2y.as_numer_denom()
        if den.free_symbols:
            exclude = solveset_real(den, y)
    ycond = solve((c*y + d)**2 >= 4*a*(b*y**2 + e*y + f))  # or Ge((c*y + d)**2, 4*a*(b*y**2 + e*y + f))
    if ycond is S.false:
        print('ycondfalse',e1,e2)
        if verbose: print('never real')
        return set()
    if u.is_number:
        assert u.equals(0) is False, "think about this..."
        return set()  # no solution possible
    # reject solutions of y that would give an imaginary x
    y2 = []
    for i in real_roots(u):
        if i in exclude:
            if verbose: print('exclude',i)
            continue
        cond = ycond.subs(y, i)
        if cond is S.false:
            if verbose: print('imagin',i)
            continue
        elif cond is not S.true:
            print('unknonwn',e1,e2)
            if verbose: print('unevaluated cond', cond)
        y2.append(i)
    sols = set([(xi.subs(y, yi), yi) for xi in x1 for yi in y2])
    if check and uu:  # we used unrad so make sure the solutions satisfy e2
        ok = []
        for s in sols:
            reps = dict(zip((x, y), s))
            # check for False since those that *do* satisfy might
            # not give a 0 with precision
            z2 = number_isnonzero(e2.subs(reps))
            if z2 is True:
                if verbose: print('reject',2)
                continue
            if z2 is None:
                if verbose: print(s, 'could not be verified in', e2)
            ok.append(s)
        sols = set(ok)
    if simplify:
        sols = set([tuple([i.simplify() for i in s]) for s in sols])
    return sols  # let user use sqrdenest if so desired


def number_isnonzero(n):
    n = S(n)
    if n.is_number:
        ri = pure_complex(n, or_real=True)
        if ri is None:
            ri = pure_complex(n.n(2), or_real=True)
            if ri is None:
                return  # undefined function with numerical args
        r, i = ri
        R = r._prec != 1
        I = i._prec != 1
        if R and r:
            return True
        if I and i:
            return True
        if R and not r:
            if I and not i:
                return False


def number_hasimaginary(n):
    n = S(n)
    if n.is_number:
        ri = pure_complex(n, or_real=True)
        if ri is None:
            ri = pure_complex(n.n(2), or_real=True)
            if ri is None:
                return  # undefined function with numerical args
        r, i = ri
        I = i._prec != 1
        if I and i:
            return True
        if not I:
            return number_isnonzero(im(n))
        return False
