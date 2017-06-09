from __future__ import print_function, division

from sympy.core import S, Symbol
from sympy.geometry.polygon import Polygon
from sympy.geometry.point import Point
from sympy.core.function import diff
from sympy.core.expr import Expr
from sympy.abc import x, y
from sympy.polys.polytools import gcd_list


def polytope_integrate(poly, expr, dims=None):
    """This is currently a basic prototype for integrating
    univariate/bivariate polynomials over 2-Polytopes.
    Parameters
    ==========
    poly : The input Polygon.
    expr : The input polynomial.
    dims : The tuple of symbols denoting axes.
    Example
    =======
    >>> from sympy.abc import x, y
    >>> from sympy.geometry.polygon import Polygon
    >>> from sympy.geometry.point import Point
    >>> from sympy.integrals.intpoly import polytope_integrate
    >>> poly = Polygon(Point(0,0), Point(0,1), Point(1,1), Point(1,0))
    >>> expr = x*y
    >>> polytope_integrate(poly, expr)
    1/4
    """
    if dims is None:
        if isinstance(expr, Expr):
            dims = tuple(expr.free_symbols)
        else:
            dims = (x, y)

    polys = decompose(expr)
    hp_params = hyperplane_parameters(poly)
    facets = poly.sides
    dim_length = len(dims)
    result = S.Zero

    for degree in polys:
        poly_contribute = S.Zero
        facet_count = 0
        for hp in hp_params:
            value_over_boundary = integration_reduction(facets, facet_count,
                                                        hp[0], hp[1],
                                                        polys[degree],
                                                        dims, degree)
            poly_contribute += value_over_boundary * (hp[1]/norm(hp[0]))
            facet_count += 1
        poly_contribute /= (dim_length + degree)
        result += poly_contribute

    return result


def integration_reduction(facets, index, a, b, expr, dims, degree):
    """Helper method for polytope_integrate.
    Returns the value of the input expression evaluated over the
    polytope facet referenced by a given index.
    Parameters
    ===========
    facets : List of facets of the polytope.
    index : Index referencing the facet to integrate the expression over.
    a : Hyperplane parameter denoting direction.
    b : Hyperplane parameter denoting distance.
    expr : The expression to integrate over the facet.
    dims : List of symbols denoting axes.
    degree : Degree of the homogeneous polynomial.
    """
    if expr == S.Zero:
        return expr

    value = S.Zero
    x0 = best_origin(a, b, facets[index], expr)
    gens = [x, y]
    m = len(facets)
    gens = [x, y]
    for i in range(0, len(dims)):
        df_i = diff(expr, gens[i])
        if df_i != 0:
            value += integration_reduction(facets, index, a, b,
                                           x0[i] * df_i, dims, degree - 1)

    for j in range(0, m):
        intersect = ()
        if j == (index - 1) % m or j == (index + 1) % m:
            intersect = intersection(facets[index], facets[j])
        if intersect:
            distance_origin = norm((intersect[0] - x0[0],
                                    intersect[1] - x0[1]))
            if is_vertex(intersect):
                if isinstance(expr, Expr):
                    value += distance_origin * \
                             expr.subs({gens[0]: intersect[0],
                                        gens[1]: intersect[1]})
                else:
                    value += distance_origin * expr
            else:
                value += integration_reduction(intersect, 0, a, b,
                                               distance_origin * expr, dims,
                                               degree)
    return value/(len(dims) + degree - 1)


def hyperplane_parameters(poly):
    """A helper function to return the hyperplane parameters
    of which the facets of the polygon are a part of.
    Currently works for only 2-Polytopes.
    Parameters
    ==========
    poly : The input Polygon
    """
    params = []
    vertices = list(poly.vertices)
    vertices.append(vertices[0])  # Close the polygon.
    for i in range(len(vertices) - 1):
        v1 = vertices[i]
        v2 = vertices[i + 1]

        a1 = v1[1] - v2[1]
        a2 = v2[0] - v1[0]
        b = v2[0] * v1[1] - v2[1] * v1[0]

        factor = gcd_list([a1, a2, b])

        b = b/factor
        a = (a1/factor, a2/factor)

        params.append((a, b))
    return params


def best_origin(a, b, lineseg, expr):
    """Helper method for polytope_integrate.
    Returns a point on the lineseg whose vector inner product with the
    divergence of `expr` yields an expression with the least maximum
    total power.
    Parameters
    ==========
    a : Hyperplane parameter denoting direction.
    b : Hyperplane parameter denoting distance.
    lineseg : Line segment on which to find the origin.
    expr : The expression which determines the best point.
    Algorithm(currently works only for 2D use case)
    ===============================================
    1 > Firstly, check for edge cases. Here that would refer to vertical
        or horizontal lines.

    2 > If input expression is a polynomial containing more than one generator
        then find out the total power of each of the generators.

        x**2 + 3 + x*y + x**4*y**5 ---> {x: 7, y: 6}

        If expression is a constant value then pick the first boundary point
        of the line segment.

    3 > First check if a point exists on the line segment where the value of
        the highest power generator becomes 0. If not check if the value of
        the next highest becomes 0. If none becomes 0 within line segment
        constraints then pick the first boundary point of the line segement.
        Actually, any point lying on the segment can be picked as best origin
        in the last case.
    Examples
    ========
    >>> from sympy.integrals.intpoly import *
    >>> from sympy.abc import x, y
    >>> from sympy.geometry.line import *
    >>> from sympy.geometry.point import *
    >>> l = Segment2D(Point(0, 3), Point(1, 1))
    >>> expr = x**3*y**7
    >>> best_origin((2, 1), 3, l, expr)
    (0, 3.0)
    """
    def x_axis_cut(ls):
        """Returns the point where the input line segment
        intersects the x-axis.
        Parameters:
        ===========
        ls : Line segment
        """
        p, q = ls.points
        if p.y == S.Zero:
            return tuple(p)
        elif q.y == S.Zero:
            return tuple(q)
        elif p.y/q.y < S.Zero:
            return p.y * (p.x - q.x)/(q.y - p.y) + p.x, S.Zero
        else:
            return ()

    def y_axis_cut(ls):
        """Returns the point where the input line segment
        intersects the y-axis.
        Parameters:
        ===========
        ls : Line segment
        """
        p, q = ls.points
        if p.x == S.Zero:
            return tuple(p)
        elif q.x == S.Zero:
            return tuple(q)
        elif p.x/q.x < S.Zero:
            return S.Zero, p.x * (p.y - q.y)/(q.x - p.x) + p.y
        else:
            return ()

    a1, b1 = lineseg.points[0]

    gens = [x, y]
    power_gens = {}

    for i in gens:
        power_gens[i] = S.Zero

    if len(gens) > 1:
        # Special case for vertical and horizontal lines
        if len(gens) == 2:
            if a[0] == S.Zero:
                if y_axis_cut(lineseg):
                    return S.Zero, b/a[1]
                else:
                    return a1, b1
            elif a[1] == S.Zero:
                if x_axis_cut(lineseg):
                    return b/a[0], S.Zero
                else:
                    return a1, b1

        if isinstance(expr, Expr):  # Find the sum total of power of each
            if expr.is_Add:         # generator and store in a dictionary.
                for monomial in expr.args:
                    if monomial.is_Pow:
                        if monomial.args[0] in gens:
                            power_gens[monomial.args[0]] += monomial.args[1]
                    else:
                        for univariate in monomial.args:
                            term_type = len(univariate.args)
                            if term_type == 0 and univariate in gens:
                                power_gens[univariate] += 1
                            elif term_type == 2 and univariate.args[0] in gens:
                                power_gens[univariate.args[0]] += univariate.args[1]
            elif expr.is_Mul:
                for term in expr.args:
                    term_type = len(term.args)
                    if term_type == 0 and term in gens:
                        power_gens[term] += 1
                    elif term_type == 2 and term.args[0] in gens:
                        power_gens[term.args[0]] += term.args[1]
            elif expr.is_Pow:
                power_gens[expr.args[0]] = expr.args[1]
            elif expr.is_Symbol:
                power_gens[expr] += 1
        else:  # If `expr` is a constant take first vertex of the line segment.
            return a1, b1

        #  TODO : This part is quite hacky. Should be made more robust with
        #  TODO : respect to symbol names and scalable w.r.t higher dimensions.
        power_gens = sorted(power_gens.items(), key=lambda k: str(k[0]))
        if power_gens[0][1] >= power_gens[1][1]:
            if y_axis_cut(lineseg):
                x0 = (S.Zero, b / a[1])
            elif x_axis_cut(lineseg):
                x0 = (b / a[0], S.Zero)
            else:
                x0 = (a1, b1)
        else:
            if x_axis_cut(lineseg):
                x0 = (b/a[0], S.Zero)
            elif y_axis_cut(lineseg):
                x0 = (S.Zero, b/a[1])
            else:
                x0 = (a1, b1)
    else:
        x0 = (b/a[0])
    return x0


def decompose(expr):
    """Decomposes an input polynomial into homogeneous ones of
    smaller or equal degree.
    Returns a dictionary with keys as the degree of the smaller
    constituting polynomials. Values are the constituting polynomials.
    Parameters
    ==========
    expr : Polynomial(SymPy expression)
    Examples
    ========
    >>> from sympy.abc import x, y
    >>> from sympy.integrals.intpoly import decompose
    >>> decompose(x**2 + x*y + x + y + x**3*y**2 + y**5)
    {1: x + y, 2: x**2 + x*y, 5: x**3*y**2 + y**5}
    """
    poly_dict = {}
    if isinstance(expr, Expr):
        if expr.is_Add:
            for monomial in expr.args:
                degree = S.Zero
                if monomial.is_Pow:
                    degree += monomial.args[1]
                else:
                    term_type = len(monomial.args)
                    if term_type == 0:
                        if monomial.is_Symbol:
                            degree += 1
                        else:
                            continue
                    for univariate in monomial.args:
                        term_type = len(univariate.args)
                        if term_type == 0 and univariate.is_Symbol:
                            degree += 1
                        elif term_type == 2:
                            degree += univariate.args[1]
                if degree in poly_dict:
                    poly_dict[degree] += monomial
                else:
                    poly_dict[degree] = monomial
        elif expr.is_Mul:
            degree = S.Zero
            for term in expr.args:
                term_type = len(term.args)
                if term_type == 0 and term.is_Symbol:
                    degree += 1
                elif term_type == 2:
                    degree += term.args[1]
            poly_dict[degree] = expr
        elif expr.is_Pow:
            poly_dict[expr.args[1]] = expr
        else:
            poly_dict[1] = expr
    else:
        poly_dict[0] = expr
    return poly_dict


def norm(point):
    """Returns the Euclidean norm of a point from origin.
    Parameters
    ==========
    point: This denotes a point in the dimensional space.
    Examples
    ========
    >>> from sympy.integrals.intpoly import norm
    >>> from sympy.geometry.point import Point
    >>> norm(Point(2, 7))
    sqrt(53)
    """
    h = 0
    half = S(1)/2
    if isinstance(point, tuple):
        h = (point[0] ** 2 + point[1] ** 2) ** half
    elif isinstance(point, Point):
        h = (point.x ** 2 + point.y ** 2) ** half
    elif isinstance(point, dict):
        s = 0
        for i in point.values():
            s += i ** 2
        h = s**half
    return h


def intersection(lineseg_1, lineseg_2):
    """Returns intersection between lines of which
    the input line segments are a part of.
    Note that this function is meant for use in integration_reduction
    and at that point in the calling function the lines denoted by the
    segments surely intersect within segment boundaries. Coincident lines
    are taken to be non-intersecting.

    Parameters
    ==========
    lineseg_1, lineseg_2: The input line segments

    Examples
    ========
    >>> from sympy.integrals.intpoly import intersection
    >>> from sympy.geometry.point import Point
    >>> from sympy.geometry.line import Segment2D
    >>> l1 = Segment2D(Point(1, 1), Point(3, 5))
    >>> l2 = Segment2D(Point(2, 0), Point(2, 5))
    >>> intersection(l1, l2)
    (2, 3)

    """
    x1, y1 = lineseg_1.points[0]
    x2, y2 = lineseg_1.points[1]
    x3, y3 = lineseg_2.points[0]
    x4, y4 = lineseg_2.points[1]

    denom = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4)

    if denom:
        t1 = x1*y2 - y1*x2
        t2 = x3*y4 - x4*y3
        return((t1*(x3 - x4) - t2*(x1 - x2))/denom,
               (t1*(y3 - y4) - t2*(y1 - y2))/denom)
    return ()


def is_vertex(ent):
    """If the input entity is a vertex return True
    Parameter
    =========
    ent : Denotes a geometric entity representing a point
    """
    if isinstance(ent, tuple):
        if len(ent) in [2, 3]:
            return True
    elif isinstance(ent, Point):
        return True
    return False


def plot_polytope(poly):
    """Plots the 2D polytope using the functions written in plotting
    module which in turn uses matplotlib backend.
    Parameter
    =========
    poly: Denotes a 2-Polytope
    """
    from sympy.plotting.plot import Plot, List2DSeries
    xl, yl = list(), list()

    xl = list(map(lambda vertex: vertex.x, poly.vertices))
    yl = list(map(lambda vertex: vertex.y, poly.vertices))

    xl.append(poly.vertices[0].x)  # Closing the polygon
    yl.append(poly.vertices[0].y)

    l2ds = List2DSeries(xl, yl)
    p = Plot(l2ds, axes='label_axes=True')
    p.show()


def plot_polynomial(expr):
    """Plots the polynomial using the functions written in
    plotting module which in turn uses matplotlib backend.
    Parameter
    =========
    expr: Denotes a polynomial(SymPy expression)
    """
    from sympy.plotting.plot import plot3d, plot
    gens = expr.free_symbols
    if len(gens) == 2:
        plot3d(expr)
    else:
        plot(expr)
