"""
Module to implement integration of uni/bivariate polynomials over
2D Polytopes(Polygons).

Uses evaluation techniques as described in Chin et al(2015)[1]

References
===========
[1] : http://dilbert.engr.ucdavis.edu/~suku/quadrature/cls-integration.pdf
"""

from __future__ import print_function, division

from functools import cmp_to_key

from sympy.core import S, diff, Expr, Symbol

from sympy.geometry import Segment2D, Polygon, Point
from sympy.abc import x, y

from sympy.polys.polytools import LC, gcd_list, degree_list


def polytope_integrate(poly, expr, **kwargs):
    """Integrates homogeneous functions over polytopes.

    This function accepts the polytope in `poly` (currently only polygons are
    implemented) and the function in `expr` (currently only
    univariate/bivariate polynomials are implemented) and returns the exact
    integral of `expr` over `poly`.
    Parameters
    ==========
    poly : The input Polygon.
    expr : The input polynomial.

    Optional Parameters:
    clockwise : Binary value to sort input points of the polygon clockwise.
    max_degree : The maximum degree of any monomial of the input polynomial.
    Examples
    ========
    >>> from sympy.abc import x, y
    >>> from sympy.geometry.polygon import Polygon
    >>> from sympy.geometry.point import Point
    >>> from sympy.integrals.intpoly import polytope_integrate
    >>> polygon = Polygon(Point(0,0), Point(0,1), Point(1,1), Point(1,0))
    >>> polys = [1, x, y, x*y, x**2*y, x*y**2]
    >>> expr = x*y
    >>> polytope_integrate(polygon, expr)
    1/4
    >>> polytope_integrate(polygon, polys, max_degree=3)
    {1: 1, x: 1/2, y: 1/2, x*y: 1/4, x*y**2: 1/6, x**2*y: 1/6}
    """
    clockwise = kwargs.get('clockwise', False)
    max_degree = kwargs.get('max_degree', None)

    if clockwise is True and isinstance(poly, Polygon):
        poly = clockwise_sort(poly)

    expr = S(expr)

    if isinstance(poly, Polygon):
        # For Vertex Representation
        hp_params = hyperplane_parameters(poly)
        facets = poly.sides
    else:
        # For Hyperplane Representation
        plen = len(poly)
        intersections = [intersection(poly[(i - 1) % plen], poly[i])
                         for i in range(0, plen)]
        hp_params = poly
        lints = len(intersections)
        facets = [Segment2D(intersections[i], intersections[(i + 1) % lints])
                  for i in range(0, lints)]

    if max_degree is not None:
        result = {}
        if not isinstance(expr, list):
            raise TypeError('Input polynomials must be list of expressions')
        result_dict = main_integrate(0, facets, hp_params, max_degree)
        for polys in expr:
            if polys not in result:
                if polys is S.Zero:
                    result[S.Zero] = S.Zero
                    continue
                integral_value = S.Zero
                monoms = decompose(polys, separate=True)
                for monom in monoms:
                    if monom.is_number:
                        integral_value += result_dict[1] * monom
                    else:
                        coeff = LC(monom)
                        integral_value += result_dict[monom / coeff] * coeff
                result[polys] = integral_value
        return result

    return main_integrate(expr, facets, hp_params)


def main_integrate(expr, facets, hp_params, max_degree=None):
    """Function to translate the problem of integrating univariate/bivariate
    polynomials over a 2-Polytope to integrating over it's boundary facets.
    This is done using Generalized Stokes Theorem and Euler Theorem.

    Parameters
    ===========
    expr : The input polynomial
    facets : Facets(Line Segments) of the 2-Polytope
    hp_params : Hyperplane Parameters of the facets

    Optional Parameters:
    max_degree : The maximum degree of any monomial of the input polynomial.

    >>> from sympy.abc import x, y
    >>> from sympy.integrals.intpoly import main_integrate,\
    hyperplane_parameters
    >>> from sympy.geometry.polygon import Polygon
    >>> from sympy.geometry.point import Point
    >>> triangle = Polygon(Point(0, 3), Point(5, 3), Point(1, 1))
    >>> facets = triangle.sides
    >>> hp_params = hyperplane_parameters(triangle)
    >>> main_integrate(x**2 + y**2, facets, hp_params)
    325/6
    """
    dims = (x, y)
    dim_length = len(dims)
    result = {}
    integral_value = S.Zero

    if max_degree:
        y_degree = max_degree
        grad_terms = [[0, 0, 0, 0]] + \
            gradient_terms(max_degree)

        for facet_count, hp in enumerate(hp_params):
            a, b = hp[0], hp[1]
            x0 = facets[facet_count].points[0]

            for i, monom in enumerate(grad_terms):
                #  Every monomial is a tuple :
                #  (term, x_degree, y_degree, value over boundary)
                m, x_d, y_d, _ = monom
                value = result.get(m, None)
                if b is S.Zero:
                    value_over_boundary = S.Zero
                else:
                    value_over_boundary = \
                        integration_reduction_dynamic(facets, facet_count, a,
                                                      b, m, dims, x_d, y_d,
                                                      y_degree, x0,
                                                      grad_terms, i)
                monom[3] = value_over_boundary
                degree = x_d + y_d
                if value is not None:
                    result[m] += value_over_boundary * \
                                        (b / norm(a)) / (dim_length + degree)
                else:
                    result[m] = value_over_boundary * \
                                (b / norm(a)) / (dim_length + degree)
        return result
    else:
        polynomials = decompose(expr)
        for deg in polynomials:
            poly_contribute = S.Zero
            facet_count = 0
            for hp in hp_params:
                value_over_boundary = integration_reduction(facets,
                                                            facet_count,
                                                            hp[0], hp[1],
                                                            polynomials[deg],
                                                            dims, deg)
                poly_contribute += value_over_boundary * (hp[1] / norm(hp[0]))
                facet_count += 1
            poly_contribute /= (dim_length + deg)
            integral_value += poly_contribute
    return integral_value


def integration_reduction(facets, index, a, b, expr, dims, degree):
    """Helper method for main_integrate. Returns the value of the input
    expression evaluated over the polytope facet referenced by a given index.

    Parameters
    ===========
    facets : List of facets of the polytope.
    index : Index referencing the facet to integrate the expression over.
    a : Hyperplane parameter denoting direction.
    b : Hyperplane parameter denoting distance.
    expr : The expression to integrate over the facet.
    dims : List of symbols denoting axes.
    degree : Degree of the homogeneous polynomial.

    >>> from sympy.abc import x, y
    >>> from sympy.integrals.intpoly import integration_reduction,\
    hyperplane_parameters
    >>> from sympy.geometry.point import Point
    >>> from sympy.geometry.polygon import Polygon
    >>> triangle = Polygon(Point(0, 3), Point(5, 3), Point(1, 1))
    >>> facets = triangle.sides
    >>> a, b = hyperplane_parameters(triangle)[0]
    >>> integration_reduction(facets, 0, a, b, 1, (x, y), 0)
    5
    """
    if expr == S.Zero:
        return expr

    a, b = (S(a[0]), S(a[1])), S(b)

    value = S.Zero
    x0 = facets[index].points[0]
    m = len(facets)
    gens = (x, y)

    inner_product = diff(expr, gens[0]) * x0[0] + diff(expr, gens[1]) * x0[1]

    if inner_product != 0:
        value += integration_reduction(facets, index, a, b,
                                       inner_product, dims, degree - 1)

    value += left_integral(m, index, facets, x0, expr, gens)

    return value/(len(dims) + degree - 1)


def left_integral(m, index, facets, x0, expr, gens):
    """Computes the left integral of Eq 10 in Chin et al.
    For the 2D case, the integral is just an evaluation of the polynomial
    at the intersection of two facets which is multiplied by the distance
    between the first point of facet and that intersection.

    Parameters
    ===========
    m : No. of hyperplanes.
    index : Index of facet to find intersections with.
    facets : List of facets(Line Segments in 2D case).
    x0 : First point on facet referenced by index.
    expr : Input polynomial
    gens : Generators which generate the polynomial

    >>> from sympy.abc import x, y
    >>> from sympy.integrals.intpoly import left_integral
    >>> from sympy.geometry.point import Point
    >>> from sympy.geometry.polygon import Polygon
    >>> triangle = Polygon(Point(0, 3), Point(5, 3), Point(1, 1))
    >>> facets = triangle.sides
    >>> left_integral(3, 0, facets, facets[0].points[0], 1, (x, y))
    5
    """
    value = S.Zero
    for j in range(0, m):
        intersect = ()
        if j == (index - 1) % m or j == (index + 1) % m:
            intersect = intersection(facets[index], facets[j])
        if intersect:
            distance_origin = norm(tuple(map(lambda x, y: x - y,
                                             intersect, x0)))
            if is_vertex(intersect):
                if isinstance(expr, Expr):
                    if len(gens) == 3:
                        expr_dict = {gens[0]: intersect[0],
                                     gens[1]: intersect[1],
                                     gens[2]: intersect[2]}
                    else:
                        expr_dict = {gens[0]: intersect[0],
                                     gens[1]: intersect[1]}
                    value += distance_origin * expr.subs(expr_dict)
                else:
                    value += distance_origin * expr
    return value


def integration_reduction_dynamic(facets, index, a, b, expr,
                                  dims, x_degree, y_degree, max_y_degree,
                                  x0, monomial_values, monom_index):
    """The same integration_reduction function which uses a dynamic
    programming approach to compute terms by using the values of the integral
    the gradient of previous terms.

    Parameters
    ===========
    facets : Facets of 2-Polytope
    index : Index of facet to find intersections with.(Used in left_integral())
    a, b : Hyperplane parameters
    expr : Input monomial
    dims : Tuple denoting axes variables
    x_degree : Exponent of 'x' in expr
    y_degree : Exponent of 'y' in expr
    max_y_degree : Maximum y-exponent of any monomial in monomial_values
    x0 : First point on facets[index]
    monomial_values : List of monomial values constituting the polynomial
    monom_index : Index of monomial whose integration is being found.

    >>> from sympy.abc import x, y
    >>> from sympy.integrals.intpoly import integration_reduction_dynamic,\
    hyperplane_parameters, gradient_terms
    >>> from sympy.geometry.point import Point
    >>> from sympy.geometry.polygon import Polygon
    >>> triangle = Polygon(Point(0, 3), Point(5, 3), Point(1, 1))
    >>> facets = triangle.sides
    >>> a, b = hyperplane_parameters(triangle)[0]
    >>> x0 = facets[0].points[0]
    >>> monomial_values = [[0, 0, 0, 0], [1, 0, 0, 5],\
                           [y, 0, 1, 15], [x, 1, 0, None]]
    >>> integration_reduction_dynamic(facets, 0, a, b, x, (x, y), 1, 0, 1, x0,\
                                      monomial_values, 3)
    25/2
    """
    expr = S(expr)
    value = S.Zero
    degree = x_degree + y_degree
    m = len(facets)
    gens = (x, y)
    if expr == S.Zero:
        return expr

    if not expr.is_number:
        a, b = (S(a[0]), S(a[1])), S(b)

        x_index = monom_index - max_y_degree + \
            x_degree - 2 if x_degree >= 1 else 0
        y_index = monom_index - 1 if y_degree >= 1 else 0

        x_value, y_value =\
            monomial_values[x_index][3], monomial_values[y_index][3]

        value += x_degree * x_value * x0[0] + y_degree * y_value * x0[1]

    value += left_integral(m, index, facets, x0, expr, gens)

    return value/(len(dims) + degree - 1)


def gradient_terms(binomial_power=0):
    """Returns a list of all the possible
    monomials between 0 and y**binomial_power

    Parameters
    ===========
    binomial_power : Power upto which terms are generated.

    Examples
    ========
    >>> from sympy.abc import x, y
    >>> from sympy.integrals.intpoly import gradient_terms
    >>> gradient_terms(2)
    [[1, 0, 0, None], [y, 0, 1, None], [y**2, 0, 2, None], [x, 1, 0, None],\
    [x*y, 1, 1, None], [x**2, 2, 0, None]]

    """
    terms = []
    for x_count in range(0, binomial_power + 1):
        for y_count in range(0, binomial_power - x_count + 1):
            terms.append([x**x_count*y**y_count,
                          x_count, y_count, None])
    return terms


def hyperplane_parameters(poly):
    """A helper function to return the hyperplane parameters
    of which the facets of the polygon are a part of.
    Currently works for only 2-Polytopes.
    Parameters
    ==========
    poly : The input Polygon

    Examples
    ========
    >>> from sympy.geometry.point import Point
    >>> from sympy.geometry.polygon import Polygon
    >>> from sympy.integrals.intpoly import hyperplane_parameters
    >>> hyperplane_parameters(Polygon(Point(0, 3), Point(5, 3), Point(1, 1)))
    [((0, 1), 3), ((1, -2), -1), ((-2, -1), -3)]
    """
    vertices = list(poly.vertices) + [poly.vertices[0]]  # Close the polygon.
    params = [None] * (len(vertices) - 1)
    for i in range(len(vertices) - 1):
        v1 = vertices[i]
        v2 = vertices[i + 1]

        a1 = v1[1] - v2[1]
        a2 = v2[0] - v1[0]
        b = v2[0] * v1[1] - v2[1] * v1[0]

        factor = gcd_list([a1, a2, b])

        b = S(b)/factor
        a = (S(a1)/factor, S(a2)/factor)
        params[i] = (a, b)

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
    >>> from sympy.integrals.intpoly import best_origin
    >>> from sympy.abc import x, y
    >>> from sympy.geometry.line import Segment2D
    >>> from sympy.geometry.point import Point
    >>> l = Segment2D(Point(0, 3), Point(1, 1))
    >>> expr = x**3*y**7
    >>> best_origin((2, 1), 3, l, expr)
    (0, 3.0)
    """
    a1, b1 = lineseg.points[0]
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

    gens = (x, y)
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
                                power_gens[univariate.args[0]] +=\
                                           univariate.args[1]
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


def decompose(expr, separate=False):
    """Decomposes an input polynomial into homogeneous ones of
    smaller or equal degree.
    Returns a dictionary with keys as the degree of the smaller
    constituting polynomials. Values are the constituting polynomials.
    Parameters
    ==========
    expr : Polynomial(SymPy expression)

    Optional Parameters :

    separate : If True then simply return a list of the constituent monomials
               If not then break up the polynomial into constituent homogeneous
               polynomials.
    Examples
    ========
    >>> from sympy.abc import x, y
    >>> from sympy.integrals.intpoly import decompose
    >>> decompose(x**2 + x*y + x + y + x**3*y**2 + y**5)
    {1: x + y, 2: x**2 + x*y, 5: x**3*y**2 + y**5}
    >>> decompose(x**2 + x*y + x + y + x**3*y**2 + y**5, True)
    [x, y, x**2, y**5, x*y, x**3*y**2]
    """
    expr = S(expr)
    poly_dict = {}

    if isinstance(expr, Expr) and not expr.is_number:
        if expr.is_Symbol:
            poly_dict[1] = expr
        elif expr.is_Add:
            symbols = expr.atoms(Symbol)
            degrees = [(sum(degree_list(monom, *symbols)), monom)
                       for monom in expr.args]
            if separate:
                return [monom[1] for monom in degrees]
            else:
                for monom in degrees:
                    degree, term = monom
                    if poly_dict.get(degree):
                        poly_dict[degree] += term
                    else:
                        poly_dict[degree] = term
        elif expr.is_Pow:
            _, degree = expr.args
            poly_dict[degree] = expr
        else:  # Now expr can only be of `Mul` type
            degree = 0
            for term in expr.args:
                term_type = len(term.args)
                if term_type == 0 and term.is_Symbol:
                    degree += 1
                elif term_type == 2:
                    degree += term.args[1]
            poly_dict[degree] = expr
    else:
        poly_dict[0] = expr

    if separate:
        return list(poly_dict.values())
    return poly_dict


def clockwise_sort(poly):
    """Returns the same polygon with points sorted in clockwise order.

    Note that it's necessary for input points to be sorted in some order
    (clockwise or anti-clockwise) for the algorithm to work. As a convention
    algorithm has been implemented keeping clockwise orientation in mind.

    Parameters
    ==========
    poly: 2-Polytope

    Examples
    ========
    >>> from sympy.integrals.intpoly import clockwise_sort
    >>> from sympy.geometry.point import Point
    >>> from sympy.geometry.polygon import Polygon
    >>> clockwise_sort(Polygon(Point(0, 0), Point(1, 0), Point(1, 1)))
    Triangle(Point2D(1, 1), Point2D(1, 0), Point2D(0, 0))

    """
    n = len(poly.vertices)
    vertices = list(poly.vertices)
    center = Point(sum(map(lambda vertex: vertex.x, poly.vertices)) / n,
                   sum(map(lambda vertex: vertex.y, poly.vertices)) / n)

    def compareTo(a, b):
        if a.x - center.x >= S.Zero and b.x - center.x < S.Zero:
            return S(-1)
        elif a.x - center.x < S.Zero and b.x - center.x >= S.Zero:
            return S(1)
        elif a.x - center.x == S.Zero and b.x - center.x == S.Zero:
            if a.y - center.y >= S.Zero or b.y - center.y >= S.Zero:
                return S(-1) if a.y > b.y else S(1)
            return S(-1) if b.y > a.y else S(1)

        det = (a.x - center.x) * (b.y - center.y) -\
              (b.x - center.x) * (a.y - center.y)
        if det < S.Zero:
            return S(-1)
        elif det > S.Zero:
            return S(1)

        first = (a.x - center.x) * (a.x - center.x) +\
                (a.y - center.y) * (a.y - center.y)
        second = (b.x - center.x) * (b.x - center.x) +\
                 (b.y - center.y) * (b.y - center.y)
        return S(-1) if first > second else S(1)

    return Polygon(*sorted(vertices, key=cmp_to_key(compareTo)))


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
    half = S(1)/2
    if isinstance(point, tuple):
        return (point[0] ** 2 + point[1] ** 2) ** half
    elif isinstance(point, Point):
        return (point.x ** 2 + point.y ** 2) ** half
    elif isinstance(point, dict):
        return sum(i**2 for i in point.values()) ** half


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
    if isinstance(lineseg_1, Segment2D):
        x1, y1 = lineseg_1.points[0]
        x2, y2 = lineseg_1.points[1]
        x3, y3 = lineseg_2.points[0]
        x4, y4 = lineseg_2.points[1]
    else:
        a1x, a1y = S(lineseg_1[0][0]), S(lineseg_1[0][1])
        a2x, a2y = S(lineseg_2[0][0]), S(lineseg_2[0][1])
        b1, b2 = S(lineseg_1[1]), S(lineseg_2[1])

        denom = a1x * a2y - a2x * a1y
        if denom:
            return (S(b1 * a2y - b2 * a1y) / denom,
                    S(b2 * a1x - b1 * a2x) / denom)
        return ()

    denom = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4)

    if denom:
        t1 = x1*y2 - y1*x2
        t2 = x3*y4 - x4*y3
        return(S(t1*(x3 - x4) - t2*(x1 - x2))/denom,
               S(t1*(y3 - y4) - t2*(y1 - y2))/denom)
    return ()


def is_vertex(ent):
    """If the input entity is a vertex return True
    Parameter
    =========
    ent : Denotes a geometric entity representing a point

    Examples
    ========
    >>> from sympy.geometry.point import Point
    >>> from sympy.integrals.intpoly import is_vertex
    >>> is_vertex((2, 3))
    True
    >>> is_vertex((2, 3, 6))
    True
    >>> is_vertex(Point(2, 3))
    True
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
