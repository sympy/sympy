from __future__ import print_function, division

from functools import cmp_to_key

from sympy.core import S, Symbol

from sympy.geometry import Segment2D
from sympy.geometry.polygon import Polygon
from sympy.geometry.point import Point
from sympy.core.expr import Expr
from sympy.abc import x, y

from sympy.polys.polytools import LC
from sympy.polys.polytools import gcd_list


def gradient_terms(x_degree, y_degree, binomial_power=None):
    terms = []
    if binomial_power:
        for x_count in range(0, binomial_power + 1):
            for y_count in range(0, binomial_power - x_count + 1):
                terms.append([x**x_count*y**y_count,
                              x_count, y_count, None])
        return terms

    for x_count in range(0, x_degree + 1):
        for y_count in range(0, y_degree + 1):
            terms.append([x**x_count*y**y_count, x_count, y_count, None])
    return terms


def polytope_integrate(poly, expr, **kwargs):
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
            if not polys in result:
                if polys is S.Zero:
                    result[S.Zero] = S.Zero
                    continue
                integral_value = S.Zero
                monoms = decompose(polys, separate=1)
                for monom in monoms:
                    term = monom[0]
                    if term.is_number:
                        integral_value += result_dict[1][0] * term
                    else:
                        coeff = LC(term)
                        integral_value += result_dict[term / coeff][0] * coeff
                result[polys] = integral_value
        return result

    return main_integrate(expr, facets, hp_params)

def main_integrate(expr, facets, hp_params, max_degree=None):
    dims = (x, y)
    if max_degree:
        expr = 0
    monoms = decompose(expr, separate=1)
    dim_length = len(dims)
    result = {}
    integral_value = S.Zero

    for term_count, term in enumerate(monoms):
        if max_degree:
            find_many = True
            y_degree = max_degree
            grad_terms = [[0, 0, 0, 0]] + \
                         gradient_terms(0, 0, max_degree)
        else:
            find_many = False
            monomial, x_degree, y_degree = term
            grad_terms = [[0, 0, 0, 0]] + \
                         gradient_terms(x_degree, y_degree)
        for facet_count, hp in enumerate(hp_params):
            a, b = hp[0], hp[1]
            x0 = facets[facet_count].points[0]

            for i, monom in enumerate(grad_terms):
                #  Every monomial is a tuple :
                #  (term, x_degree, y_degree, value over boundary)
                m, x_d, y_d, _ = monom
                value = result.get(m, None)
                value_over_boundary = \
                    integration_reduction(facets, facet_count, a, b, m,
                                          dims, x_d, y_d, y_degree, x0,
                                          grad_terms, i, find_many)
                monom[3] = value_over_boundary

                degree = x_d + y_d
                if value is not None:
                    if value[1] == term_count:
                        result[m][0] += value_over_boundary * \
                                        (b / norm(a)) / (dim_length + degree)
                else:
                    result[m] = [value_over_boundary * \
                                 (b / norm(a)) / (dim_length + degree),
                                 term_count]
    if max_degree:
        return result

    for monom in monoms:
        term = monom[0]
        if term.is_number:
            integral_value += result[1][0] * term
        else:
            coeff = LC(term)
            integral_value += result[term / coeff][0] * coeff

    return integral_value


def integration_reduction(facets, index, a, b, expr,
                          dims, x_degree, y_degree, max_y_degree,
                          x0, monomial_values, monom_index, find_many):
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
    expr = S(expr)
    value = S.Zero
    degree = x_degree + y_degree
    m = len(facets)
    gens = (x, y)
    if expr == S.Zero:
        return expr

    if not expr.is_number:
        a, b = (S(a[0]), S(a[1])), S(b)

        if find_many:
            x_index = monom_index - max_y_degree +\
                      x_degree - 2 if x_degree >= 1 else 0
        else:
            x_index = monom_index - max_y_degree - 1 if x_degree >= 1 else 0
        y_index = monom_index - 1 if y_degree >= 1 else 0

        x_value, y_value =\
            monomial_values[x_index][3], monomial_values[y_index][3]

        value += x_degree * x_value * x0[0] + y_degree * y_value * x0[1]

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

        b = S(b)/factor
        a = (S(a1)/factor, S(a2)/factor)
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


def decompose(expr, separate=0):
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
    expr = S(expr)
    if separate:
        monoms = []
    else:
        poly_dict = {}

    if isinstance(expr, Expr) and not expr.is_number:
        if expr.is_Add:
            for monomial in expr.args:
                x_degree = y_degree = 0
                if monomial.is_Pow:
                    if monomial.args[0] is x:
                        x_degree += monomial.args[1]
                    else:
                        y_degree += monomial.args[1]
                else:
                    term_type = len(monomial.args)
                    if term_type == 0:
                        if monomial.is_Symbol:
                            if monomial is x:
                                x_degree += 1
                            else:
                                y_degree += 1
                    for univariate in monomial.args:
                        term_type = len(univariate.args)
                        if term_type == 0 and univariate.is_Symbol:
                            if univariate is x:
                                x_degree += 1
                            else:
                                y_degree += 1
                        elif term_type == 2:
                            if univariate.args[0] is x:
                                x_degree += univariate.args[1]
                            else:
                                y_degree += univariate.args[1]
                if separate:
                    monoms.append((monomial, x_degree, y_degree))
                else:
                    degree = x_degree + y_degree
                    if degree in poly_dict:
                        poly_dict[degree] += monomial
                    else:
                        poly_dict[degree] = monomial
        elif expr.is_Mul:
            x_degree = y_degree = 0
            for term in expr.args:
                term_type = len(term.args)
                if term_type == 0 and term.is_Symbol:
                    if term is x:
                        x_degree += 1
                    else:
                        y_degree += 1
                elif term_type == 2:
                    if term.args[0] is x:
                        x_degree += term.args[1]
                    else:
                        y_degree += term.args[1]
            if separate:
                monoms.append((expr, x_degree, y_degree))
            else:
                poly_dict[x_degree + y_degree] = expr
        elif expr.is_Pow:
            x_degree = y_degree = 0
            if expr.args[0] is x:
                x_degree += expr.args[1]
            else:
                y_degree += expr.args[1]

            if separate:
                monoms.append((expr, x_degree, y_degree))
            else:
                poly_dict[x_degree + y_degree] = expr
        else:
            if separate:
                if expr.is_Symbol:
                    if expr is x:
                        monoms.append((expr, 1, 0))
                    else:
                        monoms.append((expr, 0, 1))
                else:
                    monoms.append((expr, 0, 0))
            else:
                poly_dict[1] = expr
    else:
        if separate:
            monoms.append((expr, 0, 0))
        else:
            poly_dict[0] = expr

    if separate:
        return monoms
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
        if a.x - center.x >= 0 and b.x - center.x < 0:
            return -1
        elif a.x - center.x < 0 and b.x - center.x >= 0:
            return 1
        elif a.x - center.x == 0 and b.x - center.x == 0:
            if a.y - center.y >= 0 or b.y - center.y >= 0:
                return -1 if a.y > b.y else 1
            return -1 if b.y > a.y else 1

        det = (a.x - center.x) * (b.y - center.y) -\
              (b.x - center.x) * (a.y - center.y)
        if det < 0:
            return -1
        elif det > 0:
            return 1

        first = (a.x - center.x) * (a.x - center.x) +\
                (a.y - center.y) * (a.y - center.y)
        second = (b.x - center.x) * (b.x - center.x) +\
                 (b.y - center.y) * (b.y - center.y)
        return -1 if first > second else 1

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
        x_num = (b1 * a2y - b2 * a1y)
        y_num = (b2 * a1x - b1 * a2x)
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
