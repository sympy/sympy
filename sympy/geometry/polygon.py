from __future__ import print_function, division

from sympy.core import Expr, S, sympify, oo, pi, Symbol
from sympy.core.compatibility import as_int, range
from sympy.functions.elementary.complexes import sign
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.elementary.trigonometric import cos, sin, tan
from sympy.geometry.exceptions import GeometryError
from sympy.logic import And
from sympy.matrices import Matrix
from sympy.simplify import simplify
from sympy.utilities import default_sort_key
from sympy.utilities.iterables import has_variety, has_dups, uniq

from .entity import GeometryEntity
from .point import Point
from .ellipse import Circle
from .line import Line, Segment
from .util import _symbol

import warnings


class Polygon(GeometryEntity):
    """A two-dimensional polygon.

    A simple polygon in space. Can be constructed from a sequence of points
    or from a center, radius, number of sides and rotation angle.

    Parameters
    ==========

    vertices : sequence of Points

    Attributes
    ==========

    area
    angles
    perimeter
    vertices
    centroid
    sides

    Raises
    ======

    GeometryError
        If all parameters are not Points.

        If the Polygon has intersecting sides.

    See Also
    ========

    sympy.geometry.point.Point, sympy.geometry.line.Segment, Triangle

    Notes
    =====

    Polygons are treated as closed paths rather than 2D areas so
    some calculations can be be negative or positive (e.g., area)
    based on the orientation of the points.

    Any consecutive identical points are reduced to a single point
    and any points collinear and between two points will be removed
    unless they are needed to define an explicit intersection (see examples).

    A Triangle, Segment or Point will be returned when there are 3 or
    fewer points provided.

    Examples
    ========

    >>> from sympy import Point, Polygon, pi
    >>> p1, p2, p3, p4, p5 = [(0, 0), (1, 0), (5, 1), (0, 1), (3, 0)]
    >>> Polygon(p1, p2, p3, p4)
    Polygon(Point(0, 0), Point(1, 0), Point(5, 1), Point(0, 1))
    >>> Polygon(p1, p2)
    Segment(Point(0, 0), Point(1, 0))
    >>> Polygon(p1, p2, p5)
    Segment(Point(0, 0), Point(3, 0))

    While the sides of a polygon are not allowed to cross implicitly, they
    can do so explicitly. For example, a polygon shaped like a Z with the top
    left connecting to the bottom right of the Z must have the point in the
    middle of the Z explicitly given:

    >>> mid = Point(1, 1)
    >>> Polygon((0, 2), (2, 2), mid, (0, 0), (2, 0), mid).area
    0
    >>> Polygon((0, 2), (2, 2), mid, (2, 0), (0, 0), mid).area
    -2

    When the the keyword `n` is used to define the number of sides of the
    Polygon then a RegularPolygon is created and the other arguments are
    interpreted as center, radius and rotation. The unrotated RegularPolygon
    will always have a vertex at Point(r, 0) where `r` is the radius of the
    circle that circumscribes the RegularPolygon. Its method `spin` can be
    used to increment that angle.

    >>> p = Polygon((0,0), 1, n=3)
    >>> p
    RegularPolygon(Point(0, 0), 1, 3, 0)
    >>> p.vertices[0]
    Point(1, 0)
    >>> p.args[0]
    Point(0, 0)
    >>> p.spin(pi/2)
    >>> p.vertices[0]
    Point(0, 1)

    """

    def __new__(cls, *args, **kwargs):
        if kwargs.get('n', 0):
            n = kwargs.pop('n')
            args = list(args)
            # return a virtual polygon with n sides
            if len(args) == 2:  # center, radius
                args.append(n)
            elif len(args) == 3:  # center, radius, rotation
                args.insert(2, n)
            return RegularPolygon(*args, **kwargs)

        vertices = [Point(a) for a in args]

        # remove consecutive duplicates
        nodup = []
        for p in vertices:
            if nodup and p == nodup[-1]:
                continue
            nodup.append(p)
        if len(nodup) > 1 and nodup[-1] == nodup[0]:
            nodup.pop()  # last point was same as first

        # remove collinear points unless they are shared points
        got = set()
        shared = set()
        for p in nodup:
            if p in got:
                shared.add(p)
            else:
                got.add(p)
        i = -3
        while i < len(nodup) - 3 and len(nodup) > 2:
            a, b, c = nodup[i], nodup[i + 1], nodup[i + 2]
            if b not in shared and Point.is_collinear(a, b, c):
                nodup.pop(i + 1)
                if a == c:
                    nodup.pop(i)
            else:
                i += 1

        vertices = list(nodup)

        if len(vertices) > 3:
            rv = GeometryEntity.__new__(cls, *vertices, **kwargs)
        elif len(vertices) == 3:
            return Triangle(*vertices, **kwargs)
        elif len(vertices) == 2:
            return Segment(*vertices, **kwargs)
        else:
            return Point(*vertices, **kwargs)

        # reject polygons that have intersecting sides unless the
        # intersection is a shared point or a generalized intersection.
        # A self-intersecting polygon is easier to detect than a
        # random set of segments since only those sides that are not
        # part of the convex hull can possibly intersect with other
        # sides of the polygon...but for now we use the n**2 algorithm
        # and check if any side intersects with any preceding side,
        # excluding the ones it is connected to
        try:
            convex = rv.is_convex()
        except ValueError:
            convex = True
        if not convex:
            sides = rv.sides
            for i, si in enumerate(sides):
                pts = si.args
                # exclude the sides connected to si
                for j in range(1 if i == len(sides) - 1 else 0, i - 1):
                    sj = sides[j]
                    if sj.p1 not in pts and sj.p2 not in pts:
                        hit = si.intersection(sj)
                        if not hit:
                            continue
                        hit = hit[0]
                        # don't complain unless the intersection is definite;
                        # if there are symbols present then the intersection
                        # might not occur; this may not be necessary since if
                        # the convex test passed, this will likely pass, too.
                        # But we are about to raise an error anyway so it
                        # won't matter too much.
                        if all(i.is_number for i in hit.args):
                            raise GeometryError(
                                "Polygon has intersecting sides.")

        return rv

    @property
    def area(self):
        """
        The area of the polygon.

        Notes
        =====

        The area calculation can be positive or negative based on the
        orientation of the points.

        See Also
        ========

        sympy.geometry.ellipse.Ellipse.area

        Examples
        ========

        >>> from sympy import Point, Polygon
        >>> p1, p2, p3, p4 = map(Point, [(0, 0), (1, 0), (5, 1), (0, 1)])
        >>> poly = Polygon(p1, p2, p3, p4)
        >>> poly.area
        3

        """
        area = 0
        args = self.args
        for i in range(len(args)):
            x1, y1 = args[i - 1].args
            x2, y2 = args[i].args
            area += x1*y2 - x2*y1
        return simplify(area) / 2

    @staticmethod
    def _isright(a, b, c):
        ba = b - a
        ca = c - a
        t_area = simplify(ba.x*ca.y - ca.x*ba.y)
        res = t_area.is_nonpositive
        if res is None:
            raise ValueError("Can't determine orientation")
        return res

    @property
    def angles(self):
        """The internal angle at each vertex.

        Returns
        =======

        angles : dict
            A dictionary where each key is a vertex and each value is the
            internal angle at that vertex. The vertices are represented as
            Points.

        See Also
        ========

        sympy.geometry.point.Point, sympy.geometry.line.LinearEntity.angle_between

        Examples
        ========

        >>> from sympy import Point, Polygon
        >>> p1, p2, p3, p4 = map(Point, [(0, 0), (1, 0), (5, 1), (0, 1)])
        >>> poly = Polygon(p1, p2, p3, p4)
        >>> poly.angles[p1]
        pi/2
        >>> poly.angles[p2]
        acos(-4*sqrt(17)/17)

        """

        # Determine orientation of points
        args = self.vertices
        cw = self._isright(args[-1], args[0], args[1])

        ret = {}
        for i in range(len(args)):
            a, b, c = args[i - 2], args[i - 1], args[i]
            ang = Line.angle_between(Line(b, a), Line(b, c))
            if cw ^ self._isright(a, b, c):
                ret[b] = 2*S.Pi - ang
            else:
                ret[b] = ang
        return ret

    @property
    def perimeter(self):
        """The perimeter of the polygon.

        Returns
        =======

        perimeter : number or Basic instance

        See Also
        ========

        sympy.geometry.line.Segment.length

        Examples
        ========

        >>> from sympy import Point, Polygon
        >>> p1, p2, p3, p4 = map(Point, [(0, 0), (1, 0), (5, 1), (0, 1)])
        >>> poly = Polygon(p1, p2, p3, p4)
        >>> poly.perimeter
        sqrt(17) + 7
        """
        p = 0
        args = self.vertices
        for i in range(len(args)):
            p += args[i - 1].distance(args[i])
        return simplify(p)

    @property
    def vertices(self):
        """The vertices of the polygon.

        Returns
        =======

        vertices : tuple of Points

        Notes
        =====

        When iterating over the vertices, it is more efficient to index self
        rather than to request the vertices and index them. Only use the
        vertices when you want to process all of them at once. This is even
        more important with RegularPolygons that calculate each vertex.

        See Also
        ========

        sympy.geometry.point.Point

        Examples
        ========

        >>> from sympy import Point, Polygon
        >>> p1, p2, p3, p4 = map(Point, [(0, 0), (1, 0), (5, 1), (0, 1)])
        >>> poly = Polygon(p1, p2, p3, p4)
        >>> poly.vertices
        (Point(0, 0), Point(1, 0), Point(5, 1), Point(0, 1))
        >>> poly.args[0]
        Point(0, 0)

        """
        return self.args

    @property
    def centroid(self):
        """The centroid of the polygon.

        Returns
        =======

        centroid : Point

        See Also
        ========

        sympy.geometry.point.Point, sympy.geometry.util.centroid

        Examples
        ========

        >>> from sympy import Point, Polygon
        >>> p1, p2, p3, p4 = map(Point, [(0, 0), (1, 0), (5, 1), (0, 1)])
        >>> poly = Polygon(p1, p2, p3, p4)
        >>> poly.centroid
        Point(31/18, 11/18)

        """
        A = 1/(6*self.area)
        cx, cy = 0, 0
        args = self.args
        for i in range(len(args)):
            x1, y1 = args[i - 1].args
            x2, y2 = args[i].args
            v = x1*y2 - x2*y1
            cx += v*(x1 + x2)
            cy += v*(y1 + y2)
        return Point(simplify(A*cx), simplify(A*cy))

    @property
    def sides(self):
        """The line segments that form the sides of the polygon.

        Returns
        =======

        sides : list of sides
            Each side is a Segment.

        Notes
        =====

        The Segments that represent the sides are an undirected
        line segment so cannot be used to tell the orientation of
        the polygon.

        See Also
        ========

        sympy.geometry.point.Point, sympy.geometry.line.Segment

        Examples
        ========

        >>> from sympy import Point, Polygon
        >>> p1, p2, p3, p4 = map(Point, [(0, 0), (1, 0), (5, 1), (0, 1)])
        >>> poly = Polygon(p1, p2, p3, p4)
        >>> poly.sides
        [Segment(Point(0, 0), Point(1, 0)),
        Segment(Point(1, 0), Point(5, 1)),
        Segment(Point(0, 1), Point(5, 1)), Segment(Point(0, 0), Point(0, 1))]

        """
        res = []
        args = self.vertices
        for i in range(-len(args), 0):
            res.append(Segment(args[i], args[i + 1]))
        return res

    def is_convex(self):
        """Is the polygon convex?

        A polygon is convex if all its interior angles are less than 180
        degrees.

        Returns
        =======

        is_convex : boolean
            True if this polygon is convex, False otherwise.

        See Also
        ========

        sympy.geometry.util.convex_hull

        Examples
        ========

        >>> from sympy import Point, Polygon
        >>> p1, p2, p3, p4 = map(Point, [(0, 0), (1, 0), (5, 1), (0, 1)])
        >>> poly = Polygon(p1, p2, p3, p4)
        >>> poly.is_convex()
        True

        """

        # Determine orientation of points
        args = self.vertices
        cw = self._isright(args[-2], args[-1], args[0])
        for i in range(1, len(args)):
            if cw ^ self._isright(args[i - 2], args[i - 1], args[i]):
                return False

        return True

    def encloses_point(self, p):
        """
        Return True if p is enclosed by (is inside of) self.

        Notes
        =====

        Being on the border of self is considered False.

        Parameters
        ==========

        p : Point

        Returns
        =======

        encloses_point : True, False or None

        See Also
        ========

        sympy.geometry.point.Point, sympy.geometry.ellipse.Ellipse.encloses_point

        Examples
        ========

        >>> from sympy import Polygon, Point
        >>> from sympy.abc import t
        >>> p = Polygon((0, 0), (4, 0), (4, 4))
        >>> p.encloses_point(Point(2, 1))
        True
        >>> p.encloses_point(Point(2, 2))
        False
        >>> p.encloses_point(Point(5, 5))
        False

        References
        ==========

        [1] http://www.ariel.com.au/a/python-point-int-poly.html

        """
        p = Point(p)
        if p in self.vertices or any(p in s for s in self.sides):
            return False

        # move to p, checking that the result is numeric
        lit = []
        for v in self.vertices:
            lit.append(v - p)  # the difference is simplified
            if lit[-1].free_symbols:
                return None

        poly = Polygon(*lit)

        # polygon closure is assumed in the following test but Polygon removes duplicate pts so
        # the last point has to be added so all sides are computed. Using Polygon.sides is
        # not good since Segments are unordered.
        args = poly.args
        indices = list(range(-len(args), 1))

        if poly.is_convex():
            orientation = None
            for i in indices:
                a = args[i]
                b = args[i + 1]
                test = ((-a.y)*(b.x - a.x) - (-a.x)*(b.y - a.y)).is_negative
                if orientation is None:
                    orientation = test
                elif test is not orientation:
                    return False
            return True

        hit_odd = False
        p1x, p1y = args[0].args
        for i in indices[1:]:
            p2x, p2y = args[i].args
            if 0 > min(p1y, p2y):
                if 0 <= max(p1y, p2y):
                    if 0 <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (-p1y)*(p2x - p1x)/(p2y - p1y) + p1x
                        if p1x == p2x or 0 <= xinters:
                            hit_odd = not hit_odd
            p1x, p1y = p2x, p2y
        return hit_odd

    def arbitrary_point(self, parameter='t'):
        """A parameterized point on the polygon.

        The parameter, varying from 0 to 1, assigns points to the position on
        the perimeter that is that fraction of the total perimeter. So the
        point evaluated at t=1/2 would return the point from the first vertex
        that is 1/2 way around the polygon.

        Parameters
        ==========

        parameter : str, optional
            Default value is 't'.

        Returns
        =======

        arbitrary_point : Point

        Raises
        ======

        ValueError
            When `parameter` already appears in the Polygon's definition.

        See Also
        ========

        sympy.geometry.point.Point

        Examples
        ========

        >>> from sympy import Polygon, S, Symbol
        >>> t = Symbol('t', real=True)
        >>> tri = Polygon((0, 0), (1, 0), (1, 1))
        >>> p = tri.arbitrary_point('t')
        >>> perimeter = tri.perimeter
        >>> s1, s2 = [s.length for s in tri.sides[:2]]
        >>> p.subs(t, (s1 + s2/2)/perimeter)
        Point(1, 1/2)

        """
        t = _symbol(parameter)
        if t.name in (f.name for f in self.free_symbols):
            raise ValueError('Symbol %s already appears in object and cannot be used as a parameter.' % t.name)
        sides = []
        perimeter = self.perimeter
        perim_fraction_start = 0
        for s in self.sides:
            side_perim_fraction = s.length/perimeter
            perim_fraction_end = perim_fraction_start + side_perim_fraction
            pt = s.arbitrary_point(parameter).subs(
                t, (t - perim_fraction_start)/side_perim_fraction)
            sides.append(
                (pt, (And(perim_fraction_start <= t, t < perim_fraction_end))))
            perim_fraction_start = perim_fraction_end
        return Piecewise(*sides)

    def plot_interval(self, parameter='t'):
        """The plot interval for the default geometric plot of the polygon.

        Parameters
        ==========

        parameter : str, optional
            Default value is 't'.

        Returns
        =======

        plot_interval : list (plot interval)
            [parameter, lower_bound, upper_bound]

        Examples
        ========

        >>> from sympy import Polygon
        >>> p = Polygon((0, 0), (1, 0), (1, 1))
        >>> p.plot_interval()
        [t, 0, 1]

        """
        t = Symbol(parameter, real=True)
        return [t, 0, 1]

    def intersection(self, o):
        """The intersection of two polygons.

        The intersection may be empty and can contain individual Points and
        complete Line Segments.

        Parameters
        ==========

        other: Polygon

        Returns
        =======

        intersection : list
            The list of Segments and Points

        See Also
        ========

        sympy.geometry.point.Point, sympy.geometry.line.Segment

        Examples
        ========

        >>> from sympy import Point, Polygon
        >>> p1, p2, p3, p4 = map(Point, [(0, 0), (1, 0), (5, 1), (0, 1)])
        >>> poly1 = Polygon(p1, p2, p3, p4)
        >>> p5, p6, p7 = map(Point, [(3, 2), (1, -1), (0, 2)])
        >>> poly2 = Polygon(p5, p6, p7)
        >>> poly1.intersection(poly2)
        [Point(2/3, 0), Point(9/5, 1/5), Point(7/3, 1), Point(1/3, 1)]

        """
        res = []
        for side in self.sides:
            inter = side.intersection(o)
            if inter is not None:
                res.extend(inter)
        return list(uniq(res))

    def distance(self, o):
        """
        Returns the shortest distance between self and o.

        If o is a point, then self does not need to be convex.
        If o is another polygon self and o must be complex.

        Examples
        ========

        >>> from sympy import Point, Polygon, RegularPolygon
        >>> p1, p2 = map(Point, [(0, 0), (7, 5)])
        >>> poly = Polygon(*RegularPolygon(p1, 1, 3).vertices)
        >>> poly.distance(p2)
        sqrt(61)
        """
        if isinstance(o, Point):
            dist = oo
            for side in self.sides:
                current = side.distance(o)
                if current == 0:
                    return S.Zero
                elif current < dist:
                    dist = current
            return dist
        elif isinstance(o, Polygon) and self.is_convex() and o.is_convex():
            return self._do_poly_distance(o)
        raise NotImplementedError()

    def _do_poly_distance(self, e2):
        """
        Calculates the least distance between the exteriors of two
        convex polygons e1 and e2. Does not check for the convexity
        of the polygons as this is checked by Polygon.distance.

        Notes
        =====

            - Prints a warning if the two polygons possibly intersect as the return
              value will not be valid in such a case. For a more through test of
              intersection use intersection().

        See Also
        ========

        sympy.geometry.point.Point.distance

        Examples
        ========

        >>> from sympy.geometry import Point, Polygon
        >>> square = Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0))
        >>> triangle = Polygon(Point(1, 2), Point(2, 2), Point(2, 1))
        >>> square._do_poly_distance(triangle)
        sqrt(2)/2

        Description of method used
        ==========================

        Method:
        [1] http://cgm.cs.mcgill.ca/~orm/mind2p.html
        Uses rotating calipers:
        [2] http://en.wikipedia.org/wiki/Rotating_calipers
        and antipodal points:
        [3] http://en.wikipedia.org/wiki/Antipodal_point
        """
        e1 = self

        '''Tests for a possible intersection between the polygons and outputs a warning'''
        e1_center = e1.centroid
        e2_center = e2.centroid
        e1_max_radius = S.Zero
        e2_max_radius = S.Zero
        for vertex in e1.vertices:
            r = Point.distance(e1_center, vertex)
            if e1_max_radius < r:
                e1_max_radius = r
        for vertex in e2.vertices:
            r = Point.distance(e2_center, vertex)
            if e2_max_radius < r:
                e2_max_radius = r
        center_dist = Point.distance(e1_center, e2_center)
        if center_dist <= e1_max_radius + e2_max_radius:
            warnings.warn("Polygons may intersect producing erroneous output")

        '''
        Find the upper rightmost vertex of e1 and the lowest leftmost vertex of e2
        '''
        e1_ymax = Point(0, -oo)
        e2_ymin = Point(0, oo)

        for vertex in e1.vertices:
            if vertex.y > e1_ymax.y or (vertex.y == e1_ymax.y and vertex.x > e1_ymax.x):
                e1_ymax = vertex
        for vertex in e2.vertices:
            if vertex.y < e2_ymin.y or (vertex.y == e2_ymin.y and vertex.x < e2_ymin.x):
                e2_ymin = vertex
        min_dist = Point.distance(e1_ymax, e2_ymin)

        '''
        Produce a dictionary with vertices of e1 as the keys and, for each vertex, the points
        to which the vertex is connected as its value. The same is then done for e2.
        '''
        e1_connections = {}
        e2_connections = {}

        for side in e1.sides:
            if side.p1 in e1_connections:
                e1_connections[side.p1].append(side.p2)
            else:
                e1_connections[side.p1] = [side.p2]

            if side.p2 in e1_connections:
                e1_connections[side.p2].append(side.p1)
            else:
                e1_connections[side.p2] = [side.p1]

        for side in e2.sides:
            if side.p1 in e2_connections:
                e2_connections[side.p1].append(side.p2)
            else:
                e2_connections[side.p1] = [side.p2]

            if side.p2 in e2_connections:
                e2_connections[side.p2].append(side.p1)
            else:
                e2_connections[side.p2] = [side.p1]

        e1_current = e1_ymax
        e2_current = e2_ymin
        support_line = Line(Point(S.Zero, S.Zero), Point(S.One, S.Zero))

        '''
        Determine which point in e1 and e2 will be selected after e2_ymin and e1_ymax,
        this information combined with the above produced dictionaries determines the
        path that will be taken around the polygons
        '''
        point1 = e1_connections[e1_ymax][0]
        point2 = e1_connections[e1_ymax][1]
        angle1 = support_line.angle_between(Line(e1_ymax, point1))
        angle2 = support_line.angle_between(Line(e1_ymax, point2))
        if angle1 < angle2:
            e1_next = point1
        elif angle2 < angle1:
            e1_next = point2
        elif Point.distance(e1_ymax, point1) > Point.distance(e1_ymax, point2):
            e1_next = point2
        else:
            e1_next = point1

        point1 = e2_connections[e2_ymin][0]
        point2 = e2_connections[e2_ymin][1]
        angle1 = support_line.angle_between(Line(e2_ymin, point1))
        angle2 = support_line.angle_between(Line(e2_ymin, point2))
        if angle1 > angle2:
            e2_next = point1
        elif angle2 > angle1:
            e2_next = point2
        elif Point.distance(e2_ymin, point1) > Point.distance(e2_ymin, point2):
            e2_next = point2
        else:
            e2_next = point1

        '''
        Loop which determins the distance between anti-podal pairs and updates the
        minimum distance accordingly. It repeats until it reaches the starting position.
        '''
        while True:
            e1_angle = support_line.angle_between(Line(e1_current, e1_next))
            e2_angle = pi - support_line.angle_between(Line(
                e2_current, e2_next))

            if (e1_angle < e2_angle) is True:
                support_line = Line(e1_current, e1_next)
                e1_segment = Segment(e1_current, e1_next)
                min_dist_current = e1_segment.distance(e2_current)

                if min_dist_current.evalf() < min_dist.evalf():
                    min_dist = min_dist_current

                if e1_connections[e1_next][0] != e1_current:
                    e1_current = e1_next
                    e1_next = e1_connections[e1_next][0]
                else:
                    e1_current = e1_next
                    e1_next = e1_connections[e1_next][1]
            elif (e1_angle > e2_angle) is True:
                support_line = Line(e2_next, e2_current)
                e2_segment = Segment(e2_current, e2_next)
                min_dist_current = e2_segment.distance(e1_current)

                if min_dist_current.evalf() < min_dist.evalf():
                    min_dist = min_dist_current

                if e2_connections[e2_next][0] != e2_current:
                    e2_current = e2_next
                    e2_next = e2_connections[e2_next][0]
                else:
                    e2_current = e2_next
                    e2_next = e2_connections[e2_next][1]
            else:
                support_line = Line(e1_current, e1_next)
                e1_segment = Segment(e1_current, e1_next)
                e2_segment = Segment(e2_current, e2_next)
                min1 = e1_segment.distance(e2_next)
                min2 = e2_segment.distance(e1_next)

                min_dist_current = min(min1, min2)
                if min_dist_current.evalf() < min_dist.evalf():
                    min_dist = min_dist_current

                if e1_connections[e1_next][0] != e1_current:
                    e1_current = e1_next
                    e1_next = e1_connections[e1_next][0]
                else:
                    e1_current = e1_next
                    e1_next = e1_connections[e1_next][1]

                if e2_connections[e2_next][0] != e2_current:
                    e2_current = e2_next
                    e2_next = e2_connections[e2_next][0]
                else:
                    e2_current = e2_next
                    e2_next = e2_connections[e2_next][1]
            if e1_current == e1_ymax and e2_current == e2_ymin:
                break
        return min_dist

    def __eq__(self, o):
        if not isinstance(o, Polygon) or len(self.args) != len(o.args):
            return False

        # See if self can ever be traversed (cw or ccw) from any of its
        # vertices to match all points of o
        args = self.args
        oargs = o.args
        n = len(args)
        o0 = oargs[0]
        for i0 in range(n):
            if args[i0] == o0:
                if all(args[(i0 + i) % n] == oargs[i] for i in range(1, n)):
                    return True
                if all(args[(i0 - i) % n] == oargs[i] for i in range(1, n)):
                    return True
        return False

    def __hash__(self):
        return super(Polygon, self).__hash__()

    def __contains__(self, o):
        """
        Return True if o is contained within the boundary lines of self.altitudes

        Parameters
        ==========

        other : GeometryEntity

        Returns
        =======

        contained in : bool
            The points (and sides, if applicable) are contained in self.

        See Also
        ========

        sympy.geometry.entity.GeometryEntity.encloses

        Examples
        ========

        >>> from sympy import Line, Segment, Point
        >>> p = Point(0, 0)
        >>> q = Point(1, 1)
        >>> s = Segment(p, q*2)
        >>> l = Line(p, q)
        >>> p in q
        False
        >>> p in s
        True
        >>> q*3 in s
        False
        >>> s in l
        True

        """

        if isinstance(o, Polygon):
            return self == o
        elif isinstance(o, Segment):
            return any(o in s for s in self.sides)
        elif isinstance(o, Point):
            if o in self.vertices:
                return True
            for side in self.sides:
                if o in side:
                    return True

        return False


class RegularPolygon(Polygon):
    """
    A regular polygon.

    Such a polygon has all internal angles equal and all sides the same length.

    Parameters
    ==========

    center : Point
    radius : number or Basic instance
        The distance from the center to a vertex
    n : int
        The number of sides

    Attributes
    ==========

    vertices
    center
    radius
    rotation
    apothem
    interior_angle
    exterior_angle
    circumcircle
    incircle
    angles

    Raises
    ======

    GeometryError
        If the `center` is not a Point, or the `radius` is not a number or Basic
        instance, or the number of sides, `n`, is less than three.

    Notes
    =====

    A RegularPolygon can be instantiated with Polygon with the kwarg n.

    Regular polygons are instantiated with a center, radius, number of sides
    and a rotation angle. Whereas the arguments of a Polygon are vertices, the
    vertices of the RegularPolygon must be obtained with the vertices method.

    See Also
    ========

    sympy.geometry.point.Point, Polygon

    Examples
    ========

    >>> from sympy.geometry import RegularPolygon, Point
    >>> r = RegularPolygon(Point(0, 0), 5, 3)
    >>> r
    RegularPolygon(Point(0, 0), 5, 3, 0)
    >>> r.vertices[0]
    Point(5, 0)

    """

    __slots__ = ['_n', '_center', '_radius', '_rot']

    def __new__(self, c, r, n, rot=0, **kwargs):
        r, n, rot = map(sympify, (r, n, rot))
        c = Point(c)
        if not isinstance(r, Expr):
            raise GeometryError("r must be an Expr object, not %s" % r)
        if n.is_Number:
            as_int(n)  # let an error raise if necessary
            if n < 3:
                raise GeometryError("n must be a >= 3, not %s" % n)

        obj = GeometryEntity.__new__(self, c, r, n, **kwargs)
        obj._n = n
        obj._center = c
        obj._radius = r
        obj._rot = rot
        return obj

    @property
    def args(self):
        """
        Returns the center point, the radius,
        the number of sides, and the orientation angle.

        Examples
        ========

        >>> from sympy import RegularPolygon, Point
        >>> r = RegularPolygon(Point(0, 0), 5, 3)
        >>> r.args
        (Point(0, 0), 5, 3, 0)
        """
        return self._center, self._radius, self._n, self._rot

    def __str__(self):
        return 'RegularPolygon(%s, %s, %s, %s)' % tuple(self.args)

    def __repr__(self):
        return 'RegularPolygon(%s, %s, %s, %s)' % tuple(self.args)

    @property
    def area(self):
        """Returns the area.

        Examples
        ========

        >>> from sympy.geometry import RegularPolygon
        >>> square = RegularPolygon((0, 0), 1, 4)
        >>> square.area
        2
        >>> _ == square.length**2
        True
        """
        c, r, n, rot = self.args
        return sign(r)*n*self.length**2/(4*tan(pi/n))

    @property
    def length(self):
        """Returns the length of the sides.

        The half-length of the side and the apothem form two legs
        of a right triangle whose hypotenuse is the radius of the
        regular polygon.

        Examples
        ========

        >>> from sympy.geometry import RegularPolygon
        >>> from sympy import sqrt
        >>> s = square_in_unit_circle = RegularPolygon((0, 0), 1, 4)
        >>> s.length
        sqrt(2)
        >>> sqrt((_/2)**2 + s.apothem**2) == s.radius
        True

        """
        return self.radius*2*sin(pi/self._n)

    @property
    def center(self):
        """The center of the RegularPolygon

        This is also the center of the circumscribing circle.

        Returns
        =======

        center : Point

        See Also
        ========

        sympy.geometry.point.Point, sympy.geometry.ellipse.Ellipse.center

        Examples
        ========

        >>> from sympy.geometry import RegularPolygon, Point
        >>> rp = RegularPolygon(Point(0, 0), 5, 4)
        >>> rp.center
        Point(0, 0)
        """
        return self._center

    centroid = center

    @property
    def circumcenter(self):
        """
        Alias for center.

        Examples
        ========

        >>> from sympy.geometry import RegularPolygon, Point
        >>> rp = RegularPolygon(Point(0, 0), 5, 4)
        >>> rp.circumcenter
        Point(0, 0)
        """
        return self.center

    @property
    def radius(self):
        """Radius of the RegularPolygon

        This is also the radius of the circumscribing circle.

        Returns
        =======

        radius : number or instance of Basic

        See Also
        ========

        sympy.geometry.line.Segment.length, sympy.geometry.ellipse.Circle.radius

        Examples
        ========

        >>> from sympy import Symbol
        >>> from sympy.geometry import RegularPolygon, Point
        >>> radius = Symbol('r')
        >>> rp = RegularPolygon(Point(0, 0), radius, 4)
        >>> rp.radius
        r

        """
        return self._radius

    @property
    def circumradius(self):
        """
        Alias for radius.

        Examples
        ========

        >>> from sympy import Symbol
        >>> from sympy.geometry import RegularPolygon, Point
        >>> radius = Symbol('r')
        >>> rp = RegularPolygon(Point(0, 0), radius, 4)
        >>> rp.circumradius
        r
        """
        return self.radius

    @property
    def rotation(self):
        """CCW angle by which the RegularPolygon is rotated

        Returns
        =======

        rotation : number or instance of Basic

        Examples
        ========

        >>> from sympy import pi
        >>> from sympy.geometry import RegularPolygon, Point
        >>> RegularPolygon(Point(0, 0), 3, 4, pi).rotation
        pi

        """
        return self._rot

    @property
    def apothem(self):
        """The inradius of the RegularPolygon.

        The apothem/inradius is the radius of the inscribed circle.

        Returns
        =======

        apothem : number or instance of Basic

        See Also
        ========

        sympy.geometry.line.Segment.length, sympy.geometry.ellipse.Circle.radius

        Examples
        ========

        >>> from sympy import Symbol
        >>> from sympy.geometry import RegularPolygon, Point
        >>> radius = Symbol('r')
        >>> rp = RegularPolygon(Point(0, 0), radius, 4)
        >>> rp.apothem
        sqrt(2)*r/2

        """
        return self.radius * cos(S.Pi/self._n)

    @property
    def inradius(self):
        """
        Alias for apothem.

        Examples
        ========

        >>> from sympy import Symbol
        >>> from sympy.geometry import RegularPolygon, Point
        >>> radius = Symbol('r')
        >>> rp = RegularPolygon(Point(0, 0), radius, 4)
        >>> rp.inradius
        sqrt(2)*r/2
        """
        return self.apothem

    @property
    def interior_angle(self):
        """Measure of the interior angles.

        Returns
        =======

        interior_angle : number

        See Also
        ========

        sympy.geometry.line.LinearEntity.angle_between

        Examples
        ========

        >>> from sympy.geometry import RegularPolygon, Point
        >>> rp = RegularPolygon(Point(0, 0), 4, 8)
        >>> rp.interior_angle
        3*pi/4

        """
        return (self._n - 2)*S.Pi/self._n

    @property
    def exterior_angle(self):
        """Measure of the exterior angles.

        Returns
        =======

        exterior_angle : number

        See Also
        ========

        sympy.geometry.line.LinearEntity.angle_between

        Examples
        ========

        >>> from sympy.geometry import RegularPolygon, Point
        >>> rp = RegularPolygon(Point(0, 0), 4, 8)
        >>> rp.exterior_angle
        pi/4

        """
        return 2*S.Pi/self._n

    @property
    def circumcircle(self):
        """The circumcircle of the RegularPolygon.

        Returns
        =======

        circumcircle : Circle

        See Also
        ========

        circumcenter, sympy.geometry.ellipse.Circle

        Examples
        ========

        >>> from sympy.geometry import RegularPolygon, Point
        >>> rp = RegularPolygon(Point(0, 0), 4, 8)
        >>> rp.circumcircle
        Circle(Point(0, 0), 4)

        """
        return Circle(self.center, self.radius)

    @property
    def incircle(self):
        """The incircle of the RegularPolygon.

        Returns
        =======

        incircle : Circle

        See Also
        ========

        inradius, sympy.geometry.ellipse.Circle

        Examples
        ========

        >>> from sympy.geometry import RegularPolygon, Point
        >>> rp = RegularPolygon(Point(0, 0), 4, 7)
        >>> rp.incircle
        Circle(Point(0, 0), 4*cos(pi/7))

        """
        return Circle(self.center, self.apothem)

    @property
    def angles(self):
        """
        Returns a dictionary with keys, the vertices of the Polygon,
        and values, the interior angle at each vertex.

        Examples
        ========

        >>> from sympy import RegularPolygon, Point
        >>> r = RegularPolygon(Point(0, 0), 5, 3)
        >>> r.angles
        {Point(-5/2, -5*sqrt(3)/2): pi/3,
         Point(-5/2, 5*sqrt(3)/2): pi/3,
         Point(5, 0): pi/3}
        """
        ret = {}
        ang = self.interior_angle
        for v in self.vertices:
            ret[v] = ang
        return ret

    def encloses_point(self, p):
        """
        Return True if p is enclosed by (is inside of) self.

        Notes
        =====

        Being on the border of self is considered False.

        The general Polygon.encloses_point method is called only if
        a point is not within or beyond the incircle or circumcircle,
        respectively.

        Parameters
        ==========

        p : Point

        Returns
        =======

        encloses_point : True, False or None

        See Also
        ========

        sympy.geometry.ellipse.Ellipse.encloses_point

        Examples
        ========

        >>> from sympy import RegularPolygon, S, Point, Symbol
        >>> p = RegularPolygon((0, 0), 3, 4)
        >>> p.encloses_point(Point(0, 0))
        True
        >>> r, R = p.inradius, p.circumradius
        >>> p.encloses_point(Point((r + R)/2, 0))
        True
        >>> p.encloses_point(Point(R/2, R/2 + (R - r)/10))
        False
        >>> t = Symbol('t', real=True)
        >>> p.encloses_point(p.arbitrary_point().subs(t, S.Half))
        False
        >>> p.encloses_point(Point(5, 5))
        False

        """

        c = self.center
        d = Segment(c, p).length
        if d >= self.radius:
            return False
        elif d < self.inradius:
            return True
        else:
            # now enumerate the RegularPolygon like a general polygon.
            return Polygon.encloses_point(self, p)

    def spin(self, angle):
        """Increment *in place* the virtual Polygon's rotation by ccw angle.

        See also: rotate method which moves the center.

        >>> from sympy import Polygon, Point, pi
        >>> r = Polygon(Point(0,0), 1, n=3)
        >>> r.vertices[0]
        Point(1, 0)
        >>> r.spin(pi/6)
        >>> r.vertices[0]
        Point(sqrt(3)/2, 1/2)

        See Also
        ========

        rotation
        rotate : Creates a copy of the RegularPolygon rotated about a Point

        """
        self._rot += angle

    def rotate(self, angle, pt=None):
        """Override GeometryEntity.rotate to first rotate the RegularPolygon
        about its center.

        >>> from sympy import Point, RegularPolygon, Polygon, pi
        >>> t = RegularPolygon(Point(1, 0), 1, 3)
        >>> t.vertices[0] # vertex on x-axis
        Point(2, 0)
        >>> t.rotate(pi/2).vertices[0] # vertex on y axis now
        Point(0, 2)

        See Also
        ========

        rotation
        spin : Rotates a RegularPolygon in place

        """

        r = type(self)(*self.args)  # need a copy or else changes are in-place
        r._rot += angle
        return GeometryEntity.rotate(r, angle, pt)

    def scale(self, x=1, y=1, pt=None):
        """Override GeometryEntity.scale since it is the radius that must be
        scaled (if x == y) or else a new Polygon must be returned.

        >>> from sympy import RegularPolygon

        Symmetric scaling returns a RegularPolygon:

        >>> RegularPolygon((0, 0), 1, 4).scale(2, 2)
        RegularPolygon(Point(0, 0), 2, 4, 0)

        Asymmetric scaling returns a kite as a Polygon:

        >>> RegularPolygon((0, 0), 1, 4).scale(2, 1)
        Polygon(Point(2, 0), Point(0, 1), Point(-2, 0), Point(0, -1))

        """
        if pt:
            pt = Point(pt)
            return self.translate(*(-pt).args).scale(x, y).translate(*pt.args)
        if x != y:
            return Polygon(*self.vertices).scale(x, y)
        c, r, n, rot = self.args
        r *= x
        return self.func(c, r, n, rot)

    def reflect(self, line):
        """Override GeometryEntity.reflect since this is not made of only
        points.

        >>> from sympy import RegularPolygon, Line

        >>> RegularPolygon((0, 0), 1, 4).reflect(Line((0, 1), slope=-2))
        RegularPolygon(Point(4/5, 2/5), -1, 4, acos(3/5))

        """
        c, r, n, rot = self.args
        cc = c.reflect(line)
        v = self.vertices[0]
        vv = v.reflect(line)
        # see how much it must get spun at the new center
        ang = Segment(cc, vv).angle_between(Segment(c, v))
        rot = (rot + ang + pi) % (2*pi/n)
        return self.func(cc, -r, n, rot)

    @property
    def vertices(self):
        """The vertices of the RegularPolygon.

        Returns
        =======

        vertices : list
            Each vertex is a Point.

        See Also
        ========

        sympy.geometry.point.Point

        Examples
        ========

        >>> from sympy.geometry import RegularPolygon, Point
        >>> rp = RegularPolygon(Point(0, 0), 5, 4)
        >>> rp.vertices
        [Point(5, 0), Point(0, 5), Point(-5, 0), Point(0, -5)]

        """
        c = self._center
        r = abs(self._radius)
        rot = self._rot
        v = 2*S.Pi/self._n

        return [Point(c.x + r*cos(k*v + rot), c.y + r*sin(k*v + rot))
                for k in range(self._n)]

    def __eq__(self, o):
        if not isinstance(o, Polygon):
            return False
        elif not isinstance(o, RegularPolygon):
            return Polygon.__eq__(o, self)
        return self.args == o.args

    def __hash__(self):
        return super(RegularPolygon, self).__hash__()


class Triangle(Polygon):
    """
    A polygon with three vertices and three sides.

    Parameters
    ==========

    points : sequence of Points
    keyword: asa, sas, or sss to specify sides/angles of the triangle

    Attributes
    ==========

    vertices
    altitudes
    orthocenter
    circumcenter
    circumradius
    circumcircle
    inradius
    incircle
    medians
    medial

    Raises
    ======

    GeometryError
        If the number of vertices is not equal to three, or one of the vertices
        is not a Point, or a valid keyword is not given.

    See Also
    ========

    sympy.geometry.point.Point, Polygon

    Examples
    ========

    >>> from sympy.geometry import Triangle, Point
    >>> Triangle(Point(0, 0), Point(4, 0), Point(4, 3))
    Triangle(Point(0, 0), Point(4, 0), Point(4, 3))

    Keywords sss, sas, or asa can be used to give the desired
    side lengths (in order) and interior angles (in degrees) that
    define the triangle:

    >>> Triangle(sss=(3, 4, 5))
    Triangle(Point(0, 0), Point(3, 0), Point(3, 4))
    >>> Triangle(asa=(30, 1, 30))
    Triangle(Point(0, 0), Point(1, 0), Point(1/2, sqrt(3)/6))
    >>> Triangle(sas=(1, 45, 2))
    Triangle(Point(0, 0), Point(2, 0), Point(sqrt(2)/2, sqrt(2)/2))

    """

    def __new__(cls, *args, **kwargs):
        if len(args) != 3:
            if 'sss' in kwargs:
                return _sss(*[simplify(a) for a in kwargs['sss']])
            if 'asa' in kwargs:
                return _asa(*[simplify(a) for a in kwargs['asa']])
            if 'sas' in kwargs:
                return _sas(*[simplify(a) for a in kwargs['sas']])
            msg = "Triangle instantiates with three points or a valid keyword."
            raise GeometryError(msg)

        vertices = [Point(a) for a in args]

        # remove consecutive duplicates
        nodup = []
        for p in vertices:
            if nodup and p == nodup[-1]:
                continue
            nodup.append(p)
        if len(nodup) > 1 and nodup[-1] == nodup[0]:
            nodup.pop()  # last point was same as first

        # remove collinear points
        i = -3
        while i < len(nodup) - 3 and len(nodup) > 2:
            a, b, c = sorted(
                [nodup[i], nodup[i + 1], nodup[i + 2]], key=default_sort_key)
            if Point.is_collinear(a, b, c):
                nodup[i] = a
                nodup[i + 1] = None
                nodup.pop(i + 1)
            i += 1

        vertices = list(filter(lambda x: x is not None, nodup))

        if len(vertices) == 3:
            return GeometryEntity.__new__(cls, *vertices, **kwargs)
        elif len(vertices) == 2:
            return Segment(*vertices, **kwargs)
        else:
            return Point(*vertices, **kwargs)

    @property
    def vertices(self):
        """The triangle's vertices

        Returns
        =======

        vertices : tuple
            Each element in the tuple is a Point

        See Also
        ========

        sympy.geometry.point.Point

        Examples
        ========

        >>> from sympy.geometry import Triangle, Point
        >>> t = Triangle(Point(0, 0), Point(4, 0), Point(4, 3))
        >>> t.vertices
        (Point(0, 0), Point(4, 0), Point(4, 3))

        """
        return self.args

    def is_similar(t1, t2):
        """Is another triangle similar to this one.

        Two triangles are similar if one can be uniformly scaled to the other.

        Parameters
        ==========

        other: Triangle

        Returns
        =======

        is_similar : boolean

        See Also
        ========

        sympy.geometry.entity.GeometryEntity.is_similar

        Examples
        ========

        >>> from sympy.geometry import Triangle, Point
        >>> t1 = Triangle(Point(0, 0), Point(4, 0), Point(4, 3))
        >>> t2 = Triangle(Point(0, 0), Point(-4, 0), Point(-4, -3))
        >>> t1.is_similar(t2)
        True

        >>> t2 = Triangle(Point(0, 0), Point(-4, 0), Point(-4, -4))
        >>> t1.is_similar(t2)
        False

        """
        if not isinstance(t2, Polygon):
            return False

        s1_1, s1_2, s1_3 = [side.length for side in t1.sides]
        s2 = [side.length for side in t2.sides]

        def _are_similar(u1, u2, u3, v1, v2, v3):
            e1 = simplify(u1/v1)
            e2 = simplify(u2/v2)
            e3 = simplify(u3/v3)
            return bool(e1 == e2) and bool(e2 == e3)

        # There's only 6 permutations, so write them out
        return _are_similar(s1_1, s1_2, s1_3, *s2) or \
            _are_similar(s1_1, s1_3, s1_2, *s2) or \
            _are_similar(s1_2, s1_1, s1_3, *s2) or \
            _are_similar(s1_2, s1_3, s1_1, *s2) or \
            _are_similar(s1_3, s1_1, s1_2, *s2) or \
            _are_similar(s1_3, s1_2, s1_1, *s2)

    def is_equilateral(self):
        """Are all the sides the same length?

        Returns
        =======

        is_equilateral : boolean

        See Also
        ========

        sympy.geometry.entity.GeometryEntity.is_similar, RegularPolygon
        is_isosceles, is_right, is_scalene

        Examples
        ========

        >>> from sympy.geometry import Triangle, Point
        >>> t1 = Triangle(Point(0, 0), Point(4, 0), Point(4, 3))
        >>> t1.is_equilateral()
        False

        >>> from sympy import sqrt
        >>> t2 = Triangle(Point(0, 0), Point(10, 0), Point(5, 5*sqrt(3)))
        >>> t2.is_equilateral()
        True

        """
        return not has_variety(s.length for s in self.sides)

    def is_isosceles(self):
        """Are two or more of the sides the same length?

        Returns
        =======

        is_isosceles : boolean

        See Also
        ========

        is_equilateral, is_right, is_scalene

        Examples
        ========

        >>> from sympy.geometry import Triangle, Point
        >>> t1 = Triangle(Point(0, 0), Point(4, 0), Point(2, 4))
        >>> t1.is_isosceles()
        True

        """
        return has_dups(s.length for s in self.sides)

    def is_scalene(self):
        """Are all the sides of the triangle of different lengths?

        Returns
        =======

        is_scalene : boolean

        See Also
        ========

        is_equilateral, is_isosceles, is_right

        Examples
        ========

        >>> from sympy.geometry import Triangle, Point
        >>> t1 = Triangle(Point(0, 0), Point(4, 0), Point(1, 4))
        >>> t1.is_scalene()
        True

        """
        return not has_dups(s.length for s in self.sides)

    def is_right(self):
        """Is the triangle right-angled.

        Returns
        =======

        is_right : boolean

        See Also
        ========

        sympy.geometry.line.LinearEntity.is_perpendicular
        is_equilateral, is_isosceles, is_scalene

        Examples
        ========

        >>> from sympy.geometry import Triangle, Point
        >>> t1 = Triangle(Point(0, 0), Point(4, 0), Point(4, 3))
        >>> t1.is_right()
        True

        """
        s = self.sides
        return Segment.is_perpendicular(s[0], s[1]) or \
            Segment.is_perpendicular(s[1], s[2]) or \
            Segment.is_perpendicular(s[0], s[2])

    @property
    def altitudes(self):
        """The altitudes of the triangle.

        An altitude of a triangle is a segment through a vertex,
        perpendicular to the opposite side, with length being the
        height of the vertex measured from the line containing the side.

        Returns
        =======

        altitudes : dict
            The dictionary consists of keys which are vertices and values
            which are Segments.

        See Also
        ========

        sympy.geometry.point.Point, sympy.geometry.line.Segment.length

        Examples
        ========

        >>> from sympy.geometry import Point, Triangle
        >>> p1, p2, p3 = Point(0, 0), Point(1, 0), Point(0, 1)
        >>> t = Triangle(p1, p2, p3)
        >>> t.altitudes[p1]
        Segment(Point(0, 0), Point(1/2, 1/2))

        """
        s = self.sides
        v = self.vertices
        return {v[0]: s[1].perpendicular_segment(v[0]),
                v[1]: s[2].perpendicular_segment(v[1]),
                v[2]: s[0].perpendicular_segment(v[2])}

    @property
    def orthocenter(self):
        """The orthocenter of the triangle.

        The orthocenter is the intersection of the altitudes of a triangle.
        It may lie inside, outside or on the triangle.

        Returns
        =======

        orthocenter : Point

        See Also
        ========

        sympy.geometry.point.Point

        Examples
        ========

        >>> from sympy.geometry import Point, Triangle
        >>> p1, p2, p3 = Point(0, 0), Point(1, 0), Point(0, 1)
        >>> t = Triangle(p1, p2, p3)
        >>> t.orthocenter
        Point(0, 0)

        """
        a = self.altitudes
        v = self.vertices
        return Line(a[v[0]]).intersection(Line(a[v[1]]))[0]

    @property
    def circumcenter(self):
        """The circumcenter of the triangle

        The circumcenter is the center of the circumcircle.

        Returns
        =======

        circumcenter : Point

        See Also
        ========

        sympy.geometry.point.Point

        Examples
        ========

        >>> from sympy.geometry import Point, Triangle
        >>> p1, p2, p3 = Point(0, 0), Point(1, 0), Point(0, 1)
        >>> t = Triangle(p1, p2, p3)
        >>> t.circumcenter
        Point(1/2, 1/2)
        """
        a, b, c = [x.perpendicular_bisector() for x in self.sides]
        return a.intersection(b)[0]

    @property
    def circumradius(self):
        """The radius of the circumcircle of the triangle.

        Returns
        =======

        circumradius : number of Basic instance

        See Also
        ========

        sympy.geometry.ellipse.Circle.radius

        Examples
        ========

        >>> from sympy import Symbol
        >>> from sympy.geometry import Point, Triangle
        >>> a = Symbol('a')
        >>> p1, p2, p3 = Point(0, 0), Point(1, 0), Point(0, a)
        >>> t = Triangle(p1, p2, p3)
        >>> t.circumradius
        sqrt(a**2/4 + 1/4)
        """
        return Point.distance(self.circumcenter, self.vertices[0])

    @property
    def circumcircle(self):
        """The circle which passes through the three vertices of the triangle.

        Returns
        =======

        circumcircle : Circle

        See Also
        ========

        sympy.geometry.ellipse.Circle

        Examples
        ========

        >>> from sympy.geometry import Point, Triangle
        >>> p1, p2, p3 = Point(0, 0), Point(1, 0), Point(0, 1)
        >>> t = Triangle(p1, p2, p3)
        >>> t.circumcircle
        Circle(Point(1/2, 1/2), sqrt(2)/2)

        """
        return Circle(self.circumcenter, self.circumradius)

    def bisectors(self):
        """The angle bisectors of the triangle.

        An angle bisector of a triangle is a straight line through a vertex
        which cuts the corresponding angle in half.

        Returns
        =======

        bisectors : dict
            Each key is a vertex (Point) and each value is the corresponding
            bisector (Segment).

        See Also
        ========

        sympy.geometry.point.Point, sympy.geometry.line.Segment

        Examples
        ========

        >>> from sympy.geometry import Point, Triangle, Segment
        >>> p1, p2, p3 = Point(0, 0), Point(1, 0), Point(0, 1)
        >>> t = Triangle(p1, p2, p3)
        >>> from sympy import sqrt
        >>> t.bisectors()[p2] == Segment(Point(0, sqrt(2) - 1), Point(1, 0))
        True

        """
        s = self.sides
        v = self.vertices
        c = self.incenter
        l1 = Segment(v[0], Line(v[0], c).intersection(s[1])[0])
        l2 = Segment(v[1], Line(v[1], c).intersection(s[2])[0])
        l3 = Segment(v[2], Line(v[2], c).intersection(s[0])[0])
        return {v[0]: l1, v[1]: l2, v[2]: l3}

    @property
    def incenter(self):
        """The center of the incircle.

        The incircle is the circle which lies inside the triangle and touches
        all three sides.

        Returns
        =======

        incenter : Point

        See Also
        ========

        incircle, sympy.geometry.point.Point

        Examples
        ========

        >>> from sympy.geometry import Point, Triangle
        >>> p1, p2, p3 = Point(0, 0), Point(1, 0), Point(0, 1)
        >>> t = Triangle(p1, p2, p3)
        >>> t.incenter
        Point(-sqrt(2)/2 + 1, -sqrt(2)/2 + 1)

        """
        s = self.sides
        l = Matrix([s[i].length for i in [1, 2, 0]])
        p = sum(l)
        v = self.vertices
        x = simplify(l.dot(Matrix([vi.x for vi in v]))/p)
        y = simplify(l.dot(Matrix([vi.y for vi in v]))/p)
        return Point(x, y)

    @property
    def inradius(self):
        """The radius of the incircle.

        Returns
        =======

        inradius : number of Basic instance

        See Also
        ========

        incircle, sympy.geometry.ellipse.Circle.radius

        Examples
        ========

        >>> from sympy.geometry import Point, Triangle
        >>> p1, p2, p3 = Point(0, 0), Point(4, 0), Point(0, 3)
        >>> t = Triangle(p1, p2, p3)
        >>> t.inradius
        1

        """
        return simplify(2 * self.area / self.perimeter)

    @property
    def incircle(self):
        """The incircle of the triangle.

        The incircle is the circle which lies inside the triangle and touches
        all three sides.

        Returns
        =======

        incircle : Circle

        See Also
        ========

        sympy.geometry.ellipse.Circle

        Examples
        ========

        >>> from sympy.geometry import Point, Triangle
        >>> p1, p2, p3 = Point(0, 0), Point(2, 0), Point(0, 2)
        >>> t = Triangle(p1, p2, p3)
        >>> t.incircle
        Circle(Point(-sqrt(2) + 2, -sqrt(2) + 2), -sqrt(2) + 2)

        """
        return Circle(self.incenter, self.inradius)

    @property
    def medians(self):
        """The medians of the triangle.

        A median of a triangle is a straight line through a vertex and the
        midpoint of the opposite side, and divides the triangle into two
        equal areas.

        Returns
        =======

        medians : dict
            Each key is a vertex (Point) and each value is the median (Segment)
            at that point.

        See Also
        ========

        sympy.geometry.point.Point.midpoint, sympy.geometry.line.Segment.midpoint

        Examples
        ========

        >>> from sympy.geometry import Point, Triangle
        >>> p1, p2, p3 = Point(0, 0), Point(1, 0), Point(0, 1)
        >>> t = Triangle(p1, p2, p3)
        >>> t.medians[p1]
        Segment(Point(0, 0), Point(1/2, 1/2))

        """
        s = self.sides
        v = self.vertices
        return {v[0]: Segment(v[0], s[1].midpoint),
                v[1]: Segment(v[1], s[2].midpoint),
                v[2]: Segment(v[2], s[0].midpoint)}

    @property
    def medial(self):
        """The medial triangle of the triangle.

        The triangle which is formed from the midpoints of the three sides.

        Returns
        =======

        medial : Triangle

        See Also
        ========

        sympy.geometry.line.Segment.midpoint

        Examples
        ========

        >>> from sympy.geometry import Point, Triangle
        >>> p1, p2, p3 = Point(0, 0), Point(1, 0), Point(0, 1)
        >>> t = Triangle(p1, p2, p3)
        >>> t.medial
        Triangle(Point(1/2, 0), Point(1/2, 1/2), Point(0, 1/2))

        """
        s = self.sides
        return Triangle(s[0].midpoint, s[1].midpoint, s[2].midpoint)

    @property
    def eulerline(self):
        """The Euler line of the triangle.

        The line which passes through circumcenter, orthocenter and the median
        of a triangle.

        Returns
        =======

        eulerline : Line if triangle is isosceles or right angled

        common center: Point if the triangle is equilateral in which case all
        those points coincide.

        Examples
        ========

        >>> from sympy.geometry import Point, Triangle
        >>> p1, p2, p3 = Point(0, 0), Point(1, 0), Point(0, 1)
        >>> t = Triangle(p1, p2, p3)
        >>> t.eulerline
        Line(Point(1/2, 1/2), Point(1/3, 1/3))

        """
        if not self.is_equilateral():
            points = list(set([self.circumcenter, self.orthocenter,
                            self.centroid]))
            return Line(points[0], points[1])
        else:
            return self.circumcenter

def rad(d):
    """Return the radian value for the given degrees (pi = 180 degrees)."""
    return d*pi/180


def deg(r):
    """Return the degree value for the given radians (pi = 180 degrees)."""
    return r/pi*180


def _slope(d):
    rv = tan(rad(d))
    return rv


def _asa(d1, l, d2):
    """Return triangle having side with length l on the x-axis."""
    xy = Line((0, 0), slope=_slope(d1)).intersection(
        Line((l, 0), slope=_slope(180 - d2)))[0]
    return Triangle((0, 0), (l, 0), xy)


def _sss(l1, l2, l3):
    """Return triangle having side of length l1 on the x-axis."""
    c1 = Circle((0, 0), l3)
    c2 = Circle((l1, 0), l2)
    inter = [a for a in c1.intersection(c2) if a.y.is_nonnegative]
    if not inter:
        return None
    pt = inter[0]
    return Triangle((0, 0), (l1, 0), pt)


def _sas(l1, d, l2):
    """Return triangle having side with length l2 on the x-axis."""
    p1 = Point(0, 0)
    p2 = Point(l2, 0)
    p3 = Point(cos(rad(d))*l1, sin(rad(d))*l1)
    return Triangle(p1, p2, p3)
