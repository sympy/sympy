from sympy.core import Basic, S, C, sympify, oo, pi
from sympy.simplify import simplify
from sympy.geometry.exceptions import GeometryError
from entity import GeometryEntity
from point import Point
from ellipse import Circle
from line import Line, Segment, Ray


class Polygon(GeometryEntity):
    """A two-dimensional polygon.

    A simple polygon in space. Can be constructed from a sequence or list
    of points.

    Parameters
    ----------
    vertices : sequence of Points

    Attributes
    ----------
    area
    angles
    perimeter
    vertices
    centroid
    sides

    Raises
    ------
    GeometryError
        If all parameters are not Points, or there are less than three
        parameters.

    See Also
    --------
    Point

    Notes
    -----
    Polygons are treated as closed paths rather than 2D areas so
    some calculations can be be negative or positive (e.g., area)
    based on the orientation of the points.

    Examples
    --------
    >>> from sympy import Point, Polygon
    >>> p1, p2, p3, p4 = map(Point, [(0, 0), (1, 0), (5, 1), (0, 1)])
    >>> Polygon(p1, p2, p3, p4)
    Polygon(Point(0, 0), Point(1, 0), Point(5, 1), Point(0, 1))

    """

    def __new__(cls, *args, **kwargs):
        vertices = GeometryEntity.extract_entities(args, remove_duplicates=False)
        if len(vertices) < 3:
            raise GeometryError("Polygon.__new__ requires at least three points")

        for p in vertices:
            if not isinstance(p, Point):
                raise GeometryError("Polygon.__new__ requires points")

        if len(vertices) == 3:
            return Triangle(*vertices, **kwargs)
        return GeometryEntity.__new__(cls, *vertices, **kwargs)

    @property
    def area(self):
        """
        The area of the polygon.

        Notes
        -----
        The area calculation can be positive or negative based on the
        orientation of the points.

        Examples
        --------
        >>> from sympy import Point, Polygon
        >>> p1, p2, p3, p4 = map(Point, [(0, 0), (1, 0), (5, 1), (0, 1)])
        >>> poly = Polygon(p1, p2, p3, p4)
        >>> poly.area
        3

        """
        area = 0
        for ind in xrange(-1, len(self.vertices) - 1):
            pi = self.vertices[ind]
            pii = self.vertices[ind + 1]
            area += pi[0]*pii[1] - pii[0]*pi[1]
        return simplify(area) / 2

    @property
    def angles(self):
        """The internal angle at each vertex.

        Returns
        -------
        angles : dict
            A dictionary where each key is a vertex and each value is the
            internal angle at that vertex. The vertices are represented as
            Points.

        See Also
        --------
        Point

        Examples
        --------
        >>> from sympy import Point, Polygon
        >>> p1, p2, p3, p4 = map(Point, [(0, 0), (1, 0), (5, 1), (0, 1)])
        >>> poly = Polygon(p1, p2, p3, p4)
        >>> poly.angles[p1]
        pi/2
        >>> poly.angles[p2]
        acos(-4*17**(1/2)/17)

        """
        def tarea(a, b, c):
            return (b[0] - a[0])*(c[1] - a[1]) - (c[0] - a[0])*(b[1] - a[1])

        def isright(a, b, c):
            return bool(tarea(a, b, c) <= 0)

        # Determine orientation of points
        cw = isright(self.vertices[-1], self.vertices[0], self.vertices[1])

        ret = {}
        for i in xrange(0, len(self.vertices)):
            a,b,c = self.vertices[i-2], self.vertices[i-1], self.vertices[i]
            ang = Line.angle_between(Line(b, a), Line(b, c))
            if cw ^ isright(a, b, c):
                ret[b] = 2*S.Pi - ang
            else:
                ret[b] = ang
        return ret

    @property
    def perimeter(self):
        """The perimeter of the polygon.

        Returns
        -------
        perimeter : number or Basic instance

        Examples
        --------
        >>> from sympy import Point, Polygon
        >>> p1, p2, p3, p4 = map(Point, [(0, 0), (1, 0), (5, 1), (0, 1)])
        >>> poly = Polygon(p1, p2, p3, p4)
        >>> poly.perimeter
        7 + 17**(1/2)
        """
        p = 0
        for ind in xrange(-1, len(self.vertices) - 1):
            p += Point.distance(self.vertices[ind], self.vertices[ind + 1])
        return simplify(p)

    @property
    def vertices(self):
        """The vertices of the polygon.

        Returns
        -------
        vertices : tuple of Points

        See Also
        --------
        Point

        Examples
        --------
        >>> from sympy import Point, Polygon
        >>> p1, p2, p3, p4 = map(Point, [(0, 0), (1, 0), (5, 1), (0, 1)])
        >>> poly = Polygon(p1, p2, p3, p4)
        >>> poly.vertices
        (Point(0, 0), Point(1, 0), Point(5, 1), Point(0, 1))

        """
        return self[:]

    @property
    def centroid(self):
        """The centroid of the polygon.

        Returns
        -------
        centroid : Point

        See Also
        --------
        Point

        Examples
        --------
        >>> from sympy import Point, Polygon
        >>> p1, p2, p3, p4 = map(Point, [(0, 0), (1, 0), (5, 1), (0, 1)])
        >>> poly = Polygon(p1, p2, p3, p4)
        >>> poly.centroid
        Point(31/18, 11/18)

        """
        A = 1 / (6*self.area)
        cx,cy = 0,0
        for ind in xrange(-1, len(self.vertices)-1):
            pi = self.vertices[ind]
            pii = self.vertices[ind+1]
            v = pi[0]*pii[1]-pii[0]*pi[1]
            cx += v*(pi[0] + pii[0])
            cy += v*(pi[1] + pii[1])
        return Point(simplify(A*cx), simplify(A*cy))

    @property
    def sides(self):
        """The line segments that form the sides of the polygon.

        Returns
        -------
        sides : list of sides
            Each side is a Segment.

        See Also
        --------
        Point
        Segment

        Examples
        --------
        >>> from sympy import Point, Polygon
        >>> p1, p2, p3, p4 = map(Point, [(0, 0), (1, 0), (5, 1), (0, 1)])
        >>> poly = Polygon(p1, p2, p3, p4)
        >>> poly.sides
        [Segment(Point(0, 0), Point(1, 0)),
        Segment(Point(1, 0), Point(5, 1)),
        Segment(Point(0, 1), Point(5, 1)), Segment(Point(0, 0), Point(0, 1))]

        """
        res = []
        for ind in xrange(0, len(self.vertices) - 1):
            res.append(Segment(self.vertices[ind], self.vertices[ind+1]))
        res.append(Segment(self.vertices[-1], self.vertices[0]))
        return res

    def is_convex(self):
        """Is the polygon convex.

        A polygon is convex if all its interior angles are less than 180
        degrees.

        Returns
        -------
        is_convex : boolean
            True if this polygon is convex, False otherwise.

        Examples
        --------
        >>> from sympy import Point, Polygon
        >>> p1, p2, p3, p4 = map(Point, [(0, 0), (1, 0), (5, 1), (0, 1)])
        >>> poly = Polygon(p1, p2, p3, p4)
        >>> poly.is_convex()
        True

        """
        def tarea(a, b, c):
            return (b[0] - a[0])*(c[1] - a[1]) - (c[0] - a[0])*(b[1] - a[1])

        def isright(a, b, c):
            return bool(tarea(a, b, c) <= 0)

        # Determine orientation of points
        cw = isright(self.vertices[-2], self.vertices[-1], self.vertices[0])

        for i in xrange(0, len(self.vertices)):
            if cw ^ isright(self.vertices[i - 2], self.vertices[i - 1], self.vertices[i]):
                return False

        return True

    def intersection(self, o):
        """The intersection of two polygons.

        The intersection may be empty and can contain individual Points and
        complete Line Segments.

        Parameters
        ----------
        other: Polygon

        Returns
        -------
        intersection : list
            The list of Segments and Points

        Examples
        --------
        >>> from sympy import Point, Polygon
        >>> p1, p2, p3, p4 = map(Point, [(0, 0), (1, 0), (5, 1), (0, 1)])
        >>> poly1 = Polygon(p1, p2, p3, p4)
        >>> p5, p6, p7, p8 = map(Point, [(3, 2), (1, -1), (0, 2), (-2, 1)])
        >>> poly2 = Polygon(p5, p6, p7, p8)
        >>> poly1.intersection(poly2)
        [Point(2/3, 0), Point(9/5, 1/5), Point(7/3, 1), Point(1/3, 1)]

        """
        res = []
        for side in self.sides:
            inter = GeometryEntity.do_intersection(side, o)
            if inter is not None:
                res.extend(inter)
        return res

    def distance(self, o):
        if isinstance(o, Point):
            dist = oo
            for side in self.sides:
                current = side.distance(o)
                if current == 0:
                    return S(0)
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
        of the polygons as it is assumed only called by Polygon.distance
        which does such checks.

        Notes:
        ======
            - Prints a warning if the two polygons possibly intersect as the return
              value will not be valid in such a case. For a more through test of
              intersection use intersection().

        Example:
        ========
            >>> from sympy.geometry import Point, Polygon
            >>> square = Polygon(Point(0, 0), Point(0, 1), Point(1, 1), Point(1, 0))
            >>> triangle = Polygon(Point(1, 2), Point(2, 2), Point(2, 1))
            >>> square._do_poly_distance(triangle)
            2**(1/2)/2

        Description of method used:
        ===========================
        Method:
            http://cgm.cs.mcgill.ca/~orm/mind2p.html
        Uses rotating calipers:
            http://en.wikipedia.org/wiki/Rotating_calipers
        and antipodal points:
            http://en.wikipedia.org/wiki/Antipodal_point
        """
        e1 = self

        '''Tests for a possible intersection between the polygons and outputs a warning'''
        e1_center = e1.centroid
        e2_center = e2.centroid
        e1_max_radius = S(0)
        e2_max_radius = S(0)
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
            print("Warning: Polygons may intersect producing erroneous output")

        '''
        Find the upper rightmost vertex of e1 and the lowest leftmost vertex of e2
        '''
        e1_ymax = (S(0), -oo)
        e2_ymin = (S(0), oo)

        for vertex in e1.vertices:
            if vertex[1] > e1_ymax[1] or (vertex[1] == e1_ymax[1] and vertex[0] > e1_ymax[0]):
                e1_ymax = vertex
        for vertex in e2.vertices:
            if vertex[1] < e2_ymin[1] or (vertex[1] == e2_ymin[1] and vertex[0] < e2_ymin[0]):
                e2_ymin = vertex
        min_dist = Point.distance(e1_ymax, e2_ymin)

        '''
        Produce a dictionary with vertices of e1 as the keys and, for each vertex, the points
        to which the vertex is connected as its value. The same is then done for e2.
        '''
        e1_connections = {}
        e2_connections = {}

        for side in e1.sides:
            if e1_connections.has_key(side.p1):
                e1_connections[side.p1].append(side.p2)
            else:
                e1_connections[side.p1] = [side.p2]

            if e1_connections.has_key(side.p2):
                e1_connections[side.p2].append(side.p1)
            else:
                e1_connections[side.p2] = [side.p1]

        for side in e2.sides:
            if e2_connections.has_key(side.p1):
                e2_connections[side.p1].append(side.p2)
            else:
                e2_connections[side.p1] = [side.p2]

            if e2_connections.has_key(side.p2):
                e2_connections[side.p2].append(side.p1)
            else:
                e2_connections[side.p2] = [side.p1]

        e1_current = e1_ymax
        e2_current = e2_ymin
        support_line = Line(Point(S(0), S(0)), Point(S(1), S(0)))

        '''
        Determine which point in e1 and e2 will be selected after e2_ymin and e1_ymax,
        this information combined with the above produced dictionaries determines the
        path that will be taken around the polygons
        '''
        point1 = e1_connections[e1_ymax][0]
        point2 = e1_connections[e1_ymax][1]
        angle1 = support_line.angle_between(Line(e1_ymax, point1))
        angle2 = support_line.angle_between(Line(e1_ymax, point2))
        if angle1 < angle2: e1_next = point1
        elif angle2 < angle1: e1_next = point2
        elif Point.distance(e1_ymax, point1) > Point.distance(e1_ymax, point2):
            e1_next = point2
        else: e1_next = point1

        point1 = e2_connections[e2_ymin][0]
        point2 = e2_connections[e2_ymin][1]
        angle1 = support_line.angle_between(Line(e2_ymin, point1))
        angle2 = support_line.angle_between(Line(e2_ymin, point2))
        if angle1 > angle2: e2_next = point1
        elif angle2 > angle1: e2_next = point2
        elif Point.distance(e2_ymin, point1) > Point.distance(e2_ymin, point2):
            e2_next = point2
        else: e2_next = point1

        '''
        Loop which determins the distance between anti-podal pairs and updates the
        minimum distance accordingly. It repeats until it reaches the starting position.
        '''
        while True:
            e1_angle = support_line.angle_between(Line(e1_current, e1_next))
            e2_angle = pi - support_line.angle_between(Line(e2_current, e2_next))

            if e1_angle < e2_angle:
                support_line = Line(e1_current, e1_next)
                e1_segment = Segment(e1_current, e1_next)
                min_dist_current = e1_segment.distance(e2_current)

                if min_dist_current.evalf() < min_dist.evalf(): min_dist = min_dist_current

                if e1_connections[e1_next][0] != e1_current:
                    e1_current = e1_next
                    e1_next = e1_connections[e1_next][0]
                else:
                    e1_current = e1_next
                    e1_next = e1_connections[e1_next][1]
            elif e1_angle > e2_angle:
                support_line = Line(e2_next, e2_current)
                e2_segment = Segment(e2_current, e2_next)
                min_dist_current = e2_segment.distance(e1_current)

                if min_dist_current.evalf() < min_dist.evalf(): min_dist = min_dist_current

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
                if min_dist_current.evalf() < min_dist.evalf(): min_dist = min_dist_current

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
            if e1_current == e1_ymax and e2_current == e2_ymin: break
        return min_dist

    def __eq__(self, o):
        if not isinstance(o, Polygon):
            return False

        # Find indices of points that are the same as the first point
        # in the other polygon
        n1,n2 = len(self.vertices), len(o.vertices)
        start_indices = []
        for ind in xrange(0, n1):
            if self.vertices[ind] == o.vertices[0]:
                start_indices.append(ind)

        if len(start_indices) == 0:
            return False

        # Check vertices clockwise and counterclockwise for equality
        imax = max(n1, n2)
        for start_ind in start_indices:
            i = start_ind

            # Check to see what orientation we should check
            dir = 0
            if self.vertices[(i + 1) % n1] == o.vertices[1]:
                dir = 1
            elif self.vertices[(i - 1) % n1] == o.vertices[1]:
                dir = -1

            # If either point to the left or right if the first point
            # is value (i.e., dir is nonzero) then check in that direction
            if dir != 0:
                areEqual = True
                for ind in xrange(2, imax):
                    if self.vertices[(i + dir*ind) % n1] != o.vertices[ind % n2]:
                        areEqual = False
                        break
                if areEqual: return True

        return False

    def __hash__(self):
        return super(Polygon, self).__hash__()

    def __contains__(self, o):
        if isinstance(o, Polygon):
            return self == o
        elif isinstance(o, Segment):
            return o in self.sides
        elif isinstance(o, Point):
            if o in self.vertices:
                return True
            for side in self.sides:
                if o in side:
                    return True
            return False
        else:
            return False


class RegularPolygon(Polygon):
    """
    A regular polygon.

    Such a polygon has all internal angles equal and all sides the same length.

    Parameters
    ----------
    center : Point
    radius : number or Basic instance
        The distance from the center to a vertex
    n : int
        The number of sides

    Attributes
    ----------
    vertices
    center
    radius
    apothem
    interior_angle
    exterior_angle
    circumcircle
    incircle
    angles

    Raises
    ------
    GeometryError
        If the `center` is not a Point, or the `radius` is not a number or Basic
        instance, or the number of sides, `n`, is less than three.

    See Also
    --------
    Point

    Examples
    --------
    >>> from sympy.geometry import RegularPolygon, Point
    >>> RegularPolygon(Point(0, 0), 5, 5)
    RegularPolygon(Point(0, 0), 5, 5)

    """

    def __new__(self, c, r, n, **kwargs):
        r = sympify(r)
        if not isinstance(c, Point):
            raise GeometryError("RegularPolygon.__new__ requires c to be a Point instance")
        if not isinstance(r, Basic):
            raise GeometryError("RegularPolygon.__new__ requires r to be a number or Basic instance")
        if n < 3:
            raise GeometryError("RegularPolygon.__new__ requires n >= 3")

        obj = GeometryEntity.__new__(self, c, r, n, **kwargs)
        return obj

    @property
    def vertices(self):
        """The vertices of the regular polygon.

        Returns
        -------
        vertices : list
            Each vertex is a Point.

        See Also
        --------
        Point

        Examples
        --------
        >>> from sympy.geometry import RegularPolygon, Point
        >>> rp = RegularPolygon(Point(0, 0), 5, 4)
        >>> rp.vertices
        [Point(5, 0), Point(0, 5), Point(-5, 0), Point(0, -5)]

        """
        points = []
        c, r, n = self
        v = 2*S.Pi/n
        for k in xrange(0, n):
            points.append(Point(c[0] + r*C.cos(k*v), c[1] + r*C.sin(k*v)))
        return points

    @property
    def center(self):
        """The center of the regular polygon

        This is also the centre of the circumscribing circle.

        Returns
        -------
        center : Point

        See Also
        --------
        Point

        Examples
        --------
        >>> from sympy.geometry import RegularPolygon, Point
        >>> rp = RegularPolygon(Point(0, 0), 5, 4)
        >>> rp.center
        Point(0, 0)

        """
        return self.__getitem__(0)

    @property
    def radius(self):
        """Radius of the regular polygon

        This is also the radius of the circumscribing circle.

        Returns
        -------
        radius : number or instance of Basic

        Examples
        --------
        >>> from sympy import Symbol
        >>> from sympy.geometry import RegularPolygon, Point
        >>> radius = Symbol('r')
        >>> rp = RegularPolygon(Point(0, 0), radius, 4)
        >>> rp.radius
        r

        """
        return self.__getitem__(1)

    @property
    def apothem(self):
        """The inradius of the regular polygon.

        The apothem/inradius is the radius of the inscribed circle.

        Returns
        -------
        apothem : number or instance of Basic

        Examples
        --------
        >>> from sympy import Symbol
        >>> from sympy.geometry import RegularPolygon, Point
        >>> radius = Symbol('r')
        >>> rp = RegularPolygon(Point(0, 0), radius, 4)
        >>> rp.apothem
        r*2**(1/2)/2

        """
        n = self.__getitem__(2)
        return self.radius * C.cos(S.Pi/n)

    @property
    def interior_angle(self):
        """Measure of the interior angles.

        Returns
        -------
        interior_angle : number

        Examples
        --------
        >>> from sympy.geometry import RegularPolygon, Point
        >>> rp = RegularPolygon(Point(0, 0), 4, 8)
        >>> rp.interior_angle
        3*pi/4

        """
        n = self.__getitem__(2)
        return (n - 2)*S.Pi/n

    @property
    def exterior_angle(self):
        """Measure of the exterior angles.

        Returns
        -------
        exterior_angle : number

        Examples
        --------
        >>> from sympy.geometry import RegularPolygon, Point
        >>> rp = RegularPolygon(Point(0, 0), 4, 8)
        >>> rp.exterior_angle
        pi/4

        """
        n = self.__getitem__(2)
        return 2*S.Pi/n

    @property
    def circumcircle(self):
        """The circumcircle of the regular polygon.

        Returns
        -------
        circumcircle : Circle

        See Also
        --------
        Circle

        Examples
        --------
        >>> from sympy.geometry import RegularPolygon, Point
        >>> rp = RegularPolygon(Point(0, 0), 4, 8)
        >>> rp.circumcircle
        Circle(Point(0, 0), 4)

        """
        return Circle(self.center, self.radius)

    @property
    def incircle(self):
        """The incircle of the regular polygon.

        Returns
        -------
        incircle : Circle

        See Also
        --------
        Circle

        Examples
        --------
        >>> from sympy.geometry import RegularPolygon, Point
        >>> rp = RegularPolygon(Point(0, 0), 4, 8)
        >>> rp.incircle
        Circle(Point(0, 0), 4*cos(pi/8))

        """
        return Circle(self.center, self.apothem)

    @property
    def angles(self):
        ret = {}
        ang = self.interior_angle
        for v in self.vertices:
            ret[v] = ang
        return ret

class Triangle(Polygon):
    """
    A polygon with three vertices and three sides.

    Parameters
    ----------
    points : sequence of Points

    Attributes
    ----------
    vertices
    altitudes
    orthocenter
    circumcenter
    circumradius
    circumcircle
    bisectors
    inradius
    incircle
    medians
    medial

    Raises
    ------
    GeometryError
        If the number of vertices is not equal to three, or one of the vertices
        is not a Point.

    See Also
    --------
    Point

    Examples
    --------
    >>> from sympy.geometry import Triangle, Point
    >>> Triangle(Point(0, 0), Point(4, 0), Point(4, 3))
    Triangle(Point(0, 0), Point(4, 0), Point(4, 3))

    """

    def __new__(cls, *args, **kwargs):
        vertices = GeometryEntity.extract_entities(args, remove_duplicates=False)
        if len(vertices) != 3:
            raise GeometryError("Triangle.__new__ requires three points")

        for p in vertices:
            if not isinstance(p, Point):
                raise GeometryError("Triangle.__new__ requires three points")

        return GeometryEntity.__new__(cls, *vertices, **kwargs)

    @property
    def vertices(self):
        """The triangle's vertices

        Returns
        -------
        vertices : tuple
            Each element in the tuple is a Point

        See Also
        --------
        Point

        Examples
        --------
        >>> from sympy.geometry import Triangle, Point
        >>> t = Triangle(Point(0, 0), Point(4, 0), Point(4, 3))
        >>> t.vertices
        (Point(0, 0), Point(4, 0), Point(4, 3))

        """
        return self[:]

    def is_similar(t1, t2):
        """Is another triangle similar to this one.

        Two triangles are similar if one can be uniformly scaled to the other.

        Parameters
        ----------
        other: Triangle

        Returns
        -------
        is_similar : boolean

        Examples
        --------
        >>> from sympy.geometry import Triangle, Point
        >>> t1 = Triangle(Point(0, 0), Point(4, 0), Point(4, 3))
        >>> t2 = Triangle(Point(0, 0), Point(-4, 0), Point(-4, -3))
        >>> t1.is_similar(t2)
        True

        >>> t2 = Triangle(Point(0, 0), Point(-4, 0), Point(-4, -4))
        >>> t1.is_similar(t2)
        False

        """
        if not isinstance(t2, Polygon) or len(t2) != 3:
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
        """Is the triangle equilateral

        Returns
        -------
        is_equilateral : boolean

        Examples
        --------
        >>> from sympy.geometry import Triangle, Point
        >>> t1 = Triangle(Point(0, 0), Point(4, 0), Point(4, 3))
        >>> t1.is_equilateral()
        False

        >>> from sympy import sqrt
        >>> t2 = Triangle(Point(0, 0), Point(10, 0), Point(5, 5*sqrt(3)))
        >>> t2.is_equilateral()
        True

        """
        s = self.sides
        return bool(s[0].length == s[1].length) and bool(s[1].length == s[2].length)

    def is_right(self):
        """Is the triangle right-angled.

        Returns
        -------
        is_right : boolean

        Examples
        --------
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

        An altitude of a triangle is a straight line through a vertex and
        perpendicular to the opposite side.

        Returns
        -------
        altitudes : dict
            The dictionary consists of keys which are vertices and values
            which are Segments.

        See Also
        --------
        Point
        Segment

        Examples
        --------
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

        Returns
        -------
        orthocenter : Point

        See Also
        --------
        Point

        Examples
        --------
        >>> from sympy.geometry import Point, Triangle
        >>> p1, p2, p3 = Point(0, 0), Point(1, 0), Point(0, 1)
        >>> t = Triangle(p1, p2, p3)

        """
        a = self.altitudes
        return GeometryEntity.intersect(a[1], a[2])[0]

    @property
    def circumcenter(self):
        """The circumcenter of the triangle

        The circumcenter is the center of the circumcircle.

        Returns
        -------
        circumcenter : Point

        See Also
        --------
        Point

        Examples
        --------
        >>> from sympy.geometry import Point, Triangle
        >>> p1, p2, p3 = Point(0, 0), Point(1, 0), Point(0, 1)
        >>> t = Triangle(p1, p2, p3)
        >>> t.circumcenter
        Point(1/2, 1/2)

        """
        a,b,c = [x.perpendicular_bisector() for x in self.sides]
        return GeometryEntity.do_intersection(a, b)[0]

    @property
    def circumradius(self):
        """The radius of the circumcircle of the triangle.

        Returns
        -------
        circumradius : number of Basic instance

        Examples
        --------
        >>> from sympy import Symbol
        >>> from sympy.geometry import Point, Triangle
        >>> a = Symbol('a')
        >>> p1, p2, p3 = Point(0, 0), Point(1, 0), Point(0, a)
        >>> t = Triangle(p1, p2, p3)
        >>> t.circumradius
        (1/4 + a**2/4)**(1/2)

        """
        return Point.distance(self.circumcenter, self.vertices[0])

    @property
    def circumcircle(self):
        """The circle which passes through the three vertices of the triangle.

        Returns
        -------
        circumcircle : Circle

        See Also
        --------
        Circle

        Examples
        --------
        >>> from sympy.geometry import Point, Triangle
        >>> p1, p2, p3 = Point(0, 0), Point(1, 0), Point(0, 1)
        >>> t = Triangle(p1, p2, p3)
        >>> t.circumcircle
        Circle(Point(1/2, 1/2), 2**(1/2)/2)

        """
        return Circle(self.circumcenter, self.circumradius)

    @property
    def bisectors(self):
        """The angle bisectors of the triangle.

        An angle bisector of a triangle is a straight line through a vertex
        which cuts the corresponding angle in half.

        Returns
        -------
        bisectors : dict
            Each key is a vertex (Point) and each value is the corresponding
            bisector (Segment).

        See Also
        --------
        Point
        Segment

        Examples
        --------
        >>> from sympy.geometry import Point, Triangle, Segment
        >>> p1, p2, p3 = Point(0, 0), Point(1, 0), Point(0, 1)
        >>> t = Triangle(p1, p2, p3)
        >>> from sympy import sqrt
        >>> t.bisectors[p2] == Segment(Point(0, sqrt(2) - 1), Point(1, 0))
        True

        """
        s = self.sides
        v = self.vertices
        c = self.incenter
        l1 = Segment(v[0], GeometryEntity.do_intersection(Line(v[0], c), s[1])[0])
        l2 = Segment(v[1], GeometryEntity.do_intersection(Line(v[1], c), s[2])[0])
        l3 = Segment(v[2], GeometryEntity.do_intersection(Line(v[2], c), s[0])[0])
        return {v[0]: l1, v[1]: l2, v[2]: l3}

    @property
    def incenter(self):
        """The center of the incircle.

        The incircle is the circle which lies inside the triangle and touches
        all three sides.

        Returns
        -------
        incenter : Point

        See Also
        --------
        incircle
        Point

        Examples
        --------
        >>> from sympy.geometry import Point, Triangle
        >>> p1, p2, p3 = Point(0, 0), Point(1, 0), Point(0, 1)
        >>> t = Triangle(p1, p2, p3)
        >>> t.incenter
        Point(1 - 2**(1/2)/2, 1 - 2**(1/2)/2)

        """
        s = self.sides
        v = self.vertices
        A,B,C = v[0],v[1],v[2]
        a,b,c = s[1].length,s[2].length,s[0].length
        x = simplify((a*A[0] + b*B[0] + c*C[0]) / (a+b+c))
        y = simplify((a*A[1] + b*B[1] + c*C[1]) / (a+b+c))
        return Point(x, y)

    @property
    def inradius(self):
        """The radius of the incircle.

        Returns
        -------
        inradius : number of Basic instance

        See Also
        --------
        incircle

        Examples
        --------
        >>> from sympy import Symbol
        >>> from sympy.geometry import Point, Triangle
        >>> a = Symbol('a')
        >>> p1, p2, p3 = Point(0, 0), Point(a, 0), Point(0, a)
        >>> t = Triangle(p1, p2, p3)
        >>> t.inradius
        (4*a**2*(a**2)**(1/2) - 2*2**(1/2)*a**2*(a**2)**(1/2))/(8*a**2)

        """
        return simplify(self.area / self.perimeter)

    @property
    def incircle(self):
        """The incircle of the triangle.

        The incircle is the circle which lies inside the triangle and touches
        all three sides.

        Returns
        -------
        incircle : Circle

        See Also
        --------
        Circle

        Examples
        --------
        >>> from sympy.geometry import Point, Triangle
        >>> p1, p2, p3 = Point(0, 0), Point(2, 0), Point(0, 2)
        >>> t = Triangle(p1, p2, p3)
        >>> t.incircle
        Circle(Point(2 - 2**(1/2), 2 - 2**(1/2)), 1 - 2**(1/2)/2)

        """
        return Circle(self.incenter, self.inradius)

    @property
    def medians(self):
        """The medians of the triangle.

        A median of a triangle is a straight line through a vertex and the
        midpoint of the opposite side, and divides the triangle into two
        equal areas.

        Returns
        -------
        medians : dict
            Each key is a vertex (Point) and each value is the median (Segment)
            at that point.

        See Also
        --------
        Point
        Segment

        Examples
        --------
        >>> from sympy.geometry import Point, Triangle
        >>> p1, p2, p3 = Point(0, 0), Point(1, 0), Point(0, 1)
        >>> t = Triangle(p1, p2, p3)
        >>> t.medians[p1]
        Segment(Point(0, 0), Point(1/2, 1/2))

        """
        s = self.sides
        v = self.vertices
        return {v[0]: Segment(s[1].midpoint, v[0]),
                v[1]: Segment(s[2].midpoint, v[1]),
                v[2]: Segment(s[0].midpoint, v[2])}

    @property
    def medial(self):
        """The medial triangle of the triangle.

        The triangle which is formed from the midpoints of the three sides.

        Returns
        -------
        medial : Triangle

        Examples
        --------
        >>> from sympy.geometry import Point, Triangle
        >>> p1, p2, p3 = Point(0, 0), Point(1, 0), Point(0, 1)
        >>> t = Triangle(p1, p2, p3)
        >>> t.medial
        Triangle(Point(1/2, 0), Point(1/2, 1/2), Point(0, 1/2))

        """
        s = self.sides
        return Triangle(s[0].midpoint, s[1].midpoint, s[2].midpoint)

    #@property
    #def excircles(self):
    #    """Returns a list of the three excircles for this triangle."""
    #    pass
