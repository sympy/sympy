from sympy.core import Basic, S, C, sympify, oo, pi
from sympy.simplify import simplify
from sympy.geometry.exceptions import GeometryError
from entity import GeometryEntity
from point import Point
from ellipse import Circle
from line import Line, Segment, Ray


class Polygon(GeometryEntity):
    """
    A simple polygon in space. Can be constructed from a sequence or list
    of points.

    Notes:
    ======
        - Polygons are treated as closed paths rather than 2D areas so
          some calculations can be be negative or positive (e.g., area)
          based on the orientation of the points.
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
        Returns the area of the polygon.

        Notes:
        ======
            - The area calculation can be positive or negative based on
              the orientation of the points.
        """
        area = 0
        for ind in xrange(-1, len(self.vertices)-1):
            pi = self.vertices[ind]
            pii = self.vertices[ind+1]
            area += pi[0]*pii[1]-pii[0]*pi[1]
        return simplify(area) / 2

    @property
    def angles(self):
        """
        Returns a dictionary of {point: angle} entries containing the
        measure of all the internal angles of this polygon formed at
        each vertex.

        Examples:
        ======
            >>> from sympy import Point, Polygon
            >>> p1,p2,p3,p4 = map(Point, [(0,0), (1,0), (5,1), (0,1)])
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
        """Returns the perimeter of the polygon."""
        p = 0
        for ind in xrange(-1, len(self.vertices)-1):
            p += Point.distance(self.vertices[ind], self.vertices[ind+1])
        return simplify(p)

    @property
    def vertices(self):
        """The vertices that define the polygon."""
        return self[:]

    @property
    def centroid(self):
        """The centroid of the polygon."""
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
        """A list of the segments that form the sides of the polygon."""
        res = []
        for ind in xrange(0, len(self.vertices)-1):
            res.append( Segment(self.vertices[ind], self.vertices[ind+1]) )
        res.append( Segment(self.vertices[-1], self.vertices[0]) )
        return res

    def is_convex(self):
        """Returns True if this polygon is convex, False otherwise."""
        def tarea(a, b, c):
            return (b[0] - a[0])*(c[1] - a[1]) - (c[0] - a[0])*(b[1] - a[1])

        def isright(a, b, c):
            return bool(tarea(a, b, c) <= 0)

        # Determine orientation of points
        cw = isright(self.vertices[-2], self.vertices[-1], self.vertices[0])



        for i in xrange(0, len(self.vertices)):
            if cw ^ isright(self.vertices[i-2], self.vertices[i-1], self.vertices[i]):
                return False

        return True

    def intersection(self, o):
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
        e1_max_radius = 0
        e2_max_radius = 0
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
        e1_ymax = (0, -oo)
        e2_ymin = (0, oo)

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
        support_line = Line(Point(0, 0), Point(1, 0))

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
    A polygon which is regular. Such a polygon has all internal angles
    equal and all sides the same length.

    Usage:
    ======
        The polygon is constructed from a center point, a radius value,
        and the number of sides. The radius value is the distance from
        the center to a vertex.

    Example:
    ========
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
        points = []
        c, r, n = self
        v = 2*S.Pi/n
        for k in xrange(0, n):
            points.append( Point(c[0] + r*C.cos(k*v), c[1] + r*C.sin(k*v)) )
        return points

    @property
    def center(self):
        """
        Returns the center of the regular polygon (i.e., the center of the
        circumscribing circle).
        """
        return self.__getitem__(0)

    @property
    def radius(self):
        """
        Returns the radius of the regular polygon (i.e., the radius of the
        circumscribing circle).
        """
        return self.__getitem__(1)

    @property
    def apothem(self):
        """
        Returns the apothem/inradius of the regular polygon (i.e., the
        radius of the inscribed circle).
        """
        n = self.__getitem__(2)
        return self.radius * C.cos(S.Pi/n)

    @property
    def interior_angle(self):
        """Returns the measure of the interior angles."""
        n = self.__getitem__(2)
        return (n-2)*S.Pi/n

    @property
    def exterior_angle(self):
        """Returns the measure of the exterior angles."""
        n = self.__getitem__(2)
        return 2*S.Pi/n

    @property
    def circumcircle(self):
        """Returns a Circle instance describing the circumcircle."""
        return Circle(self.center, self.radius)

    @property
    def incircle(self):
        """Returns a Circle instance describing the inscribed circle."""
        return Circle(self.center, self.apothem)

    @property
    def angles(self):
        ret = {}
        ang = self.interior_angle
        for v in self.vertices:
            ret[v] = ang
        return ret

class Triangle(Polygon):
    """Any 3-sided polygon."""
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
        return self[:]

    def is_similar(t1, t2):
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
        """Returns True if the triangle is equilateral, False otherwise."""
        s = self.sides
        return bool(s[0].length == s[1].length) and bool(s[1].length == s[2].length)

    def is_right(self):
        """Returns True if the triangle is right-angled, False otherwise."""
        s = self.sides
        return Segment.is_perpendicular(s[0], s[1]) or \
               Segment.is_perpendicular(s[1], s[2]) or \
               Segment.is_perpendicular(s[0], s[2])

    @property
    def altitudes(self):
        """
        The altitudes of the triangle in a dictionary where the key
        is the vertex and the value is the altitude at that point.

        Example:
        ========
            >>> from sympy.geometry import Point, Triangle

            >>> p1,p2,p3 = Point(0,0), Point(1,0), Point(0,1)
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
        """The orthocenter of the triangle."""
        a = self.altitudes
        return GeometryEntity.intersect(a[1], a[2])[0]

    @property
    def circumcenter(self):
        """The circumcenter of the triangle."""
        a,b,c = [x.perpendicular_bisector() for x in self.sides]
        return GeometryEntity.do_intersection(a, b)[0]

    @property
    def circumradius(self):
        """The radius of the circumcircle of the triangle."""
        return Point.distance(self.circumcenter, self.vertices[0])

    @property
    def circumcircle(self):
        """The circumcircle of the triangle."""
        return Circle(self.circumcenter, self.circumradius)

    @property
    def bisectors(self):
        """
        The angle bisectors of the triangle in a dictionary where the
        key is the vertex and the value is the bisector at that point.

        Example:
        ========
            >>> from sympy.geometry import Point, Triangle, Segment

            >>> p1,p2,p3 = Point(0,0), Point(1,0), Point(0,1)
            >>> t = Triangle(p1, p2, p3)

            >>> from sympy import sqrt
            >>> t.bisectors[p2] == Segment(Point(0, sqrt(2)-1), Point(1, 0))
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
        """The incenter of the triangle."""
        s = self.sides
        v = self.vertices
        A,B,C = v[0],v[1],v[2]
        a,b,c = s[1].length,s[2].length,s[0].length
        x = simplify( (a*A[0] + b*B[0] + c*C[0]) / (a+b+c) )
        y = simplify( (a*A[1] + b*B[1] + c*C[1]) / (a+b+c) )
        return Point(x, y)

    @property
    def inradius(self):
        """The inradius of the triangle."""
        return simplify(self.area / self.perimeter)

    @property
    def incircle(self):
        """The incircle of the triangle."""
        return Circle(self.incenter, self.inradius)

    @property
    def medians(self):
        """
        The medians of the triangle in a dictionary where the key
        is the vertex and the value is the median at that point.

        Example:
        ========
            >>> from sympy.geometry import Point, Triangle

            >>> p1,p2,p3 = Point(0,0), Point(1,0), Point(0,1)
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
        """
        The medial triangle of the triangle, the triangle which is formed
        from the midpoints of the three sides.
        """
        s = self.sides
        return Triangle(s[0].midpoint, s[1].midpoint, s[2].midpoint)

    #@property
    #def excircles(self):
    #    """Returns a list of the three excircles for this triangle."""
    #    pass
