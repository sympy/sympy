from sympy.core.basic import Basic, S, C, sympify
from sympy.simplify import simplify
from sympy.geometry.exceptions import GeometryError
from entity import GeometryEntity
from point import Point
from ellipse import Circle
from line import Line, Segment


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
            >>> from sympy.geometry import *
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
        """A list of the segments that form the sides fo the polygon."""
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
