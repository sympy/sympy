def intersection(*entities):
    """
    Finds the intersection between a list GeometryEntity instances. Returns a
    list of all the intersections, Will raise a NotImplementedError exception
    if unable to calculate the intersection.

    Examples:
    =========
        >>> from sympy.geometry import Point, Line, Circle, intersection
        >>> p1,p2,p3 = Point(0,0), Point(1,1), Point(-1, 5)
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
        [Point(1 - 5**(1/2)/5, 1 + 2*5**(1/2)/5), Point(1 + 5**(1/2)/5, 1 - 2*5**(1/2)/5)]

    Notes:
    ======
        - The intersection of any geometrical entity with itself should return
          a list with one item: the entity in question.
        - An intersection requires two or more entities. If only a single
          entity is given then one will receive an empty intersection list.
        - It is possible for intersection() to miss intersections that one
          knows exists because the required quantities were not fully
          simplified internally.
        - Reals should be converted to Rationals, e.g. Rational(str(real_num))
          or else failures due to floating point issues may result.
    """
    from entity import GeometryEntity

    entities = GeometryEntity.extract_entities(entities, False)
    if len(entities) <= 1:
        return []

    res = GeometryEntity.do_intersection(entities[0], entities[1])
    for entity in entities[2:]:
        newres = []
        for x in res:
            newres.extend( GeometryEntity.do_intersection(x, entity) )
        res = newres
    return res


def convex_hull(*args):
    """
    Returns a Polygon representing the convex hull of a set of 2D points.

    Notes:
    ======
        This can only be performed on a set of non-symbolic points.

    Example:
    ========
        >>> from sympy.geometry import Point, convex_hull
        >>> points = [ Point(x) for x in [(1,1), (1,2), (3,1), (-5,2), (15,4)] ]
        >>> convex_hull(points)
        Polygon(Point(-5, 2), Point(1, 1), Point(3, 1), Point(15, 4))

    Description of method used:
    ===========================
        See http://en.wikipedia.org/wiki/Graham_scan.
    """
    from point import Point
    from line import Segment
    from polygon import Polygon

    def uniquify(a):
        # not order preserving
        return list(set(a))

    p = args[0]
    if isinstance(p, Point):
        p = uniquify(args)

    if len(p) == 1:
        return p[0]
    elif len(p) == 2:
        return Segment(p[0], p[1])

    def orientation(p, q, r):
        '''Return positive if p-q-r are clockwise, neg if ccw, zero if colinear.'''
        return (q[1]-p[1])*(r[0]-p[0]) - (q[0]-p[0])*(r[1]-p[1])

    '''scan to find upper and lower convex hulls of a set of 2d points.'''
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
    return Polygon(convexHull)


def are_similar(e1, e2):
    """
    Returns True if e1 and e2 are similar (one can be uniformly scaled to
    the other) or False otherwise.

    Notes:
    ======
        - If the two objects are equal then they are always similar.
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


def poly_poly_distance(e1, e2):
    """
    Calculates the least distance between the exteriors of two
    convex polygons e1 and e2.

    Notes:
    ======
        - Prints a warning if the two polygons possibly intersect as the return
          value will not be valid in such a case. For a more through test of
          intersection use intersection().

    Example:
    ========
        >>> from sympy.geometry import Point, Polygon, poly_poly_distance
        >>> square = Polygon(Point(0,0), Point(0,1), Point(1,1), Point(1,0))
        >>> triangle = Polygon(Point(1,2), Point(2,2), Point(2,1))
        >>> poly_poly_distance(square, triangle)
        2**(1/2)/2

    Description of method used:
    ===========================
    http://cgm.cs.mcgill.ca/~orm/mind2p.html
    """
    from sympy import oo, pi
    from sympy.geometry import Polygon, Ray, Point, Line, Segment, GeometryError

    '''Calculates the distance between a point and a line segment'''
    def pt_seg_dist(pt, seg):
        seg_vector = Point(seg.p2[0] - seg.p1[0], seg.p2[1] - seg.p1[1])
        pt_vector = Point(pt[0] - seg.p1[0], pt[1] - seg.p1[1])
        t = (seg_vector[0]*pt_vector[0] + seg_vector[1]*pt_vector[1])/seg.length**2
        if t>=1:
            distance = Point.distance(seg.p2, pt)
        elif t<=0:
            distance = Point.distance(seg.p1, pt)
        else:
            distance = Point.distance((seg.p1+Point(t*seg_vector[0], t*seg_vector[1])), pt)
        return distance

    '''Error checking'''
    if not (isinstance(e1, Polygon) and isinstance(e2, Polygon)):
        n1 = e1.__class__.__name__
        n2 = e2.__class__.__name__
        raise GeometryError("Cannot calculate distance between %s and %s" % (n1, n2))
    elif not (e1.is_convex() and e2.is_convex()):
        raise GeometryError("poly_poly_distance requires convex polygons as arguments")

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
    if center_dist<=e1_max_radius+e2_max_radius:
        print("Warning: Polygons may intersect producing erroneous output")

    '''
    Find the upper rightmost vertex of e1 and the lowest leftmost vertex of e2
    '''
    e1_ymax = (0,-oo)
    e2_ymin = (0,oo)

    for vertex in e1.vertices:
        if vertex[1] > e1_ymax[1] or (vertex[1]==e1_ymax[1] and vertex[0]>e1_ymax[0]):
            e1_ymax = vertex
    for vertex in e2.vertices:
        if vertex[1] < e2_ymin[1] or (vertex[1]==e2_ymin[1] and vertex[0]<e2_ymin[0]):
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
    support_line = Line(Point(0,0),Point(1,0))

    '''
    Determine which point in e1 and e2 will be selected after e2_ymin and e1_ymax,
    this information combined with the above produced dictionaries determines the
    path that will be taken around the polygons
    '''
    point1 = e1_connections[e1_ymax][0]
    point2 = e1_connections[e1_ymax][1]
    angle1 = support_line.angle_between(Line(e1_ymax, point1))
    angle2 = support_line.angle_between(Line(e1_ymax, point2))
    if angle1<angle2: e1_next = point1
    elif angle2<angle1: e1_next = point2
    elif Point.distance(e1_ymax, point1) > Point.distance(e1_ymax, point2):
        e1_next = point2
    else: e1_next = point1

    point1 = e2_connections[e2_ymin][0]
    point2 = e2_connections[e2_ymin][1]
    angle1 = support_line.angle_between(Line(e2_ymin, point1))
    angle2 = support_line.angle_between(Line(e2_ymin, point2))
    if angle1>angle2: e2_next = point1
    elif angle2>angle1: e2_next = point2
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
            min_dist_current = pt_seg_dist(e2_current, e1_segment)

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
            min_dist_current = pt_seg_dist(e1_current, e2_segment)

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
            min1 = pt_seg_dist(e2_next, e1_segment)
            min2 = pt_seg_dist(e1_next, e2_segment)

            min_dist_current = min(min1,min2)
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
