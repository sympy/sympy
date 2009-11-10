
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
    if len(entities) <= 1: return []

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
        ret = {}
        # not order preserving
        map(ret.__setitem__, a, [])
        return ret.keys()

    p = args[0]
    if isinstance(p, Point):
        p = uniquify(args)

    if len(p) is 1:
        return p[0]
    else:
        if len(p) is 2:
            return Segment(p[0],p[1])



    def orientation(p,q,r):
        '''Return positive if p-q-r are clockwise, neg if ccw, zero if colinear.'''
        return (q[1]-p[1])*(r[0]-p[0]) - (q[0]-p[0])*(r[1]-p[1])


    '''scan to find upper and lower convex hulls of a set of 2d points.'''
    U = []
    L = []
    #print "[graham] points:",Points
    p.sort()
    for p_i in p:
        while len(U) > 1 and orientation(U[-2],U[-1],p_i) <= 0: U.pop()
        while len(L) > 1 and orientation(L[-2],L[-1],p_i) >= 0: L.pop()
        U.append(p_i)
        L.append(p_i)        #print "[graham] result:", U,L
    U.reverse()
    convexHull = tuple(L+U[1:-1])

    #convexHull = uniquify(convexHull)
    #print "U:",U.reverse()
    #print "L:",L
    #print "ch(",len(convexHull),"):",convexHull
    if len(convexHull) is 2:
        return Segment(convexHull[0],convexHull[1])
    return Polygon(convexHull)





def are_similar(e1, e2):
    """
    Returns True if e1 and e2 are similar (one can be uniformly scaled to
    the other) or False otherwise.

    Notes:
    ======
        - If the two objects are equal then they are always similar.
    """
    if e1 == e2: return True
    try:
        return e1.is_similar(e2)
    except AttributeError:
        try:
            return e2.is_similar(e1)
        except AttributeError:
            n1 = e1.__class__.__name__
            n2 = e2.__class__.__name__
            raise GeometryError("Cannot test similarity between %s and %s" % (n1, n2))
