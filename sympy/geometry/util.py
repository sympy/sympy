
def intersection(*entities):
    """
    Finds the intersection between a list GeometryEntity instances. Returns a
    list of all the intersections, Will raise a NotImplementedError exception
    if unable to calculate the intersection.

    Examples:
    =========
        >>> from sympy.geometry import *
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
        [Point(1 - 1/5*5**(1/2), 1 + (2/5)*5**(1/2)), Point(1 + (1/5)*5**(1/2), 1 - 2/5*5**(1/2))]

    Notes:
    ======
        - The intersection of any geometrical entity with itself should return
          a list with one item: the entity in question.
        - An intersection requires two or more entities. If only a single
          entity is given then one will receive an empty intersection list.
        - It is possible for intersection() to miss intersections that one
          knows exists because the required quantities were not fully
          simplified internally.
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
        >>> from sympy.geometry import Point
        >>> points = [ Point(x) for x in [(1,1), (1,2), (3,1), (-5,2), (15,4)] ]
        >>> convex_hull(points)
        Polygon(Point(3, 1), Point(15, 4), Point(-5, 2), Point(1, 1))

    Description of method used:
    ===========================
        See http://en.wikipedia.org/wiki/Graham_scan.
    """
    from point import Point
    from line import Segment
    from polygon import Polygon

    p = args[0]
    if isinstance(p, Point):
        p = args

    # Basic checks
    if len(p) == 1:
        return p[0]
    elif len(p) == 2:
        return Segment(p[0], p[1])

    # Find lowest+rightmost point
    m = 0
    for i in xrange(1, len(p)):
        if (p[i][1] < p[m][1]) or ((p[i][1] == p[m][1]) and (p[i][0] > p[m][0])):
            m = i
    p[0], p[m] = p[m], p[0]

    def tarea(a, b, c):
        return (b[0] - a[0])*(c[1] - a[1]) - (c[0] - a[0])*(b[1] - a[1])

    # Radial sort of points with respect to p[0] (our pivot)
    destroy = {}
    p0 = p[0]
    def pcompare(p1, p2):
        a = tarea(p0, p1, p2)
        if a > 0:
            return -1
        elif a < 0:
            return 1
        else:
            x = abs(p1[0] - p0[0]) - abs(p2[0] - p0[0])
            y = abs(p1[1] - p0[1]) - abs(p2[1] - p0[1])
            if (x < 0) or (y < 0):
                destroy[p1] = True
                return -1
            elif (x > 0) or (y > 0):
                destroy[p2] = True
                return 1
            else:
                destroy[p1] = True
                return 0
    p = p[1:]
    p.sort(pcompare)
    p.insert(0, p0)

    # Destroy points as found by sorting
    for i in xrange(len(p)-1, -1, -1):
        if p[i] in destroy:
            del p[i]

    # Graham scan
    def isleft(a, b, c):
        return (tarea(a, b, c) > 0)

    top = [p[0], p[1]]
    i = 2
    while i < len(p):
        p1 = top[-2]
        p2 = top[-1]
        if isleft(p1, p2, p[i]):
            top.append(p[i])
            i += 1
        else:
            top.pop()
    return Polygon(top)

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