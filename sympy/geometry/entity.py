# How entities are ordered; used by __cmp__ in GeometryEntity
ordering_of_classes = [
    "Point",
    "Segment",
    "Ray",
    "Line",
    "Triangle",
    "RegularPolygon",
    "Polygon",
    "Circle",
    "Ellipse",
    "Curve"
]

class GeometryEntity(tuple):
    """The base class for any geometrical entity."""

    def __new__(cls, *args, **kwargs):
        return tuple.__new__(cls, args)

    def __getnewargs__(self):
        return tuple(self)

    @staticmethod
    def do_intersection(e1, e2):
        """
        Determines the intersection between two geometrical entities. Returns
        a list of all of the intersections.
        """
        try:
            return e1.intersection(e2)
        except Exception:
            pass

        try:
            return e2.intersection(e1)
        except NotImplementedError:
            n1,n2 = type(e1).__name__, type(e2).__name__
            raise NotImplementedError("Unable to determine intersection between '%s' and '%s'" % (n1, n2))

    def is_similar(self, other):
        """
        Return True if self and other are similar. Two entities are similar
        if a uniform scaling (enlarging or shrinking) of one of the entities
        will allow one to obtain the other.

        Notes:
        ======
            - This method is not intended to be used directly but rather
              through the are_similar() method found in util.py.
            - An entity is not required to implement this method.
            - If two different types of entities can be similar, it is only
              required that one of them be able to determine this.
        """
        raise NotImplementedError()

    def intersection(self, o):
        """
        Returns a list of all of the intersections of this entity and another
        entity.

        Notes:
        ======
            - This method is not intended to be used directly but rather
              through the intersection() method found in util.py.
            - An entity is not required to implement this method.
            - If two different types of entities can intersect, it is only
              required that one of them be able to determine this.
        """
        raise NotImplementedError()


    @staticmethod
    def extract_entities(args, remove_duplicates=True):
        """
        Takes a set of arguments and extracts all of the GeometryEntity
        instances (recursively). Returns a tuple of all of the instances
        found.

        Notes:
        ======
            - Duplicates of entities are removed if the remove_duplicates
              argument is set to True, otherwise duplicates remain.
              The default is True.
            - Anything that is not a GeometryEntity instance is simply
              ignored.
            - Ordering of arguments is always maintained. If duplicates
              are removed then the entry with the lowest index is kept.
        """
        ret = list()
        for arg in args:
            if isinstance(arg, GeometryEntity):
                ret.append(arg)
            elif isinstance(arg, (list, tuple, set)):
                ret.extend(GeometryEntity.extract_entities(arg))
        if remove_duplicates:
            temp = set(ret)
            ind, n = 0, len(ret)
            for counter in xrange(0, n):
                x = ret[ind]
                if x in temp:
                    temp.remove(x)
                    ind += 1
                else:
                    del ret[ind]
        return tuple(ret)

    def __ne__(self, o):
        return not self.__eq__(o)

    def __radd__(self, a):
        return a.__add__(self)

    def __rsub__(self, a):
        return a.__sub__(self)

    def __rmul__(self, a):
        return a.__mul__(self)

    def __rdiv__(self, a):
        return a.__div__(self)

    def __str__(self):
        from sympy.printing import sstr
        return type(self).__name__ + sstr (tuple(self))

    def __repr__(self):
        from sympy.printing import srepr
        return type(self).__name__ + srepr(tuple(self))

    def __cmp__(self, other):
        n1 = self.__class__.__name__
        n2 = other.__class__.__name__
        c = cmp(n1, n2)
        if not c: return 0

        i1 = -1
        for cls in self.__class__.__mro__:
            try:
                i1 = ordering_of_classes.index(cls.__name__)
                break
            except ValueError:
                i1 = -1
        if i1 == -1: return c

        i2 = -1
        for cls in other.__class__.__mro__:
            try:
                i2 = ordering_of_classes.index(cls.__name__)
                break
            except ValueError:
                i2 = -1
        if i2 == -1: return c

        return cmp(i1, i2)
