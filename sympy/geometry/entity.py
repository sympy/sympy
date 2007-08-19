class GeometryEntity(object):
    """The base class for any geometrical entity."""

    def __new__(cls, *args, **kwargs):
        obj = object.__new__(cls)
        obj._args = tuple(args)
        return obj

    @staticmethod
    def do_intersection(e1, e2):
        """
        Determines the intersection between two geometrical entities. Returns
        a list of all of the intersections.
        """
        try:
            return e1._intersection(e2)
        except Exception:
            pass

        try:
            return e2._intersection(e1)
        except NotImplementedError:
            n1,n2 = type(e1).__name__, type(e2).__name__
            raise NotImplementedError("Unable to determine intersection between '%s' and '%s'" % (n1, n2))

    @staticmethod
    def _normalize_args(args):
        """Removes duplicates of arguments."""
        try:
            if not isinstance(args[0], GeometryEntity):
                args = args[0]
            return list(set(args))
        except:
            return args

    def _intersection(self, o):
        """
        Returns a list of all of the intersections of this entity and another
        entity.
        """
        raise NotImplementedError()

    def __ne__(self, o):
        return not self.__eq__(o)

    def __hash__(self):
        return hash(self._args)

    def __radd__(self, a):
        return a.__add__(self)

    def __rsub__(self, a):
        return a.__sub__(self)

    def __rmul__(self, a):
        return a.__mul__(self)

    def __rdiv__(self, a):
        return a.__div__(self)

    def __repr__(self):
        return str(self)
