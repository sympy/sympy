from sympy.core.symbol import Symbol


class BaseScalar(Symbol):
    """
    A coordinate symbol/base scalar.

    Ideally, users should not instantiate this class.

    """

    def __new__(cls, name, index, system):
        from sympy.vector.coordsysrect import CoordSysRect
        obj = super(BaseScalar, cls).__new__(cls, name)
        if not isinstance(system, CoordSysRect):
            raise TypeError("system should be a CoordSysRect")
        if index not in range(0, 3):
            raise ValueError("Invalid index specified.")
        #The _id is used for equating purposes, and for hashing
        obj._id = (index, system)
        obj._name = name
        obj._system = system

        return obj

    @property
    def system(self):
        return self._system

    def __eq__(self, other):
        #Check if the other object is a BaseScalar of same index
        #and coordinate system
        if isinstance(other, BaseScalar):
            if other._id == self._id:
                return True
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return self._id.__hash__()

    def __str__(self, printer=None):
        return self._name

    __repr__ = __str__
    _sympystr = __str__
