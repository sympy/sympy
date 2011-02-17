def array(data):
    """
    Constructs an Array object.

    If NumPy is present, it uses numpy.array, otherwise it uses sympy's Array.

    Example::

    >>> from sympy import array
    >>> array([1, 2, 3])
    array([1, 2, 3])

    """
    try:
        import numpy
        numpy_ok = True
    except ImportError:
        numpy_ok = False

    if numpy_ok:
        return numpy.array(data)
    else:
        return Array(data)

class Array(object):
    """
    Implements a NumPy ndarray object in pure Python.

    This class should be 100% compatible with NumPy's ndarray, and you can use
    it in cases when NumPy is not available. The name of this class is
    deliberately chosen to be different than ndarray, to prevent confusion
    which implementation one is using, but otherwise it should behave exactly
    the same.

    Currently it only implements 1D arrays.

    Example:

    >>> from sympy.utilities.python_numpy import Array
    >>> a = Array([1, 2, 3])
    >>> a * 3
    array([3, 6, 9])

    """

    def __init__(self, data):
        try:
            self._data = list(data)
        except TypeError:
            self._data = [data]

    def __str__(self):
        return "array(%s)" % (str(self._data))

    def __repr__(self):
        return "array(%s)" % (str(self._data))

    def __pow__(self, i):
        return Array([x**i for x in self._data])

    def __mul__(self, i):
        return Array([x*i for x in self._data])

    def __rmul__(self, i):
        return Array([i*x for x in self._data])

    def __div__(self, i):
        return Array([x/i for x in self._data])

    def __rdiv__(self, i):
        return Array([i/x for x in self._data])

    def __add__(self, i):
        return Array([x+i for x in self._data])

    def __radd__(self, i):
        return Array([i+x for x in self._data])

    def __sub__(self, i):
        return Array([x-i for x in self._data])

    def __rsub__(self, i):
        return Array([i-x for x in self._data])

    def __pow__(self, i):
        return Array([x**i for x in self._data])

    def __rpow__(self, i):
        return Array([i**x for x in self._data])

    def __eq__(self, o):
        if isinstance(o, Array):
            return Array([x == y for x, y in zip(self._data, o._data)])
        if len(self._data) == 1:
            return self._data[0] == o
        raise ValueError("Unsupported operation.")

    def __ne__(self, o):
        if isinstance(o, Array):
            return Array([x != y for x, y in zip(self._data, o._data)])
        if len(self._data) == 1:
            return self._data[0] != o
        raise ValueError("Unsupported operation.")

    def __iter__(self):
        return self._data.__iter__()

    def __len__(self):
        return len(self._data)

    def __getitem__(self, key):
        return self._data[key]

    def all(self):
        """
        Returns True if all elements equal to True.

        Example::

            >>> from sympy import array
            >>> array([True, True, True]).all()
            True
            >>> array([True, False, True]).all()
            False
            >>> array([False, False, False]).any()
            False

        """
        for x in self._data:
            if not x:
                return False
        return True

    def any(self):
        """
        Returns True if any element equals to True.

        Example::

            >>> from sympy import array
            >>> array([True, True, True]).any()
            True
            >>> array([True, False, True]).any()
            True
            >>> array([False, False, False]).any()
            False

        """
        for x in self._data:
            if x:
                return True
        return False

    def __nonzero__(self):
        if len(self._data) == 1:
            return self._data[0] == True
        raise ValueError("ValueError: The truth value of an array with more than one element is ambiguous. Use a.any() or a.all()")
