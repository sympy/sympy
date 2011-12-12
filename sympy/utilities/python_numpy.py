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

    def __new__(cls, data=None, _empty=False, shape=None):
        """
        Constructs a new Array.

        If _empty == True, then it constructs an empty array of the shape
            'shape'.
        """
        if _empty:
            obj = object.__new__(cls)
            if not isinstance(shape, (tuple, list)):
                shape = (shape,)
            s = 1
            for x in shape:
                s *= x
            obj._shape = shape
            obj._data = [0]*s
            return obj
        else:
            assert data is not None
        if hasattr(data, "__array__"):
            a = empty(data.shape)
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    a[i, j] = data[i, j]
            return a
        elif hasattr(data, "__iter__"):
            assert len(data) > 0
            if hasattr(data[0], "__iter__"):
                return vstack([Array(x) for x in data])
            else:
                obj = object.__new__(cls)
                obj._data = data
                obj._shape = (len(data),)
                return obj
        else:
            obj = object.__new__(cls)
            obj._data = [data]
            obj._shape = ()
            return obj

    @property
    def shape(self):
        return self._shape

    def ravel(self):
        return Array(self._data)

    def tolist(self):
        if len(self.shape) == 1:
            return self._data
        elif len(self.shape) == 2:
            d = []
            for n in range(self.shape[0]):
                d.append(self._data[n*self.shape[1]:(n+1)*self.shape[1]])
            return d
        else:
            raise NotImplementedError()

    def __array__(self):
        return self

    def __str__(self):
        return "array(%s)" % (str(self._data))

    def __repr__(self):
        return "array(%s)" % (str(self._data))

    def _op(self, i, op):
        a = empty(self._shape)
        if isinstance(i, Array):
            assert self._shape == i.shape
            a._data = [op(x, y) for x, y in zip(self._data, i._data)]
        elif hasattr(i, "__array__"):
            i = i.__array__()
            assert self._shape == i.shape
            a._data = [op(x, y) for x, y in zip(self._data, i.ravel().tolist())]
        else:
            a._data = [op(x, i) for x in self._data]
        return a

    def __mul__(self, i):
        return self._op(i, lambda x, o: x*o)

    def __rmul__(self, i):
        return self._op(i, lambda x, o: o*x)

    def __div__(self, i):
        return self._op(i, lambda x, o: x/o)

    def __rdiv__(self, i):
        return self._op(i, lambda x, o: o/x)

    def __add__(self, i):
        return self._op(i, lambda x, o: x+o)

    def __radd__(self, i):
        return self._op(i, lambda x, o: o+x)

    def __sub__(self, i):
        return self._op(i, lambda x, o: x-o)

    def __rsub__(self, i):
        return self._op(i, lambda x, o: o-x)

    def __pow__(self, i):
        return self._op(i, lambda x, o: x**i)

    def __rpow__(self, i):
        return self._op(i, lambda x, o: o**x)

    def __eq__(self, o):
        if isinstance(o, Array):
            if self.shape == o.shape:
                return Array([x == y for x, y in zip(self._data, o._data)])
            else:
                return False
        if len(self._data) == 1:
            return self._data[0] == o
        raise ValueError("Unsupported operation.")

    def __ne__(self, o):
        if isinstance(o, Array):
            if self.shape == o.shape:
                return Array([x != y for x, y in zip(self._data, o._data)])
            else:
                return True
        if len(self._data) == 1:
            return self._data[0] != o
        raise ValueError("Unsupported operation.")

    def __iter__(self):
        return self._data.__iter__()

    def __len__(self):
        return len(self._data)

    def _key2index(self, key):
        if not isinstance(key, (list, tuple)):
            key = [key]
        assert len(key) == len(self._shape)
        for d, k, s in zip(range(len(key)), key, self._shape):
            if not (k >= 0 and k < s):
                msg = "index (%d) out of range (0<=index<%d) in dimension %d" \
                        % (k, s, d)
                raise IndexError(msg)
        if len(key) == 3:
            return key[0] * self._shape[1] + key[1] * self._shape[2] + key[2]
        elif len(key) == 2:
            return key[0] * self._shape[1] + key[1]
        elif len(key) == 1:
            return key[0]
        else:
            raise NotImplementedError("Not implemented yet")

    def __getitem__(self, key):
        return self._data[self._key2index(key)]

    def __setitem__(self, key, val):
        self._data[self._key2index(key)] = val

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

def empty(shape, dtype=None):
    return Array(_empty=True, shape=shape)

def vstack(tup):
    h = len(tup)
    assert h > 0
    assert len(tup[0].shape) == 1
    w = tup[0].shape[0]
    a = empty((h, w))
    for i, x in enumerate(tup):
        assert len(x.shape) == 1
        assert x.shape[0] == w
        for j in range(w):
            a[i, j] = x[j]
    return a
