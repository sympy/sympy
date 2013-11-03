# -*- coding:utf-8 -*-

"""
Definition of physical dimensions.

Unit systems will be constructed on top of these dimensions.

Most of the examples below used MKS system and are presented from the computer
point of view: from a human point, adding length to time is not legal in MKS
but it is in natural system; for a computer in natural system there is no time
dimension (but a velocity dimension instead) so the question of adding time
to length has no meaning.
"""

from __future__ import division
from copy import copy
import numbers

from sympy.core.containers import Dict, Tuple
from sympy import Number, sympify


class Dimension(Dict):
    """
    This class represent the dimension of a physical quantities.

    The dimensions may have a name and a symbol. All other
    arguments are dimensional powers. They represent a characteristic of a
    quantity, giving an interpretation to it: for example (in classical
    mechanics) we know that time is different from temperature, and dimensions
    make this difference (but they do not provide any measure of these
    quantites).

        >>> from sympy.physics.unitsystems.dimensions import Dimension
        >>> length = Dimension(length=1)
        >>> length
        {length: 1}
        >>> time = Dimension(time=1)
    """

    def __new__(cls, *args, **kwargs):
        """
        Create a new dimension.

        Possibilities are (examples given with list/tuple work also with
        tuple/list):

            >>> from sympy.physics.unitsystems.dimensions import Dimension
            >>> Dimension(length=1)
            {length: 1}
            >>> Dimension({"length": 1})
            {length: 1}
            >>> Dimension([("length", 1), ("time", -1)])
            {length: 1, time: -1}
        """

        # before setting the dict, check if a name and/or a symbol are defined
        # if so, remove them from the dict
        name = kwargs.pop('name', None)
        symbol = kwargs.pop('symbol', None)

        # pairs of (dimension, power)
        pairs = []

        # add first items from args to the pairs
        for arg in args:
            # construction with {"length": 1}
            if isinstance(arg, dict):
                arg = copy(arg)
                pairs.extend(arg.items())
            elif isinstance(arg, (Tuple, tuple, list)):
                #TODO: add construction with ("length", 1); not trivial because
                #      e.g. [("length", 1), ("time", -1)] has also length = 2

                for p in arg:
                    #TODO: check that p is a tuple
                    if len(p) != 2:
                        raise ValueError("Length of iterable has to be 2; "
                                         "'%d' found" % len(p))

                # construction with [("length", 1), ...]
                pairs.extend(arg)
            else:
            # error if the arg is not of previous types
                raise TypeError("Positional arguments can only be: "
                                "dict, tuple, list; '%s' found" % type(arg))

        pairs.extend(kwargs.items())

        # check validity of dimension key and power
        for pair in pairs:
            #if not isinstance(p[0], str):
            #    raise TypeError("key %s is not a string." % p[0])
            if not isinstance(pair[1], (numbers.Real, Number)):
                raise TypeError("Power corresponding to '%s' is not a number"
                                % pair[0])

        # filter dimensions set to zero; this avoid the following odd result:
        # Dimension(length=1) == Dimension(length=1, mass=0) => False
        pairs = [pair for pair in pairs if pair[1] != 0]

        new = Dict.__new__(cls, *pairs)
        new.name = name
        new.symbol = symbol

        return new

    def __str__(self):
        """
        Display the string representation of the dimension.

        Usually one will always use a symbol to denote the dimension. If no
        symbol is defined then it uses the name or, if there is no name, the
        default dict representation.
        """

        if self.symbol is not None:
            return self.symbol
        elif self.name is not None:
            return self.name
        else:
            return repr(self)

    @property
    def is_dimensionless(self):
        """
        Check if the dimension object really has a dimension.

        A dimension should have at least one component with non-zero power.
        """

        for key in self:
            if self[key] != 0:
                return False
        else:
            return True

    def __neg__(self):
        return self

    def __add__(self, other):
        """
        Define the addition for Dimension.

        Addition of dimension has a sense only if the second object is the same
        dimension (we don't add length to time).
        """

        if not isinstance(other, Dimension):
            raise TypeError("Only dimension can be added; '%s' is not valid"
                            % type(other))
        elif isinstance(other, Dimension) and self != other:
            raise ValueError("Only dimension which are equal can be added; "
                             "'%s' and '%s' are different" % (self, other))

        return self

    def __sub__(self, other):
        # there is no notion of ordering (or magnitude) among dimension,
        # subtraction is equivalent to addition when the operation is legal
        return self + other

    def __pow__(self, other):
        #TODO: be sure that it works with rational numbers (e.g. when dealing
        #      with dimension under a fraction)

        other = sympify(other)
        if isinstance(other, (numbers.Real, Number)):
            return Dimension([(x, y*other) for x, y in self.items()])
        else:
            raise TypeError("Dimensions can be exponentiated only with "
                            "numbers; '%s' is not valid" % type(other))

    def __mul__(self, other):
        if not isinstance(other, Dimension):
            #TODO: improve to not raise error: 2*L could be a legal operation
            #      (the same comment apply for __div__)
            raise TypeError("Only dimension can be multiplied; '%s' is not "
                            "valid" % type(other))

        d = dict(self)
        for key in other:
            try:
                d[key] += other[key]
            except KeyError:
                d[key] = other[key]
        d = Dimension(d)

        # if all dimensions are zero, then return 1 so that there is no more
        # dimensions
        if d.is_dimensionless:
            return 1
        else:
            return d

    def __div__(self, other):
        if not isinstance(other, Dimension):
            raise TypeError("Only dimension can be divided; '%s' is not valid"
                            % type(other))

        d = dict(self)
        for key in other:
            try:
                d[key] -= other[key]
            except KeyError:
                d[key] = -other[key]
        d = Dimension(d)

        if d.is_dimensionless:
            return 1
        else:
            return d

    __truediv__ = __div__

    def __rdiv__(self, other):
        return other * pow(self, -1)

    __rtruediv__ = __rdiv__


class DimensionSystem(object):
    """
    DimensionSystem represents a coherent set of dimensions.

    In a system dimensions are of two types:

    - base dimensions;
    - derived dimensions: these are defined in terms of the base dimensions
      (for example velocity is defined from the division of length by time).
    """

    def __init__(self, base, dims=(), name="", descr=""):

        self.name = name
        self.descr = descr

        self._base_dims = self._sort_dims(base)
        self._dims = self._sort_dims(list(dims) + [d for d in base
                                                   if d not in dims])

    @staticmethod
    def _sort_dims(dims):
        """
        Sort dimensions given in argument using their str function.

        This function will ensure that we get always the same tuple for a given
        set of dimensions.
        """

        return tuple(sorted(dims, key=str))
