# -*- coding:utf-8 -*-

"""
Definition of physical dimensions.

Unit systems will be constructed on top of these dimensions.
"""

from __future__ import division
from copy import copy
import numbers

from sympy.core.containers import Dict, Tuple
from sympy import Number


class Dimension(Dict):
    """
    This class represent the dimension of a physical quantities.

    The dimensions may have a name and a symbol. All other
    arguments are dimensional powers. They represent a characteristic of a
    quantity, giving an interpretation to it: for example (in classical
    mechanics) we know that time is different from temperature, and dimensions
    make this difference (but they do not provide any measure of these
    quantites).
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
                    if len(p) != 2:
                        raise ValueError("length of iterable has to be 2.")

                # construction with [("length", 1), ...]
                pairs.extend(arg)
            else:
            # error if the arg is not of previous types
                raise TypeError("positional arguments can only be: "
                                "dict, tuple, list.")

        pairs.extend(kwargs.items())

        # check validity of dimension key and power
        for pair in pairs:
            #if not isinstance(p[0], str):
            #    raise TypeError("key %s is not a string." % p[0])
            if not isinstance(pair[1], (numbers.Real, Number)):
                raise TypeError("power corresponding to %s is not a number."
                                % pair[0])

        # filter dimensions set to zero; this avoid the following odd result:
        # Dimension(length=1) == Dimension(length=1, mass=0) => False
        pairs = [pair for pair in pairs if pair[1] != 0]

        new = Dict.__new__(cls, *pairs)
        new.name = name
        new.symbol = symbol

        return new
