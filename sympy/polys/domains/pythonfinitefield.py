"""Implementation of :class:`PythonFiniteField` class. """

from __future__ import print_function, division

from sympy.polys.domains.finitefield import FiniteField
from sympy.polys.domains.pythonintegerring import PythonIntegerRing

from sympy.utilities import public

@public
class PythonFiniteField(FiniteField):
    """Finite field based on Python's integers. """

    alias = 'FF_python'

    def __new__(cls, mod, symmetric=True):
        return super(PythonFiniteField, cls).__new__(
            cls, mod, PythonIntegerRing(), symmetric)

    def __getnewargs__(self):
        return self.args[0], self.args[2]
