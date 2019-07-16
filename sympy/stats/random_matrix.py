from __future__ import print_function, division

from sympy import Basic, MatrixSymbol
from sympy.stats.rv import PSpace, _symbol_converter

class RandomMatrixPSpace(PSpace):
    """
    Represents probability space for
    random matrices. It contains the mechanics
    for handling the API calls for random matrices.
    """
    def __new__(cls, sym, model=None):
        sym = _symbol_converter(sym)
        return Basic.__new__(cls, sym, model)

    model = property(lambda self: self.args[1])
