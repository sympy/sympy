from sympy.core import Expr, S
from sympy.core.basic import Basic
from sympy.core.symbol import Str
from sympy.core.sympify import _sympify
from sympy.printing.pretty.stringpict import prettyForm
from sympy.printing.precedence import PRECEDENCE


class BaseScalar(Expr):
    """
    A coordinate symbol/base scalar.

    Ideally, users should not instantiate this class.

    Unicode pretty forms in Python 2 should use the `u` prefix.

    """

    def __new__(cls, index, system, pretty_str=None, latex_str=None):
        from sympy.vector.coordsysrect import CoordSys3D
        if pretty_str is None:
            pretty_str = "x{}".format(index)
        if latex_str is None:
            latex_str = "x_{}".format(index)
        if not isinstance(pretty_str, Basic):
            pretty_str = Str(pretty_str)
        if not isinstance(latex_str, Basic):
            latex_str = Str(latex_str)

        index = _sympify(index)
        system = _sympify(system)
        obj = super().__new__(cls, index, system, pretty_str, latex_str)
        if not isinstance(system, CoordSys3D):
            raise TypeError("system should be a CoordSys3D")
        if index not in range(0, 3):
            raise ValueError("Invalid index specified.")
        # The _id is used for equating purposes, and for hashing
        obj._id = (index, system)
        obj._name = obj.name = system._name + '.' + system._variable_names[index]
        obj._pretty_form = '' + pretty_str.name
        obj._latex_form = latex_str.name
        obj._system = system

        return obj

    is_commutative = True
    is_symbol = True

    @property
    def free_symbols(self):
        return {self}

    _diff_wrt = True

    def _eval_derivative(self, s):
        if self == s:
            return S.One
        return S.Zero

    def _latex(self, printer=None):
        return self._latex_form

    def _pretty(self, printer=None):
        return prettyForm(self._pretty_form)

    precedence = PRECEDENCE['Atom']

    @property
    def system(self):
        return self._system

    def _sympystr(self, printer):
        return self._name
