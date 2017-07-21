from sympy.core import Expr, Symbol, S
from sympy.core.sympify import _sympify
from sympy.core.compatibility import range
from sympy.printing.pretty.stringpict import prettyForm
from sympy.printing.precedence import PRECEDENCE


class BaseScalar(Expr):
    """
    A coordinate symbol/base scalar.

    Ideally, users should not instantiate this class.

    Unicode pretty forms in Python 2 should use the `u` prefix.

    """

    def __new__(cls, name, index, system, pretty_str, latex_str):
        from sympy.vector.coordsysrect import CoordSys3D
        if isinstance(name, Symbol):
            name = name.name
        if isinstance(pretty_str, Symbol):
            pretty_str = pretty_str.name
        if isinstance(latex_str, Symbol):
            latex_str = latex_str.name

        index = _sympify(index)
        system = _sympify(system)
        obj = super(BaseScalar, cls).__new__(cls, Symbol(name), index, system,
                                             Symbol(pretty_str),
                                             Symbol(latex_str))
        if not isinstance(system, CoordSys3D):
            raise TypeError("system should be a CoordSys3D")
        if index not in range(0, 3):
            raise ValueError("Invalid index specified.")
        # The _id is used for equating purposes, and for hashing
        obj._id = (index, system)
        obj._name = obj.name = name
        obj._pretty_form = u'' + pretty_str
        obj._latex_form = latex_str
        obj._system = system

        return obj

    is_commutative = True

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

    def __str__(self, printer=None):
        return self._name

    __repr__ = __str__
    _sympystr = __str__
