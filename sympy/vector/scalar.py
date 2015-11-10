from sympy.core import Expr, Symbol
from sympy.core.sympify import _sympify
from sympy.core.compatibility import u, range
from sympy.printing.pretty.stringpict import prettyForm


class BaseScalar(Expr):
    """
    A coordinate symbol/base scalar.

    Ideally, users should not instantiate this class.

    """

    def __new__(cls, name, index, system, pretty_str, latex_str):
        from sympy.vector.coordsysrect import CoordSysCartesian
        if isinstance(name, Symbol):
            name = name.name
        if isinstance(pretty_str, Symbol):
            pretty_str = pretty_str.name
        if isinstance(latex_str, Symbol):
            latex_str = latex_str.name

        index = _sympify(index)
        system = _sympify(system)
        obj = super(BaseScalar, cls).__new__(cls, Symbol(name), index, system, Symbol(pretty_str), Symbol(latex_str))
        if not isinstance(system, CoordSysCartesian):
            raise TypeError("system should be a CoordSysCartesian")
        if index not in range(0, 3):
            raise ValueError("Invalid index specified.")
        # The _id is used for equating purposes, and for hashing
        obj._id = (index, system)
        obj._name = name
        obj._pretty_form = u(pretty_str)
        obj._latex_form = latex_str
        obj._system = system

        return obj

    @property
    def free_symbols(self):
        return set([self])

    _diff_wrt = True

    def _latex(self, printer=None):
        return self._latex_form

    def _pretty(self, printer=None):
        return prettyForm(self._pretty_form)

    @property
    def system(self):
        return self._system

    def __str__(self, printer=None):
        return self._name

    __repr__ = __str__
    _sympystr = __str__
