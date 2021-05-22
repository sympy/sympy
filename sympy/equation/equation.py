from sympy.core import Basic
from sympy.core.sympify import _sympify

from .sideproxy import SideProxy


class SymbolicRelation(Basic):
    """
    Base class for all symbolic binary relations.

    Explanation
    ===========

    Unlike boolean relation, symbolic relation behaves as a container for
    the arguments. Its truth value is never evaluated, and all features are
    aimed for symbolic manipulation of the arguments.

    See the docstring of :obj:`~.Equation` for the examples.

    See Also
    ========

    sympy.core.relational.Relational : Boolean relation

    """

    def __new__(cls, lhs, rhs, **kwargs):
        lhs = _sympify(lhs)
        rhs = _sympify(rhs)
        return super().__new__(cls, lhs, rhs)

    @property
    def lhs(self):
        """The left-hand side of the relation."""
        return self.args[0]

    @property
    def rhs(self):
        """The right-hand side of the relation."""
        return self.args[1]

    @property
    def apply(self):
        """Proxy object to apply operation on both sides."""
        return SideProxy(self, "both")

    @property
    def applylhs(self):
        """Proxy object to apply operation on left hand sides."""
        return SideProxy(self, "lhs")

    @property
    def applyrhs(self):
        """Proxy object to apply operation on right hand sides."""
        return SideProxy(self, "rhs")


class Equation(SymbolicRelation):
    """
    Symbolic equation.

    Examples
    ========

    Symbolic equation is not reduced to boolean value.

    >>> from sympy import Eqn
    >>> Eqn(1, 1).simplify()
    Eqn(1, 1)

    Arguments can be manipulated by ``applylhs``, ``applyrhs``, or ``apply``
    properties. Note that this is purely structural manipulation, and does not
    guarantee mathematical correctness.

    >>> from sympy import cos, sin, gamma, trigsimp
    >>> from sympy.abc import x
    >>> eqn = Eqn(sin(x)**2 + cos(x)**2, gamma(x)/gamma(x-2))
    >>> eqn.apply.simplify()    # apply simplify method on both sides
    Eqn(1, (x - 2)*(x - 1))
    >>> eqn.applyrhs.simplify()     # apply simplify method on right hand side
    Eqn(sin(x)**2 + cos(x)**2, (x - 2)*(x - 1))
    >>> eqn.applylhs(trigsimp)      # apply trigsimp function on left hand side
    Eqn(1, gamma(x)/gamma(x - 2))

    Equation is a valid argument for ``subs()`` method.

    >>> (x**2).subs(Eqn(x, 2))
    4

    See Also
    ========

    sympy.core.relational.Equality : Boolean equality

    """
    rel_op = "=="

    @property
    def reversed(self):
        """
        Return the equation with sides reversed.

        Examples
        ========

        >>> from sympy import Eqn
        >>> from sympy.abc import x
        >>> Eqn(x, 1)
        Eqn(x, 1)
        >>> _.reversed
        Eqn(1, x)

        """
        return self.func(self.rhs, self.lhs)


Eqn = Equation
