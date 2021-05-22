"""
Module to implement proxy object for symbolic manipulation of relations.
"""

from functools import partial


class SideProxy:
    """
    Proxy object to apply methods on relation.

    This object passes the expression manipulation down to the arguments
    of applied relation object, and returns a new one. It is implemented
    to distinguish the operation on relation itself, and operation on
    the arguments.

    Applying attributes on proxy object calls the attributes of relation's
    arguments and assembles a new relation. Calling the proxy object and
    passing the function as first argument applies the function on
    relation's arguments and assembles a new one.

    This class is intended to be accessed by ``.apply``, ``.applylhs``,
    and ``.applyrhs`` properties of ``SymbolicRelation`` object. Do not
    construct this object directly.

    Parameters
    ==========

    rel : SymbolicRelation

    side : "both", "lhs" or "rhs"

    Examples
    ========

    ``SideProxy`` can be accessed via ``apply`` attribute of
    ``SymbolicRelation`` object.

    >>> from sympy import exp, sin, Function, Eqn
    >>> from sympy.abc import x
    >>> f = Function('f')
    >>> eqn = Eqn(f(x), sin(x))
    >>> eqn.apply
    <Side proxy for 'Eqn(f(x), sin(x))' on both side>
    >>> eqn.applyrhs
    <Side proxy for 'Eqn(f(x), sin(x))' on rhs side>

    Applying methods on ``SideProxy`` applies the operation down to the
    corresponding argument(s).

    >>> eqn.apply.diff(x)
    Eqn(Derivative(f(x), x), cos(x))
    >>> eqn.applyrhs.rewrite(exp)
    Eqn(f(x), -I*(exp(I*x) - exp(-I*x))/2)

    Passing a function or a class to ``SideProxy`` applies it down to the
    corresponding argument(s). Note that keyword arguments can be passed, but
    the function or the class must be unary.

    >>> Eqn(x, 2).apply(lambda x: x**2)
    Eqn(x**2, 4)
    >>> Eqn(x, 0).apply(sin, evaluate=False)
    Eqn(sin(x), sin(0))
    """
    def __init__(self, rel, side):
        if side not in ("both", "lhs", "rhs"):
            raise ValueError("Only 'both', 'lhs' or 'rhs' is allowed for side.")
        self.rel = rel
        self.side = side

    def __repr__(self):
        return "<Side proxy for '%s' on %s side>" % (self.rel, self.side)

    def _apply_func(self, func, *args, **kwargs):
        lhs, rhs = self.rel.args
        if self.side in ("both", "lhs"):
            lhs = func(lhs, *args, **kwargs)
        if self.side in ("both", "rhs"):
            rhs = func(rhs, *args, **kwargs)
        return self.rel.func(lhs, rhs)

    def _apply_attr(self, attrname):
        lhs, rhs = self.rel.args
        if self.side in ("both", "lhs"):
            lhs = getattr(lhs, attrname)
        if self.side in ("both", "rhs"):
            rhs = getattr(rhs, attrname)
        return self.rel.func(lhs, rhs)

    def _apply_method(self, methodname, *args, **kwargs):
        lhs, rhs = self.rel.args
        if self.side in ("both", "lhs"):
            lhs = getattr(lhs, methodname)(*args, **kwargs)
        if self.side in ("both", "rhs"):
            rhs = getattr(rhs, methodname)(*args, **kwargs)
        return self.rel.func(lhs, rhs)

    def __getattr__(self, attrname):
        lhs, rhs = self.rel.args

        attrs = []
        if self.side in ("both", "lhs"):
            attrs.append(getattr(lhs, attrname))
        if self.side in ("both", "rhs"):
            attrs.append(getattr(rhs, attrname))

        if not all(callable(attr) for attr in attrs):
            return partial(self._apply_attr, attrname)
        elif all(callable(attr) for attr in attrs):
            return partial(self._apply_method, attrname)
        else:
            raise AttributeError("Inconsistent methods are called.")

    def __call__(self, func, **kwargs):
        lhs, rhs = self.rel.args
        if self.side in ("both", "lhs"):
            lhs = func(lhs, **kwargs)
        if self.side in ("both", "rhs"):
            rhs = func(rhs, **kwargs)
        return self.rel.func(lhs, rhs)
