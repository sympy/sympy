"""
Module to implement proxy object for symbolic manipulation of relations.
"""

from functools import partial

class SideProxy:
    """
    Proxy object to apply methods on relation.

    Parameters
    ==========

    rel : AppliedBinaryRelation

    side : "both", "lhs" or "rhs"

    """
    def __init__(self, rel, side):
        if side not in ("both", "lhs", "rhs"):
            raise ValueError("Only 'both', 'lhs' or 'rhs' is allowed for side.")
        self.rel = rel
        self.side = side

    def __repr__(self):
        return "<Side proxy for '%s' on %s side.>" % (self.rel, self.side)

    def _apply_func(self, func, *args, **kwargs):
        lhs, rhs = self.rel.arguments
        if self.side in ("both", "lhs"):
            lhs = func(lhs, *args, **kwargs)
        if self.side in ("both", "rhs"):
            rhs = func(rhs, *args, **kwargs)
        return self.rel.function(lhs, rhs)

    def _apply_attr(self, attrname):
        lhs, rhs = self.rel.arguments
        if self.side in ("both", "lhs"):
            lhs = getattr(lhs, attrname)
        if self.side in ("both", "rhs"):
            rhs = getattr(rhs, attrname)
        return self.rel.function(lhs, rhs)

    def _apply_method(self, methodname, *args, **kwargs):
        lhs, rhs = self.rel.arguments
        if self.side in ("both", "lhs"):
            lhs = getattr(lhs, methodname)(*args, **kwargs)
        if self.side in ("both", "rhs"):
            rhs = getattr(rhs, methodname)(*args, **kwargs)
        return self.rel.function(lhs, rhs)

    def __getattr__(self, attrname):
        lhs, rhs = self.rel.arguments

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

    def __call__(self, func, *args, **kwargs):
        lhs, rhs = self.rel.arguments
        if self.side in ("both", "lhs"):
            lhs = func(lhs, *args, **kwargs)
        if self.side in ("both", "rhs"):
            rhs = func(rhs, *args, **kwargs)
        return self.rel.function(lhs, rhs)
