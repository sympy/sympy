"""
Module to allow operation between binary relations by multipledispatch.
This is useful for symbolic manipulation of equation.
"""
from collections import defaultdict

from sympy.multipledispatch.dispatcher import (Dispatcher, str_signature,
    RaiseNotImplementedError)
from .binrel import AppliedBinaryRelation


class RelOpDispatcher:
    """
    Multiple dispatcher to define the operation between ``AppliedBinaryRelation``
    with different relations.

    Multipledispatcher registers the function over the class of arguments.
    Since ``AppliedBinaryRelation`` is best classified by its ``.function``, not
    class, we use this dispatcher to define the operation between these.

    Parameters
    ==========

    name : str
    doc : str, optional

    Examples
    ========

    We define a operation on two equalities, which adds the both sides
    and divide them by 2.

    >>> from sympy import Q
    >>> from sympy.relation import Equal
    >>> from sympy.relation.relop import RelOpDispatcher
    >>> from sympy.abc import x,y
    >>> my_relop = RelOpDispatcher("my_relop")
    >>> @my_relop.register(Equal, Equal)
    ... def _(eq1, eq2):
    ...     lhs = (eq1.lhs + eq2.lhs)/2
    ...     rhs = (eq1.rhs + eq2.rhs)/2
    ...     return Q.eq(lhs, rhs)
    >>> my_relop(Q.eq(x,1), Q.eq(y,3))
    x/2 + y/2 = 2
    """
    def __init__(self, name, doc=None):
        self.name = name
        self.doc = doc
        self._dispatcher = Dispatcher(name)

    def __repr__(self):
        return "<dispatched %s>" % self.name

    def register(self, typ1, typ2, **kwargs):
        """
        Register the binary handler for two binary relation classes.
        """
        return self._dispatcher.register(typ1, typ2, **kwargs)

    def __call__(self, arg1, arg2, **kwargs):
        if isinstance(arg1, AppliedBinaryRelation):
            typ1 = type(arg1.function)
        else:
            typ1 = type(arg1)
        if isinstance(arg2, AppliedBinaryRelation):
            typ2 = type(arg2.function)
        else:
            typ2 = type(arg2)
        func = self._dispatcher.dispatch(typ1, typ2)

        if func is None:
            raise NotImplementedError(
            "%s and %s are not dispatched on %s" % (typ1, typ2, self))

        return func(arg1, arg2, **kwargs)

    @property
    def __doc__(self):
        docs = [
            "Dispatcher for operation between relations : %s" % self.name,
        ]

        if self.doc:
            docs.append(self.doc)

        s = "Registered relation classes\n"
        s += '=' * len(s)
        docs.append(s)

        amb_sigs = []

        typ_sigs = defaultdict(list)
        for sigs in self._dispatcher.ordering[::-1]:
            key = self._dispatcher.funcs[sigs]
            typ_sigs[key].append(sigs)

        for func, sigs in typ_sigs.items():

            sigs_str = ', '.join('<%s>' % str_signature(sig) for sig in sigs)

            if isinstance(func, RaiseNotImplementedError):
                amb_sigs.append(sigs_str)
                continue

            s = 'Inputs: %s\n' % sigs_str
            s += '-' * len(s) + '\n'
            if func.__doc__:
                s += func.__doc__.strip()
            else:
                s += func.__name__
            docs.append(s)

        return '\n\n'.join(docs)


relop_add = RelOpDispatcher('add')

relop_mul = RelOpDispatcher('mul')

relop_pow = RelOpDispatcher('pow')
