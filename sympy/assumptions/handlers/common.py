"""
This module defines base class for handlers, and some core handlers:
commutative and is_true.
"""

from sympy.core import Symbol, Basic, Number
from sympy.core.logic import _fuzzy_group
from sympy.core.numbers import NaN
from sympy.logic.boolalg import (
    conjuncts, BooleanTrue, BooleanFalse, Not, Or, And, Implies, Equivalent
)
from sympy.assumptions import Q, ask, AppliedPredicate
from sympy.multipledispatch.dispatcher import (
    Dispatcher, MDNotImplementedError
)


class AskHandlerClass(Dispatcher):
    """
    Base class that all Ask Handlers must inherit.

    This class dispatches various types, and return fuzzy boolean values
    when called.
    When the dispatched function return ``None``, the next function for
    the signature is queried until some other value is returned, or no
    more registered function is left.

    """

    def __call__(self, *args, **kwargs):
        # If `None` is returned, search for next function
        # If all signature return `None`, return `None`
        types = tuple([type(arg) for arg in args])

        try:
            func = self._cache[types]
        except KeyError:
            func = self.dispatch(*types)
            if not func:
                return None
            self._cache[types] = func

        try:
            # MDNotImplementedError may be raised here...
            result = func(*args, **kwargs)
            if result is None:
                # ...or here
                raise MDNotImplementedError
            return result

        except MDNotImplementedError:
            funcs = self.dispatch_iter(*types)
            next(funcs)  # burn first
            for func in funcs:
                try:
                    result = func(*args, **kwargs)
                    if result is None:
                        raise MDNotImplementedError
                    return result
                except MDNotImplementedError:
                    pass
            return None

    def copy(self, name=None, doc=None):
        """
        Create a new handler with new name and new document, with
        all dispatched functions preserved.

        Examples
        ========

        >>> from sympy import Symbol, Add
        >>> from sympy.assumptions.handlers import AskHandlerClass
        >>> from sympy.abc import x
        >>> MyHandler = AskHandlerClass('myhandler')
        >>> @MyHandler.register(Symbol)
        ... def _(expr, assumptions):
        ...     return True
        >>> @MyHandler.register(Add)
        ... def _(expr, assumptions):
        ...     return False
        >>> print(MyHandler(x, assumptions=True))
        True
        >>> print(MyHandler(x+1, assumptions=True))
        False
        >>> print(MyHandler(2*x, assumptions=True))
        None

        Generate ``MyHandler2`` from ``MyHandler``, and override the
        registered function for ``Add``.

        >>> MyHandler2 = MyHandler.copy('myhandler2')
        >>> @MyHandler2.register(Add)
        ... def _(expr, assumptions):
        ...     return True
        >>> print(MyHandler2(x, assumptions=True))
        True
        >>> print(MyHandler2(x+1, assumptions=True))
        True
        >>> print(MyHandler2(2*x, assumptions=True))
        None

        """
        if name is None:
            name = self.name
        new = self.__class__(name, doc)
        new.funcs.update(self.funcs)
        new._cache.update(self._cache)
        new.ordering.extend(self.ordering)
        return new

    @staticmethod
    def AlwaysTrue(expr, assumptions):
        return True

    @staticmethod
    def AlwaysFalse(expr, assumptions):
        return False

    @staticmethod
    def AlwaysNone(expr, assumptions):
        return None

CommonHandler = AskHandlerClass('CommonHandler')

CommonHandler.register(NaN)(CommonHandler.AlwaysFalse)


### AskCommutativeHandler ###

AskCommutativeHandler = CommonHandler.copy(
    'AskCommutativeHandler',
    doc="""Handler for key 'commutative'"""
)

for sig in (NaN, Number):
    AskCommutativeHandler.register(sig)(AskCommutativeHandler.AlwaysTrue)

@AskCommutativeHandler.register(Basic)
def _(expr, assumptions):
    for arg in expr.args:
        if not ask(Q.commutative(arg), assumptions):
            return False
    return True

@AskCommutativeHandler.register(Symbol)
def _(expr, assumptions):
    """Objects are expected to be commutative unless otherwise stated"""
    assumps = conjuncts(assumptions)
    if expr.is_commutative is not None:
        return expr.is_commutative and not ~Q.commutative(expr) in assumps
    if Q.commutative(expr) in assumps:
        return True
    elif ~Q.commutative(expr) in assumps:
        return False
    return True


### TautologicalHandler ###

TautologicalHandler = AskHandlerClass(
    'TautologicalHandler',
    doc = """Wrapper allowing to query the truth value of a boolean expression."""
)

TautologicalHandler.register(BooleanTrue)(TautologicalHandler.AlwaysTrue)
TautologicalHandler.register(BooleanFalse)(TautologicalHandler.AlwaysFalse)

@TautologicalHandler.register(bool)
def _(expr, assumptions):
    return expr

@TautologicalHandler.register(AppliedPredicate)
def _(expr, assumptions):
    return ask(expr, assumptions)

@TautologicalHandler.register(Not)
def _(expr, assumptions):
    value = ask(expr.args[0], assumptions=assumptions)
    if value in (True, False):
        return not value
    else:
        return None

@TautologicalHandler.register(Or)
def _(expr, assumptions):
    result = False
    for arg in expr.args:
        p = ask(arg, assumptions=assumptions)
        if p is True:
            return True
        if p is None:
            result = None
    return result

@TautologicalHandler.register(And)
def _(expr, assumptions):
    result = True
    for arg in expr.args:
        p = ask(arg, assumptions=assumptions)
        if p is False:
            return False
        if p is None:
            result = None
    return result

@TautologicalHandler.register(Implies)
def _(expr, assumptions):
    p, q = expr.args
    return ask(~p | q, assumptions=assumptions)

@TautologicalHandler.register(Equivalent)
def _(expr, assumptions):
    p, q = expr.args
    pt = ask(p, assumptions=assumptions)
    if pt is None:
        return None
    qt = ask(q, assumptions=assumptions)
    if qt is None:
        return None
    return pt == qt

#### Helper methods
def test_closed_group(expr, assumptions, key):
    """
    Test for membership in a group with respect
    to the current operation
    """
    return _fuzzy_group(
        (ask(key(a), assumptions) for a in expr.args), quick_exit=True)
