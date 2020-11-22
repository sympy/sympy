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
from sympy.multipledispatch.conflict import ordering


class AskHandlerClass(Dispatcher):
    """
    Class for handler which evaluates the ``AppliedPredicate``.

    Explanation
    ===========

    This class dispatches various types, and return fuzzy boolean values
    when called.
    When the dispatched function return ``None``, the next function for
    the signature is queried until some other value is returned, or no
    more registered function is left.

    Parameters
    ==========

    name : str

    doc : str, optional

    base : tuple of AskHandlerClass, optional
        Base dispatchers, which inherites the dispatched functions.

    Examples
    ========

    Generate a base dispatcher ``MyHandler``.

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
    >>> MyHandler(x, assumptions=True)
    True
    >>> MyHandler(x+1, assumptions=True)
    False
    >>> print(MyHandler(2*x, assumptions=True)) # unregistred type
    None

    Generate new dispatcher ``MyHandler2``, which uses ``MyHandler`` as
    base dispatcher

    >>> MyHandler2 = AskHandlerClass('myhandler2', base=MyHandler)
    >>> @MyHandler2.register(Add)
    ... def _(expr, assumptions):
    ...     return True
    >>> MyHandler2(x, assumptions=True) # refer to base
    True
    >>> MyHandler2(x+1, assumptions=True)
    True

    """

    def __init__(self, name, doc=None, base=()):
        super().__init__(name, doc)
        if not type(base) is tuple:
            base = (base,)
        for b in reversed(base):
            self.funcs.update(b.funcs)
        self.ordering = ordering(self.funcs)
        self.base = base

    def __call__(self, *args, **kwargs):
        # If `None` is returned, search for next function
        # If all signature return `None`, return `None`
        types = tuple([type(arg) for arg in args])

        try:
            func = self._cache[types]
        except KeyError:
            func = self.dispatch(*types)
            self._cache[types] = func

        if not func:
            return None

        try:
            # MDNotImplementedError may be raised here...
            result = func(*args, **kwargs)
            if result is None:
                # ...or here
                raise MDNotImplementedError
            return result

        except MDNotImplementedError:
            # Search for next function which does not return None
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

    @staticmethod
    def AlwaysTrue(expr, assumptions):
        return True

    @staticmethod
    def AlwaysFalse(expr, assumptions):
        return False

    @staticmethod
    def AlwaysNone(expr, assumptions):
        return None

    def __getstate__(self):
        return {'name': self.name,
                'funcs': self.funcs,
                'base': self.base}

    def __setstate__(self, d):
        super().__setstate__(d)
        self.base = d['base']

CommonHandler = AskHandlerClass('CommonHandler')

CommonHandler.register(NaN)(CommonHandler.AlwaysFalse)


### AskCommutativeHandler ###

AskCommutativeHandler = AskHandlerClass(
    'AskCommutativeHandler',
    doc="""Handler for key 'commutative'""",
    base=CommonHandler
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
