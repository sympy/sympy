"""A module which implements predicates and assumption context."""

from sympy.core.cache import cacheit
from sympy.core.singleton import S
from sympy.core.sympify import _sympify
from sympy.logic.boolalg import Boolean
from sympy.utilities.source import get_class
from contextlib import contextmanager


class AssumptionsContext(set):
    """
    Set representing assumptions.

    Explanation
    ===========

    This is used to represent global assumptions, but you can also use this
    class to create your own local assumptions contexts. It is basically a thin
    wrapper to Python's set, so see its documentation for advanced usage.

    Examples
    ========

    >>> from sympy import Q
    >>> from sympy.assumptions.assume import global_assumptions
    >>> from sympy.abc import x
    >>> global_assumptions
    AssumptionsContext()

    Adding applied predicate to global assumptions:

    >>> global_assumptions.add(Q.real(x))
    >>> global_assumptions
    AssumptionsContext({Q.real(x)})

    Removing applied predicatae from global assumptions:

    >>> global_assumptions.remove(Q.real(x))
    >>> global_assumptions
    AssumptionsContext()

    Clearing the global assumptions:

    >>> global_assumptions.clear()

    """

    def add(self, *assumptions):
        """Add an assumption."""
        for a in assumptions:
            super().add(a)

    def _sympystr(self, printer):
        if not self:
            return "%s()" % self.__class__.__name__
        return "{}({})".format(self.__class__.__name__, printer._print_set(self))

global_assumptions = AssumptionsContext()


class AppliedPredicate(Boolean):
    """
    The class of expressions resulting from applying ``Predicate`` to
    its arguments. ``AppliedPredicate`` merely wraps its argument and
    remain unevaluated. To evaluate it, use ``ask`` function.

    Examples
    ========

    >>> from sympy import Q, ask
    >>> Q.integer(1)
    Q.integer(1)

    ``func`` attribute returns the predicate, and ``args`` attribute
    returns the argument.

    >>> type(Q.integer(1))
    <class 'sympy.assumptions.assume.AppliedPredicate'>
    >>> Q.integer(1).func
    Q.integer
    >>> Q.integer(1).args
    (1,)

    Applied predicate can be evaluated to boolean value.

    >>> ask(Q.integer(1))
    True

    """
    __slots__ = ()

    is_Atom = True  # do not attempt to decompose this

    def __new__(cls, predicate, arg):
        arg = _sympify(arg)
        return Boolean.__new__(cls, predicate, arg)

    def __hash__(self):
        return super().__hash__()

    @property
    def arg(self):
        """
        Return the expression used by this assumption.

        Examples
        ========

        >>> from sympy import Q, Symbol
        >>> x = Symbol('x')
        >>> a = Q.integer(x + 1)
        >>> a.arg
        x + 1

        """
        return self._args[1]

    @property
    def args(self):
        return self._args[1:]

    @property
    def func(self):
        return self._args[0]

    @cacheit
    def sort_key(self, order=None):
        return (self.class_key(), (2, (self.func.name, self.arg.sort_key())),
                S.One.sort_key(), S.One)

    def __eq__(self, other):
        if type(other) is AppliedPredicate:
            return self._args == other._args
        return False

    def _eval_ask(self, assumptions):
        return self.func.eval(self.arg, assumptions)

    @property
    def binary_symbols(self):
        from sympy.core.relational import Eq, Ne
        if self.func.name in ['is_true', 'is_false']:
            i = self.arg
            if i.is_Boolean or i.is_Symbol or isinstance(i, (Eq, Ne)):
                return i.binary_symbols
        return set()


class Predicate(Boolean):
    """
    A predicate is a function that returns a boolean value.

    Explanation
    ===========

    Predicate consists of the name and multipledispatch handler.

    >>> from sympy import Q
    >>> type(Q.prime)
    <class 'sympy.assumptions.assume.Predicate'>
    >>> Q.prime.name
    'prime'
    >>> Q.prime.handler
    <dispatched AskPrimeHandler>

    When applied to its argument, ``AppliedPredicate`` instance is
    returned. This merely wrap the argument and remain unevaluated:

    >>> expr = Q.prime(7)
    >>> type(expr)
    <class 'sympy.assumptions.assume.AppliedPredicate'>
    >>> expr.func
    Q.prime
    >>> expr.args
    (7,)

    To obtain the truth value of an expression containing predicates,
    use the function ``ask``:

    >>> from sympy import ask
    >>> ask(Q.prime(7))
    True

    The tautological predicate ``Q.is_true`` can be used to wrap other objects:
    >>> from sympy.abc import x
    >>> Q.is_true(x > 1)
    Q.is_true(x > 1)

    Currently, only unary predicate is supported. Polyadic predicate
    will be implemented in future.

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Predicate_(mathematical_logic)

    """
    __slots__ = ('name', '_handler')

    is_Atom = True

    def __new__(cls, name, handler=None):
        obj = Boolean.__new__(cls)
        obj.name = name
        obj._handler = handler
        return obj

    @property
    def handler(self):
        """
        Multipledispatch handler instance of *self*.

        Since the handlers in ``assumptions.handlers module`` cannot be
        imported due to circular import problem, we store the path in
        ``self._handler`` and get it in runtime by this property.
        """
        from .handlers import CommonHandler
        if self._handler is None:
            name = ''.join(['Ask', self.name.capitalize(), 'Handler'])
            _handler = CommonHandler.copy(name)
            self._handler = _handler
        return get_class(self._handler)

    def _hashable_content(self):
        return (self.name, self.handler)

    def __getnewargs__(self):
        return (self.name, self.handler)

    def __call__(self, expr):
        return AppliedPredicate(self, expr)

    @cacheit
    def sort_key(self, order=None):
        return self.class_key(), (1, (self.name,)), S.One.sort_key(), S.One

    def eval(self, expr, assumptions=True):
        """
        Evaluate self(expr) under the given assumptions.

        This uses only direct resolution methods, not logical inference.
        """
        return self.handler(expr, assumptions=assumptions)


@contextmanager
def assuming(*assumptions):
    """
    Context manager for assumptions.

    Examples
    ========

    >>> from sympy.assumptions import assuming, Q, ask
    >>> from sympy.abc import x, y
    >>> print(ask(Q.integer(x + y)))
    None
    >>> with assuming(Q.integer(x), Q.integer(y)):
    ...     print(ask(Q.integer(x + y)))
    True
    """
    old_global_assumptions = global_assumptions.copy()
    global_assumptions.update(assumptions)
    try:
        yield
    finally:
        global_assumptions.clear()
        global_assumptions.update(old_global_assumptions)
