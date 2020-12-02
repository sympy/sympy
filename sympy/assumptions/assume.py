"""A module which implements predicates and assumption context."""

from contextlib import contextmanager
import inspect
from sympy.core.cache import cacheit
from sympy.core.singleton import S
from sympy.core.symbol import Str
from sympy.core.sympify import _sympify
from sympy.logic.boolalg import Boolean
from sympy.multipledispatch.dispatcher import (
    Dispatcher, MDNotImplementedError
)
from sympy.utilities.source import get_class


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

    Default assumption context is ``global_assumptions``, which is empty
    by default.

    >>> from sympy import Q
    >>> from sympy.assumptions import global_assumptions
    >>> global_assumptions
    AssumptionsContext()

    You can add default assumption.

    >>> from sympy.abc import x
    >>> global_assumptions.add(Q.real(x))
    >>> global_assumptions
    AssumptionsContext({Q.real(x)})
    >>> global_assumptions.add(Q.positive(x))
    >>> global_assumptions
    AssumptionsContext({Q.positive(x), Q.real(x)})

    And you can remove it.

    >>> global_assumptions.remove(Q.real(x))
    >>> global_assumptions
    AssumptionsContext({Q.positive(x)})

    ``clear()`` method removes every assumptions.

    >>> global_assumptions.clear()
    >>> global_assumptions
    AssumptionsContext()

    See Also
    ========

    assuming

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
    the arguments. ``AppliedPredicate`` merely wraps its argument and
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
        return (self.class_key(),
                (2, (self.func.sort_key(), self.arg.sort_key())),
                S.One.sort_key(), S.One)

    def __eq__(self, other):
        if type(other) is AppliedPredicate:
            return self._args == other._args
        return False

    def __hash__(self):
        return super().__hash__()

    def _eval_ask(self, assumptions):
        return self.func.eval(self.arg, assumptions)

    @property
    def binary_symbols(self):
        from sympy.core.relational import Eq, Ne
        if self.func.name.name in ['is_true', 'is_false']:
            i = self.arg
            if i.is_Boolean or i.is_Symbol or isinstance(i, (Eq, Ne)):
                return i.binary_symbols
        return set()


class Predicate(Boolean):
    """
    A predicate is a function that returns a boolean value [1].

    Explanation
    ===========

    When a predicate is applied to arguments, ``AppliedPredicate``
    instance is returned. This merely wraps the argument and remain
    unevaluated. To obtain the truth value of applied predicate, use the
    function ``ask``.

    Every predicate in SymPy can be accessed via the property of ``Q``.
    For example, ``Q.even`` returns the predicate which checks if the
    argument is even number.

    Currently, only unary predicate is supported. Polyadic predicate
    will be implemented in future.

    Examples
    ========

    Structure of predicate:

    >>> from sympy import Q, Predicate
    >>> isinstance(Q.prime, Predicate)
    True
    >>> Q.prime.name
    Str('prime')

    Applying and evaluating to boolean value:

    >>> from sympy import ask
    >>> expr = Q.prime(7)
    >>> ask(expr)
    True

    The tautological predicate ``Q.is_true`` can be used to wrap other objects:
    >>> from sympy.abc import x
    >>> Q.is_true(x > 1)
    Q.is_true(x > 1)

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Predicate_(mathematical_logic)

    """

    is_Atom = True
    _handler = None

    def __new__(cls, name, handlers=None):
        if cls is Predicate:
            return UndefinedPredicate(name, handlers)
        if not isinstance(name, Str):
            name = Str(name)
        obj = super().__new__(cls, name)
        return obj

    @classmethod
    def get_handler(cls):
        if cls._handler is None:
            name = ''.join(["Ask", cls.__name__.capitalize(), "Handler"])
            handler = Dispatcher(name, doc="Handler for key %s" % name)
            cls._handler = handler
        return cls._handler

    @property
    def handler(self):
        return self.get_handler()

    @property
    def name(self):
        return self.args[0]

    def __call__(self, expr):
        return AppliedPredicate(self, expr)

    def eval(self, expr, assumptions=True):
        """
        Evaluate self(expr) under the given assumptions.

        This uses only direct resolution methods, not logical inference.
        """
        result = None
        for func in self.handler.dispatch_iter(type(expr)):
            try:
                result = func(expr, assumptions)
            except MDNotImplementedError:
                continue
            else:
                if result is not None:
                    return result
        return result


class UndefinedPredicate(Predicate):
    """
    Predicate without handler.

    Explanation
    ===========

    This predicate is generated by using ``Predicate`` directly for
    construction. It does not have a handler, and evaluating this with
    arguments is done by SAT solver.

    """

    def __new__(cls, name, handlers=None):
        obj = super().__new__(cls, name, handlers)
        # support old design
        obj.handlers = handlers or []
        return obj

    @property
    def handler(self):
        return None

    def _hashable_content(self):
        return (self.name,)

    def __getnewargs__(self):
        return (self.name,)

    def __call__(self, expr):
        return AppliedPredicate(self, expr)

    def add_handler(self, handler):
        # Will be deprecated
        self.handlers.append(handler)

    def remove_handler(self, handler):
        # Will be deprecated
        self.handlers.remove(handler)

    def eval(self, expr, assumptions=True):
        # Support for deprecated design
        # When old design is removed, this will always return None
        res, _res = None, None
        mro = inspect.getmro(type(expr))
        for handler in self.handlers:
            cls = get_class(handler)
            for subclass in mro:
                eval_ = getattr(cls, subclass.__name__, None)
                if eval_ is None:
                    continue
                res = eval_(expr, assumptions)
                # Do not stop if value returned is None
                # Try to check for higher classes
                if res is None:
                    continue
                if _res is None:
                    _res = res
                elif res is None:
                    # since first resolutor was conclusive, we keep that value
                    res = _res
                else:
                    # only check consistency if both resolutors have concluded
                    if _res != res:
                        raise ValueError('incompatible resolutors')
                break
        return res


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
