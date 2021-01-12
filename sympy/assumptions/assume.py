"""A module which implements predicates and assumption context."""

from contextlib import contextmanager
import inspect
from sympy.core.assumptions import ManagedProperties
from sympy.core.symbol import Str
from sympy.core.sympify import _sympify
from sympy.logic.boolalg import Boolean
from sympy.multipledispatch.dispatcher import (
    Dispatcher, MDNotImplementedError
)
from sympy.utilities.source import get_class


class AssumptionsContext(set):
    """
    Set containing default assumptions which are applied to the ``ask()``
    function.

    Explanation
    ===========

    This is used to represent global assumptions, but you can also use this
    class to create your own local assumptions contexts. It is basically a thin
    wrapper to Python's set, so see its documentation for advanced usage.

    Examples
    ========

    The default assumption context is ``global_assumptions``, which is initially empty:

    >>> from sympy import ask, Q
    >>> from sympy.assumptions import global_assumptions
    >>> global_assumptions
    AssumptionsContext()

    You can add default assumptions:

    >>> from sympy.abc import x
    >>> global_assumptions.add(Q.real(x))
    >>> global_assumptions
    AssumptionsContext({Q.real(x)})
    >>> ask(Q.real(x))
    True

    And remove them:

    >>> global_assumptions.remove(Q.real(x))
    >>> print(ask(Q.real(x)))
    None

    The ``clear()`` method removes every assumption:

    >>> global_assumptions.add(Q.positive(x))
    >>> global_assumptions
    AssumptionsContext({Q.positive(x)})
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
    remain unevaluated. To evaluate it, use the ``ask()`` function.

    Examples
    ========

    >>> from sympy import Q, ask
    >>> Q.integer(1)
    Q.integer(1)

    The ``function`` attribute returns the predicate, and the ``arguments``
    attribute returns the tuple of arguments.

    >>> type(Q.integer(1))
    <class 'sympy.assumptions.assume.AppliedPredicate'>
    >>> Q.integer(1).function
    Q.integer
    >>> Q.integer(1).arguments
    (1,)

    Applied predicates can be evaluated to a boolean value with ``ask``:

    >>> ask(Q.integer(1))
    True

    """
    __slots__ = ()

    is_Atom = True  # do not attempt to decompose this

    def __new__(cls, predicate, *args):
        if not isinstance(predicate, Predicate):
            raise TypeError("%s is not a Predicate." % predicate)
        args = map(_sympify, args)
        return super().__new__(cls, predicate, *args)

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
        # Will be deprecated
        args = self._args
        if len(args) == 2:
            # backwards compatibility
            return args[1]
        raise TypeError("'arg' property is allowed only for unary predicates.")

    @property
    def args(self):
        # Will be deprecated and return normal Basic.func
        return self._args[1:]

    @property
    def func(self):
        # Will be deprecated and return normal Basic.func
        return self._args[0]

    @property
    def function(self):
        """
        Return the predicate.
        """
        # Will be changed to self.args[0] after args overridding is removed
        return self._args[0]

    @property
    def arguments(self):
        """
        Return the arguments which are applied to the predicate.
        """
        # Will be changed to self.args[1:] after args overridding is removed
        return self._args[1:]

    def _eval_ask(self, assumptions):
        return self.function.eval(self.arguments, assumptions)

    @property
    def binary_symbols(self):
        from sympy.core.relational import Eq, Ne
        from .ask import Q
        if self.function == Q.is_true:
            i = self.arguments[0]
            if i.is_Boolean or i.is_Symbol or isinstance(i, (Eq, Ne)):
                return i.binary_symbols
        return set()


class PredicateMeta(ManagedProperties):
    """
    Metaclass for ``Predicate``

    If class attribute ``handler`` is not defined, assigns empty Dispatcher
    to it.
    """
    def __new__(cls, clsname, bases, dct):
        if "handler" not in dct:
            name = f"Ask{clsname.capitalize()}Handler"
            handler = Dispatcher(name, doc="Handler for key %s" % name)
            dct["handler"] = handler
        return super().__new__(cls, clsname, bases, dct)


class Predicate(Boolean, metaclass=PredicateMeta):
    """
    Base class for mathematical predicates. It also serves as a
    constructor for undefined predicate objects.

    Explanation
    ===========

    Predicate is a function that returns a boolean value [1].

    Predicate function is object, and it is instance of predicate class.
    When a predicate is applied to arguments, ``AppliedPredicate``
    instance is returned. This merely wraps the argument and remain
    unevaluated. To obtain the truth value of applied predicate, use the
    function ``ask``.

    Evaluation of predicate is done by multiple dispatching. You can
    register new handler to the predicate to support new types.

    Every predicate in SymPy can be accessed via the property of ``Q``.
    For example, ``Q.even`` returns the predicate which checks if the
    argument is even number.

    To define a predicate which can be evaluated, you must subclass this
    class, make an instance of it, and register it to ``Q``. After then,
    dispatch the handler by argument types.

    If you directly construct predicate using this class, you will get
    ``UndefinedPredicate`` which cannot be dispatched. This is useful
    when you are building boolean expressions which do not need to be
    evaluated.

    Examples
    ========

    Applying and evaluating to boolean value:

    >>> from sympy import Q, ask
    >>> from sympy.abc import x
    >>> ask(Q.prime(7))
    True

    You can define a new predicate by subclassing and dispatching. Here,
    we define a predicate for sexy primes [2] as an example.

    >>> from sympy import Predicate, Integer
    >>> class SexyPrimePredicate(Predicate):
    ...     name = "sexyprime"
    >>> Q.sexyprime = SexyPrimePredicate()
    >>> @Q.sexyprime.register(Integer, Integer)
    ... def _(int1, int2, assumptions):
    ...     args = sorted([int1, int2])
    ...     if not all(ask(Q.prime(a), assumptions) for a in args):
    ...         return False
    ...     return args[1] - args[0] == 6
    >>> ask(Q.sexyprime(5, 11))
    True

    Direct constructing returns ``UndefinedPredicate``, which can be
    applied but cannot be dispatched.

    >>> from sympy import Predicate, Integer
    >>> Q.P = Predicate("P")
    >>> type(Q.P)
    <class 'sympy.assumptions.assume.UndefinedPredicate'>
    >>> Q.P(1)
    Q.P(1)
    >>> Q.P.register(Integer)(lambda expr, assump: True)
    Traceback (most recent call last):
      ...
    TypeError: <class 'sympy.assumptions.assume.UndefinedPredicate'> cannot be dispatched.

    The tautological predicate ``Q.is_true`` can be used to wrap other objects:
    >>> Q.is_true(x > 1)
    Q.is_true(x > 1)

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Predicate_(mathematical_logic)
    .. [2] https://en.wikipedia.org/wiki/Sexy_prime

    """

    is_Atom = True

    def __new__(cls, *args, **kwargs):
        if cls is Predicate:
            return UndefinedPredicate(*args, **kwargs)
        obj = super().__new__(cls, *args)
        return obj

    @property
    def name(self):
        # May be overridden
        return type(self).__name__

    def register(self, *types, **kwargs):
        if self.handler is None:
            # condition for UndefinedPredicate
            raise TypeError("%s cannot be dispatched." % type(self))
        return self.handler.register(*types, **kwargs)

    def __call__(self, *args):
        return AppliedPredicate(self, *args)

    def eval(self, args, assumptions=True):
        """
        Evaluate ``self(*args)`` under the given assumptions.

        This uses only direct resolution methods, not logical inference.
        """
        types = tuple(type(a) for a in args)
        result = None
        for func in self.handler.dispatch_iter(*types):
            try:
                result = func(*args, assumptions)
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

    Examples
    ========

    >>> from sympy import Predicate, Q
    >>> Q.P = Predicate('P')
    >>> Q.P.func
    <class 'sympy.assumptions.assume.UndefinedPredicate'>
    >>> Q.P.name
    Str('P')

    """

    def __new__(cls, name, handlers=None):
        # "handlers" parameter supports old design
        if not isinstance(name, Str):
            name = Str(name)
        obj = super(Boolean, cls).__new__(cls, name)
        obj.handlers = handlers or []
        return obj

    @property
    def name(self):
        return self.args[0]

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

    def eval(self, args, assumptions=True):
        # Support for deprecated design
        # When old design is removed, this will always return None
        expr, = args
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
