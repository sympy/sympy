# doctests are disabled because of issue #1521
from sympy.logic.boolalg import Boolean, Not

class AssumptionsContext(set):
    """Set representing assumptions.

    This is used to represent global assumptions, but you can also use this
    class to create your own local assumptions contexts. It is basically a thin
    wrapper to Python's set, so see its documentation for advanced usage.

    Examples:
        >>> from sympy import global_assumptions, Assume, Q
        >>> global_assumptions
        AssumptionsContext()
        >>> from sympy.abc import x
        >>> global_assumptions.add(Assume(x, Q.real))
        >>> global_assumptions
        AssumptionsContext([Assume(x, Q.real)])
        >>> global_assumptions.remove(Assume(x, Q.real))
        >>> global_assumptions
        AssumptionsContext()
        >>> global_assumptions.clear()

    """

    def add(self, *assumptions):
        """Add an assumption."""
        for a in assumptions:
            assert isinstance(a, Assume), 'can only store instances of Assume'
            super(AssumptionsContext, self).add(a)

global_assumptions = AssumptionsContext()

class Assume(Boolean):
    """New-style assumptions.

    >>> from sympy import Assume, Q
    >>> from sympy.abc import x
    >>> Assume(x, Q.integer)
    Assume(x, Q.integer)
    >>> Assume(x, Q.integer, False)
    Not(Assume(x, Q.integer))
    >>> Assume( x > 1 )
    Assume(1 < x, Q.is_true)

    """
    def __new__(cls, expr, predicate=None, value=True):
        from sympy import Q
        if predicate is None:
            predicate = Q.is_true
        elif not isinstance(predicate, Predicate):
            key = str(predicate)
            try:
                predicate = getattr(Q, key)
            except AttributeError:
                predicate = Predicate(key)
        if value:
            return Boolean.__new__(cls, expr, predicate)
        else:
            return Not(Boolean.__new__(cls, expr, predicate))

    is_Atom = True # do not attempt to decompose this

    @property
    def expr(self):
        """
        Return the expression used by this assumption.

        Examples:
            >>> from sympy import Assume, Q
            >>> from sympy.abc import x
            >>> a = Assume(x+1, Q.integer)
            >>> a.expr
            1 + x

        """
        return self._args[0]

    @property
    def key(self):
        """
        Return the key used by this assumption.
        It is a string, e.g. 'integer', 'rational', etc.

        Examples:
            >>> from sympy import Assume, Q
            >>> from sympy.abc import x
            >>> a = Assume(x, Q.integer)
            >>> a.key
            Q.integer

        """
        return self._args[1]

    def __eq__(self, other):
        if type(other) == Assume:
            return self._args == other._args
        return False

    def __hash__(self):
        return super(Assume, self).__hash__()

def eliminate_assume(expr, symbol=None):
    """
    Convert an expression with assumptions to an equivalent with all assumptions
    replaced by symbols.

    Assume(x, integer=True) --> integer
    Assume(x, integer=False) --> ~integer

    Examples:
        >>> from sympy.assumptions.assume import eliminate_assume
        >>> from sympy import Assume, Q
        >>> from sympy.abc import x
        >>> eliminate_assume(Assume(x, Q.positive))
        Q.positive
        >>> eliminate_assume(Assume(x, Q.positive, False))
        Not(Q.positive)

    """
    if expr.func is Assume:
        if symbol is not None:
            if not expr.expr.has(symbol):
                return
        return expr.key

    if expr.func is Predicate:
        return expr

    # To be used when Python 2.4 compatability is no longer required
    #return expr.func(*[x for x in (eliminate_assume(arg, symbol) for arg in expr.args) if x is not None])
    return expr.func(*filter(lambda x: x is not None, [eliminate_assume(arg, symbol) for arg in expr.args]))

class Predicate(Boolean):

    is_Atom = True

    def __new__(cls, name, handlers=None):
        obj = Boolean.__new__(cls)
        obj.name = name
        obj.handlers = handlers or []
        return obj

    def _hashable_content(self):
        return (self.name,)

    def __getnewargs__(self):
        return (self.name,)

    def __call__(self, expr):
        return Assume(expr, self.name)

    def add_handler(self, handler):
        self.handlers.append(handler)

    def remove_handler(self, handler):
        self.handlers.remove(handler)
