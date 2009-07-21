# doctests are disabled because of issue #1521
from sympy.core import Basic, Symbol
from sympy.core.relational import Relational

__global_assumptions = []

def register_global_assumptions(*assump):
    """Register an assumption as global

    Examples:
        >>> from sympy import *
        >>> list_global_assumptions()
        []
        >>> x = Symbol('x')
        >>> register_global_assumptions(Assume(x, Q.real))
        >>> list_global_assumptions()
        [Assume(x, 'real', True)]
        >>> clean_global_assumptions()
    """
    __global_assumptions.extend(assump)

def list_global_assumptions():
    """List all global assumptions

    Examples:
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> list_global_assumptions()
        []
        >>> register_global_assumptions(Assume(x>0))
        >>> list_global_assumptions()
        [Assume(0 < x, 'relational', True)]
        >>> clean_global_assumptions()
    """
    return __global_assumptions[:] # make a copy

def remove_global_assumptions(*assump):
    """Remove a global assumption. If argument is not
    a global assumption, it will raise  ValueError

    Examples:
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> register_global_assumptions(Assume(x, Q.even))
        >>> list_global_assumptions()
        [Assume(x, 'even', True)]
        >>> remove_global_assumptions(Assume(x, Q.even))
        >>> list_global_assumptions()
        []
    """
    for assumption in assump:
        __global_assumptions.remove(assumption)

def clean_global_assumptions():
    """
    Remove all global assumptions

    Examples:
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> register_global_assumptions(Assume(x, Q.integer))
        >>> list_global_assumptions()
        [Assume(x, 'integer', True)]
        >>> clean_global_assumptions()
        >>> list_global_assumptions()
        []
    """
    global __global_assumptions
    __global_assumptions = []

class Assume(Basic):
    """New-style assumptions

    >>> from sympy import *
    >>> x = Symbol('x')
    >>> Assume(x, Q.integer)
    Assume(x, 'integer', True)
    >>> Assume(x, Q.integer, False)
    Assume(x, 'integer', False)
    >>> Assume( x > 1 )
    Assume(1 < x, 'relational', True)
    """
    def __init__(self, expr, key='relational', value=True):
        self._args = (expr, key, value)

    is_Atom = True # do not attempt to decompose this

    @property
    def expr(self):
        """
        Returns the expression used by this assumption

        Examples:
            >>> from sympy import *
            >>> x = Symbol('x')
            >>> a = Assume(x+1, Q.integer)
            >>> a.expr
            1 + x
        """
        return self._args[0]

    @property
    def key(self):
        """
        Returns the key used by this assumption.
        It is a string, e.g. 'integer', 'rational', etc.

        Examples:
            >>> from sympy import *
            >>> x = Symbol('x')
            >>> a = Assume(x, Q.integer)
            >>> a.key
            'integer'
        """
        return self._args[1]

    @property
    def value(self):
        """
        Returns the value stored by this assumptions.
        It's a boolean. True means that the assumption
        holds always, and False means the assumption
        does not hold

        Examples:
            >>> from sympy import *
            >>> x = Symbol('x')
            >>> a = Assume(x, Q.integer)
            >>> a.value
            True
            >>> b = Assume(x, Q.integer, False)
            >>> b.value
            False
        """
        return self._args[2]

    def __eq__(self, other):
        if type(other) == Assume:
            return self._args == other._args
        return False

def eliminate_assume(expr, symbol=None):
    """
    Will convert an expression with assumptions to an equivalent with all assumptions
    replaced by symbols
    Assume(x, integer=True) --> integer
    Assume(x, integer=False) --> ~integer

    Examples:
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> eliminate_assume(Assume(x, Q.positive))
        positive
        >>> eliminate_assume(Assume(x, Q.positive, False))
        Not(positive)
    """
    if type(expr) == Assume:
        if symbol is not None:
            if not expr.expr.has(symbol): return
        if expr.value: return Symbol(expr.key)
        return ~Symbol(expr.key)
    args = []
    for a in expr.args:
        args.append(eliminate_assume(a))
    return type(expr)(*args)

