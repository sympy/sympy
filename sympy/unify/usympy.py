""" SymPy interface to Unificaiton engine

See sympy.unify for module level docstring
See sympy.unify.core for algorithmic docstring """


from sympy.core import Basic, Expr, Tuple, Add, Mul, Pow, FiniteSet
from sympy.core import Wild as ExprWild
from sympy.matrices import MatAdd, MatMul
from sympy.core.sets import Union, Intersection, FiniteSet
from sympy.core.operations import AssocOp
from sympy.unify.core import Compound, Variable
from sympy.unify import core

class Wild(object):
    def __init__(self, arg):
        self.arg = arg
    def __eq__(self, other):
        return type(self) == type(other) and self.arg == other.arg

def sympy_associative(op):
    assoc_ops = (AssocOp, MatAdd, MatMul, Union, Intersection, FiniteSet)
    return any(issubclass(op, aop) for aop in assoc_ops)

def sympy_commutative(op):
    comm_ops = (Add, MatAdd, Union, Intersection, FiniteSet)
    return any(issubclass(op, cop) for cop in comm_ops)

def is_associative(x):
    return isinstance(x, Compound) and sympy_associative(x.op)

def is_commutative(x):
    if not isinstance(x, Compound):
        return False
    if sympy_commutative(x.op):
        return True
    if isinstance(x.op, Expr):
        return _build(x).is_commutative

def patternify(expr, *wilds):
    """ Create a matching pattern from an expression

    Example
    =======

    >>> from sympy import symbols, sin, cos
    >>> from sympy.unify.usympy import patternify
    >>> a, b, c, x, y = symbols('a b c x y')

    >>> # Search for anything of the form sin(foo)**2 + cos(foo)**2
    >>> pattern = patternify(sin(x)**2 + cos(x)**2, x)

    >>> # Search for any two things added to c. Note that here c is not a wild
    >>> a, b, c = symbols('a,b,c')
    >>> pattern = patternify(a + b + c, a, b)
    """
    from sympy.rules.tools import subs
    return subs(dict(zip(wilds, map(Wild, wilds))))(expr)

def deconstruct(s):
    """ Turn a SymPy object into a Compound """
    if isinstance(s, ExprWild):
        return Variable(s)
    if isinstance(s, Wild):
        return Variable(s.arg)
    if not isinstance(s, Basic) or s.is_Atom:
        return s
    return Compound(s.__class__, tuple(map(deconstruct, s.args)))

def construct(t):
    """ Turn a Compound into a SymPy object """
    if isinstance(t, Variable):
        return t.arg
    if not isinstance(t, Compound):
        return t
    if t.op in (Add, Mul, Pow, FiniteSet):
        return t.op(*map(construct, t.args), **{'evaluate': False})
    else:
        return Basic.__new__(t.op, *map(construct, t.args))

def rebuild(s):
    """ Rebuild a SymPy expression

    This removes harm caused by Expr-Rules interactions
    """
    return construct(deconstruct(s))

def unify(x, y, s=None, **kwargs):
    """ Structural unification of two expressions/patterns

    Examples
    ========

    >>> from sympy.unify.usympy import unify
    >>> from sympy import Wild
    >>> from sympy.abc import x, y, z
    >>> from sympy.core.compatibility import next
    >>> expr = 2*x + y + z
    >>> pattern = 2*Wild('p') + Wild('q')
    >>> next(unify(expr, pattern, {}))
    {p_: x, q_: y + z}

    >>> expr = x + y + z
    >>> pattern = Wild('p') + Wild('q')
    >>> len(list(unify(expr, pattern, {})))
    12
    """
    s = s or {}
    s = dict((deconstruct(k), deconstruct(v)) for k, v in s.items())

    ds = core.unify(deconstruct(x), deconstruct(y), s,
                                                is_associative=is_associative,
                                                is_commutative=is_commutative,
                                                **kwargs)
    for d in ds:
        yield dict((construct(k), construct(v)) for k, v in d.items())
