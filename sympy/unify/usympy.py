from sympy import Basic, Expr, Tuple, Add, Mul, Pow
from sympy import Wild as ExprWild

from core import Compound, Variable
import core

class Wild(object):
    def __init__(self, arg):
        self.arg = arg
    def __eq__(self, other):
        return type(self) == type(other) and self.arg == other.arg

def sympy_associative(op):
    from sympy import MatAdd, MatMul, Union, Intersection
    from sympy.core.operations import AssocOp
    assoc_ops = (AssocOp, MatAdd, MatMul, Union, Intersection)
    return any(issubclass(op, aop) for aop in assoc_ops)

def sympy_commutative(op):
    from sympy import Add, MatAdd, Union, Intersection
    comm_ops = (Add, MatAdd, Union, Intersection)
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
    if t.op in (Add, Mul, Pow):
        return t.op(*map(construct, t.args), **{'evaluate': False})
    else:
        return Basic.__new__(t.op, *map(construct, t.args))

def rebuild(s):
    """ Rebuild a SymPy expression

    This removes harm caused by Expr-Rules interactions
    """
    return construct(deconstruct(s))

def unify(x, y, s={}, **kwargs):
    """ Structural unification of two expressions/patterns

    Examples
    ========

    >>> from sympy.unify.usympy import unify
    >>> from sympy import Wild
    >>> from sympy.abc import x, y, z
    >>> expr = 2*x + y + z
    >>> pattern = 2*Wild('p') + Wild('q')
    >>> next(unify(expr, pattern, {}))
    {p_: x, q_: y + z}

    >>> expr = x + y + z
    >>> pattern = Wild('p') + Wild('q')
    >>> len(list(unify(expr, pattern, {})))
    12
    """

    ds = core.unify(deconstruct(x), deconstruct(y), {},
                                                is_associative=is_associative,
                                                is_commutative=is_commutative,
                                                **kwargs)
    for d in ds:
        yield dict((construct(k), construct(v)) for k, v in d.items())
