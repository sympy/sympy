from sympy import Basic, Wild
from sympy.core.operations import AssocOp
from unify import Compound, Variable, _unify
from unify import *

def is_associative(x):
    return isinstance(x, Compound) and issubclass(x.op, AssocOp)

def is_commutative(x):
    # return isinstance(x, Compound) and construct(x).is_commutative
    return isinstance(x, Compound) and _build(x).is_commutative

def destruct(s):
    """ Turn a SymPy object into a Compound Tuple """
    if isinstance(s, Wild):
        return Variable(s)
    if not isinstance(s, Basic) or s.is_Atom:
        return s
    return Compound(s.__class__, tuple(map(destruct, s.args)))

def rebuild(s):
    """ Rebuild a SymPy expression using auto-evaluation """
    return _build(destruct(s))

def _build(t):
    # This function is the same as construct but builds the op in an
    # autoevaluated way. We do this only to get is_commutative. This function
    # should not be used otherwise.
    if isinstance(t, Variable):
        return t.arg
    if not isinstance(t, Compound):
        return t
    # This does auto-evaluation. Watch out!
    return t.op(*map(_build, t.args))

def construct(t):
    """ Turn a Compound Tuple into a SymPy object """
    if isinstance(t, Variable):
        return t.arg
    if not isinstance(t, Compound):
        return t
    return Basic.__new__(t.op, *map(construct, t.args))

def unify(x, y, s, **kwargs):
    """ Structural unification of two expressions possibly containing Wilds

    >>> from unify_sympy import unify
    >>> from sympy import Wild
    >>> from sympy.abc import x, y, z
    >>> expr = 2*x + y + z
    >>> pattern = 2*Wild('p') + Wild('q')
    >>> list(unify(expr, pattern, {}))
    [{p_: x, q_: y + z}]

    >>> expr = x + y + z
    >>> pattern = Wild('p') + Wild('q')
    >>> list(unify(expr, pattern, {}))
    [{p_: z, q_: x + y}, {p_: y + z, q_: x}]
    """

    ds = _unify(destruct(x), destruct(y), {}, is_associative=is_associative,
                                              is_commutative=is_commutative,
                                              **kwargs)
    for d in ds:
        yield {construct(k): construct(v) for k, v in d.items()}

def rewriterule(p1, p2):
    from sympy.rules.tools import subs
    from sympy import Expr
    def rewrite_rl(expr):
        match = unify(p1, expr, {})
        for m in match:
            expr2 = subs(m)(p2)
            if isinstance(expr2, Expr):
                expr2 = rebuild(expr2)
            yield expr2
    return rewrite_rl
