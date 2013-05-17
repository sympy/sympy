from sympy.core import S
from sympy.core.sympify import sympify
from sympy.core.basic import Basic
from sympy.core.numbers import Integer

from sympy.combinatorics import Permutation
from sympy.matrices.expressions import MatMul, MatrixExpr, Trace

from sympy.physics.matrices import mgamma

from sympy.rules.strat_pure import condition, do_one
from sympy.rules.traverse import top_down # TODO top_down_early_stop

class DiracMatrix(MatrixExpr):
    """
    Examples
    ========

    >>> from sympy import symbols
    >>> from sympy.physics.diracalgebra import DiracMatrix, dirac_simplify
    >>> from sympy.matrices.expressions import Trace
    >>> ga, gb, gc, gd, ge = map(DiracMatrix, symbols('a b c d e'))
    >>> g0, g1, g2, g3, g5 = map(DiracMatrix, range(4)+[5])
    >>> dirac_simplify(Trace(ga))
    0
    >>> dirac_simplify(Trace(gb*g5*ga))
    0

    """
    is_commutative = False
    shape = 2, 2

    def __new__(cls, index):
        obj = Basic.__new__(cls, sympify(index))
        return obj

    def _entry(self, i, j):
        mgamma(self.args[0])[i, j]

    @property
    def free_symbols(self):
        return set((self.args[0],))


def cond_mul(expr):
    return isinstance(expr, MatMul) and all(isinstance(a, DiracMatrix) for a in expr.args)


def cond_trace(expr):
    return isinstance(expr, Trace) and (isinstance(expr.args[0], DiracMatrix)
                                        or cond_mul(expr.args[0]))


def get_args(trace_mul_args):
    return trace_mul_args.args[0].args


def tr_zero(expr):
    """Tr(even nb of gammas * gamma5) -> 0, etc."""
    g5s = sum(a == DiracMatrix(5) for a in get_args(expr))
    gs  = sum(a != DiracMatrix(5) for a in get_args(expr))
    if (g5s%2 == 1 and gs%2 == 0) or (g5s%2 == 0 and gs%2 == 1):
        return 0
    return expr


def canonical_order(expr):
    """product of gammas -> numeric gammas * alphabetic gammas * g5)"""
    def order_args(args):
        perm_id = range(len(args))
        to_be_sorted = zip(args, perm_id)
        def sort_key(k):
            k = k[0].args[0]
            if k in [S.One*i for i in range(4)]:
                return int(k)
            elif k == S.One*5:
                return 5
            else:
                return -abs(hash(k.name))
        ordered, perm = zip(*sorted(to_be_sorted, key=sort_key))
        return ordered, Permutation(perm).signature()

    ordered, sign = order_args(expr.args)
    return sign*MatMul(*ordered)


dirac_simplify = top_down(do_one(condition(cond_trace, tr_zero),
                                 condition(cond_mul,   canonical_order)))
