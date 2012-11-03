from itertools import groupby

from sympy.core import S, sympify, Mul, Basic, Expr, Dummy, Symbol, Integer
from sympy.combinatorics import Permutation

from sympy.matrices.expressions import MatMul, MatrixExpr, Trace
from sympy.matrices import MatrixBase, ImmutableMatrix

from sympy.physics.matrices import mgamma

from sympy.functions.special.tensor_functions import LeviCivita

from sympy.diffgeom.einstein_notation import MetricTangentSpace, U, D

from sympy.rules.strat_pure import condition, do_one, chain
from sympy.rules.traverse import top_down, bottom_up # TODO top_down_early_stop

from sympy.printing.pretty.stringpict import prettyForm
from sympy.printing.pretty.pretty_symbology import pretty_symbol


#  A globally used MetricTangentSpace
minkowski_metric_matrix = ImmutableMatrix([[-1,0,0,0],
                                           [ 0,1,0,0],
                                           [ 0,0,1,0],
                                           [ 0,0,0,1]])
qft_minkowski = MetricTangentSpace('minkowski', minkowski_metric_matrix)


class DiracMatrix(MatrixExpr):
    # XXX This should **not** subclass `TensorComponent` because we do not want
    # the usual rules for tensors to catch it.
    #  This is so because we implement the "components of the tensor being
    # operators" abuse of notation, instead of implementing dirac matrices as
    # tensors living on the product of Lorentz and Dirac spaces.
    #  We choose this way because the other alternative, namely moving all the
    # quantum-specific Hilbert space implementations to a more independent part
    # of SymPy, would be too costly.
    #  The drawback is that now the devs should understand both the
    # `einstein_notation` module and the `quantum` module.
    r"""
    If you want to create tensors in the same MetricTangentSpace use
    ``space = DiracMatrix.tangent_space`` or ``space =
    diracalgebra.qft_minkowski_space``.

    Examples
    ========

    >>> from sympy import symbols, pprint
    >>> from sympy.physics.diracalgebra import *
    >>> from sympy.matrices.expressions import Trace

    >>> a, b, c, d, e = symbols('a b c d e')
    >>> ga, gb, gc, gd, ge = [DiracMatrix(U(i)) for i in symbols('a b c d e')]
    >>> g_a, g_b, g_c, g_d, g_e = [DiracMatrix(D(i)) for i in symbols('a b c d e')]
    >>> g5 = DiracMatrix(U(5))
    >>> g0, g1, g2, g3 = [DiracMatrix(U(i)) for i in range(4)]
    >>> g_0, g_1, g_2, g_3 = [DiracMatrix(D(i)) for i in range(4)]

    >>> pprint( qft_minkowski.metric(U(Symbol('a')), U(Symbol('b'))))
      ab
     g

    >>> pprint( g_0*g_2)
    gamma *gamma
         0      2

    >>> pprint( Trace(DiracMatrix(U(a))*DiracMatrix(U(5))))
         /     a      5\
    Trace\gamma *gamma /

    >>> dirac_simplify(Trace(ga))
    0
    >>> dirac_simplify(Trace(g5))
    0
    >>> dirac_simplify(Trace(g5*ga))
    0
    >>> dirac_simplify(Trace(gb*g5*ga))
    0

    >>> pprint( dirac_simplify(Trace(gb*gb*ga*gc)))
    0
    >>> pprint( dirac_simplify(Trace(gb*gd*ga*gc)))
    0
    """
    #__new__(self, index):

    tangent_space = qft_minkowski
    is_commutative = False
    shape = 4*S.One, 4*S.One # TODO these should just be 4, 4
    @property
    def indices(self):
        return self.args

    def _entry(self, i, j):
        raise NotImplementedError

    @property
    def free_symbols(self):
        return set((self.args[0],))

    def _pretty(self, printer, *args):
        index = self.indices[0]
        index_pretty = printer._print(index.args[0])
        gamma = prettyForm(pretty_symbol('gamma'))
        if isinstance(index, U):
            return gamma**index_pretty
        top = prettyForm(*index_pretty.left(' '*gamma.width()))
        bot = prettyForm(*gamma.right(' '*index_pretty.width()))
        return prettyForm(binding=prettyForm.POW, *bot.below(top))


##############################################################################
# Rules.
##############################################################################

# The general idea:
#  1. Move all gamma matrix indices up
#  2. Do the simplifications
#      - traces that are zero
#      - canonical ordering
#      - squares of matrices that give g*Id
#      - traces that give the metric or LeviCivita          TODO
#  3. Contact with the metric, so the indices are back down TODO
#      - contactions that give Id                           TODO
#      - contractions that just lower the index             TODO
#


# Simplification rules unrelated to the formalism here.
# TODO move them to a better location (inside the concerned objects)
def trace_remove_scalars(a_trace): #TODO should not be here
    """Tr(scalar*matrix) -> scalar*Tr(matrix)"""
    if isinstance(a_trace, Trace) and isinstance(a_trace.args[0], MatMul):
        s = [s for s in a_trace.args[0].args if not isinstance(s, (MatrixBase, MatrixExpr))]
        m = [m for m in a_trace.args[0].args if isinstance(m, (MatrixBase, MatrixExpr))]
        return Mul(*s)*Trace(MatMul(*m)) if m else Mul(*s) #XXX Why is Trace a MatExpr
    return a_trace


# Conditions used for matching nodes on which to work.
def is_gamma_index_down(expr):
    """is a dirac matrix with lower index"""
    return isinstance(expr, DiracMatrix) and isinstance(expr.indices[0], D)

def is_mul_gamma(expr):
    """is a product of gamma matrices"""
    return isinstance(expr, MatMul) and all(isinstance(a, DiracMatrix) for a in expr.args)

def is_trace_gamma(expr):
    """is a trace of a gamma matrix or a product of gamma matrices"""
    return isinstance(expr, Trace) and (isinstance(expr.args[0], DiracMatrix)
                                        or is_mul_gamma(expr.args[0]))

# 1. Move all gamma matrix indices up
def rise_index(expr):
    """gamma_a -> g_ab*gamma^b"""
    d = Dummy()
    orig_index = expr.indices[0].args[0]
    if isinstance(orig_index, Symbol):
        return qft_minkowski.metric(D(orig_index), D(d))*DiracMatrix(U(d))
    else:
        return qft_minkowski.metric_matrix[orig_index, orig_index]*DiracMatrix(U(orig_index))


# 2. Do te simplifications
def tr_zero(expr):
    """Tr(odd nb of gammas * gamma5) -> 0 and similar"""
    g5s = sum(a == DiracMatrix(U(5)) for a in expr.args[0].args)
    gs  = sum(a != DiracMatrix(U(5)) for a in expr.args[0].args)
    if gs % 2 == 1 or (g5s%2 == 1 and gs in (0, 2)):
        return 0
    return expr


def canonical_order(expr):
    """product of gammas -> numeric gammas * alphabetic gammas * g5"""
    def order_args(args):
        perm_id = range(len(args))
        to_be_sorted = zip(args, perm_id)
        def sort_key(k):
            k = k[0].args[0].args[0] # DiracMatrix.U_or_D.Symbol
            if k in [S.One*i for i in range(4)]:
                return int(k)
            elif k == S.One*5:
                return 5
            else:
                return -abs(hash(k.name))#TODO hash collisions
        ordered, perm = zip(*sorted(to_be_sorted, key=sort_key))
        return ordered, Permutation(perm).signature()

    ordered, sign = order_args(expr.args)
    return sign*MatMul(*ordered)


def squares(expr):
    """gamma^a*gamma^a -> g^{aa}*Id for canonically ordered products"""
    args = []
    margs = []
    for k, g in groupby(expr.args):
        print k
        g = list(g)
        squares = len(g)//2
        margs.extend([k]*(len(g)%2))
        if k != DiracMatrix(U(5)):
            args.extend([qft_minkowski.metric(*k.indices*2)]*squares)
    return Mul(*args)*MatMul(*margs)


# Very naive implementation
dirac_simplify = chain(
    top_down( condition( is_gamma_index_down,
                         rise_index)),
    top_down( trace_remove_scalars),
    top_down( condition( is_trace_gamma,
                         tr_zero)),
    top_down( condition( is_mul_gamma,
                         chain( canonical_order,
                                squares))),
                                #))),
    top_down( trace_remove_scalars))
