"""
Classes and functions useful for rewriting expressions for optimized code
generation. Some languages (or standards thereof), e.g. C99, offer specialized
math functions for better performance and/or precision.

Using the ``optimize`` function in this module, together with a collection of
rules (represented as instances of ``Optimization``), one can rewrite the
expressions for this purpose::

    >>> from sympy import Symbol, exp, log
    >>> from sympy.codegen.rewriting import optimize, optims_c99
    >>> x = Symbol('x')
    >>> optimize(3*exp(2*x) - 3, optims_c99)
    3*expm1(2*x)
    >>> optimize(exp(2*x) - 3, optims_c99)
    exp(2*x) - 3
    >>> optimize(log(3*x + 3), optims_c99)
    log1p(x) + log(3)
    >>> optimize(log(2*x + 3), optims_c99)
    log(2*x + 3)

The ``optims_c99`` imported above is tuple containing the following instances
(which may be imported from ``sympy.codegen.rewriting``):

- ``expm1_opt``
- ``log1p_opt``
- ``exp2_opt``
- ``log2_opt``
- ``log2const_opt``


"""
from itertools import chain
from sympy import cos, exp, log, Max, Min, Wild, expand_log, Dummy, sin, sinc
from sympy.assumptions import Q, ask
from sympy.codegen.cfunctions import log1p, log2, exp2, expm1
from sympy.codegen.matrix_nodes import MatrixSolve
from sympy.core.expr import UnevaluatedExpr
from sympy.core.power import Pow
from sympy.codegen.numpy_nodes import logaddexp, logaddexp2
from sympy.codegen.scipy_nodes import cosm1
from sympy.core.mul import Mul
from sympy.matrices.expressions.matexpr import MatrixSymbol
from sympy.utilities.iterables import sift


class Optimization:
    """ Abstract base class for rewriting optimization.

    Subclasses should implement ``__call__`` taking an expression
    as argument.

    Parameters
    ==========
    cost_function : callable returning number
    priority : number

    """
    def __init__(self, cost_function=None, priority=1):
        self.cost_function = cost_function
        self.priority=priority


class ReplaceOptim(Optimization):
    """ Rewriting optimization calling replace on expressions.

    Explanation
    ===========

    The instance can be used as a function on expressions for which
    it will apply the ``replace`` method (see
    :meth:`sympy.core.basic.Basic.replace`).

    Parameters
    ==========

    query :
        First argument passed to replace.
    value :
        Second argument passed to replace.

    Examples
    ========

    >>> from sympy import Symbol
    >>> from sympy.codegen.rewriting import ReplaceOptim
    >>> from sympy.codegen.cfunctions import exp2
    >>> x = Symbol('x')
    >>> exp2_opt = ReplaceOptim(lambda p: p.is_Pow and p.base == 2,
    ...     lambda p: exp2(p.exp))
    >>> exp2_opt(2**x)
    exp2(x)

    """

    def __init__(self, query, value, **kwargs):
        super().__init__(**kwargs)
        self.query = query
        self.value = value

    def __call__(self, expr):
        return expr.replace(self.query, self.value)


def optimize(expr, optimizations):
    """ Apply optimizations to an expression.

    Parameters
    ==========

    expr : expression
    optimizations : iterable of ``Optimization`` instances
        The optimizations will be sorted with respect to ``priority`` (highest first).

    Examples
    ========

    >>> from sympy import log, Symbol
    >>> from sympy.codegen.rewriting import optims_c99, optimize
    >>> x = Symbol('x')
    >>> optimize(log(x+3)/log(2) + log(x**2 + 1), optims_c99)
    log1p(x**2) + log2(x + 3)

    """

    for optim in sorted(optimizations, key=lambda opt: opt.priority, reverse=True):
        new_expr = optim(expr)
        if optim.cost_function is None:
            expr = new_expr
        else:
            before, after = map(lambda x: optim.cost_function(x), (expr, new_expr))
            if before > after:
                expr = new_expr
    return expr


exp2_opt = ReplaceOptim(
    lambda p: p.is_Pow and p.base == 2,
    lambda p: exp2(p.exp)
)


_d = Wild('d', properties=[lambda x: x.is_Dummy])
_u = Wild('u', properties=[lambda x: not x.is_number and not x.is_Add])
_v = Wild('v')
_w = Wild('w')
_n = Wild('n', properties=[lambda x: x.is_number])

sinc_opt1 = ReplaceOptim(
    sin(_w)/_w, sinc(_w)
)
sinc_opt2 = ReplaceOptim(
    sin(_n*_w)/_w, _n*sinc(_n*_w)
)
sinc_opts = (sinc_opt1, sinc_opt2)

log2_opt = ReplaceOptim(_v*log(_w)/log(2), _v*log2(_w), cost_function=lambda expr: expr.count(
    lambda e: (  # division & eval of transcendentals are expensive floating point operations...
        e.is_Pow and e.exp.is_negative  # division
        or (isinstance(e, (log, log2)) and not e.args[0].is_number))  # transcendental
    )
)

log2const_opt = ReplaceOptim(log(2)*log2(_w), log(_w))

logsumexp_2terms_opt = ReplaceOptim(
    lambda l: (isinstance(l, log)
               and l.args[0].is_Add
               and len(l.args[0].args) == 2
               and all(isinstance(t, exp) for t in l.args[0].args)),
    lambda l: (
        Max(*[e.args[0] for e in l.args[0].args]) +
        log1p(exp(Min(*[e.args[0] for e in l.args[0].args])))
    )
)


class _FuncMinusOne:
    def __init__(self, func, func_m_1):
        self.func = func
        self.func_m_1 = func_m_1

    def _try_func_m_1(self, expr):
        protected, old_new =  expr.replace(self.func, lambda arg: Dummy(), map=True)
        factored = protected.factor()
        new_old = {v: k for k, v in old_new.items()}
        return factored.replace(_d - 1, lambda d: self.func_m_1(new_old[d].args[0])).xreplace(new_old)

    def __call__(self, e):
        numbers, non_num = sift(e.args, lambda arg: arg.is_number, binary=True)
        non_num_func, non_num_other = sift(non_num, lambda arg: arg.has(self.func),
            binary=True)
        numsum = sum(numbers)
        new_func_terms, done = [], False
        for func_term in non_num_func:
            if done:
                new_func_terms.append(func_term)
            else:
                looking_at = func_term + numsum
                attempt = self._try_func_m_1(looking_at)
                if looking_at == attempt:
                    new_func_terms.append(func_term)
                else:
                    done = True
                    new_func_terms.append(attempt)
        if not done:
            new_func_terms.append(numsum)
        return e.func(*chain(new_func_terms, non_num_other))


expm1_opt = ReplaceOptim(lambda e: e.is_Add, _FuncMinusOne(exp, expm1))
cosm1_opt = ReplaceOptim(lambda e: e.is_Add, _FuncMinusOne(cos, cosm1))

log1p_opt = ReplaceOptim(
    lambda e: isinstance(e, log),
    lambda l: expand_log(l.replace(
        log, lambda arg: log(arg.factor())
    )).replace(log(_u+1), log1p(_u))
)

def create_expand_pow_optimization(limit, *, base_req=lambda b: b.is_symbol):
    """ Creates an instance of :class:`ReplaceOptim` for expanding ``Pow``.

    Explanation
    ===========

    The requirements for expansions are that the base needs to be a symbol
    and the exponent needs to be an Integer (and be less than or equal to
    ``limit``).

    Parameters
    ==========

    limit : int
         The highest power which is expanded into multiplication.
    base_req : function returning bool
         Requirement on base for expansion to happen, default is to return
         the ``is_symbol`` attribute of the base.

    Examples
    ========

    >>> from sympy import Symbol, sin
    >>> from sympy.codegen.rewriting import create_expand_pow_optimization
    >>> x = Symbol('x')
    >>> expand_opt = create_expand_pow_optimization(3)
    >>> expand_opt(x**5 + x**3)
    x**5 + x*x*x
    >>> expand_opt(x**5 + x**3 + sin(x)**3)
    x**5 + sin(x)**3 + x*x*x
    >>> opt2 = create_expand_pow_optimization(3 , base_req=lambda b: not b.is_Function)
    >>> opt2((x+1)**2 + sin(x)**2)
    sin(x)**2 + (x + 1)*(x + 1)

    """
    return ReplaceOptim(
        lambda e: e.is_Pow and base_req(e.base) and e.exp.is_Integer and abs(e.exp) <= limit,
        lambda p: (
            UnevaluatedExpr(Mul(*([p.base]*+p.exp), evaluate=False)) if p.exp > 0 else
            1/UnevaluatedExpr(Mul(*([p.base]*-p.exp), evaluate=False))
        ))

# Optimization procedures for turning A**(-1) * x into MatrixSolve(A, x)
def _matinv_predicate(expr):
    # TODO: We should be able to support more than 2 elements
    if expr.is_MatMul and len(expr.args) == 2:
        left, right = expr.args
        if left.is_Inverse and right.shape[1] == 1:
            inv_arg = left.arg
            if isinstance(inv_arg, MatrixSymbol):
                return bool(ask(Q.fullrank(left.arg)))

    return False

def _matinv_transform(expr):
    left, right = expr.args
    inv_arg = left.arg
    return MatrixSolve(inv_arg, right)


matinv_opt = ReplaceOptim(_matinv_predicate, _matinv_transform)


logaddexp_opt = ReplaceOptim(log(exp(_v)+exp(_w)), logaddexp(_v, _w))
logaddexp2_opt = ReplaceOptim(log(Pow(2, _v)+Pow(2, _w)), logaddexp2(_v, _w)*log(2))

# Collections of optimizations:
optims_c99 = (expm1_opt, log1p_opt, exp2_opt, log2_opt, log2const_opt)

optims_numpy = optims_c99 + (logaddexp_opt, logaddexp2_opt,) + sinc_opts

optims_scipy = (cosm1_opt,)
