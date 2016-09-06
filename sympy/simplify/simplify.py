from __future__ import print_function, division

from collections import defaultdict

from sympy.core import (Basic, S, Add, Mul, Pow,
    Symbol, sympify, expand_mul, expand_func,
    Function, Dummy, Expr, factor_terms,
    symbols, expand_power_exp)
from sympy.core.compatibility import (iterable,
    ordered, range, as_int)
from sympy.core.numbers import Float, I, pi, Rational, Integer
from sympy.core.function import expand_log, count_ops, _mexpand, _coeff_isneg
from sympy.core.rules import Transform
from sympy.core.evaluate import global_evaluate
from sympy.functions import (
    gamma, exp, sqrt, log, exp_polar, piecewise_fold)
from sympy.functions.elementary.exponential import ExpBase
from sympy.functions.elementary.hyperbolic import HyperbolicFunction
from sympy.functions.elementary.integers import ceiling
from sympy.functions.elementary.complexes import unpolarify
from sympy.functions.elementary.trigonometric import TrigonometricFunction
from sympy.functions.combinatorial.factorials import CombinatorialFunction
from sympy.functions.special.bessel import besselj, besseli, besselk, jn, bessely

from sympy.utilities.iterables import has_variety

from sympy.simplify.radsimp import radsimp, fraction
from sympy.simplify.trigsimp import trigsimp, exptrigsimp
from sympy.simplify.powsimp import powsimp
from sympy.simplify.cse_opts import sub_pre, sub_post
from sympy.simplify.sqrtdenest import sqrtdenest
from sympy.simplify.combsimp import combsimp

from sympy.polys import (together, cancel, factor)


import mpmath



def separatevars(expr, symbols=[], dict=False, force=False):
    """
    Separates variables in an expression, if possible.  By
    default, it separates with respect to all symbols in an
    expression and collects constant coefficients that are
    independent of symbols.

    If dict=True then the separated terms will be returned
    in a dictionary keyed to their corresponding symbols.
    By default, all symbols in the expression will appear as
    keys; if symbols are provided, then all those symbols will
    be used as keys, and any terms in the expression containing
    other symbols or non-symbols will be returned keyed to the
    string 'coeff'. (Passing None for symbols will return the
    expression in a dictionary keyed to 'coeff'.)

    If force=True, then bases of powers will be separated regardless
    of assumptions on the symbols involved.

    Notes
    =====
    The order of the factors is determined by Mul, so that the
    separated expressions may not necessarily be grouped together.

    Although factoring is necessary to separate variables in some
    expressions, it is not necessary in all cases, so one should not
    count on the returned factors being factored.

    Examples
    ========

    >>> from sympy.abc import x, y, z, alpha
    >>> from sympy import separatevars, sin
    >>> separatevars((x*y)**y)
    (x*y)**y
    >>> separatevars((x*y)**y, force=True)
    x**y*y**y

    >>> e = 2*x**2*z*sin(y)+2*z*x**2
    >>> separatevars(e)
    2*x**2*z*(sin(y) + 1)
    >>> separatevars(e, symbols=(x, y), dict=True)
    {'coeff': 2*z, x: x**2, y: sin(y) + 1}
    >>> separatevars(e, [x, y, alpha], dict=True)
    {'coeff': 2*z, alpha: 1, x: x**2, y: sin(y) + 1}

    If the expression is not really separable, or is only partially
    separable, separatevars will do the best it can to separate it
    by using factoring.

    >>> separatevars(x + x*y - 3*x**2)
    -x*(3*x - y - 1)

    If the expression is not separable then expr is returned unchanged
    or (if dict=True) then None is returned.

    >>> eq = 2*x + y*sin(x)
    >>> separatevars(eq) == eq
    True
    >>> separatevars(2*x + y*sin(x), symbols=(x, y), dict=True) == None
    True

    """
    expr = sympify(expr)
    if dict:
        return _separatevars_dict(_separatevars(expr, force), symbols)
    else:
        return _separatevars(expr, force)


def _separatevars(expr, force):
    if len(expr.free_symbols) == 1:
        return expr
    # don't destroy a Mul since much of the work may already be done
    if expr.is_Mul:
        args = list(expr.args)
        changed = False
        for i, a in enumerate(args):
            args[i] = separatevars(a, force)
            changed = changed or args[i] != a
        if changed:
            expr = expr.func(*args)
        return expr

    # get a Pow ready for expansion
    if expr.is_Pow:
        expr = Pow(separatevars(expr.base, force=force), expr.exp)

    # First try other expansion methods
    expr = expr.expand(mul=False, multinomial=False, force=force)

    _expr, reps = posify(expr) if force else (expr, {})
    expr = factor(_expr).subs(reps)

    if not expr.is_Add:
        return expr

    # Find any common coefficients to pull out
    args = list(expr.args)
    commonc = args[0].args_cnc(cset=True, warn=False)[0]
    for i in args[1:]:
        commonc &= i.args_cnc(cset=True, warn=False)[0]
    commonc = Mul(*commonc)
    commonc = commonc.as_coeff_Mul()[1]  # ignore constants
    commonc_set = commonc.args_cnc(cset=True, warn=False)[0]

    # remove them
    for i, a in enumerate(args):
        c, nc = a.args_cnc(cset=True, warn=False)
        c = c - commonc_set
        args[i] = Mul(*c)*Mul(*nc)
    nonsepar = Add(*args)

    if len(nonsepar.free_symbols) > 1:
        _expr = nonsepar
        _expr, reps = posify(_expr) if force else (_expr, {})
        _expr = (factor(_expr)).subs(reps)

        if not _expr.is_Add:
            nonsepar = _expr

    return commonc*nonsepar


def _separatevars_dict(expr, symbols):
    if symbols:
        if not all((t.is_Atom for t in symbols)):
            raise ValueError("symbols must be Atoms.")
        symbols = list(symbols)
    elif symbols is None:
        return {'coeff': expr}
    else:
        symbols = list(expr.free_symbols)
        if not symbols:
            return None

    ret = dict(((i, []) for i in symbols + ['coeff']))

    for i in Mul.make_args(expr):
        expsym = i.free_symbols
        intersection = set(symbols).intersection(expsym)
        if len(intersection) > 1:
            return None
        if len(intersection) == 0:
            # There are no symbols, so it is part of the coefficient
            ret['coeff'].append(i)
        else:
            ret[intersection.pop()].append(i)

    # rebuild
    for k, v in ret.items():
        ret[k] = Mul(*v)

    return ret


def _is_sum_surds(p):
    args = p.args if p.is_Add else [p]
    for y in args:
        if not ((y**2).is_Rational and y.is_real):
            return False
    return True


def posify(eq):
    """Return eq (with generic symbols made positive) and a
    dictionary containing the mapping between the old and new
    symbols.

    Any symbol that has positive=None will be replaced with a positive dummy
    symbol having the same name. This replacement will allow more symbolic
    processing of expressions, especially those involving powers and
    logarithms.

    A dictionary that can be sent to subs to restore eq to its original
    symbols is also returned.

    >>> from sympy import posify, Symbol, log, solve
    >>> from sympy.abc import x
    >>> posify(x + Symbol('p', positive=True) + Symbol('n', negative=True))
    (_x + n + p, {_x: x})

    >>> eq = 1/x
    >>> log(eq).expand()
    log(1/x)
    >>> log(posify(eq)[0]).expand()
    -log(_x)
    >>> p, rep = posify(eq)
    >>> log(p).expand().subs(rep)
    -log(x)

    It is possible to apply the same transformations to an iterable
    of expressions:

    >>> eq = x**2 - 4
    >>> solve(eq, x)
    [-2, 2]
    >>> eq_x, reps = posify([eq, x]); eq_x
    [_x**2 - 4, _x]
    >>> solve(*eq_x)
    [2]
    """
    eq = sympify(eq)
    if iterable(eq):
        f = type(eq)
        eq = list(eq)
        syms = set()
        for e in eq:
            syms = syms.union(e.atoms(Symbol))
        reps = {}
        for s in syms:
            reps.update(dict((v, k) for k, v in posify(s)[1].items()))
        for i, e in enumerate(eq):
            eq[i] = e.subs(reps)
        return f(eq), {r: s for s, r in reps.items()}

    reps = dict([(s, Dummy(s.name, positive=True))
                 for s in eq.free_symbols if s.is_positive is None])
    eq = eq.subs(reps)
    return eq, {r: s for s, r in reps.items()}


def hypersimp(f, k):
    """Given combinatorial term f(k) simplify its consecutive term ratio
       i.e. f(k+1)/f(k).  The input term can be composed of functions and
       integer sequences which have equivalent representation in terms
       of gamma special function.

       The algorithm performs three basic steps:

       1. Rewrite all functions in terms of gamma, if possible.

       2. Rewrite all occurrences of gamma in terms of products
          of gamma and rising factorial with integer,  absolute
          constant exponent.

       3. Perform simplification of nested fractions, powers
          and if the resulting expression is a quotient of
          polynomials, reduce their total degree.

       If f(k) is hypergeometric then as result we arrive with a
       quotient of polynomials of minimal degree. Otherwise None
       is returned.

       For more information on the implemented algorithm refer to:

       1. W. Koepf, Algorithms for m-fold Hypergeometric Summation,
          Journal of Symbolic Computation (1995) 20, 399-417
    """
    f = sympify(f)

    g = f.subs(k, k + 1) / f

    g = g.rewrite(gamma)
    g = expand_func(g)
    g = powsimp(g, deep=True, combine='exp')

    if g.is_rational_function(k):
        return simplify(g, ratio=S.Infinity)
    else:
        return None


def hypersimilar(f, g, k):
    """Returns True if 'f' and 'g' are hyper-similar.

       Similarity in hypergeometric sense means that a quotient of
       f(k) and g(k) is a rational function in k.  This procedure
       is useful in solving recurrence relations.

       For more information see hypersimp().

    """
    f, g = list(map(sympify, (f, g)))

    h = (f/g).rewrite(gamma)
    h = h.expand(func=True, basic=False)

    return h.is_rational_function(k)


def signsimp(expr, evaluate=None):
    """Make all Add sub-expressions canonical wrt sign.

    If an Add subexpression, ``a``, can have a sign extracted,
    as determined by could_extract_minus_sign, it is replaced
    with Mul(-1, a, evaluate=False). This allows signs to be
    extracted from powers and products.

    Examples
    ========

    >>> from sympy import signsimp, exp, symbols
    >>> from sympy.abc import x, y
    >>> i = symbols('i', odd=True)
    >>> n = -1 + 1/x
    >>> n/x/(-n)**2 - 1/n/x
    (-1 + 1/x)/(x*(1 - 1/x)**2) - 1/(x*(-1 + 1/x))
    >>> signsimp(_)
    0
    >>> x*n + x*-n
    x*(-1 + 1/x) + x*(1 - 1/x)
    >>> signsimp(_)
    0

    Since powers automatically handle leading signs

    >>> (-2)**i
    -2**i

    signsimp can be used to put the base of a power with an integer
    exponent into canonical form:

    >>> n**i
    (-1 + 1/x)**i

    By default, signsimp doesn't leave behind any hollow simplification:
    if making an Add canonical wrt sign didn't change the expression, the
    original Add is restored. If this is not desired then the keyword
    ``evaluate`` can be set to False:

    >>> e = exp(y - x)
    >>> signsimp(e) == e
    True
    >>> signsimp(e, evaluate=False)
    exp(-(x - y))

    """
    if evaluate is None:
        evaluate = global_evaluate[0]
    expr = sympify(expr)
    if not isinstance(expr, Expr) or expr.is_Atom:
        return expr
    e = sub_post(sub_pre(expr))
    if not isinstance(e, Expr) or e.is_Atom:
        return e
    if e.is_Add:
        return e.func(*[signsimp(a) for a in e.args])
    if evaluate:
        e = e.xreplace({m: -(-m) for m in e.atoms(Mul) if -(-m) != m})
    return e


def simplify(expr, ratio=1.7, measure=count_ops, fu=False):
    """
    Simplifies the given expression.

    Simplification is not a well defined term and the exact strategies
    this function tries can change in the future versions of SymPy. If
    your algorithm relies on "simplification" (whatever it is), try to
    determine what you need exactly  -  is it powsimp()?, radsimp()?,
    together()?, logcombine()?, or something else? And use this particular
    function directly, because those are well defined and thus your algorithm
    will be robust.

    Nonetheless, especially for interactive use, or when you don't know
    anything about the structure of the expression, simplify() tries to apply
    intelligent heuristics to make the input expression "simpler".  For
    example:

    >>> from sympy import simplify, cos, sin
    >>> from sympy.abc import x, y
    >>> a = (x + x**2)/(x*sin(y)**2 + x*cos(y)**2)
    >>> a
    (x**2 + x)/(x*sin(y)**2 + x*cos(y)**2)
    >>> simplify(a)
    x + 1

    Note that we could have obtained the same result by using specific
    simplification functions:

    >>> from sympy import trigsimp, cancel
    >>> trigsimp(a)
    (x**2 + x)/x
    >>> cancel(_)
    x + 1

    In some cases, applying :func:`simplify` may actually result in some more
    complicated expression. The default ``ratio=1.7`` prevents more extreme
    cases: if (result length)/(input length) > ratio, then input is returned
    unmodified.  The ``measure`` parameter lets you specify the function used
    to determine how complex an expression is.  The function should take a
    single argument as an expression and return a number such that if
    expression ``a`` is more complex than expression ``b``, then
    ``measure(a) > measure(b)``.  The default measure function is
    :func:`count_ops`, which returns the total number of operations in the
    expression.

    For example, if ``ratio=1``, ``simplify`` output can't be longer
    than input.

    ::

        >>> from sympy import sqrt, simplify, count_ops, oo
        >>> root = 1/(sqrt(2)+3)

    Since ``simplify(root)`` would result in a slightly longer expression,
    root is returned unchanged instead::

       >>> simplify(root, ratio=1) == root
       True

    If ``ratio=oo``, simplify will be applied anyway::

        >>> count_ops(simplify(root, ratio=oo)) > count_ops(root)
        True

    Note that the shortest expression is not necessary the simplest, so
    setting ``ratio`` to 1 may not be a good idea.
    Heuristically, the default value ``ratio=1.7`` seems like a reasonable
    choice.

    You can easily define your own measure function based on what you feel
    should represent the "size" or "complexity" of the input expression.  Note
    that some choices, such as ``lambda expr: len(str(expr))`` may appear to be
    good metrics, but have other problems (in this case, the measure function
    may slow down simplify too much for very large expressions).  If you don't
    know what a good metric would be, the default, ``count_ops``, is a good
    one.

    For example:

    >>> from sympy import symbols, log
    >>> a, b = symbols('a b', positive=True)
    >>> g = log(a) + log(b) + log(a)*log(1/b)
    >>> h = simplify(g)
    >>> h
    log(a*b**(-log(a) + 1))
    >>> count_ops(g)
    8
    >>> count_ops(h)
    5

    So you can see that ``h`` is simpler than ``g`` using the count_ops metric.
    However, we may not like how ``simplify`` (in this case, using
    ``logcombine``) has created the ``b**(log(1/a) + 1)`` term.  A simple way
    to reduce this would be to give more weight to powers as operations in
    ``count_ops``.  We can do this by using the ``visual=True`` option:

    >>> print(count_ops(g, visual=True))
    2*ADD + DIV + 4*LOG + MUL
    >>> print(count_ops(h, visual=True))
    2*LOG + MUL + POW + SUB

    >>> from sympy import Symbol, S
    >>> def my_measure(expr):
    ...     POW = Symbol('POW')
    ...     # Discourage powers by giving POW a weight of 10
    ...     count = count_ops(expr, visual=True).subs(POW, 10)
    ...     # Every other operation gets a weight of 1 (the default)
    ...     count = count.replace(Symbol, type(S.One))
    ...     return count
    >>> my_measure(g)
    8
    >>> my_measure(h)
    14
    >>> 15./8 > 1.7 # 1.7 is the default ratio
    True
    >>> simplify(g, measure=my_measure)
    -log(a)*log(b) + log(a) + log(b)

    Note that because ``simplify()`` internally tries many different
    simplification strategies and then compares them using the measure
    function, we get a completely different result that is still different
    from the input expression by doing this.
    """
    expr = sympify(expr)

    try:
        return expr._eval_simplify(ratio=ratio, measure=measure)
    except AttributeError:
        pass

    original_expr = expr = signsimp(expr)

    from sympy.simplify.hyperexpand import hyperexpand
    from sympy.functions.special.bessel import BesselBase
    from sympy import Sum, Product

    if not isinstance(expr, Basic) or not expr.args:  # XXX: temporary hack
        return expr

    if not isinstance(expr, (Add, Mul, Pow, ExpBase)):
        if isinstance(expr, Function) and hasattr(expr, "inverse"):
            if len(expr.args) == 1 and len(expr.args[0].args) == 1 and \
               isinstance(expr.args[0], expr.inverse(argindex=1)):
                return simplify(expr.args[0].args[0], ratio=ratio,
                                measure=measure, fu=fu)
        return expr.func(*[simplify(x, ratio=ratio, measure=measure, fu=fu)
                         for x in expr.args])

    # TODO: Apply different strategies, considering expression pattern:
    # is it a purely rational function? Is there any trigonometric function?...
    # See also https://github.com/sympy/sympy/pull/185.

    def shorter(*choices):
        '''Return the choice that has the fewest ops. In case of a tie,
        the expression listed first is selected.'''
        if not has_variety(choices):
            return choices[0]
        return min(choices, key=measure)

    expr = bottom_up(expr, lambda w: w.normal())
    expr = Mul(*powsimp(expr).as_content_primitive())
    _e = cancel(expr)
    expr1 = shorter(_e, _mexpand(_e).cancel())  # issue 6829
    expr2 = shorter(together(expr, deep=True), together(expr1, deep=True))

    if ratio is S.Infinity:
        expr = expr2
    else:
        expr = shorter(expr2, expr1, expr)
    if not isinstance(expr, Basic):  # XXX: temporary hack
        return expr

    expr = factor_terms(expr, sign=False)

    # hyperexpand automatically only works on hypergeometric terms
    expr = hyperexpand(expr)

    expr = piecewise_fold(expr)

    if expr.has(BesselBase):
        expr = besselsimp(expr)

    if expr.has(TrigonometricFunction) and not fu or expr.has(
            HyperbolicFunction):
        expr = trigsimp(expr, deep=True)

    if expr.has(log):
        expr = shorter(expand_log(expr, deep=True), logcombine(expr))

    if expr.has(CombinatorialFunction, gamma):
        expr = combsimp(expr)

    if expr.has(Sum):
        expr = sum_simplify(expr)

    if expr.has(Product):
        expr = product_simplify(expr)

    short = shorter(powsimp(expr, combine='exp', deep=True), powsimp(expr), expr)
    short = shorter(short, factor_terms(short), expand_power_exp(expand_mul(short)))
    if short.has(TrigonometricFunction, HyperbolicFunction, ExpBase):
        short = exptrigsimp(short, simplify=False)

    # get rid of hollow 2-arg Mul factorization
    hollow_mul = Transform(
        lambda x: Mul(*x.args),
        lambda x:
        x.is_Mul and
        len(x.args) == 2 and
        x.args[0].is_Number and
        x.args[1].is_Add and
        x.is_commutative)
    expr = short.xreplace(hollow_mul)

    numer, denom = expr.as_numer_denom()
    if denom.is_Add:
        n, d = fraction(radsimp(1/denom, symbolic=False, max_terms=1))
        if n is not S.One:
            expr = (numer*n).expand()/d

    if expr.could_extract_minus_sign():
        n, d = fraction(expr)
        if d != 0:
            expr = signsimp(-n/(-d))

    if measure(expr) > ratio*measure(original_expr):
        expr = original_expr

    return expr


def sum_simplify(s):
    """Main function for Sum simplification"""
    from sympy.concrete.summations import Sum
    from sympy.core.function import expand

    terms = Add.make_args(expand(s))
    s_t = [] # Sum Terms
    o_t = [] # Other Terms

    for term in terms:
        if isinstance(term, Mul):
            other = 1
            sum_terms = []

            if not term.has(Sum):
                o_t.append(term)
                continue

            mul_terms = Mul.make_args(term)
            for mul_term in mul_terms:
                if isinstance(mul_term, Sum):
                    r = mul_term._eval_simplify()
                    sum_terms.extend(Add.make_args(r))
                else:
                    other = other * mul_term
            if len(sum_terms):
                #some simplification may have happened
                #use if so
                s_t.append(Mul(*sum_terms) * other)
            else:
                o_t.append(other)
        elif isinstance(term, Sum):
            #as above, we need to turn this into an add list
            r = term._eval_simplify()
            s_t.extend(Add.make_args(r))
        else:
            o_t.append(term)


    result = Add(sum_combine(s_t), *o_t)

    return result

def sum_combine(s_t):
    """Helper function for Sum simplification

       Attempts to simplify a list of sums, by combining limits / sum function's
       returns the simplified sum
    """
    from sympy.concrete.summations import Sum


    used = [False] * len(s_t)

    for method in range(2):
        for i, s_term1 in enumerate(s_t):
            if not used[i]:
                for j, s_term2 in enumerate(s_t):
                    if not used[j] and i != j:
                        temp = sum_add(s_term1, s_term2, method)
                        if isinstance(temp, Sum) or isinstance(temp, Mul):
                            s_t[i] = temp
                            s_term1 = s_t[i]
                            used[j] = True

    result = S.Zero
    for i, s_term in enumerate(s_t):
        if not used[i]:
            result = Add(result, s_term)

    return result

def factor_sum(self, limits=None):
    """Helper function for Sum simplification

       if limits is specified, "self" is the inner part of a sum

       Returns the sum with constant factors brought outside
    """
    from sympy.core.exprtools import factor_terms
    from sympy.concrete.summations import Sum

    result = self.function if limits is None else self
    limits = self.limits if limits is None else limits
    #avoid any confusion w/ as_independent
    if result == 0:
        return S.Zero

    #get the summation variables
    sum_vars = set([limit.args[0] for limit in limits])

    #finally we try to factor out any common terms
    #and remove the from the sum if independent
    retv = factor_terms(result)
    #avoid doing anything bad
    if not result.is_commutative:
        return Sum(result, *limits)

    i, d = retv.as_independent(*sum_vars)
    if isinstance(retv, Add):
        return i * Sum(1, *limits) + Sum(d, *limits)
    else:
        return i * Sum(d, *limits)

def sum_add(self, other, method=0):
    """Helper function for Sum simplification"""
    from sympy.concrete.summations import Sum
    from sympy import Mul

    #we know this is something in terms of a constant * a sum
    #so we temporarily put the constants inside for simplification
    #then simplify the result
    def __refactor(val):
        args = Mul.make_args(val)
        sumv = next(x for x in args if isinstance(x, Sum))
        constant = Mul(*[x for x in args if x != sumv])
        return Sum(constant * sumv.function, *sumv.limits)

    if isinstance(self, Mul):
        rself = __refactor(self)
    else:
        rself = self

    if isinstance(other, Mul):
        rother = __refactor(other)
    else:
        rother = other

    if type(rself) == type(rother):
        if method == 0:
            if rself.limits == rother.limits:
                return factor_sum(Sum(rself.function + rother.function, *rself.limits))
        elif method == 1:
            if simplify(rself.function - rother.function) == 0:
                if len(rself.limits) == len(rother.limits) == 1:
                    i = rself.limits[0][0]
                    x1 = rself.limits[0][1]
                    y1 = rself.limits[0][2]
                    j = rother.limits[0][0]
                    x2 = rother.limits[0][1]
                    y2 = rother.limits[0][2]

                    if i == j:
                        if x2 == y1 + 1:
                            return factor_sum(Sum(rself.function, (i, x1, y2)))
                        elif x1 == y2 + 1:
                            return factor_sum(Sum(rself.function, (i, x2, y1)))

    return Add(self, other)


def product_simplify(s):
    """Main function for Product simplification"""
    from sympy.concrete.products import Product

    terms = Mul.make_args(s)
    p_t = [] # Product Terms
    o_t = [] # Other Terms

    for term in terms:
        if isinstance(term, Product):
            p_t.append(term)
        else:
            o_t.append(term)

    used = [False] * len(p_t)

    for method in range(2):
        for i, p_term1 in enumerate(p_t):
            if not used[i]:
                for j, p_term2 in enumerate(p_t):
                    if not used[j] and i != j:
                        if isinstance(product_mul(p_term1, p_term2, method), Product):
                            p_t[i] = product_mul(p_term1, p_term2, method)
                            used[j] = True

    result = Mul(*o_t)

    for i, p_term in enumerate(p_t):
        if not used[i]:
            result = Mul(result, p_term)

    return result


def product_mul(self, other, method=0):
    """Helper function for Product simplification"""
    from sympy.concrete.products import Product

    if type(self) == type(other):
        if method == 0:
            if self.limits == other.limits:
                return Product(self.function * other.function, *self.limits)
        elif method == 1:
            if simplify(self.function - other.function) == 0:
                if len(self.limits) == len(other.limits) == 1:
                    i = self.limits[0][0]
                    x1 = self.limits[0][1]
                    y1 = self.limits[0][2]
                    j = other.limits[0][0]
                    x2 = other.limits[0][1]
                    y2 = other.limits[0][2]

                    if i == j:
                        if x2 == y1 + 1:
                            return Product(self.function, (i, x1, y2))
                        elif x1 == y2 + 1:
                            return Product(self.function, (i, x2, y1))

    return Mul(self, other)


def _nthroot_solve(p, n, prec):
    """
     helper function for ``nthroot``
     It denests ``p**Rational(1, n)`` using its minimal polynomial
    """
    from sympy.polys.numberfields import _minimal_polynomial_sq
    from sympy.solvers import solve
    while n % 2 == 0:
        p = sqrtdenest(sqrt(p))
        n = n // 2
    if n == 1:
        return p
    pn = p**Rational(1, n)
    x = Symbol('x')
    f = _minimal_polynomial_sq(p, n, x)
    if f is None:
        return None
    sols = solve(f, x)
    for sol in sols:
        if abs(sol - pn).n() < 1./10**prec:
            sol = sqrtdenest(sol)
            if _mexpand(sol**n) == p:
                return sol


def logcombine(expr, force=False):
    """
    Takes logarithms and combines them using the following rules:

    - log(x) + log(y) == log(x*y) if both are not negative
    - a*log(x) == log(x**a) if x is positive and a is real

    If ``force`` is True then the assumptions above will be assumed to hold if
    there is no assumption already in place on a quantity. For example, if
    ``a`` is imaginary or the argument negative, force will not perform a
    combination but if ``a`` is a symbol with no assumptions the change will
    take place.

    Examples
    ========

    >>> from sympy import Symbol, symbols, log, logcombine, I
    >>> from sympy.abc import a, x, y, z
    >>> logcombine(a*log(x) + log(y) - log(z))
    a*log(x) + log(y) - log(z)
    >>> logcombine(a*log(x) + log(y) - log(z), force=True)
    log(x**a*y/z)
    >>> x,y,z = symbols('x,y,z', positive=True)
    >>> a = Symbol('a', real=True)
    >>> logcombine(a*log(x) + log(y) - log(z))
    log(x**a*y/z)

    The transformation is limited to factors and/or terms that
    contain logs, so the result depends on the initial state of
    expansion:

    >>> eq = (2 + 3*I)*log(x)
    >>> logcombine(eq, force=True) == eq
    True
    >>> logcombine(eq.expand(), force=True)
    log(x**2) + I*log(x**3)

    See Also
    ========
    posify: replace all symbols with symbols having positive assumptions

    """

    def f(rv):
        if not (rv.is_Add or rv.is_Mul):
            return rv

        def gooda(a):
            # bool to tell whether the leading ``a`` in ``a*log(x)``
            # could appear as log(x**a)
            return (a is not S.NegativeOne and  # -1 *could* go, but we disallow
                (a.is_real or force and a.is_real is not False))

        def goodlog(l):
            # bool to tell whether log ``l``'s argument can combine with others
            a = l.args[0]
            return a.is_positive or force and a.is_nonpositive is not False

        other = []
        logs = []
        log1 = defaultdict(list)
        for a in Add.make_args(rv):
            if a.func is log and goodlog(a):
                log1[()].append(([], a))
            elif not a.is_Mul:
                other.append(a)
            else:
                ot = []
                co = []
                lo = []
                for ai in a.args:
                    if ai.is_Rational and ai < 0:
                        ot.append(S.NegativeOne)
                        co.append(-ai)
                    elif ai.func is log and goodlog(ai):
                        lo.append(ai)
                    elif gooda(ai):
                        co.append(ai)
                    else:
                        ot.append(ai)
                if len(lo) > 1:
                    logs.append((ot, co, lo))
                elif lo:
                    log1[tuple(ot)].append((co, lo[0]))
                else:
                    other.append(a)

        # if there is only one log at each coefficient and none have
        # an exponent to place inside the log then there is nothing to do
        if not logs and all(len(log1[k]) == 1 and log1[k][0] == [] for k in log1):
            return rv

        # collapse multi-logs as far as possible in a canonical way
        # TODO: see if x*log(a)+x*log(a)*log(b) -> x*log(a)*(1+log(b))?
        # -- in this case, it's unambiguous, but if it were were a log(c) in
        # each term then it's arbitrary whether they are grouped by log(a) or
        # by log(c). So for now, just leave this alone; it's probably better to
        # let the user decide
        for o, e, l in logs:
            l = list(ordered(l))
            e = log(l.pop(0).args[0]**Mul(*e))
            while l:
                li = l.pop(0)
                e = log(li.args[0]**e)
            c, l = Mul(*o), e
            if l.func is log:  # it should be, but check to be sure
                log1[(c,)].append(([], l))
            else:
                other.append(c*l)

        # logs that have the same coefficient can multiply
        for k in list(log1.keys()):
            log1[Mul(*k)] = log(logcombine(Mul(*[
                l.args[0]**Mul(*c) for c, l in log1.pop(k)]),
                force=force))

        # logs that have oppositely signed coefficients can divide
        for k in ordered(list(log1.keys())):
            if not k in log1:  # already popped as -k
                continue
            if -k in log1:
                # figure out which has the minus sign; the one with
                # more op counts should be the one
                num, den = k, -k
                if num.count_ops() > den.count_ops():
                    num, den = den, num
                other.append(num*log(log1.pop(num).args[0]/log1.pop(den).args[0]))
            else:
                other.append(k*log1.pop(k))

        return Add(*other)

    return bottom_up(expr, f)


def bottom_up(rv, F, atoms=False, nonbasic=False):
    """Apply ``F`` to all expressions in an expression tree from the
    bottom up. If ``atoms`` is True, apply ``F`` even if there are no args;
    if ``nonbasic`` is True, try to apply ``F`` to non-Basic objects.
    """
    try:
        if rv.args:
            args = tuple([bottom_up(a, F, atoms, nonbasic)
                for a in rv.args])
            if args != rv.args:
                rv = rv.func(*args)
            rv = F(rv)
        elif atoms:
            rv = F(rv)
    except AttributeError:
        if nonbasic:
            try:
                rv = F(rv)
            except TypeError:
                pass

    return rv


def besselsimp(expr):
    """
    Simplify bessel-type functions.

    This routine tries to simplify bessel-type functions. Currently it only
    works on the Bessel J and I functions, however. It works by looking at all
    such functions in turn, and eliminating factors of "I" and "-1" (actually
    their polar equivalents) in front of the argument. Then, functions of
    half-integer order are rewritten using strigonometric functions and
    functions of integer order (> 1) are rewritten using functions
    of low order.  Finally, if the expression was changed, compute
    factorization of the result with factor().

    >>> from sympy import besselj, besseli, besselsimp, polar_lift, I, S
    >>> from sympy.abc import z, nu
    >>> besselsimp(besselj(nu, z*polar_lift(-1)))
    exp(I*pi*nu)*besselj(nu, z)
    >>> besselsimp(besseli(nu, z*polar_lift(-I)))
    exp(-I*pi*nu/2)*besselj(nu, z)
    >>> besselsimp(besseli(S(-1)/2, z))
    sqrt(2)*cosh(z)/(sqrt(pi)*sqrt(z))
    >>> besselsimp(z*besseli(0, z) + z*(besseli(2, z))/2 + besseli(1, z))
    3*z*besseli(0, z)/2
    """
    # TODO
    # - better algorithm?
    # - simplify (cos(pi*b)*besselj(b,z) - besselj(-b,z))/sin(pi*b) ...
    # - use contiguity relations?

    def replacer(fro, to, factors):
        factors = set(factors)

        def repl(nu, z):
            if factors.intersection(Mul.make_args(z)):
                return to(nu, z)
            return fro(nu, z)
        return repl

    def torewrite(fro, to):
        def tofunc(nu, z):
            return fro(nu, z).rewrite(to)
        return tofunc

    def tominus(fro):
        def tofunc(nu, z):
            return exp(I*pi*nu)*fro(nu, exp_polar(-I*pi)*z)
        return tofunc

    orig_expr = expr

    ifactors = [I, exp_polar(I*pi/2), exp_polar(-I*pi/2)]
    expr = expr.replace(
        besselj, replacer(besselj,
        torewrite(besselj, besseli), ifactors))
    expr = expr.replace(
        besseli, replacer(besseli,
        torewrite(besseli, besselj), ifactors))

    minusfactors = [-1, exp_polar(I*pi)]
    expr = expr.replace(
        besselj, replacer(besselj, tominus(besselj), minusfactors))
    expr = expr.replace(
        besseli, replacer(besseli, tominus(besseli), minusfactors))

    z0 = Dummy('z')

    def expander(fro):
        def repl(nu, z):
            if (nu % 1) == S(1)/2:
                return exptrigsimp(trigsimp(unpolarify(
                        fro(nu, z0).rewrite(besselj).rewrite(jn).expand(
                            func=True)).subs(z0, z)))
            elif nu.is_Integer and nu > 1:
                return fro(nu, z).expand(func=True)
            return fro(nu, z)
        return repl

    expr = expr.replace(besselj, expander(besselj))
    expr = expr.replace(bessely, expander(bessely))
    expr = expr.replace(besseli, expander(besseli))
    expr = expr.replace(besselk, expander(besselk))

    if expr != orig_expr:
        expr = expr.factor()

    return expr


def nthroot(expr, n, max_len=4, prec=15):
    """
    compute a real nth-root of a sum of surds

    Parameters
    ==========

    expr : sum of surds
    n : integer
    max_len : maximum number of surds passed as constants to ``nsimplify``

    Algorithm
    =========

    First ``nsimplify`` is used to get a candidate root; if it is not a
    root the minimal polynomial is computed; the answer is one of its
    roots.

    Examples
    ========

    >>> from sympy.simplify.simplify import nthroot
    >>> from sympy import Rational, sqrt
    >>> nthroot(90 + 34*sqrt(7), 3)
    sqrt(7) + 3

    """
    expr = sympify(expr)
    n = sympify(n)
    p = expr**Rational(1, n)
    if not n.is_integer:
        return p
    if not _is_sum_surds(expr):
        return p
    surds = []
    coeff_muls = [x.as_coeff_Mul() for x in expr.args]
    for x, y in coeff_muls:
        if not x.is_rational:
            return p
        if y is S.One:
            continue
        if not (y.is_Pow and y.exp == S.Half and y.base.is_integer):
            return p
        surds.append(y)
    surds.sort()
    surds = surds[:max_len]
    if expr < 0 and n % 2 == 1:
        p = (-expr)**Rational(1, n)
        a = nsimplify(p, constants=surds)
        res = a if _mexpand(a**n) == _mexpand(-expr) else p
        return -res
    a = nsimplify(p, constants=surds)
    if _mexpand(a) is not _mexpand(p) and _mexpand(a**n) == _mexpand(expr):
        return _mexpand(a)
    expr = _nthroot_solve(expr, n, prec)
    if expr is None:
        return p
    return expr


def nsimplify(expr, constants=[], tolerance=None, full=False, rational=None):
    """
    Find a simple representation for a number or, if there are free symbols or
    if rational=True, then replace Floats with their Rational equivalents. If
    no change is made and rational is not False then Floats will at least be
    converted to Rationals.

    For numerical expressions, a simple formula that numerically matches the
    given numerical expression is sought (and the input should be possible
    to evalf to a precision of at least 30 digits).

    Optionally, a list of (rationally independent) constants to
    include in the formula may be given.

    A lower tolerance may be set to find less exact matches. If no tolerance
    is given then the least precise value will set the tolerance (e.g. Floats
    default to 15 digits of precision, so would be tolerance=10**-15).

    With full=True, a more extensive search is performed
    (this is useful to find simpler numbers when the tolerance
    is set low).

    Examples
    ========

    >>> from sympy import nsimplify, sqrt, GoldenRatio, exp, I, exp, pi
    >>> nsimplify(4/(1+sqrt(5)), [GoldenRatio])
    -2 + 2*GoldenRatio
    >>> nsimplify((1/(exp(3*pi*I/5)+1)))
    1/2 - I*sqrt(sqrt(5)/10 + 1/4)
    >>> nsimplify(I**I, [pi])
    exp(-pi/2)
    >>> nsimplify(pi, tolerance=0.01)
    22/7

    See Also
    ========
    sympy.core.function.nfloat

    """
    try:
        return sympify(as_int(expr))
    except (TypeError, ValueError):
        pass
    expr = sympify(expr).xreplace({
        Float('inf'): S.Infinity,
        Float('-inf'): S.NegativeInfinity,
        })
    if expr is S.Infinity or expr is S.NegativeInfinity:
        return expr
    if rational or expr.free_symbols:
        return _real_to_rational(expr, tolerance)

    # SymPy's default tolerance for Rationals is 15; other numbers may have
    # lower tolerances set, so use them to pick the largest tolerance if None
    # was given
    if tolerance is None:
        tolerance = 10**-min([15] +
             [mpmath.libmp.libmpf.prec_to_dps(n._prec)
             for n in expr.atoms(Float)])
    # XXX should prec be set independent of tolerance or should it be computed
    # from tolerance?
    prec = 30
    bprec = int(prec*3.33)

    constants_dict = {}
    for constant in constants:
        constant = sympify(constant)
        v = constant.evalf(prec)
        if not v.is_Float:
            raise ValueError("constants must be real-valued")
        constants_dict[str(constant)] = v._to_mpmath(bprec)

    exprval = expr.evalf(prec, chop=True)
    re, im = exprval.as_real_imag()

    # safety check to make sure that this evaluated to a number
    if not (re.is_Number and im.is_Number):
        return expr

    def nsimplify_real(x):
        orig = mpmath.mp.dps
        xv = x._to_mpmath(bprec)
        try:
            # We'll be happy with low precision if a simple fraction
            if not (tolerance or full):
                mpmath.mp.dps = 15
                rat = mpmath.pslq([xv, 1])
                if rat is not None:
                    return Rational(-int(rat[1]), int(rat[0]))
            mpmath.mp.dps = prec
            newexpr = mpmath.identify(xv, constants=constants_dict,
                tol=tolerance, full=full)
            if not newexpr:
                raise ValueError
            if full:
                newexpr = newexpr[0]
            expr = sympify(newexpr)
            if x and not expr:  # don't let x become 0
                raise ValueError
            if expr.is_finite is False and not xv in [mpmath.inf, mpmath.ninf]:
                raise ValueError
            return expr
        finally:
            # even though there are returns above, this is executed
            # before leaving
            mpmath.mp.dps = orig
    try:
        if re:
            re = nsimplify_real(re)
        if im:
            im = nsimplify_real(im)
    except ValueError:
        if rational is None:
            return _real_to_rational(expr)
        return expr

    rv = re + im*S.ImaginaryUnit
    # if there was a change or rational is explicitly not wanted
    # return the value, else return the Rational representation
    if rv != expr or rational is False:
        return rv
    return _real_to_rational(expr)


def _real_to_rational(expr, tolerance=None):
    """
    Replace all reals in expr with rationals.

    >>> from sympy import nsimplify
    >>> from sympy.abc import x

    >>> nsimplify(.76 + .1*x**.5, rational=True)
    sqrt(x)/10 + 19/25

    """
    inf = Float('inf')
    p = expr
    reps = {}
    reduce_num = None
    if tolerance is not None and tolerance < 1:
        reduce_num = ceiling(1/tolerance)
    for float in p.atoms(Float):
        key = float
        if reduce_num is not None:
            r = Rational(float).limit_denominator(reduce_num)
        elif (tolerance is not None and tolerance >= 1 and
                float.is_Integer is False):
            r = Rational(tolerance*round(float/tolerance)
                ).limit_denominator(int(tolerance))
        else:
            r = nsimplify(float, rational=False)
            # e.g. log(3).n() -> log(3) instead of a Rational
            if float and not r:
                r = Rational(float)
            elif not r.is_Rational:
                if float == inf or float == -inf:
                    r = S.ComplexInfinity
                elif float < 0:
                    float = -float
                    d = Pow(10, int((mpmath.log(float)/mpmath.log(10))))
                    r = -Rational(str(float/d))*d
                elif float > 0:
                    d = Pow(10, int((mpmath.log(float)/mpmath.log(10))))
                    r = Rational(str(float/d))*d
                else:
                    r = Integer(0)
        reps[key] = r
    return p.subs(reps, simultaneous=True)


def clear_coefficients(expr, rhs=S.Zero):
    """Return `p, r` where `p` is the expression obtained when Rational
    additive and multiplicative coefficients of `expr` have been stripped
    away in a naive fashion (i.e. without simplification). The operations
    needed to remove the coefficients will be applied to `rhs` and returned
    as `r`.

    Examples
    ========

    >>> from sympy.simplify.simplify import clear_coefficients
    >>> from sympy.abc import x, y
    >>> from sympy import Dummy
    >>> expr = 4*y*(6*x + 3)
    >>> clear_coefficients(expr - 2)
    (y*(2*x + 1), 1/6)

    When solving 2 or more expressions like `expr = a`,
    `expr = b`, etc..., it is advantageous to provide a Dummy symbol
    for `rhs` and  simply replace it with `a`, `b`, etc... in `r`.

    >>> rhs = Dummy('rhs')
    >>> clear_coefficients(expr, rhs)
    (y*(2*x + 1), _rhs/12)
    >>> _[1].subs(rhs, 2)
    1/6
    """
    was = None
    free = expr.free_symbols
    if expr.is_Rational:
        return (S.Zero, rhs - expr)
    while expr and was != expr:
        was = expr
        m, expr = (
            expr.as_content_primitive()
            if free else
            factor_terms(expr).as_coeff_Mul(rational=True))
        rhs /= m
        c, expr = expr.as_coeff_Add(rational=True)
        rhs -= c
    expr = signsimp(expr, evaluate = False)
    if _coeff_isneg(expr):
        expr = -expr
        rhs = -rhs
    return expr, rhs
