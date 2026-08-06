"""
The Risch Algorithm for transcendental function integration.

The core algorithms for the Risch algorithm are here.  The subproblem
algorithms are in the rde.py and prde.py files for the Risch
Differential Equation solver and the parametric problems solvers,
respectively.  All important information concerning the differential extension
for an integrand is stored in a DifferentialExtension object, which in the code
is usually called DE.  Throughout the code and Inside the DifferentialExtension
object, the conventions/attribute names are that the base domain is QQ and each
differential extension is x, t0, t1, ..., tn-1 = DE.t. DE.x is the variable of
integration (Dx == 1), DE.D is a list of the derivatives of
x, t1, t2, ..., tn-1 = t, DE.T is the list [x, t1, t2, ..., tn-1], DE.t is the
outer-most variable of the differential extension at the given level (the level
can be adjusted using DE.increment_level() and DE.decrement_level()),
k is the field C(x, t0, ..., tn-2), where C is the constant field.  The
numerator of a fraction is denoted by a and the denominator by
d.  If the fraction is named f, fa == numer(f) and fd == denom(f).
Fractions are returned as tuples (fa, fd).  DE.d and DE.t are used to
represent the topmost derivation and extension variable, respectively.
The docstring of a function signifies whether an argument is in k[t], in
which case it will just return a Poly in t, or in k(t), in which case it
will return the fraction (fa, fd). Other variable names probably come
from the names used in Bronstein's book.
"""
from __future__ import annotations
from types import GeneratorType
from functools import reduce

from sympy.core.add import Add
from sympy.core.function import Lambda
from sympy.core.mul import Mul
from sympy.core.intfunc import ilcm
from sympy.core.numbers import I, pi
from sympy.core.power import Pow
from sympy.core.relational import Ne
from sympy.core.singleton import S
from sympy.core.sorting import ordered, default_sort_key
from sympy.core.symbol import Dummy, Symbol
from sympy.functions.combinatorial.factorials import binomial
from sympy.functions.elementary.exponential import log, exp
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.special.error_functions import (Ei, erf, erfi, li)
from sympy.functions.elementary.hyperbolic import (cosh, coth, sinh,
    tanh)
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.elementary.trigonometric import (atan, sin, cos,
    tan, acot, cot, asin, acos, sec, csc)
from .integrals import integrate, Integral
from .heurisch import _symbols
from sympy.polys.polyerrors import (DomainError, ExactQuotientFailed,
    GeneratorsError, PolynomialError)
from sympy.polys.polytools import (real_roots, cancel, Poly, gcd,
    reduced)
from sympy.polys.rootoftools import RootSum
from sympy.utilities.iterables import numbered_symbols


def integer_powers(exprs):
    """
    Rewrites a list of expressions as integer multiples of each other.

    Explanation
    ===========

    For example, if you have [x, x/2, x**2 + 1, 2*x/3], then you can rewrite
    this as [(x/6) * 6, (x/6) * 3, (x**2 + 1) * 1, (x/6) * 4]. This is useful
    in the Risch integration algorithm, where we must write exp(x) + exp(x/2)
    as (exp(x/2))**2 + exp(x/2), but not as exp(x) + sqrt(exp(x)) (this is
    because only the transcendental case is implemented and we therefore cannot
    integrate algebraic extensions). The integer multiples returned by this
    function for each term are the smallest possible (their content equals 1).

    Returns a list of tuples where the first element is the base term and the
    second element is a list of `(item, factor)` terms, where `factor` is the
    integer multiplicative factor that must multiply the base term to obtain
    the original item.

    The easiest way to understand this is to look at an example:

    >>> from sympy.abc import x
    >>> from sympy.integrals.risch import integer_powers
    >>> integer_powers([x, x/2, x**2 + 1, 2*x/3])
    [(x/6, [(x, 6), (x/2, 3), (2*x/3, 4)]), (x**2 + 1, [(x**2 + 1, 1)])]

    We can see how this relates to the example at the beginning of the
    docstring.  It chose x/6 as the first base term.  Then, x can be written as
    (x/2) * 2, so we get (0, 2), and so on. Now only element (x**2 + 1)
    remains, and there are no other terms that can be written as a rational
    multiple of that, so we get that it can be written as (x**2 + 1) * 1.

    """
    # Here is the strategy:

    # First, go through each term and determine if it can be rewritten as a
    # rational multiple of any of the terms gathered so far.
    # cancel(a/b).is_Rational is sufficient for this.  If it is a multiple, we
    # add its multiple to the dictionary.

    terms = {}
    for term in exprs:
        for trm, trm_list in terms.items():
            a = cancel(term/trm)
            if a.is_Rational:
                trm_list.append((term, a))
                break
        else:
            terms[term] = [(term, S.One)]

    # After we have done this, we have all the like terms together, so we just
    # need to find a common denominator so that we can get the base term and
    # integer multiples such that each term can be written as an integer
    # multiple of the base term, and the content of the integers is 1.

    newterms = {}
    for term, term_list in terms.items():
        common_denom = reduce(ilcm, [i.as_numer_denom()[1] for _, i in
            term_list])
        newterm = term/common_denom
        newmults = [(i, j*common_denom) for i, j in term_list]
        newterms[newterm] = newmults

    return sorted(iter(newterms.items()), key=lambda item: item[0].sort_key())


def _tan_of_multiple(n, t):
    """
    Rewrite tan(n*a) as a rational function of ``t`` == tan(a).

    Explanation
    ===========

    For an integer n, tan(n*a) == Im((1 + I*t)**n)/Re((1 + I*t)**n), that
    is::

        tan(n*a) == (binomial(n, 1)*t - binomial(n, 3)*t**3 + ...)/
            (1 - binomial(n, 2)*t**2 + ...)

    This is used to rewrite the tangents of integer multiples of a common
    argument in terms of a single tangent, like exp of integer multiples
    of a common argument is rewritten as powers of a single exponential.
    """
    n = int(n)
    if n < 0:
        return -_tan_of_multiple(-n, t)
    num = Add(*[(-1)**k*binomial(n, 2*k + 1)*t**(2*k + 1)
        for k in range((n + 1)//2)])
    den = Add(*[(-1)**k*binomial(n, 2*k)*t**(2*k)
        for k in range(n//2 + 1)])
    return num/den


def _tan_of_combination(terms):
    """
    Rewrite tan(Sum(n_i*a_i)) as a rational function of t_i == tan(a_i).

    ``terms`` is a list of (n_i, t_i) pairs with integer n_i; the multiple
    angle formula and the addition formula
    tan(a + b) == (tan(a) + tan(b))/(1 - tan(a)*tan(b)) are applied.
    """
    vals = [_tan_of_multiple(n, t) for n, t in terms if n]
    if not vals:
        return S.Zero
    while len(vals) > 1:
        ta, tb = vals.pop(), vals.pop()
        vals.append(cancel((ta + tb)/(1 - ta*tb)))
    return vals[0]


class DifferentialExtension:
    """
    A container for all the information relating to a differential extension.

    Explanation
    ===========

    The attributes of this object are (see also the docstring of __init__):

    - f: The original (Expr) integrand.
    - x: The variable of integration.
    - T: List of variables in the extension.
    - D: List of derivations in the extension; corresponds to the elements of T.
    - fa: Poly of the numerator of the integrand.
    - fd: Poly of the denominator of the integrand.
    - Tfuncs: Lambda() representations of each element of T (except for x).
      For back-substitution after integration.
    - backsubs: A (possibly empty) list of further substitutions to be made on
      the final integral to make it look more like the integrand.
    - exts:
    - extargs:
    - cases: List of string representations of the cases of T.
    - t: The top level extension variable, as defined by the current level
      (see level below).
    - d: The top level extension derivation, as defined by the current
      derivation (see level below).
    - case: The string representation of the case of self.d.
    (Note that self.T and self.D will always contain the complete extension,
    regardless of the level.  Therefore, you should ALWAYS use DE.t and DE.d
    instead of DE.T[-1] and DE.D[-1].  If you want to have a list of the
    derivations or variables only up to the current level, use
    DE.D[:len(DE.D) + DE.level + 1] and DE.T[:len(DE.T) + DE.level + 1].  Note
    that, in particular, the derivation() function does this.)

    The following are also attributes, but will probably not be useful other
    than in internal use:
    - newf: Expr form of fa/fd.
    - level: The number (between -1 and -len(self.T)) such that
      self.T[self.level] == self.t and self.D[self.level] == self.d.
      Use the methods self.increment_level() and self.decrement_level() to change
      the current level.
    """
    # __slots__ is defined mainly so we can iterate over all the attributes
    # of the class easily (the memory use doesn't matter too much, since we
    # only create one DifferentialExtension per integration).  Also, it's nice
    # to have a safeguard when debugging.
    __slots__ = ('f', 'x', 'T', 'D', 'fa', 'fd', 'Tfuncs', 'backsubs',
        'exts', 'extargs', 'cases', 'case', 't', 'd', 'newf', 'level',
        'ts', 'dummy')

    def __init__(self, f=None, x=None, handle_first='log', dummy=False,
            extension=None, rewrite_complex=None, trig=True):
        """
        Tries to build a transcendental extension tower from ``f`` with respect to ``x``.

        Explanation
        ===========

        If it is successful, creates a DifferentialExtension object with, among
        others, the attributes fa, fd, D, T, Tfuncs, and backsubs such that
        fa and fd are Polys in T[-1] with rational coefficients in T[:-1],
        fa/fd == f, and D[i] is a Poly in T[i] with rational coefficients in
        T[:i] representing the derivative of T[i] for each i from 1 to len(T).
        Tfuncs is a list of Lambda objects for back replacing the functions
        after integrating.  Lambda() is only used (instead of lambda) to make
        them easier to test and debug. Note that Tfuncs corresponds to the
        elements of T, except for T[0] == x, but they should be back-substituted
        in reverse order.  backsubs is a (possibly empty) back-substitution list
        that should be applied on the completed integral to make it look more
        like the original integrand.

        If it is unsuccessful, it raises NotImplementedError.

        If ``trig=True`` (the default), tangent extensions are built for
        the trigonometric functions (with sin, cos, sec and csc rewritten
        as rational functions of tan of the half angle) and inverse
        tangent extensions for atan and acot; ``trig=False`` restores the
        old behavior of raising NotImplementedError for trigonometric
        integrands.  In both cases, if the integrand contains I or
        ``rewrite_complex=True``, trigonometric functions are rewritten in
        terms of complex exponentials and logarithms instead.

        You can also create an object by manually setting the attributes as a
        dictionary to the extension keyword argument.  You must include at least
        D.  Warning, any attribute that is not given will be set to None. The
        attributes T, t, d, cases, case, x, and level are set automatically and
        do not need to be given.  The functions in the Risch Algorithm will NOT
        check to see if an attribute is None before using it.  This also does not
        check to see if the extension is valid (non-algebraic) or even if it is
        self-consistent.  Therefore, this should only be used for
        testing/debugging purposes.
        """
        # XXX: If you need to debug this function, set the break point here

        if extension:
            if 'D' not in extension:
                raise ValueError("At least the key D must be included with "
                    "the extension flag to DifferentialExtension.")
            for attr in extension:
                setattr(self, attr, extension[attr])

            self._auto_attrs()

            return
        elif f is None or x is None:
            raise ValueError("Either both f and x or a manual extension must "
            "be given.")

        if handle_first not in ('log', 'exp'):
            raise ValueError("handle_first must be 'log' or 'exp', not %s." %
                str(handle_first))

        # f will be the original function, self.f might change if we reset
        # (e.g., we pull out a constant from an exponential)
        self.f = f
        self.x = x
        # setting the default value 'dummy'
        self.dummy = dummy
        self.reset()
        exp_new_extension, log_new_extension = True, True
        tan_new_extension = bool(trig)

        # case of 'automatic' choosing
        if rewrite_complex is None:
            rewrite_complex = I in self.f.atoms()

        if rewrite_complex:
            rewritables = {
                (sin, cos, cot, tan, sinh, cosh, coth, tanh): exp,
                (asin, acos, acot, atan): log,
            }
            # rewrite the trigonometric components
            for candidates, rule in rewritables.items():
                self.newf = self.newf.rewrite(candidates, rule)
            self.newf = cancel(self.newf)
        elif not trig:
            if any(i.has(x) for i in self.f.atoms(sin, cos, cot, tan, sinh,
                    cosh, coth, tanh, asin, acos, acot, atan)):
                raise NotImplementedError("Trigonometric and hyperbolic "
                    "extensions are not supported (yet!).  Try rewriting in "
                    "terms of exp and log, using rewrite_complex=True, or "
                    "using trig=True for the real trigonometric functions.")
        else:
            if any(i.has(x) for i in self.f.atoms(asin, acos)):
                raise NotImplementedError("The inverse trigonometric "
                "extensions asin and acos are algebraic and are not "
                "supported.")

        exps = set()
        pows = set()
        numpows = set()
        sympows = set()
        logs = set()
        symlogs = set()
        tans = set()
        atans = set()

        while True:
            if self.newf.is_rational_function(*self.T):
                break

            if not exp_new_extension and not log_new_extension and \
                    not tan_new_extension:
                # We couldn't find a new extension on the last pass, so I guess
                # we can't do it.
                raise NotImplementedError("Couldn't find an elementary "
                    "transcendental extension for %s.  Try using a " % str(f) +
                    "manual extension with the extension flag.")

            exps, pows, numpows, sympows, log_new_extension = \
                    self._rewrite_exps_pows(exps, pows, numpows, sympows, log_new_extension)

            logs, symlogs = self._rewrite_logs(logs, symlogs)

            if trig:
                tans = self._rewrite_tans(tans)
                # acot(b) == atan(1/b) (with sympy's conventions)
                self.newf = self.newf.xreplace({i: i.rewrite(atan) for i in
                    self.newf.atoms(acot)
                    if i.args[0].is_rational_function(*self.T) and
                    i.args[0].has(*self.T)})
                atans = update_sets(atans, self.newf.atoms(atan),
                    lambda i: i.args[0].is_rational_function(*self.T) and
                    i.args[0].has(*self.T))

            if handle_first == 'exp' or not log_new_extension:
                exp_new_extension = self._exp_part(exps)
                if exp_new_extension is None:
                    # reset and restart
                    self.f = self.newf
                    self.reset()
                    exp_new_extension = True
                    continue

            if handle_first == 'log' or not exp_new_extension:
                log_new_extension = self._log_part(logs)

            if trig:
                if exp_new_extension or log_new_extension:
                    # Prefer building the exponential and logarithmic
                    # extensions lower in the tower, so that the tangents
                    # end up on top, where integrate_hypertangent() can
                    # handle them.
                    tan_new_extension = True
                else:
                    tan_new_extension = self._atan_part(atans)
                    tan_new_extension = self._tan_part(tans) or \
                        tan_new_extension

        self.fa, self.fd = frac_in(self.newf, self.t)
        self._auto_attrs()

        return

    def __getattr__(self, attr):
        # Avoid AttributeErrors when debugging
        if attr not in self.__slots__:
            raise AttributeError("%s has no attribute %s" % (repr(self), repr(attr)))
        return None

    def _rewrite_exps_pows(self, exps, pows, numpows,
            sympows, log_new_extension):
        """
        Rewrite exps/pows for better processing.
        """
        from .prde import is_deriv_k

        # Pre-preparsing.
        #################
        # Get all exp arguments, so we can avoid ahead of time doing
        # something like t1 = exp(x), t2 = exp(x/2) == sqrt(t1).

        # Things like sqrt(exp(x)) do not automatically simplify to
        # exp(x/2), so they will be viewed as algebraic.  The easiest way
        # to handle this is to convert all instances of exp(a)**Rational
        # to exp(Rational*a) before doing anything else.  Note that the
        # _exp_part code can generate terms of this form, so we do need to
        # do this at each pass (or else modify it to not do that).

        ratpows = [i for i in self.newf.atoms(Pow)
                   if (isinstance(i.base, exp) and i.exp.is_Rational)]

        ratpows_repl = [
            (i, i.base.base**(i.exp*i.base.exp)) for i in ratpows]
        self.backsubs += [(j, i) for i, j in ratpows_repl]
        self.newf = self.newf.xreplace(dict(ratpows_repl))

        # To make the process deterministic, the args are sorted
        # so that functions with smaller op-counts are processed first.
        # Ties are broken with the default_sort_key.

        # XXX Although the method is deterministic no additional work
        # has been done to guarantee that the simplest solution is
        # returned and that it would be affected be using different
        # variables. Though it is possible that this is the case
        # one should know that it has not been done intentionally, so
        # further improvements may be possible.

        # TODO: This probably doesn't need to be completely recomputed at
        # each pass.
        exps = update_sets(exps, self.newf.atoms(exp),
            lambda i: i.exp.is_rational_function(*self.T) and
            i.exp.has(*self.T))
        pows = update_sets(pows, self.newf.atoms(Pow),
            lambda i: i.exp.is_rational_function(*self.T) and
            i.exp.has(*self.T))
        numpows = update_sets(numpows, set(pows),
            lambda i: not i.base.has(*self.T))
        sympows = update_sets(sympows, set(pows) - set(numpows),
            lambda i: i.base.is_rational_function(*self.T) and
            not i.exp.is_Integer)

        # The easiest way to deal with non-base E powers is to convert them
        # into base E, integrate, and then convert back.
        for i in ordered(pows):
            old = i
            new = exp(i.exp*log(i.base))
            # If exp is ever changed to automatically reduce exp(x*log(2))
            # to 2**x, then this will break.  The solution is to not change
            # exp to do that :)
            if i in sympows:
                if i.exp.is_Rational:
                    raise NotImplementedError("Algebraic extensions are "
                        "not supported (%s)." % str(i))
                # We can add a**b only if log(a) in the extension, because
                # a**b == exp(b*log(a)).
                basea, based = frac_in(i.base, self.t)
                A = is_deriv_k(basea, based, self)
                if A is None:
                    # Nonelementary monomial (so far)

                    # TODO: Would there ever be any benefit from just
                    # adding log(base) as a new monomial?
                    # ANSWER: Yes, otherwise we can't integrate x**x (or
                    # rather prove that it has no elementary integral)
                    # without first manually rewriting it as exp(x*log(x))
                    self.newf = self.newf.xreplace({old: new})
                    self.backsubs += [(new, old)]
                    log_new_extension = self._log_part([log(i.base)])
                    exps = update_sets(exps, self.newf.atoms(exp), lambda i:
                        i.exp.is_rational_function(*self.T) and i.exp.has(*self.T))
                    continue
                ans, u, const = A
                newterm = exp(i.exp*(log(const) + u))
                # Under the current implementation, exp kills terms
                # only if they are of the form a*log(x), where a is a
                # Number.  This case should have already been killed by the
                # above tests.  Again, if this changes to kill more than
                # that, this will break, which maybe is a sign that you
                # shouldn't be changing that.  Actually, if anything, this
                # auto-simplification should be removed.  See
                # https://groups.google.com/group/sympy/browse_thread/thread/a61d48235f16867f

                self.newf = self.newf.xreplace({i: newterm})

            elif i not in numpows:
                continue
            else:
                # i in numpows
                newterm = new
            # TODO: Just put it in self.Tfuncs
            self.backsubs.append((new, old))
            self.newf = self.newf.xreplace({old: newterm})
            exps.append(newterm)

        return exps, pows, numpows, sympows, log_new_extension

    def _rewrite_tans(self, tans):
        """
        Rewrite trigonometric functions in terms of tan for better processing.

        Explanation
        ===========

        sin, cos, sec and csc of ``a`` are rational functions of tan(a/2)
        (the Weierstrass substitution), and cot(a) == 1/tan(a), so they are
        first rewritten in terms of tangents.  Then, like in the
        exponential case, the arguments of the tangents are grouped with
        integer_powers(), and tan of an integer multiple of the base
        argument is rewritten as a rational function of tan of the base
        argument, so that e.g. tan(x) and sin(x) can both be expressed in
        the field generated by tan(x/2).
        """
        trigs = [i for i in self.newf.atoms(sin, cos, sec, csc, cot)
            if i.args[0].is_rational_function(*self.T) and
            i.args[0].has(*self.T)]
        for i in ordered(trigs):
            self.newf = self.newf.xreplace({i: i.rewrite(tan)})

        atoms = self.newf.atoms(tan)
        tans = update_sets(tans, atoms,
            lambda i: i.args[0].is_rational_function(*self.T) and
            i.args[0].has(*self.T))

        ip = integer_powers([i.args[0] for i in tans])
        for base, others in ip:
            self.newf = self.newf.xreplace({tan(arg):
                _tan_of_multiple(n, tan(base)) for arg, n in others if n != 1})

        # Recollect, since the tangents of the multiples were rewritten.
        tans = update_sets([], self.newf.atoms(tan),
            lambda i: i.args[0].is_rational_function(*self.T) and
            i.args[0].has(*self.T))

        return tans

    def _rewrite_logs(self, logs, symlogs):
        """
        Rewrite logs for better processing.
        """
        atoms = self.newf.atoms(log)
        logs = update_sets(logs, atoms,
            lambda i: i.args[0].is_rational_function(*self.T) and
            i.args[0].has(*self.T))
        symlogs = update_sets(symlogs, atoms,
            lambda i: i.has(*self.T) and i.args[0].is_Pow and
            i.args[0].base.is_rational_function(*self.T) and
            not i.args[0].exp.is_Integer)

        # We can handle things like log(x**y) by converting it to y*log(x)
        # This will fix not only symbolic exponents of the argument, but any
        # non-Integer exponent, like log(sqrt(x)).  The exponent can also
        # depend on x, like log(x**x).
        for i in ordered(symlogs):
            # Unlike in the exponential case above, we do not ever
            # potentially add new monomials (above we had to add log(a)).
            # Therefore, there is no need to run any is_deriv functions
            # here.  Just convert log(a**b) to b*log(a) and let
            # log_new_extension() handle it from there.
            lbase = log(i.args[0].base)
            logs.append(lbase)
            new = i.args[0].exp*lbase
            self.newf = self.newf.xreplace({i: new})
            self.backsubs.append((new, i))

        # remove any duplicates
        logs = sorted(set(logs), key=default_sort_key)

        return logs, symlogs

    def _auto_attrs(self):
        """
        Set attributes that are generated automatically.
        """
        if not self.T:
            # i.e., when using the extension flag and T isn't given
            self.T = [i.gen for i in self.D]
        if not self.x:
            self.x = self.T[0]
        self.cases = [get_case(d, t) for d, t in zip(self.D, self.T)]
        self.level = -1
        self.t = self.T[self.level]
        self.d = self.D[self.level]
        self.case = self.cases[self.level]

    def _exp_part(self, exps):
        """
        Try to build an exponential extension.

        Returns
        =======

        Returns True if there was a new extension, False if there was no new
        extension but it was able to rewrite the given exponentials in terms
        of the existing extension, and None if the entire extension building
        process should be restarted.  If the process fails because there is no
        way around an algebraic extension (e.g., exp(log(x)/2)), it will raise
        NotImplementedError.
        """
        from .prde import is_log_deriv_k_t_radical
        new_extension = False
        restart = False
        expargs = [i.exp for i in exps]
        ip = integer_powers(expargs)
        for arg, others in ip:
            # Minimize potential problems with algebraic substitution
            others.sort(key=lambda i: i[1])

            arga, argd = frac_in(arg, self.t)
            A = is_log_deriv_k_t_radical(arga, argd, self)

            if A is not None:
                ans, u, n, const = A
                # if n is 1 or -1, it's algebraic, but we can handle it
                if n == -1:
                    # This probably will never happen, because
                    # Rational.as_numer_denom() returns the negative term in
                    # the numerator.  But in case that changes, reduce it to
                    # n == 1.
                    n = 1
                    u **= -1
                    const *= -1
                    ans = [(i, -j) for i, j in ans]

                if n == 1:
                    # Example: exp(x + x**2) over QQ(x, exp(x), exp(x**2))
                    self.newf = self.newf.xreplace({exp(arg): exp(const)*Mul(*[
                        u**power for u, power in ans])})
                    self.newf = self.newf.xreplace({exp(p*exparg):
                        exp(const*p) * Mul(*[u**power for u, power in ans])
                        for exparg, p in others})
                    # TODO: Add something to backsubs to put exp(const*p)
                    # back together.

                    continue

                else:
                    # Bad news: we have an algebraic radical.  But maybe we
                    # could still avoid it by choosing a different extension.
                    # For example, integer_powers() won't handle exp(x/2 + 1)
                    # over QQ(x, exp(x)), but if we pull out the exp(1), it
                    # will.  Or maybe we have exp(x + x**2/2), over
                    # QQ(x, exp(x), exp(x**2)), which is exp(x)*sqrt(exp(x**2)),
                    # but if we use QQ(x, exp(x), exp(x**2/2)), then they will
                    # all work.
                    #
                    # So here is what we do: If there is a non-zero const, pull
                    # it out and retry.  Also, if len(ans) > 1, then rewrite
                    # exp(arg) as the product of exponentials from ans, and
                    # retry that.  If const == 0 and len(ans) == 1, then we
                    # assume that it would have been handled by either
                    # integer_powers() or n == 1 above if it could be handled,
                    # so we give up at that point.  For example, you can never
                    # handle exp(log(x)/2) because it equals sqrt(x).

                    if const or len(ans) > 1:
                        rad = Mul(*[term**(power/n) for term, power in ans])
                        self.newf = self.newf.xreplace({exp(p*exparg):
                            exp(const*p)*rad for exparg, p in others})
                        self.newf = self.newf.xreplace(dict(list(zip(reversed(self.T),
                            reversed([f(self.x) for f in self.Tfuncs])))))
                        restart = True
                        break
                    else:
                        # TODO: give algebraic dependence in error string
                        raise NotImplementedError("Cannot integrate over "
                            "algebraic extensions.")

            else:
                arga, argd = frac_in(arg, self.t)
                darga = (argd*derivation(Poly(arga, self.t), self) -
                    arga*derivation(Poly(argd, self.t), self))
                dargd = argd**2
                darga, dargd = darga.cancel(dargd, include=True)
                darg = darga.as_expr()/dargd.as_expr()
                self.t = next(self.ts)
                self.T.append(self.t)
                self.extargs.append(arg)
                self.exts.append('exp')
                self.D.append(darg.as_poly(self.t, expand=False)*Poly(self.t,
                    self.t, expand=False))
                if self.dummy:
                    i = Dummy("i")
                else:
                    i = Symbol('i')
                self.Tfuncs += [Lambda(i, exp(arg.subs(self.x, i)))]
                self.newf = self.newf.xreplace(
                        {exp(exparg): self.t**p for exparg, p in others})
                new_extension = True

        if restart:
            return None
        return new_extension

    def _log_part(self, logs):
        """
        Try to build a logarithmic extension.

        Returns
        =======

        Returns True if there was a new extension and False if there was no new
        extension but it was able to rewrite the given logarithms in terms
        of the existing extension.  Unlike with exponential extensions, there
        is no way that a logarithm is not transcendental over and cannot be
        rewritten in terms of an already existing extension in a non-algebraic
        way, so this function does not ever return None or raise
        NotImplementedError.
        """
        from .prde import is_deriv_k
        new_extension = False
        logargs = [i.args[0] for i in logs]
        for arg in ordered(logargs):
            # The log case is easier, because whenever a logarithm is algebraic
            # over the base field, it is of the form a1*t1 + ... an*tn + c,
            # which is a polynomial, so we can just replace it with that.
            # In other words, we don't have to worry about radicals.
            arga, argd = frac_in(arg, self.t)
            A = is_deriv_k(arga, argd, self)
            if A is not None:
                ans, u, const = A
                newterm = log(const) + u
                self.newf = self.newf.xreplace({log(arg): newterm})
                continue

            else:
                arga, argd = frac_in(arg, self.t)
                darga = (argd*derivation(Poly(arga, self.t), self) -
                    arga*derivation(Poly(argd, self.t), self))
                dargd = argd**2
                darg = darga.as_expr()/dargd.as_expr()
                self.t = next(self.ts)
                self.T.append(self.t)
                self.extargs.append(arg)
                self.exts.append('log')
                self.D.append(cancel(darg.as_expr()/arg).as_poly(self.t,
                    expand=False))
                if self.dummy:
                    i = Dummy("i")
                else:
                    i = Symbol('i')
                self.Tfuncs += [Lambda(i, log(arg.subs(self.x, i)))]
                self.newf = self.newf.xreplace({log(arg): self.t})
                new_extension = True

        return new_extension

    def _atan_part(self, atans):
        """
        Try to build an inverse tangent extension.

        Explanation
        ===========

        atan(b) is a primitive over k: D(atan(b)) == Db/(1 + b**2) is in
        k.  It is a monomial over k unless Db/(1 + b**2) == Dv +
        Sum(c_i*Dt_i) for v in k, constants c_i and the existing inverse
        tangent extensions t_i (the arguments of the tangent extensions
        are already in k, so they are absorbed by Dv), which is checked
        with limited_integrate().  Since rewriting a dependent inverse
        tangent in terms of the existing extensions requires determining
        piecewise constants (like atan(x) + atan(1/x) == +-pi/2 or
        atan(2*x/(1 - x**2)) == 2*atan(x) +- pi), dependent inverse
        tangents raise NotImplementedError.

        Returns True if there was a new extension and False if there was
        no new extension.
        """
        from .prde import limited_integrate
        new_extension = False
        for arg in ordered([i.args[0] for i in atans]):
            # limited_integrate() needs the level attributes (case, d, ...)
            # that are otherwise only set after the extension is built.
            self._auto_attrs()

            arga, argd = frac_in(arg, self.t)
            darga = (argd*derivation(Poly(arga, self.t), self) -
                arga*derivation(Poly(argd, self.t), self))
            dargd = argd**2
            darg = darga.as_expr()/dargd.as_expr()
            dcand = cancel(darg/(1 + arg**2))

            G = [frac_in(self.D[i].as_expr(), self.t) for i, ext in
                enumerate(self.exts) if ext == 'atan']
            dca, dcd = frac_in(dcand, self.t)
            try:
                A = limited_integrate(dca, dcd, G, self)
            except NonElementaryIntegralException:
                A = None
            if A is not None:
                raise NotImplementedError("%s is algebraically dependent "
                    "on the existing extension tower, which is not "
                    "supported." % atan(arg))

            self.t = next(self.ts)
            self.T.append(self.t)
            self.extargs.append(arg)
            self.exts.append('atan')
            self.D.append(dcand.as_poly(self.t, expand=False))
            if self.dummy:
                i = Dummy("i")
            else:
                i = Symbol('i')
            self.Tfuncs += [Lambda(i, atan(arg.subs(self.x, i)))]
            self.newf = self.newf.xreplace({atan(arg): self.t})
            new_extension = True

        return new_extension

    def _tan_part(self, tans):
        """
        Try to build a tangent (hypertangent) extension.

        Returns
        =======

        Returns True if there was a new extension and False if there was
        no new extension.  A tangent whose argument is an integer linear
        combination of the arguments of the existing tangent extensions is
        rewritten as a rational function of them; other algebraically
        dependent tangents (whose argument differs from a rational linear
        combination of the existing arguments by a constant, like tan(x)
        and tan(x + 1), or tan(x/2) appearing after tan(x)) raise
        NotImplementedError.
        """
        new_extension = False
        for arg in ordered([i.args[0] for i in tans]):
            # The tangent of a multiple of an inverse tangent extension is
            # not transcendental: tan(n*atan(b)) is a rational function of
            # b for an integer n (and algebraic otherwise).
            atanvars = [(self.T[i], self.extargs[i]) for i, ext in
                enumerate(self.exts) if ext == 'atan']
            if any(arg.has(v) for v, _ in atanvars):
                for v, b in atanvars:
                    n = cancel(arg/v)
                    if n.is_Integer:
                        self.newf = self.newf.xreplace({tan(arg):
                            _tan_of_multiple(n, b)})
                        break
                else:
                    raise NotImplementedError("The tangent of an "
                        "expression involving an inverse tangent, like "
                        "%s, is not supported." % tan(arg))
                continue

            # tan(arg) is algebraic over the existing tower if and only if
            # arg differs from a rational linear combination of the
            # arguments of the existing tangent extensions by a constant
            # (tan of a sum is a rational function of the tangents of the
            # summands).  Solving Darg == Sum(q_i*Da_i) for constant q_i
            # finds all such relations at once, including ones spanning
            # several tangents, like tan(x + x**2) with tan(x) and
            # tan(x**2), which a pairwise check would miss.
            darg = derivation(arg, self, basic=True)
            tanidx = [i for i, ext in enumerate(self.exts) if ext == 'tan']
            if tanidx:
                from .prde import constant_system
                from sympy.polys.polymatrix import PolyMatrix
                dargs = [derivation(self.extargs[i], self, basic=True)
                    for i in tanidx]
                dum = Dummy()
                A, u = constant_system(PolyMatrix([dargs], dum),
                    PolyMatrix([darg], dum), self)
                u = u.to_Matrix()
                if A and all(derivation(qi, self, basic=True).is_zero
                        for qi in u) and cancel(darg -
                        Add(*[qi*dai for qi, dai in zip(u, dargs)])) == 0:
                    # arg == Sum(q_i*a_i) + const
                    if all(qi.is_Integer for qi in u) and cancel(arg -
                            Add(*[qi*self.extargs[i] for qi, i in
                                zip(u, tanidx)])) == 0:
                        self.newf = self.newf.xreplace({tan(arg):
                            _tan_of_combination([(qi, self.T[i])
                                for qi, i in zip(u, tanidx)])})
                        continue
                    raise NotImplementedError("%s is algebraically "
                        "dependent on the tangents of %s, which is not "
                        "supported." % (tan(arg), ', '.join(
                            str(self.extargs[i]) for i in tanidx)))

            arga, argd = frac_in(arg, self.t)
            darga = (argd*derivation(Poly(arga, self.t), self) -
                arga*derivation(Poly(argd, self.t), self))
            dargd = argd**2
            darga, dargd = darga.cancel(dargd, include=True)
            darg = darga.as_expr()/dargd.as_expr()
            self.t = next(self.ts)
            self.T.append(self.t)
            self.extargs.append(arg)
            self.exts.append('tan')
            self.D.append(darg.as_poly(self.t, expand=False)*Poly(self.t**2 +
                1, self.t, expand=False))
            if self.dummy:
                i = Dummy("i")
            else:
                i = Symbol('i')
            self.Tfuncs += [Lambda(i, tan(arg.subs(self.x, i)))]
            self.newf = self.newf.xreplace({tan(arg): self.t})
            # The polynomial part of the integral of a hypertangent
            # element contains terms c*log(tan(arg)**2 + 1); rewrite them
            # in terms of cos to make the answer look nicer.  The pattern
            # is expressed in terms of the original functions, since the
            # back-substitution is applied after the T variables have been
            # replaced by them.
            args = arg.subs(list(zip(reversed(self.T[1:]),
                reversed([f(self.x) for f in self.Tfuncs]))))
            self.backsubs += [(log(tan(args)**2 + 1), -2*log(cos(args)))]
            new_extension = True

        return new_extension

    @property
    def _important_attrs(self):
        """
        Returns some of the more important attributes of self.

        Explanation
        ===========

        Used for testing and debugging purposes.

        The attributes are (fa, fd, D, T, Tfuncs, backsubs,
        exts, extargs).
        """
        return (self.fa, self.fd, self.D, self.T, self.Tfuncs,
            self.backsubs, self.exts, self.extargs)

    # NOTE: this printing doesn't follow the Python's standard
    # eval(repr(DE)) == DE, where DE is the DifferentialExtension object,
    # also this printing is supposed to contain all the important
    # attributes of a DifferentialExtension object
    def __repr__(self):
        # no need to have GeneratorType object printed in it
        r = [(attr, getattr(self, attr)) for attr in self.__slots__
                if not isinstance(getattr(self, attr), GeneratorType)]
        return self.__class__.__name__ + '(dict(%r))' % (r)

    # fancy printing of DifferentialExtension object
    def __str__(self):
        return (self.__class__.__name__ + '({fa=%s, fd=%s, D=%s})' %
                (self.fa, self.fd, self.D))

    # should only be used for debugging purposes, internally
    # f1 = f2 = log(x) at different places in code execution
    # may return D1 != D2 as True, since 'level' or other attribute
    # may differ
    def __eq__(self, other):
        for attr in self.__class__.__slots__:
            d1, d2 = getattr(self, attr), getattr(other, attr)
            if not (isinstance(d1, GeneratorType) or d1 == d2):
                return False
        return True

    def reset(self):
        """
        Reset self to an initial state.  Used by __init__.
        """
        self.t = self.x
        self.T = [self.x]
        self.D = [Poly(1, self.x)]
        self.level = -1
        self.exts = [None]
        self.extargs = [None]
        if self.dummy:
            self.ts = numbered_symbols('t', cls=Dummy)
        else:
            # For testing
            self.ts = numbered_symbols('t')
        # For various things that we change to make things work that we need to
        # change back when we are done.
        self.backsubs = []
        self.Tfuncs = []
        self.newf = self.f

    def indices(self, extension):
        """
        Parameters
        ==========

        extension : str
            Represents a valid extension type.

        Returns
        =======

        list: A list of indices of 'exts' where extension of
            type 'extension' is present.

        Examples
        ========

        >>> from sympy.integrals.risch import DifferentialExtension
        >>> from sympy import log, exp
        >>> from sympy.abc import x
        >>> DE = DifferentialExtension(log(x) + exp(x), x, handle_first='exp')
        >>> DE.indices('log')
        [2]
        >>> DE.indices('exp')
        [1]

        """
        return [i for i, ext in enumerate(self.exts) if ext == extension]

    def increment_level(self):
        """
        Increment the level of self.

        Explanation
        ===========

        This makes the working differential extension larger.  self.level is
        given relative to the end of the list (-1, -2, etc.), so we do not need
        do worry about it when building the extension.
        """
        if self.level >= -1:
            raise ValueError("The level of the differential extension cannot "
                "be incremented any further.")

        self.level += 1
        self.t = self.T[self.level]
        self.d = self.D[self.level]
        self.case = self.cases[self.level]
        return None

    def decrement_level(self):
        """
        Decrease the level of self.

        Explanation
        ===========

        This makes the working differential extension smaller.  self.level is
        given relative to the end of the list (-1, -2, etc.), so we do not need
        do worry about it when building the extension.
        """
        if self.level <= -len(self.T):
            raise ValueError("The level of the differential extension cannot "
                "be decremented any further.")

        self.level -= 1
        self.t = self.T[self.level]
        self.d = self.D[self.level]
        self.case = self.cases[self.level]
        return None


def update_sets(seq, atoms, func):
    s = set(seq)
    s = atoms.intersection(s)
    new = atoms - s
    s.update(list(filter(func, new)))
    return list(s)


class DecrementLevel:
    """
    A context manager for decrementing the level of a DifferentialExtension.
    """
    __slots__ = ('DE',)

    def __init__(self, DE):
        self.DE = DE
        return

    def __enter__(self):
        self.DE.decrement_level()

    def __exit__(self, exc_type, exc_value, traceback):
        self.DE.increment_level()


class NonElementaryIntegralException(Exception):
    """
    Exception used by subroutines within the Risch algorithm to indicate to one
    another that the function being integrated does not have an elementary
    integral in the given differential field.
    """
    # TODO: Rewrite algorithms below to use this (?)

    # TODO: Pass through information about why the integral was nonelementary,
    # and store that in the resulting NonElementaryIntegral somehow.
    pass


def gcdex_diophantine(a, b, c):
    """
    Extended Euclidean Algorithm, Diophantine version.

    Explanation
    ===========

    Given ``a``, ``b`` in K[x] and ``c`` in (a, b), the ideal generated by ``a`` and
    ``b``, return (s, t) such that s*a + t*b == c and either s == 0 or s.degree()
    < b.degree().

    This is the Diophantine version of ``ExtendedEuclidean`` from Section
    1.3 of Bronstein's book.
    """
    # TODO: This should go in densetools.py.
    # XXX: Better name?

    s, g = a.half_gcdex(b)
    s *= c.exquo(g)  # Inexact division means c is not in (a, b)
    if s and s.degree() >= b.degree():
        _, s = s.div(b)
    t = (c - s*a).exquo(b)
    return (s, t)


def frac_in(f, t, *, cancel=False, **kwargs):
    """
    Returns the tuple (fa, fd), where fa and fd are Polys in t.

    Explanation
    ===========

    This is a common idiom in the Risch Algorithm functions, so we abstract
    it out here. ``f`` should be a basic expression, a Poly, or a tuple (fa, fd),
    where fa and fd are either basic expressions or Polys, and f == fa/fd.
    **kwargs are applied to Poly.
    """
    if isinstance(f, tuple):
        fa, fd = f
        f = fa.as_expr()/fd.as_expr()
    fa, fd = f.as_expr().as_numer_denom()
    fa, fd = fa.as_poly(t, **kwargs), fd.as_poly(t, **kwargs)
    if cancel:
        fa, fd = fa.cancel(fd, include=True)
    if fa is None or fd is None:
        raise ValueError("Could not turn %s into a fraction in %s." % (f, t))
    return (fa, fd)


def as_poly_1t(p, t, z):
    """
    (Hackish) way to convert an element ``p`` of K[t, 1/t] to K[t, z].

    In other words, ``z == 1/t`` will be a dummy variable that Poly can handle
    better.

    See issue 5131.

    Examples
    ========

    >>> from sympy import random_poly
    >>> from sympy.integrals.risch import as_poly_1t
    >>> from sympy.abc import x, z

    >>> p1 = random_poly(x, 10, -10, 10)
    >>> p2 = random_poly(x, 10, -10, 10)
    >>> p = p1 + p2.subs(x, 1/x)
    >>> as_poly_1t(p, x, z).as_expr().subs(z, 1/x) == p
    True
    """
    # TODO: Use this on the final result.  That way, we can avoid answers like
    # (...)*exp(-x).
    pa, pd = frac_in(p, t, cancel=True)
    if not pd.is_monomial:
        # XXX: Is there a better Poly exception that we could raise here?
        # Either way, if you see this (from the Risch Algorithm) it indicates
        # a bug.
        raise PolynomialError("%s is not an element of K[%s, 1/%s]." % (p, t, t))

    t_part, remainder = pa.div(pd)

    ans = t_part.as_poly(t, z, expand=False)

    if remainder:
        one = remainder.one
        tp = t*one
        r = pd.degree() - remainder.degree()
        z_part = remainder.transform(one, tp) * tp**r
        z_part = z_part.replace(t, z).to_field().quo_ground(pd.LC())
        ans += z_part.as_poly(t, z, expand=False)

    return ans


def derivation(p, DE, coefficientD=False, basic=False):
    """
    Computes Dp.

    Explanation
    ===========

    Given the derivation D with D = d/dx and p is a polynomial in t over
    K(x), return Dp.

    If coefficientD is True, it computes the derivation kD
    (kappaD), which is defined as kD(sum(ai*Xi**i, (i, 0, n))) ==
    sum(Dai*Xi**i, (i, 1, n)) (Definition 3.2.2).  X in this case is
    T[-1], so coefficientD computes the derivative just with respect to T[:-1],
    with T[-1] treated as a constant.

    If ``basic=True``, the returns a Basic expression.  Elements of D can still be
    instances of Poly.

    See Definition 3.2.2 in Section 3.2 of Bronstein's book.
    """
    if basic:
        r = 0
    else:
        r = Poly(0, DE.t)

    t = DE.t
    if coefficientD:
        if DE.level <= -len(DE.T):
            # 'base' case, the answer is 0.
            return r
        DE.decrement_level()

    D = DE.D[:len(DE.D) + DE.level + 1]
    T = DE.T[:len(DE.T) + DE.level + 1]

    for d, v in zip(D, T):
        pv = p.as_poly(v)
        if pv is None or basic:
            pv = p.as_expr()

        if basic:
            r += d.as_expr()*pv.diff(v)
        else:
            r += (d.as_expr()*pv.diff(v).as_expr()).as_poly(t)

    if basic:
        r = cancel(r)
    if coefficientD:
        DE.increment_level()

    return r


def get_case(d, t):
    """
    Returns the type of the derivation d.

    Returns one of {'exp', 'tan', 'base', 'primitive', 'other_linear',
    'other_nonlinear'}.
    """
    if not d.expr.has(t):
        if d.is_one:
            return 'base'
        return 'primitive'
    if d.rem(Poly(t, t)).is_zero:
        return 'exp'
    if d.rem(Poly(1 + t**2, t)).is_zero:
        return 'tan'
    if d.degree(t) > 1:
        return 'other_nonlinear'
    return 'other_linear'


def splitfactor(p, DE, coefficientD=False, z=None):
    """
    Splitting factorization.

    Explanation
    ===========

    Given a derivation D on k[t] and ``p`` in k[t], return (p_n, p_s) in
    k[t] x k[t] such that p = p_n*p_s, p_s is special, and each square
    factor of p_n is normal.

    This is ``SplitFactor`` from Section 3.5 of Bronstein's book.
    """
    kinv = [1/x for x in DE.T[:DE.level]]
    if z:
        kinv.append(z)

    One = Poly(1, DE.t, domain=p.get_domain())
    Dp = derivation(p, DE, coefficientD=coefficientD)
    # XXX: Is this right?
    if p.is_zero:
        return (p, One)

    if not p.expr.has(DE.t):
        s = p.as_poly(*kinv).gcd(Dp.as_poly(*kinv)).as_poly(DE.t)
        n = p.exquo(s)
        return (n, s)

    if not Dp.is_zero:
        h = p.gcd(Dp).to_field()
        g = p.gcd(p.diff(DE.t)).to_field()
        s = h.exquo(g)

        if s.degree(DE.t) == 0:
            return (p, One)

        q_split = splitfactor(p.exquo(s), DE, coefficientD=coefficientD)

        return (q_split[0], q_split[1]*s)
    else:
        return (p, One)


def splitfactor_sqf(p, DE, coefficientD=False, z=None, basic=False):
    """
    Splitting Square-free Factorization.

    Explanation
    ===========

    Given a derivation D on k[t] and ``p`` in k[t], returns (N1, ..., Nm)
    and (S1, ..., Sm) in k[t]^m such that p =
    (N1*N2**2*...*Nm**m)*(S1*S2**2*...*Sm**m) is a splitting
    factorization of ``p`` and the Ni and Si are square-free and coprime.

    This is ``SplitSquarefreeFactor`` from Section 3.5 of Bronstein's book.
    """
    # TODO: This algorithm appears to be faster in every case
    # TODO: Verify this and splitfactor() for multiple extensions
    kkinv = [1/x for x in DE.T[:DE.level]] + DE.T[:DE.level]
    if z:
        kkinv = [z]

    S = []
    N = []
    p_sqf = p.sqf_list_include()
    if p.is_zero:
        return (((p, 1),), ())

    for pj, i in p_sqf:
        Si = pj.as_poly(*kkinv).gcd(derivation(pj, DE,
            coefficientD=coefficientD,basic=basic).as_poly(*kkinv)).as_poly(DE.t)
        pi = Poly(pj, DE.t)
        Ni = pi.exquo(Si)
        if not Si.is_one:
            S.append((Si, i))
        if not Ni.is_one:
            N.append((Ni, i))

    return (tuple(N), tuple(S))


def canonical_representation(a, d, DE):
    """
    Canonical Representation.

    Explanation
    ===========

    Given a derivation D on k[t] and f = a/d in k(t), return (f_p, f_s,
    f_n) in k[t] x k(t) x k(t) such that f = f_p + f_s + f_n is the
    canonical representation of f (f_p is a polynomial, f_s is reduced
    (has a special denominator), and f_n is simple (has a normal
    denominator).

    This is ``CanonicalRepresentation`` from Section 3.5 of Bronstein's
    book.
    """
    # Make d monic
    l = Poly(1/d.LC(), DE.t)
    a, d = a.mul(l), d.mul(l)

    q, r = a.div(d)
    dn, ds = splitfactor(d, DE)

    b, c = gcdex_diophantine(dn.as_poly(DE.t), ds.as_poly(DE.t), r.as_poly(DE.t))
    b, c = b.as_poly(DE.t), c.as_poly(DE.t)

    return (q, (b, ds), (c, dn))


def hermite_reduce(a, d, DE):
    """
    Hermite Reduction - Mack's Linear Version.

    Given a derivation D on k(t) and f = a/d in k(t), returns g, h, r in
    k(t) such that f = Dg + h + r, h is simple, and r is reduced.

    This is ``HermiteReduce`` (Mack's linear version) from Section 5.3 of
    Bronstein's book.
    """
    # Make d monic
    l = Poly(1/d.LC(), DE.t)
    a, d = a.mul(l), d.mul(l)

    fp, fs, fn = canonical_representation(a, d, DE)
    a, d = fn
    l = Poly(1/d.LC(), DE.t)
    a, d = a.mul(l), d.mul(l)

    ga = Poly(0, DE.t)
    gd = Poly(1, DE.t)

    dd = derivation(d, DE)
    dm = gcd(d.to_field(), dd.to_field()).as_poly(DE.t)
    ds, _ = d.div(dm)

    while dm.degree(DE.t) > 0:

        ddm = derivation(dm, DE)
        dm2 = gcd(dm.to_field(), ddm.to_field())
        dms, _ = dm.div(dm2)
        ds_ddm = ds.mul(ddm)
        ds_ddm_dm, _ = ds_ddm.div(dm)

        b, c = gcdex_diophantine(-ds_ddm_dm.as_poly(DE.t),
            dms.as_poly(DE.t), a.as_poly(DE.t))
        b, c = b.as_poly(DE.t), c.as_poly(DE.t)

        db = derivation(b, DE).as_poly(DE.t)
        ds_dms, _ = ds.div(dms)
        a = c.as_poly(DE.t) - db.mul(ds_dms).as_poly(DE.t)

        ga = ga*dm + b*gd
        gd = gd*dm
        ga, gd = ga.cancel(gd, include=True)
        dm = dm2

    q, r = a.div(ds)
    ga, gd = ga.cancel(gd, include=True)

    r, d = r.cancel(ds, include=True)
    rra = q*fs[1] + fp*fs[1] + fs[0]
    rrd = fs[1]
    rra, rrd = rra.cancel(rrd, include=True)

    return ((ga, gd), (r, d), (rra, rrd))


def polynomial_reduce(p, DE):
    """
    Polynomial Reduction.

    Explanation
    ===========

    Given a derivation D on k(t) and p in k[t] where t is a nonlinear
    monomial over k, return q, r in k[t] such that p = Dq  + r, and
    deg(r) < deg_t(Dt).

    This is ``PolynomialReduce`` from Section 5.4 of Bronstein's book.
    """
    q = Poly(0, DE.t)
    while p.degree(DE.t) >= DE.d.degree(DE.t):
        m = p.degree(DE.t) - DE.d.degree(DE.t) + 1
        q0 = Poly(DE.t**m, DE.t).mul(Poly(p.as_poly(DE.t).LC()/
            (m*DE.d.LC()), DE.t))
        q += q0
        p = p - derivation(q0, DE)

    return (q, p)


def laurent_series(a, d, F, n, DE):
    """
    Contribution of ``F`` to the full partial fraction decomposition of A/D.

    Explanation
    ===========

    Given a field K of characteristic 0 and ``A``,``D``,``F`` in K[x] with D monic,
    nonzero, coprime with A, and ``F`` the factor of multiplicity n in the square-
    free factorization of D, return the principal parts of the Laurent series of
    A/D at all the zeros of ``F``.

    This is ``LaurentSeries`` from Section 2.7 of Bronstein's book.
    """
    if F.degree()==0:
        return 0
    Z = _symbols('z', n)
    z = Symbol('z')
    Z.insert(0, z)
    delta_a = Poly(0, DE.t)
    delta_d = Poly(1, DE.t)

    E = d.quo(F**n)
    ha, hd = (a, E*Poly(z**n, DE.t))
    dF = derivation(F,DE)
    B, _ = gcdex_diophantine(E, F, Poly(1,DE.t))
    C, _ = gcdex_diophantine(dF, F, Poly(1,DE.t))

    # initialization
    F_store = F
    V, DE_D_list, H_list= [], [], []

    for j in range(0, n):
    # jth derivative of z would be substituted with dfnth/(j+1) where dfnth =(d^n)f/(dx)^n
        F_store = derivation(F_store, DE)
        v = (F_store.as_expr())/(j + 1)
        V.append(v)
        DE_D_list.append(Poly(Z[j + 1],Z[j]))

    DE_new = DifferentialExtension(extension = {'D': DE_D_list}) #a differential indeterminate
    for j in range(0, n):
        zEha = Poly(z**(n + j), DE.t)*E**(j + 1)*ha
        zEhd = hd
        Pa, Pd = cancel((zEha, zEhd))[1], cancel((zEha, zEhd))[2]
        Q = Pa.quo(Pd)
        for i in range(0, j + 1):
            Q = Q.subs(Z[i], V[i])
        # ha/hd == h**(j)/j!; differentiate once more (in both x and the
        # differential indeterminate z) and divide by j + 1, using the
        # quotient rule: D(ha/hd) == (hd*D(ha) - ha*D(hd))/hd**2.
        Dha = (hd*derivation(ha, DE, basic=True).as_poly(DE.t)
             - ha*derivation(hd, DE, basic=True).as_poly(DE.t)
             + hd*derivation(ha, DE_new, basic=True).as_poly(DE.t)
             - ha*derivation(hd, DE_new, basic=True).as_poly(DE.t))
        Dhd = Poly(j + 1, DE.t)*hd**2
        ha, hd = Dha, Dhd

        QBC = Poly(Q, DE.t)*B**(1 + j)*C**(n + j)
        H_list.append(QBC)

        Ff, _ = F.div(gcd(F, Q))
        F_stara, F_stard = frac_in(Ff, DE.t)
        if F_stara.degree(DE.t) - F_stard.degree(DE.t) > 0:
            # XXX: Only the principal parts at the real zeros of F are
            # included; contributions from complex zeros are dropped (the
            # H_list return value covers all zeros, since it is not
            # evaluated at the roots).
            H = (QBC*F_stard).rem(F_stara)
            alphas = real_roots(F_stara)
            for alpha in list(alphas):
                # delta += H(alpha)/(t - alpha)**(n - j)
                pa = Poly((DE.t - alpha)**(n - j), DE.t)
                delta_a = delta_a*pa + Poly(H.eval(alpha), DE.t)*delta_d
                delta_d = delta_d*pa
    delta_a, delta_d = delta_a.cancel(delta_d, include=True)
    return (delta_a, delta_d, H_list)


def recognize_derivative(a, d, DE, z=None):
    """
    Compute the squarefree factorization of the denominator of f
    and for each Di the polynomial H in K[x] (see Theorem 2.7.1), using the
    LaurentSeries algorithm. Write Di = GiEi where Gj = gcd(Hn, Di) and
    gcd(Ei,Hn) = 1. Since the residues of f at the roots of Gj are all 0, and
    the residue of f at a root alpha of Ei is Hi(a) != 0, f is the derivative of a
    rational function if and only if Ei = 1 for each i, which is equivalent to
    Di | H[-1] for each i.

    There is no named counterpart in Bronstein's book; this is based on
    Theorem 2.7.1 (Section 2.7) via ``LaurentSeries``.
    """
    flag = True
    a, d = a.cancel(d, include=True)
    _, r = a.div(d)
    Np, Sp = splitfactor_sqf(d, DE, coefficientD=True, z=z)

    # Degree-zero (content) factors are not poles
    Np = [(s, n) for s, n in Np if s.degree(DE.t) > 0]
    if any(n == 1 for s, n in Np):
        # A simple pole at a normal prime always has a nonzero residue
        # (nu_p(Dv) == nu_p(v) - 1 <= -2 for any pole of Dv at a normal p),
        # so f is not the derivative of a rational function.
        return False
    if Np:
        # The Laurent series machinery below is only justified for factors
        # with constant roots (Theorem 2.7.1 is stated for K[x] with the
        # roots algebraic over the constant field K); deciding whether the
        # residues vanish at nonconstant poles of order > 1 is not yet
        # implemented.
        raise NotImplementedError("recognize_derivative() cannot decide "
            "residue vanishing at nonconstant poles of order greater "
            "than 1.")

    for s, n in Sp:
        delta_a, delta_d, H = laurent_series(r, d, s, n, DE)
        if not H[-1].rem(s.as_poly(DE.t)).is_zero:  # Di does not divide H[-1]
            flag = False
            break
    return flag


def recognize_log_derivative(a, d, DE, z=None):
    """
    There exists a v in K(x)* such that f = dv/v
    where f a rational function if and only if f can be written as f = A/D
    where D is squarefree,deg(A) < deg(D), gcd(A, D) = 1,
    and all the roots of the Rothstein-Trager resultant are integers. In that case,
    any of the Rothstein-Trager, Lazard-Rioboo-Trager or Czichowski algorithm
    produces u in K(x) such that du/dx = uf.

    This is the in-field logarithmic derivative problem from Section 5.12
    of Bronstein's book (the "Recognizing Logarithmic Derivatives"
    subsection; neither edition gives it as named pseudocode).
    """

    z = z or Dummy('z')
    a, d = a.cancel(d, include=True)
    _, a = a.div(d)

    pz = Poly(z, DE.t)
    Dd = derivation(d, DE)
    q = a - pz*Dd
    r, _ = d.resultant(q, includePRS=True)
    r = Poly(r, z)
    Np, Sp = splitfactor_sqf(r, DE, coefficientD=True, z=z)

    if any(s.as_poly(z).degree() > 0 for s, _ in Np):
        # The normal part of the splitting factorization contains the
        # factors of the resultant whose roots are not constants; such
        # roots cannot be integers, so f is not the logarithmic derivative
        # of a k(t)-radical.
        return False

    for s, _ in Sp:
        # All roots of the resultant must be integers, which holds iff the
        # rational roots of s, counted with multiplicity, account for its
        # full degree (i.e., s splits into linear factors over QQ) and are
        # all integers.  This covers complex and irrational real roots with
        # the same degree count, and is much faster than isolating the real
        # roots.
        s = s.as_poly(z)
        r = s.ground_roots()

        if sum(r.values()) != s.degree() or not all(j.is_Integer for j in r):
            return False
    return True

def residue_reduce(a, d, DE, z=None, invert=True):
    """
    Lazard-Rioboo-Rothstein-Trager resultant reduction.

    Explanation
    ===========

    Given a derivation ``D`` on k(t) and f in k(t) simple, return g
    elementary over k(t) and a Boolean b in {True, False} such that f -
    Dg in k[t] if b == True or f + h and f + h - Dg do not have an
    elementary integral over k(t) for any h in k<t> (reduced) if b ==
    False.

    Returns (G, b), where G is a tuple of tuples of the form (s_i, S_i),
    such that g = Add(*[RootSum(s_i, lambda z: z*log(S_i(z, t))) for
    S_i, s_i in G]). f - Dg is the remaining integral, which is elementary
    only if b == True, and hence the integral of f is elementary only if
    b == True.

    f - Dg is not calculated in this function because that would require
    explicitly calculating the RootSum.  Use residue_reduce_derivation().

    This is ``ResidueReduce`` (the Lazard-Rioboo-Trager version) from
    Section 5.6 of Bronstein's book.
    """
    # TODO: Use log_to_atan() from rationaltools.py
    # If r = residue_reduce(...), then the logarithmic part is given by:
    # sum([RootSum(a[0].as_poly(z), lambda i: i*log(a[1].as_expr()).subs(z,
    # i)).subs(t, log(x)) for a in r[0]])

    z = z or Dummy('z')
    a, d = a.cancel(d, include=True)
    a, d = a.to_field().mul_ground(1/d.LC()), d.to_field().mul_ground(1/d.LC())
    kkinv = [1/x for x in DE.T[:DE.level]] + DE.T[:DE.level]

    if a.is_zero:
        return ([], True)
    _, a = a.div(d)

    pz = Poly(z, DE.t)

    Dd = derivation(d, DE)
    q = a - pz*Dd

    if Dd.degree(DE.t) <= d.degree(DE.t):
        r, R = d.resultant(q, includePRS=True)
    else:
        r, R = q.resultant(d, includePRS=True)

    R_map, H = {}, []
    for i in R:
        R_map[i.degree()] = i

    r = Poly(r, z)
    Np, Sp = splitfactor_sqf(r, DE, coefficientD=True, z=z)

    for s, i in Sp:
        if i == d.degree(DE.t):
            s = Poly(s, z).monic()
            H.append((s, d))
        else:
            h = R_map.get(i)
            if h is None:
                continue
            h_lc = Poly(h.as_poly(DE.t).LC(), DE.t, field=True)

            h_lc_sqf = h_lc.sqf_list_include(all=True)

            for a, j in h_lc_sqf:
                h = Poly(h, DE.t, field=True).exquo(Poly(gcd(a, s**j, *kkinv),
                    DE.t))

            s = Poly(s, z).monic()

            if invert:
                h_lc = Poly(h.as_poly(DE.t).LC(), DE.t, field=True, expand=False)
                inv, coeffs = h_lc.as_poly(z, field=True).invert(s), [S.One]

                for coeff in h.coeffs()[1:]:
                    L = reduced(inv*coeff.as_poly(inv.gens), [s])[1]
                    coeffs.append(L.as_expr())

                h = Poly(dict(list(zip(h.monoms(), coeffs))), DE.t)

            H.append((s, h))

    b = not any(cancel(i.as_expr()).has(DE.t, z) for i, _ in Np)

    return (H, b)


def residue_reduce_to_basic(H, DE, z):
    """
    Converts the tuple returned by residue_reduce() into a Basic expression.
    """
    # TODO: check what Lambda does with RootOf
    i = Dummy('i')
    s = list(zip(reversed(DE.T), reversed([f(DE.x) for f in DE.Tfuncs])))

    return sum(RootSum(a[0].as_poly(z), Lambda(i, i*log(a[1].as_expr()).subs(
        {z: i}).subs(s))) for a in H)


def residue_reduce_derivation(H, DE, z):
    """
    Computes the derivation of an expression returned by residue_reduce().

    In general, this is a rational function in t, so this returns an
    as_expr() result.
    """
    # TODO: verify that this is correct for multiple extensions
    i = Dummy('i')
    return S(sum(RootSum(a[0].as_poly(z), Lambda(i, i*derivation(a[1],
        DE).as_expr().subs(z, i)/a[1].as_expr().subs(z, i))) for a in H))


def integrate_primitive_polynomial(p, DE):
    """
    Integration of primitive polynomials.

    Explanation
    ===========

    Given a primitive monomial t over k, and ``p`` in k[t], return q in k[t],
    r in k, and a bool b in {True, False} such that r = p - Dq is in k if b is
    True, or r = p - Dq does not have an elementary integral over k(t) if b is
    False.

    This is ``IntegratePrimitivePolynomial`` from Section 5.8 of
    Bronstein's book.
    """
    Zero = Poly(0, DE.t)
    q = Poly(0, DE.t)

    if not p.expr.has(DE.t):
        return (Zero, p, True)

    from .prde import limited_integrate
    while True:
        if not p.expr.has(DE.t):
            return (q, p, True)

        Dta, Dtb = frac_in(DE.d, DE.T[DE.level - 1])

        with DecrementLevel(DE):  # We had better be integrating the lowest extension (x)
                                  # with ratint().
            a = p.LC()
            aa, ad = frac_in(a, DE.t)

            try:
                rv = limited_integrate(aa, ad, [(Dta, Dtb)], DE)
                if rv is None:
                    raise NonElementaryIntegralException
                (ba, bd), c = rv
            except NonElementaryIntegralException:
                return (q, p, False)

        m = p.degree(DE.t)
        q0 = c[0].as_poly(DE.t)*Poly(DE.t**(m + 1)/(m + 1), DE.t) + \
            (ba.as_expr()/bd.as_expr()).as_poly(DE.t)*Poly(DE.t**m, DE.t)

        p = p - derivation(q0, DE)
        q = q + q0


def integrate_primitive_special(a, d, DE):
    """
    Integration of logarithmic functions in terms of li.

    Explanation
    ===========

    Given t == log(w) over k and f == a/d in k(t), try to write
    f == Dv + Sum(r_j*Dw/(t - c_j)) for v in k(t) and constants r_j and
    c_j, using limited_integrate(); since
    D(e**c*li(w*e**(-c))) == Dw/(log(w) - c), the integral is then
    v + Sum(r_j*e**(c_j)*li(w*e**(-c_j))).  Returns the integral as an
    Expr in terms of the original functions, or None.
    """
    from .prde import limited_integrate
    w = DE.extargs[DE.level]
    s = list(zip(reversed(DE.T), reversed([f(DE.x) for f in DE.Tfuncs])))
    dw = derivation(w, DE, basic=True)

    gens = []
    terms = []
    for c in _constant_residue_candidates(DE.t, d, DE):
        # D(e**(n*c)*li(w**n*e**(-n*c))) == w**(n - 1)*Dw/(t - c); the
        # relevant integer multiples n of the logarithm are read off
        # from the numerator of f at the pole t == c, which must be a
        # combination of w**n*Dw/w with constant coefficients
        ns = {1}
        p = Poly(DE.t - c, DE.t)
        from .rde import order_at
        k = order_at(d, p, DE.t)
        try:
            ar = (a*d.to_field().exquo(p**k).invert(p)).rem(p)
        except (DomainError, NotImplementedError, PolynomialError):
            ar = None
        if ar is not None and not ar.is_zero:
            beta = cancel(ar.as_expr()*w/dw)
            try:
                bp = Poly(beta, w)
                if all(derivation(co, DE, basic=True).is_zero
                        for co in bp.coeffs()):
                    ns.update(m[0] for m in bp.monoms() if m[0] != 0)
            except (PolynomialError, GeneratorsError):
                pass
        for n in sorted(ns):
            if not n:
                continue
            gens.append(frac_in(cancel(w**(n - 1)*dw/(DE.t - c)), DE.t))
            terms.append(exp(n*c)*li(w**n*exp(-n*c)))
    if not gens:
        return None
    try:
        A = limited_integrate(a, d, gens, DE, unsure=[])
    except (NonElementaryIntegralException, NotImplementedError):
        # NotImplementedError means the candidates cannot be decided;
        # since li is only a best-effort addition on top of a part of the
        # integrand that was already proven nonelementary, give up on it.
        return None
    if A is None:
        return None
    (va, vd), R = A
    return (va.as_expr()/vd.as_expr() +
        Add(*[R[j].as_expr()*terms[j] for j in range(len(terms))])).subs(s)


def integrate_primitive(a, d, DE, z=None, special=False):
    """
    Integration of primitive functions.

    Explanation
    ===========

    Given a primitive monomial t over k and f in k(t), return g elementary over
    k(t), i in k(t), and b in {True, False} such that i = f - Dg is in k if b
    is True or i = f - Dg does not have an elementary integral over k(t) if b
    is False.

    This function returns a Basic expression for the first argument.  If b is
    True, the second argument is Basic expression in k to recursively integrate.
    If b is False, the second argument is an unevaluated Integral, which has
    been proven to be nonelementary.

    If ``special=True`` and t is a logarithm, integrals with no elementary
    form are expressed in terms of li when possible (see
    integrate_primitive_special()).

    This is ``IntegratePrimitive`` from Section 5.8 of Bronstein's book.
    """
    # XXX: a and d must be canceled, or this might return incorrect results
    z = z or Dummy("z")
    s = list(zip(reversed(DE.T), reversed([f(DE.x) for f in DE.Tfuncs])))

    g1, h, r = hermite_reduce(a, d, DE)
    g2, b = residue_reduce(h[0], h[1], DE, z=z)
    if not b:
        if special and DE.exts and DE.exts[DE.level] == 'log' and \
                DE.extargs[DE.level] is not None:
            A = integrate_primitive_special(a, d, DE)
            if A is not None:
                return (A, S.Zero, True)
        i = cancel(a.as_expr()/d.as_expr() - (g1[1]*derivation(g1[0], DE) -
            g1[0]*derivation(g1[1], DE)).as_expr()/(g1[1]**2).as_expr() -
            residue_reduce_derivation(g2, DE, z))
        i = NonElementaryIntegral(cancel(i).subs(s), DE.x)
        return ((g1[0].as_expr()/g1[1].as_expr()).subs(s) +
            residue_reduce_to_basic(g2, DE, z), i, b)

    # h - Dg2 + r
    p = cancel(h[0].as_expr()/h[1].as_expr() - residue_reduce_derivation(g2,
        DE, z) + r[0].as_expr()/r[1].as_expr())
    p = p.as_poly(DE.t)

    q, i, b = integrate_primitive_polynomial(p, DE)

    ret = ((g1[0].as_expr()/g1[1].as_expr() + q.as_expr()).subs(s) +
        residue_reduce_to_basic(g2, DE, z))
    if not b:
        # i == p - Dq from integrate_primitive_polynomial(), which has been
        # proven to have no elementary integral over k(t), and ret contains
        # the partial q, so f == D(ret) + i holds for the returned values.
        i = NonElementaryIntegral(cancel(i.as_expr()).subs(s), DE.x)
    else:
        i = cancel(i.as_expr())

    return (ret, i, b)


def integrate_hyperexponential_polynomial(p, DE, z, unsure=None):
    """
    Integration of hyperexponential polynomials.

    Explanation
    ===========

    Given a hyperexponential monomial t over k and ``p`` in k[t, 1/t], return q in
    k[t, 1/t] and a bool b in {True, False} such that p - Dq in k if b is True,
    or p - Dq does not have an elementary integral over k(t) if b is False.

    If ``unsure`` is a list, powers of t whose Risch differential equation
    cannot be decided (NotImplementedError) are appended to it and treated
    like unsolved ones, instead of propagating the error; b == False then
    only proves nonelementarity if ``unsure`` remains empty.

    This is ``IntegrateHyperexponentialPolynomial`` from Section 5.9 of
    Bronstein's book.
    """
    t1 = DE.t
    dtt = DE.d.exquo(Poly(DE.t, DE.t))
    qa = Poly(0, DE.t)
    qd = Poly(1, DE.t)
    b = True

    if p.is_zero:
        return(qa, qd, b)

    from sympy.integrals.rde import rischDE

    with DecrementLevel(DE):
        for i in range(-p.degree(z), p.degree(t1) + 1):
            if not i:
                continue
            elif i < 0:
                # If you get AttributeError: 'NoneType' object has no attribute 'nth'
                # then this should really not have expand=False
                # But it shouldn't happen because p is already a Poly in t and z
                a = p.as_poly(z, expand=False).nth(-i)
            else:
                # If you get AttributeError: 'NoneType' object has no attribute 'nth'
                # then this should really not have expand=False
                a = p.as_poly(t1, expand=False).nth(i)

            aa, ad = frac_in(a, DE.t, field=True)
            aa, ad = aa.cancel(ad, include=True)
            iDt = Poly(i, t1)*dtt
            iDta, iDtd = frac_in(iDt, DE.t, field=True)
            try:
                va, vd = rischDE(iDta, iDtd, Poly(aa, DE.t), Poly(ad, DE.t), DE)
                va, vd = frac_in((va, vd), t1, cancel=True)
            except NonElementaryIntegralException:
                b = False
            except NotImplementedError:
                if unsure is None:
                    raise
                # It could not be decided whether this term has an
                # elementary integral; record it, so that the caller does
                # not claim a nonelementarity proof for it.
                unsure.append(i)
                b = False
            else:
                # q += v*t**i
                if i > 0:
                    ti = Poly(t1**i, t1)
                else:
                    ti = Poly(z**-i, z)

                qa = qa*vd + va*ti*qd
                qd *= vd

    return (qa, qd, b)


def _parametric_term_rde(fa, fd, ga, gd, gens, DE):
    """
    Solve Dy + f*y == g + Sum(r_j*gens_j) parametrically.

    Explanation
    ===========

    Given f, g in k(t) and gens a list of (num, den) pairs in k(t),
    either return None, in which case the equation has no solution with
    y in k(t) and constants r_j, or return (y, R), where y is an Expr and
    R the list of the constants r_j.  This is used to express integrals
    in terms of special functions: the gens are the densities of the
    derivatives of the special function candidates.
    """
    from .prde import param_rischDE
    try:
        # A solution found with unsure bounds is still a solution, and no
        # nonelementarity claims are derived from a failure here, so unsure
        # special denominator bounds are acceptable.
        h, A = param_rischDE(fa, fd, [(ga, gd)] + gens, DE, unsure=[])
    except NotImplementedError:
        # The candidates cannot be decided; since the special functions
        # are only a best-effort addition on top of a part of the
        # integrand that was already proven nonelementary, treat this
        # like the absence of a solution.
        return None
    V = [v for v in A.nullspace() if v[0] != 0]
    if not V:
        return None
    # normalize so that the coefficient of g is 1; then y == Sum(d_i*h_i)
    # solves Dy + f*y == g + Sum(c_j*gens_j) for [1, c_1, ..., c_m,
    # d_1, ..., d_r] the normalized nullspace vector
    v = V[0]/V[0][0]
    v = [i.as_expr() if isinstance(i, Poly) else S(i) for i in v]
    r = len(h)
    m = len(v) - r - 1
    y = sum(v[m + 1 + j]*h[j][0].as_expr()/h[j][1].as_expr()
        for j in range(r))
    return (cancel(y), v[1:m + 1])


def _complete_square_candidates(s, DE):
    """
    Ways of writing s in k as u**2 + c with u in k and c a constant.

    Returns a list of (u, c) pairs, one for each variable of k in which
    ``s`` is a quadratic polynomial with constant leading coefficient.
    """
    cands = []
    for v in DE.T[:len(DE.T) + DE.level + 1]:
        try:
            p = Poly(s, v)
        except PolynomialError:
            continue
        if p.degree(v) != 2:
            continue
        A, B = p.nth(2), p.nth(1)
        if not derivation(A, DE, basic=True).is_zero:
            continue
        c = cancel(p.nth(0) - B**2/(4*A))
        if not derivation(c, DE, basic=True).is_zero:
            continue
        cands.append((sqrt(A)*v + B/(2*sqrt(A)), c))
    return cands


def _constant_residue_candidates(w, d, DE):
    """
    Constants c such that w - c vanishes on a factor of d.

    Explanation
    ===========

    Given w in k and a Poly d in DE.t over k, return the constants c for
    which some irreducible factor p of d divides the numerator of w - c.
    These give the candidate arguments w - c of Ei in the integration of
    e**w-type integrands with denominator d (and similarly the poles
    t == c for li in the logarithmic case).
    """
    from sympy.polys.polytools import factor_list
    T = DE.T[:len(DE.T) + DE.level + 1]
    try:
        _, factors = factor_list(d.as_expr())
    except (DomainError, PolynomialError, NotImplementedError,
            GeneratorsError):
        _, factors = d.sqf_list()
        factors = [(p.as_expr(), m) for p, m in factors]
    cands = []
    for p, _ in factors:
        # the main variable of the factor is the topmost variable of the
        # tower appearing in it
        vs = [v for v in reversed(T) if p.has(v)]
        if not vs:
            continue
        v = vs[0]
        try:
            pp = Poly(p, v)
            wa, wd = frac_in(w, v)
        except PolynomialError:
            continue
        if pp.degree(v) < 1 or pp.gcd(wd).degree(v) > 0:
            continue
        try:
            c = (wa*wd.invert(pp)).rem(pp)
        except (DomainError, NotImplementedError):
            continue
        if c.degree(v) > 0:
            continue
        c = c.as_expr()
        if derivation(c, DE, basic=True).is_zero and \
                all(cancel(c - i) != 0 for i in cands):
            cands.append(c)
    return cands


def integrate_hyperexponential_special(pp, DE, z, unsure=None):
    """
    Integration of hyperexponential polynomials in terms of special
    functions.

    Explanation
    ===========

    Given t == exp(g) hyperexponential over k and p in k[t, 1/t] whose
    nonzero powers have no elementary integral, try to integrate each
    term a*t**i in terms of erf, erfi and Ei: with u**2 + c == -i*g for u
    in k and constant c, D(erf(u)) is a constant multiple of Du*t**i, and
    with constant c such that i*g - c vanishes on a factor of the
    denominator of a, D(Ei(i*g - c)) == e**(-c)*D(i*g)*t**i/(i*g - c).
    The remaining part y*t**i is found by solving
    Dy + i*Dg*y == a - Sum(r_j*gens_j) parametrically.

    Returns (q, sp, r, b) such that q is an Expr in k(t), sp is an Expr
    in terms of the original functions of the special part, r in k[t, 1/t]
    is the untreated part, and p == Dq + D(sp untranslated) + r, with
    b == True if all the nonzero powers were integrated.
    """
    t1 = DE.t
    g = DE.extargs[DE.level]
    s = list(zip(reversed(DE.T), reversed([f(DE.x) for f in DE.Tfuncs])))
    dg = (DE.d.exquo(Poly(DE.t, DE.t))).as_expr()

    q = S.Zero
    sp = S.Zero
    b = True

    if pp.is_zero:
        return (q, sp, S.Zero, b)
    r = pp.nth(0, 0)

    for i in range(-pp.degree(z), pp.degree(t1) + 1):
        if not i:
            continue
        if i < 0:
            a = pp.as_poly(z, expand=False).nth(-i)
        else:
            a = pp.as_poly(t1, expand=False).nth(i)
        if a == 0:
            continue

        with DecrementLevel(DE):
            gens = []
            terms = []
            for u, c in _complete_square_candidates(cancel(-i*g), DE):
                # e**(i*g) == e**(-c)*e**(-u**2), and D(erf(u)) is
                # 2/sqrt(pi)*Du*e**(-u**2); for a negative leading
                # coefficient u is imaginary, e**(-u**2) == e**(w**2)
                # with u == I*w, and erfi is used instead
                w = cancel(u/I)
                if not w.has(I):
                    if w.could_extract_minus_sign():
                        w = -w
                    gens.append(frac_in(derivation(w, DE, basic=True),
                        DE.t))
                    terms.append(sqrt(pi)/2*exp(-c)*erfi(w))
                else:
                    if u.could_extract_minus_sign():
                        u = -u
                    gens.append(frac_in(derivation(u, DE, basic=True),
                        DE.t))
                    terms.append(sqrt(pi)/2*exp(-c)*erf(u))
            aa, ad = frac_in(a, DE.t)
            w = cancel(i*g)
            dw = cancel(i*dg)
            for c in _constant_residue_candidates(w, ad, DE):
                gens.append(frac_in(cancel(dw/(w - c)), DE.t))
                terms.append(exp(c)*Ei(w - c))

            A = None
            if gens:
                fa, fd = frac_in(dw, DE.t)
                A = _parametric_term_rde(fa, fd, aa, ad, gens, DE)

        if A is None:
            b = False
            r += a*t1**i if i > 0 else a*z**-i
            continue
        y, R = A
        if unsure and i in unsure:
            # this power is solved after all
            unsure.remove(i)
        q += y*t1**i
        sp -= Add(*[R[j]*terms[j] for j in range(len(terms))])

    return (q, sp.subs(s), r, b)


def integrate_hyperexponential(a, d, DE, z=None, conds='piecewise',
        special=False):
    """
    Integration of hyperexponential functions.

    Explanation
    ===========

    Given a hyperexponential monomial t over k and f in k(t), return g
    elementary over k(t), i in k(t), and a bool b in {True, False} such that
    i = f - Dg is in k if b is True or i = f - Dg does not have an elementary
    integral over k(t) if b is False.

    This function returns a Basic expression for the first argument.  If b is
    True, the second argument is Basic expression in k to recursively integrate.
    If b is False, the second argument is an unevaluated Integral, which has
    been proven to be nonelementary.

    If ``special=True``, parts of the integral without an elementary form
    are expressed in terms of erf, erfi and Ei when possible (see
    integrate_hyperexponential_special()).

    This is ``IntegrateHyperexponential`` from Section 5.9 of Bronstein's
    book.
    """
    # XXX: a and d must be canceled, or this might return incorrect results
    z = z or Dummy("z")
    s = [(z, DE.t**-1)] + list(zip(reversed(DE.T), reversed([f(DE.x) for f in DE.Tfuncs])))

    g1, h, r = hermite_reduce(a, d, DE)
    g2, b = residue_reduce(h[0], h[1], DE, z=z)
    if not b:
        i = cancel(a.as_expr()/d.as_expr() - (g1[1]*derivation(g1[0], DE) -
            g1[0]*derivation(g1[1], DE)).as_expr()/(g1[1]**2).as_expr() -
            residue_reduce_derivation(g2, DE, z))
        i = NonElementaryIntegral(cancel(i.subs(s)), DE.x)
        return ((g1[0].as_expr()/g1[1].as_expr()).subs(s) +
            residue_reduce_to_basic(g2, DE, z), i, b)

    # p should be a polynomial in t and 1/t, because Sirr == k[t, 1/t]
    # h - Dg2 + r
    p = cancel(h[0].as_expr()/h[1].as_expr() - residue_reduce_derivation(g2,
        DE, z) + r[0].as_expr()/r[1].as_expr())
    pp = as_poly_1t(p, DE.t, z)

    # With special=True, terms whose Risch differential equation cannot be
    # decided are collected in unsure and attempted with the special
    # function candidates below; they are only fatal if they remain
    # unsolved there.
    unsure = [] if special else None
    qa, qd, b = integrate_hyperexponential_polynomial(pp, DE, z,
        unsure=unsure)

    i = pp.nth(0, 0)

    ret = ((g1[0].as_expr()/g1[1].as_expr()).subs(s) \
        + residue_reduce_to_basic(g2, DE, z))

    qas = qa.as_expr().subs(s)
    qds = qd.as_expr().subs(s)
    if conds == 'piecewise' and DE.x not in qds.free_symbols:
        # We have to be careful if the exponent is S.Zero!

        # XXX: Does qd = 0 always necessarily correspond to the exponential
        # equaling 1?
        ret += Piecewise(
                (qas/qds, Ne(qds, 0)),
                (integrate((p - i).subs(DE.t, 1).subs(s), DE.x), True)
            )
    else:
        ret += qas/qds

    if not b and special and DE.exts and DE.exts[DE.level] == 'exp' and \
            DE.extargs[DE.level] is not None:
        rem = cancel(p - (qd*derivation(qa, DE) -
            qa*derivation(qd, DE)).as_expr()/(qd**2).as_expr())
        q2, sp, rem2, b = integrate_hyperexponential_special(
            as_poly_1t(rem, DE.t, z), DE, z, unsure=unsure)
        ret += q2.subs(s) + sp
        if not b:
            if unsure:
                raise NotImplementedError("Cannot decide whether the "
                    "remaining part of the integral is elementary.")
            i = NonElementaryIntegral(cancel(rem2).subs(s), DE.x)
            return (ret, i, b)

    if not b:
        if unsure:
            raise NotImplementedError("Cannot decide whether the remaining "
                "part of the integral is elementary.")
        i = p - (qd*derivation(qa, DE) - qa*derivation(qd, DE)).as_expr()/\
            (qd**2).as_expr()
        i = NonElementaryIntegral(cancel(i).subs(s), DE.x)
    return (ret, i, b)


def integrate_hypertangent_polynomial(p, DE):
    """
    Integration of hypertangent polynomials.

    Explanation
    ===========

    Given a differential field k such that sqrt(-1) is not in k, a
    hypertangent monomial t over k, and p in k[t], return q in k[t] and
    c in k such that p - Dq - c*D(t**2 + 1)/(t**1 + 1) is in k and p -
    Dq does not have an elementary integral over k(t) if Dc != 0.

    This is ``IntegrateHypertangentPolynomial`` from Section 5.10 of
    Bronstein's book.
    """
    # XXX: Make sure that sqrt(-1) is not in k.
    q, r = polynomial_reduce(p, DE)
    a = DE.d.exquo(Poly(DE.t**2 + 1, DE.t))
    c = Poly(r.nth(1)/(2*a.as_expr()), DE.t)
    return (q, c)


def integrate_hypertangent_reduced(a, d, DE):
    """
    Integration of hypertangent reduced elements.

    Explanation
    ===========

    Given a differential field k such that sqrt(-1) is not in k, a
    hypertangent monomial t over k, and f == a/d in k(t) with d special
    (a power of t**2 + 1 times an element of k), return qa, qd in k[t]
    and b in {True, False} such that f - D(qa/qd) is in k[t] if b is
    True, or f - D(qa/qd) does not have an elementary integral over k(t)
    if b is False.

    The terms of f with denominator (t**2 + 1)**m are eliminated by
    solving a coupled differential system over k for the coefficients of
    the ansatz (u*t + v)/(t**2 + 1)**m, for decreasing m.
    """
    from .rde import coupled_DE_system, order_at
    p = Poly(DE.t**2 + 1, DE.t)
    eta = DE.d.exquo(p).as_expr()
    qa, qd = Poly(0, DE.t), Poly(1, DE.t)

    while True:
        a, d = a.cancel(d, include=True)
        m = order_at(d, p, DE.t)
        if m <= 0:
            return (qa, qd, True)

        # f == (c1*t + c0)/(t**2 + 1)**m + lower order terms, and
        # D((u*t + v)/(t**2 + 1)**m) == ((Du - 2*m*eta*v)*t + Dv +
        # 2*m*eta*u)/(t**2 + 1)**m + lower order terms.
        w = d.to_field().exquo(p**m)  # in k
        if w.degree(DE.t) > 0:
            raise ValueError("The denominator of %s/%s is not special." %
                (a, d))
        w = w.TC()
        crem = a.rem(p)
        c0, c1 = cancel(crem.nth(0)/w), cancel(crem.nth(1)/w)
        with DecrementLevel(DE):
            try:
                u, v = coupled_DE_system(S.Zero, 2*m*eta, c1, c0, DE)
            except NonElementaryIntegralException:
                return (qa, qd, False)
        q0a = Poly(u*DE.t + v, DE.t)
        q0d = p**m

        qa, qd = (qa*q0d + q0a*qd).cancel(qd*q0d, include=True)
        # f == f - D(q0a/q0d)
        a, d = (a*q0d**2 - (derivation(q0a, DE)*q0d -
            q0a*derivation(q0d, DE))*d).cancel(d*q0d**2, include=True)


def integrate_hypertangent(a, d, DE, z=None):
    """
    Integration of hypertangent functions.

    Explanation
    ===========

    Given a differential field k such that sqrt(-1) is not in k, a
    hypertangent monomial t over k and f in k(t), return g elementary
    over k(t), i in k(t), and b in {True, False} such that i == f - Dg
    is in k if b is True or i == f - Dg does not have an elementary
    integral over k(t) if b is False.

    This function returns a Basic expression for the first argument.  If
    b is True, the second argument is Basic expression in k to
    recursively integrate.  If b is False, the second argument is an
    unevaluated Integral, which has been proven to be nonelementary.
    """
    # XXX: a and d must be canceled, or this might return incorrect results
    z = z or Dummy("z")
    s = list(zip(reversed(DE.T), reversed([f(DE.x) for f in DE.Tfuncs])))

    g1, h, r = hermite_reduce(a, d, DE)
    g2, b = residue_reduce(h[0], h[1], DE, z=z)

    if not b:
        i = cancel(a.as_expr()/d.as_expr() - (g1[1]*derivation(g1[0], DE) -
            g1[0]*derivation(g1[1], DE)).as_expr()/(g1[1]**2).as_expr() -
            residue_reduce_derivation(g2, DE, z))
        i = NonElementaryIntegral(cancel(i).subs(s), DE.x)
        return ((g1[0].as_expr()/g1[1].as_expr()).subs(s) +
            residue_reduce_to_basic(g2, DE, z), i, b)

    # p == h - Dg2 + r is in k[t] plus a reduced element whose denominator
    # is a power of the special polynomial t**2 + 1.
    p = cancel(h[0].as_expr()/h[1].as_expr() - residue_reduce_derivation(g2,
        DE, z) + r[0].as_expr()/r[1].as_expr())
    pa, pd = frac_in(p, DE.t, cancel=True)

    q1a, q1d, b = integrate_hypertangent_reduced(pa, pd, DE)
    # pp == p - D(q1a/q1d)
    ppa, ppd = (pa*q1d**2 - (derivation(q1a, DE)*q1d -
        q1a*derivation(q1d, DE))*pd).cancel(pd*q1d**2, include=True)

    ret = (g1[0].as_expr()/g1[1].as_expr() +
        q1a.as_expr()/q1d.as_expr()).subs(s) + \
        residue_reduce_to_basic(g2, DE, z)
    if not b:
        i = NonElementaryIntegral(cancel(ppa.as_expr()/ppd.as_expr()).subs(s),
            DE.x)
        return (ret, i, b)

    try:
        pp = ppa.to_field().exquo(ppd).as_poly(DE.t)
    except ExactQuotientFailed:
        raise ValueError("The part of %s/%s remaining after the reduced "
            "part should be a polynomial in %s; the given extension is "
            "probably not a hypertangent monomial." % (a, d, DE.t))
    q2, c = integrate_hypertangent_polynomial(pp, DE)
    if derivation(c, DE).is_zero:
        # The polynomial part contributes c*log(t**2 + 1), since
        # D(log(t**2 + 1)) == D(t**2 + 1)/(t**2 + 1) == 2*eta*t.
        eta = DE.d.exquo(Poly(DE.t**2 + 1, DE.t)).as_expr()
        i = pp - derivation(q2, DE) - Poly(c.as_expr()*2*eta*DE.t, DE.t)
        ret += (q2.as_expr() + c.as_expr()*log(DE.t**2 + 1)).subs(s)
        return (ret, cancel(i.as_expr()), True)
    else:
        i = pp - derivation(q2, DE)
        ret += q2.as_expr().subs(s)
        i = NonElementaryIntegral(cancel(i.as_expr()).subs(s), DE.x)
        return (ret, i, False)


def integrate_nonlinear_no_specials(a, d, DE, z=None):
    """
    Integration of nonlinear monomials with no specials.

    Explanation
    ===========

    Given a nonlinear monomial t over k such that Sirr ({p in k[t] | p is
    special, monic, and irreducible}) is empty, and f in k(t), returns g
    elementary over k(t) and a Boolean b in {True, False} such that f - Dg is
    in k if b == True, or f - Dg does not have an elementary integral over k(t)
    if b == False.

    This function is applicable to all nonlinear extensions, but in the case
    where it returns b == False, it will only have proven that the integral of
    f - Dg is nonelementary if Sirr is empty.

    This function returns a Basic expression.

    This is ``IntegrateNonLinearNoSpecial`` from Section 5.11 of
    Bronstein's book.
    """
    # TODO: Integral from k?
    # TODO: split out nonelementary integral
    # XXX: a and d must be canceled, or this might not return correct results
    z = z or Dummy("z")
    s = list(zip(reversed(DE.T), reversed([f(DE.x) for f in DE.Tfuncs])))

    g1, h, r = hermite_reduce(a, d, DE)
    g2, b = residue_reduce(h[0], h[1], DE, z=z)
    if not b:
        return ((g1[0].as_expr()/g1[1].as_expr()).subs(s) +
            residue_reduce_to_basic(g2, DE, z), b)

    # Because f has no specials, this should be a polynomial in t, or else
    # there is a bug.
    p = cancel(h[0].as_expr()/h[1].as_expr() - residue_reduce_derivation(g2,
        DE, z).as_expr() + r[0].as_expr()/r[1].as_expr()).as_poly(DE.t)
    q1, q2 = polynomial_reduce(p, DE)

    if q2.expr.has(DE.t):
        b = False
    else:
        b = True

    ret = (cancel(g1[0].as_expr()/g1[1].as_expr() + q1.as_expr()).subs(s) +
        residue_reduce_to_basic(g2, DE, z))
    return (ret, b)


class NonElementaryIntegral(Integral):
    """
    Represents a nonelementary Integral.

    Explanation
    ===========

    If the result of integrate() is an instance of this class, it is
    guaranteed to be nonelementary.  Note that integrate() by default will try
    to find any closed-form solution, even in terms of special functions which
    may themselves not be elementary.  To make integrate() only give
    elementary solutions, or, in the cases where it can prove the integral to
    be nonelementary, instances of this class, use integrate(risch=True).
    In this case, integrate() may raise NotImplementedError if it cannot make
    such a determination.

    integrate() uses the deterministic Risch algorithm to integrate elementary
    functions or prove that they have no elementary integral.  In some cases,
    this algorithm can split an integral into an elementary and nonelementary
    part, so that the result of integrate will be the sum of an elementary
    expression and a NonElementaryIntegral.

    Examples
    ========

    >>> from sympy import integrate, exp, log, Integral
    >>> from sympy.abc import x

    >>> a = integrate(exp(-x**2), x, risch=True)
    >>> print(a)
    Integral(exp(-x**2), x)
    >>> type(a)
    <class 'sympy.integrals.risch.NonElementaryIntegral'>

    >>> expr = (2*log(x)**2 - log(x) - x**2)/(log(x)**3 - x**2*log(x))
    >>> b = integrate(expr, x, risch=True)
    >>> print(b)
    -log(-x + log(x))/2 + log(x + log(x))/2 + Integral(1/log(x), x)
    >>> type(b.atoms(Integral).pop())
    <class 'sympy.integrals.risch.NonElementaryIntegral'>

    """
    # TODO: This is useful in and of itself, because isinstance(result,
    # NonElementaryIntegral) will tell if the integral has been proven to be
    # elementary. But should we do more?  Perhaps a no-op .doit() if
    # elementary=True?  Or maybe some information on why the integral is
    # nonelementary.
    pass


def risch_integrate(f, x, extension=None, handle_first='log',
                    separate_integral=False, rewrite_complex=None,
                    conds='piecewise', trig=True, special=False):
    r"""
    The Risch Integration Algorithm.

    Explanation
    ===========

    Only transcendental functions are supported.  Currently, exponentials,
    logarithms, real trigonometric functions (which are handled as
    rational functions of tangents, using the half angle substitution for
    sin, cos, sec and csc) and inverse tangents (which are primitive
    extensions like logarithms) are supported.  If ``trig=False``,
    trigonometric functions raise NotImplementedError instead (unless the
    integrand contains I, or ``rewrite_complex=True``, in which case they
    are rewritten in terms of complex exponentials and logarithms, as
    before).  asin and acos are algebraic over the tangents and are not
    supported.

    If this function returns an unevaluated Integral in the result, it means
    that it has proven that integral to be nonelementary.  Any errors will
    result in raising NotImplementedError.  The unevaluated Integral will be
    an instance of NonElementaryIntegral, a subclass of Integral.

    handle_first may be either 'exp' or 'log'.  This changes the order in
    which the extension is built, and may result in a different (but
    equivalent) solution (for an example of this, see issue 5109).  It is also
    possible that the integral may be computed with one but not the other,
    because not all cases have been implemented yet.  It defaults to 'log' so
    that the outer extension is exponential when possible, because more of the
    exponential case has been implemented.

    If ``separate_integral`` is ``True``, the result is returned as a tuple (ans, i),
    where the integral is ans + i, ans is elementary, and i is either a
    NonElementaryIntegral or 0.  This useful if you want to try further
    integrating the NonElementaryIntegral part using other algorithms to
    possibly get a solution in terms of special functions.  It is False by
    default.

    If ``special`` is ``True``, parts of the integral that have no elementary
    antiderivative are expressed in terms of the special functions erf,
    erfi, Ei and li when they have such a form (see the examples below).
    An unevaluated Integral remaining in the result still means that the
    integrand has been proven to have no elementary antiderivative.  It
    is False by default.

    Examples
    ========

    >>> from sympy.integrals.risch import risch_integrate
    >>> from sympy import exp, log, pprint, tan
    >>> from sympy.abc import x

    First, we try integrating exp(-x**2). Except for a constant factor of
    2/sqrt(pi), this is the famous error function.

    >>> pprint(risch_integrate(exp(-x**2), x))
      /
     |
     |    2
     |  -x
     | e    dx
     |
    /

    The unevaluated Integral in the result means that risch_integrate() has
    proven that exp(-x**2) does not have an elementary anti-derivative.

    In many cases, risch_integrate() can split out the elementary
    anti-derivative part from the nonelementary anti-derivative part.
    For example,

    >>> pprint(risch_integrate((2*log(x)**2 - log(x) - x**2)/(log(x)**3 -
    ... x**2*log(x)), x))
                                             /
                                            |
      log(-x + log(x))   log(x + log(x))    |   1
    - ---------------- + --------------- +  | ------ dx
             2                  2           | log(x)
                                            |
                                           /

    This means that it has proven that the integral of 1/log(x) is
    nonelementary.  This function is also known as the logarithmic integral,
    and is often denoted as Li(x).

    risch_integrate() currently only accepts purely transcendental functions
    with exponentials, logarithms and trigonometric functions, though note
    that this can include nested exponentials and logarithms, as well as
    exponentials with bases other than E.

    >>> pprint(risch_integrate(exp(x)*exp(exp(x)), x))
     / x\
     \e /
    e
    >>> pprint(risch_integrate(exp(exp(x)), x))
      /
     |
     |  / x\
     |  \e /
     | e     dx
     |
    /

    Trigonometric functions are handled as rational functions of tangents.

    >>> risch_integrate(tan(x), x)
    -log(cos(x))
    >>> pprint(risch_integrate(tan(x)/x, x))
      /
     |
     | tan(x)
     | ------ dx
     |   x
     |
    /

    With ``special=True``, nonelementary parts are expressed in terms of
    special functions when possible.

    >>> risch_integrate(exp(-x**2), x, special=True)
    sqrt(pi)*erf(x)/2
    >>> risch_integrate(exp(x)/x + 1/log(x), x, special=True)
    Ei(x) + li(x)

    >>> pprint(risch_integrate(x*x**x*log(x) + x**x + x*x**x, x))
       x
    x*x
    >>> pprint(risch_integrate(x**x, x))
      /
     |
     |  x
     | x  dx
     |
    /

    >>> pprint(risch_integrate(-1/(x*log(x)*log(log(x))**2), x))
         1
    -----------
    log(log(x))

    This implements the integration procedure outlined in Section 5.2 of
    Bronstein's book (there is no single corresponding pseudocode
    function).
    """
    f = S(f)

    DE = extension or DifferentialExtension(f, x, handle_first=handle_first,
            dummy=True, rewrite_complex=rewrite_complex, trig=trig)
    fa, fd = DE.fa, DE.fd

    result = S.Zero
    for case in reversed(DE.cases):
        if not fa.expr.has(DE.t) and not fd.expr.has(DE.t) and not case == 'base':
            DE.decrement_level()
            fa, fd = frac_in((fa, fd), DE.t)
            continue

        fa, fd = fa.cancel(fd, include=True)
        if case == 'exp':
            ans, i, b = integrate_hyperexponential(fa, fd, DE, conds=conds,
                special=special)
        elif case == 'primitive':
            ans, i, b = integrate_primitive(fa, fd, DE, special=special)
        elif case == 'tan':
            ans, i, b = integrate_hypertangent(fa, fd, DE)
        elif case == 'base':
            # XXX: We can't call ratint() directly here because it doesn't
            # handle polynomials correctly.
            ans = integrate(fa.as_expr()/fd.as_expr(), DE.x, risch=False)
            b = False
            i = S.Zero
        else:
            raise NotImplementedError("Only exponential, logarithmic and "
            "tangent extensions are currently supported.")

        result += ans
        if b:
            DE.decrement_level()
            fa, fd = frac_in(i, DE.t)
        else:
            result = result.subs(DE.backsubs)
            if not i.is_zero:
                i = NonElementaryIntegral(i.function.subs(DE.backsubs),i.limits)
            if not separate_integral:
                result += i
                return result
            else:

                if isinstance(i, NonElementaryIntegral):
                    return (result, i)
                else:
                    return (result, 0)
