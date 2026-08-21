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
from sympy.core.function import Function, Lambda
from sympy.core.mul import Mul
from sympy.core.intfunc import ilcm
from sympy.core.numbers import Float, I, Rational, oo, zoo
from sympy.core.power import Pow
from sympy.core.relational import Ne
from sympy.core.singleton import S
from sympy.core.sorting import ordered, default_sort_key
from sympy.core.symbol import Dummy, Symbol
from sympy.functions.elementary.exponential import log, exp
from sympy.functions.elementary.hyperbolic import (cosh, coth, sinh,
    tanh)
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.elementary.trigonometric import (atan, sin, cos,
    tan, acot, cot, asin, acos)
from .integrals import integrate, Integral
from .heurisch import _symbols
from .rationaltools import log_to_real
from sympy.polys.rationaltools import together
from sympy.polys.polyerrors import NotInvertible, PolynomialError
from sympy.polys.polytools import (real_roots, cancel, Poly, gcd, rem, factor,
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
    - exts: The type ('exp' or 'log') of each extension; exts[i]
      describes T[i + 1] (T[0] == x is not an extension).
    - extargs: The argument of the exp or log of each extension, indexed
      like exts.
    - cases: List of string representations of the cases of T.
    - t: The top level extension variable, as defined by the current level
      (see level below).
    - d: The top level extension derivation, as defined by the current
      derivation (see level below).
    - case: The string representation of the case of self.d.
    - transcendental: Whether all the extensions in the tower are transcendental
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
    __slots__ = ('f', 'origf', 'x', 'T', 'D', 'fa', 'fd', 'Tfuncs',
            'backsubs', 'exts', 'extargs', 'cases', 'case', 'transcendental',
            't', 'd', 'newf', 'level', 'ts', 'dummy', 'algebraic',
            'sign_consts')

    def __init__(self, f=None, x=None, handle_first='log',
            dummy=False, extension=None, rewrite_complex=None,
            algebraic=False):
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

        self.algebraic = algebraic

        if extension:
            if 'D' not in extension:
                raise ValueError("At least the key D must be included with "
                                 "the extension flag to DifferentialExtension.")
            for attr in extension:
                setattr(self, attr, extension[attr])
            if 'transcendental' not in extension:
                # A manual extension may declare transcendental=False to
                # mark an algebraic tower (e.g. t representing x**(S(1)/2)
                # as exp(log(x)/2)); nonelementary conclusions are then
                # invalid and are degraded to unevaluated Integrals.
                self.transcendental = True

            self._auto_attrs()

            return
        elif f is None or x is None:
            raise ValueError("Either both f and x or a manual extension must "
                             "be given.")

        if handle_first not in ('log', 'exp'):
            raise ValueError("handle_first must be 'log' or 'exp', not %s." %
                             str(handle_first))

        # f will be the original function, self.f might change if we reset
        # (e.g., we pull out a constant from an exponential); origf always
        # stays the expression the user handed in
        self.f = f
        self.origf = f
        self.x = x
        # setting the default value 'dummy'
        self.dummy = dummy
        self.reset()
        exp_new_extension, log_new_extension = True, True

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
        else:
            if any(i.has(x) for i in self.f.atoms(sin, cos, cot, tan, sinh,
                    cosh, coth, tanh, asin, acos, acot, atan)):
                raise NotImplementedError("Trigonometric and hyperbolic "
                    "extensions are not supported (yet!).  Try rewriting in "
                    "terms of exp and log, or using rewrite_complex=True.")

        exps = set()
        pows = set()
        numpows = set()
        sympows = set()
        logs = set()
        symlogs = set()

        while True:
            if self.newf.is_rational_function(*self.T):
                break

            if not exp_new_extension and not log_new_extension:
                # We couldn't find a new extension on the last pass, so I guess
                # we can't do it.
                raise NotImplementedError("Couldn't find an elementary "
                                          "transcendental extension for %s.  Try using a " % str(f) +
                                          "manual extension with the extension flag.")

            exps, pows, numpows, sympows, log_new_extension = \
                self._rewrite_exps_pows(exps, pows, numpows, sympows, log_new_extension)

            logs, symlogs = self._rewrite_logs(logs, symlogs)

            if handle_first == 'exp' or not log_new_extension:
                exp_new_extension = self._exp_part(exps)
                if exp_new_extension is None:
                    # Reset and restart.  reset() clears backsubs, so any
                    # rewrite recorded there -- a branch-constant Dummy
                    # from _log_part(), a user radical's notation -- must
                    # first be folded back into newf, or the Dummy leaks
                    # into the final result as a free symbol and the
                    # original notation is never re-encountered by the
                    # rebuild.  The branch-ratio constants are the
                    # exception: their Dummys stay in newf as opaque
                    # constants, so their mapping must instead survive
                    # the reset and be re-recorded.
                    sign_consts = self.sign_consts
                    ratio_dummies = {s for s, _, _ in sign_consts}
                    self.newf = self.newf.subs(
                        [(a, b) for a, b in self.backsubs
                         if a not in ratio_dummies])
                    self.f = self.newf
                    self.reset()
                    self.sign_consts = sign_consts
                    self.backsubs += [(s, R) for s, R, _ in sign_consts]
                    if sign_consts:
                        self.transcendental = False
                    exp_new_extension = True
                    continue

            if handle_first == 'log' or not exp_new_extension:
                log_new_extension = self._log_part(logs)

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
        # exp(u)**q == exp(q*u) is exact for integer q, but for
        # fractional q the left side is a principal root that differs
        # from exp(q*u) by a locally constant root of unity off the real
        # line.  Radicals of exponentials also arise internally (the
        # radical case of _exp_part substitutes t**(p/n) and restarts);
        # those denote a root the construction is free to choose, so
        # folding them to exp(q*u) is exact, and rewriting the answer
        # back into principal-root notation is what turned correctly
        # integrated results into non-antiderivatives at complex points.
        # A fractional power the user wrote, however, means the
        # principal root specifically: fold it as ratio*exp(q*u) with an
        # opaque locally constant ratio (a root of unity on each
        # component), restored exactly on backsubstitution -- which also
        # keeps integrands mixing exp(u)**q with exp(q*u) itself
        # pointwise correct.
        subs_map = {}
        for i, j in ratpows_repl:
            if i.exp.is_Integer:
                self.backsubs.append((j, i))
                subs_map[i] = j
            elif self.origf is not None and self.origf.has(i):
                ratio = Dummy('exp_branch')
                self.backsubs.append((ratio, i/j))
                subs_map[i] = ratio*j
            else:
                subs_map[i] = j
        self.newf = self.newf.xreplace(subs_map)

        if self.algebraic:
            # Factor radicands and distribute the fractional power over
            # the factors, so that content like a perfect square does
            # not become a spurious tower generator.  Distributing a
            # fractional power is not branch-faithful (sqrt(u*v) ==
            # sqrt(u)*sqrt(v) fails when u and v are both negative), so
            # unless every factor is provably nonnegative the split is
            # multiplied by a fresh constant s standing for the exact
            # branch ratio original/split.  The ratio is locally
            # constant away from branch points (its logarithmic
            # derivative vanishes identically) and satisfies s**q == 1,
            # so s really is a constant of the extension; substituting
            # it back after integration restores the principal branch
            # everywhere.
            reps = {}
            for i in self.newf.atoms(Pow):
                if (i.exp.is_Rational and not i.exp.is_Integer and
                        not i.base.is_Symbol and not isinstance(i.base, exp)):
                    bf = factor(i.base)
                    pd = bf.as_powers_dict()
                    if len(pd) > 1 or any(e != 1 for e in pd.values()):
                        # factor() is not canonical about factor signs
                        # (they can depend on the symbol names), which
                        # would make structurally equal radicands split
                        # differently; normalize each factor to a
                        # positive leading coefficient in x.  The branch
                        # ratio absorbs any sign choice, so this is
                        # purely a canonicalization.
                        norm = {}
                        negexp = S.Zero
                        for b, e in pd.items():
                            if b.is_Add:
                                try:
                                    lc = Poly(b, self.x).LC()
                                except PolynomialError:
                                    lc = None
                                if lc is not None and lc.is_negative:
                                    b = -b
                                    negexp += e
                            norm[b] = norm.get(b, S.Zero) + e
                        if negexp:
                            norm[S.NegativeOne] = norm.get(S.NegativeOne,
                                S.Zero) + negexp
                        pd = norm
                        split = Mul(*[Pow(b, e*i.exp) for b, e in
                                      pd.items()])
                        if all(b.is_nonnegative for b in pd):
                            self.backsubs.append((split, i))
                            reps[i] = split
                        else:
                            s = Dummy('s%d' % len(self.sign_consts))
                            self.sign_consts.append((s, i/split, i.exp.q))
                            self.backsubs.append((s, i/split))
                            reps[i] = s*split
                            # results are only generic in s (it is not a
                            # transcendental constant but a root of
                            # unity), so even if every radical collapses
                            # to integer powers the candidates need the
                            # acceptance filter and nonelementary
                            # conclusions are invalid
                            self.transcendental = False
            if reps:
                self.newf = self.newf.xreplace(reps)

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
                           (i.exp.has(*self.T) or (self.algebraic and
                            i.exp.is_Rational and not i.exp.is_Integer)))
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
                    if not self.algebraic:
                        raise NotImplementedError("Algebraic extensions are "
                            "not supported (%s).  Rerun with algebraic=True "
                            "to represent it as an exp-log tower "
                            "(experimental)." % str(i))
                    self.transcendental = False
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
                    if new != old:
                        # new == old happens when exp auto-simplifies the
                        # rewritten form straight back (e.g.
                        # exp(log(x)/2) -> sqrt(x) on the first pass,
                        # before log(x) is a tower symbol); an identity
                        # pair would just pollute backsubs.
                        self.backsubs += [(new, old)]
                    log_new_extension = self._log_part([log(i.base)])
                    exps = update_sets(exps, self.newf.atoms(exp), lambda i:
                                       i.exp.is_rational_function(*self.T) and i.exp.has(*self.T))
                    continue
                ans, u, const = A
                newterm = exp(i.exp*(log(const) + u))
                if i.exp.is_Rational and not i.exp.is_Integer:
                    # log(base) can differ from log(const) + u by
                    # 2*pi*I times an integer winding (log(-1) above is
                    # the principal choice), so this rewrite needs the
                    # same branch-ratio correction as the radicand
                    # split: e.g. it maps sqrt(-w) to I*sqrt(w), which
                    # has the wrong sign where w < 0.  Both i and
                    # newterm live in tower symbols here, and Tfuncs
                    # images can reference lower generators, so the
                    # substitutions must be sequential, highest first.
                    tower_subs = list(zip(reversed(self.T),
                        reversed([f(self.x) for f in self.Tfuncs])))
                    ratio = (i/newterm).subs(tower_subs)
                    if ratio != 1:
                        s = Dummy('s%d' % len(self.sign_consts))
                        self.sign_consts.append((s, ratio, i.exp.q))
                        self.backsubs.append((s, ratio))
                        self.transcendental = False
                        newterm = s*newterm
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
            if new != old:
                self.backsubs.append((new, old))
            self.newf = self.newf.xreplace({old: newterm})
            if isinstance(newterm, exp):
                exps.append(newterm)
            elif isinstance(newterm, Mul):
                # a branch-ratio constant may ride along; only the exp
                # factor belongs in the worklist
                exps.extend(a for a in newterm.args if isinstance(a, exp))

        return exps, pows, numpows, sympows, log_new_extension

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
                        # If we're building a non-transcendental extension,
                        # don't restart - just continue to avoid infinite loops
                        if self.transcendental:
                            rad = Mul(*[term**(power/n) for term, power in ans])
                            self.newf = self.newf.xreplace({exp(p*exparg):
                                                            exp(const*p)*rad for exparg, p in others})
                            self.newf = self.newf.xreplace(dict(list(zip(reversed(self.T),
                                                                         reversed([f(self.x) for f in self.Tfuncs])))))
                            restart = True
                            break
                        # else: fall through to add extension normally
                    elif not self.algebraic:
                        # TODO: give algebraic dependence in error string
                        raise NotImplementedError("Cannot integrate over "
                            "algebraic extensions.")

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
                # u + log(const) equals log(arg) only up to a locally
                # constant branch term: log(const) is the principal-branch
                # choice, valid where every argument involved is positive
                # (e.g. log(-x**2) == 2*log(x) + I*pi only holds for
                # x > 0).  The one case where it is exact everywhere is a
                # single tower logarithm with coefficient one and a
                # positive constant: log(c*w) == log(c) + log(w) for
                # c > 0 on the whole complex plane.
                if (u is not self.x and u in self.T and const.is_positive and
                        self.exts[self.T.index(u) - 1] == 'log'):
                    self.newf = self.newf.xreplace({log(arg): log(const) + u})
                    continue
                # Otherwise the difference log(arg) - u is locally
                # constant wherever the functions are defined, so
                # integrating with an opaque constant in its place and
                # restoring the exact difference on backsubstitution keeps
                # the answer, now expressed through the original
                # logarithm, correct on every connected component.
                branch_const = Dummy('log_branch')
                # arg and u may involve tower symbols; express the
                # difference through the concrete functions (and through
                # any originals recorded in backsubs so far, so that
                # e.g. x**x reappears as x**x rather than exp(x*log(x))).
                concrete = list(zip(reversed(self.T),
                    reversed([f(self.x) for f in self.Tfuncs])))
                diff_expr = (log(arg.subs(concrete)) -
                    u.subs(concrete)).subs(self.backsubs)
                self.backsubs.append((branch_const, diff_expr))
                self.newf = self.newf.xreplace({log(arg): branch_const + u})
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
        self.transcendental = True
        self.level = -1
        self.exts = []
        self.extargs = []
        if self.dummy:
            self.ts = numbered_symbols('t', cls=Dummy)
        else:
            # For testing
            self.ts = numbered_symbols('t')
            # For various things that we change to make things work that we need to
            # change back when we are done.
        self.backsubs = []
        self.sign_consts = []
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

        list: A list of indices into T (and D) of the extensions of
            type 'extension'.  Note that self.exts[i] describes the
            extension T[i + 1], since T[0] == x is not an extension.

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
        return [i for i, ext in enumerate(self.exts, 1) if ext == extension]

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

    ans = Poly(t_part, t, z, expand=False)

    if remainder:
        one = remainder.one
        tp = t*one
        r = pd.degree() - remainder.degree()
        z_part = remainder.transform(one, tp) * tp**r
        z_part = z_part.replace(t, z).to_field().quo_ground(pd.LC())
        ans += Poly(z_part, t, z, expand=False)

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

    for pi, i in p_sqf:
        Si = pi.as_poly(*kkinv).gcd(derivation(pi, DE,
            coefficientD=coefficientD,basic=basic).as_poly(*kkinv)).as_poly(DE.t)
        pi = Poly(pi, DE.t)
        Si = Poly(Si, DE.t)
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
    dm = Poly(gcd(d.to_field(), dd.to_field()), DE.t)
    ds, _ = d.div(dm)

    while dm.degree(DE.t) > 0:

        ddm = derivation(dm, DE)
        dm2 = gcd(dm.to_field(), ddm.to_field())
        dms, _ = dm.div(dm2)
        ds_ddm = ds.mul(ddm)
        ds_ddm_dm, _ = ds_ddm.div(dm)

        b, c = gcdex_diophantine(-Poly(ds_ddm_dm, DE.t),
            Poly(dms, DE.t), Poly(a, DE.t))
        b, c = Poly(b, DE.t), Poly(c, DE.t)

        db = Poly(derivation(b, DE), DE.t)
        ds_dms, _ = ds.div(dms)
        a = Poly(c, DE.t) - Poly(db.mul(ds_dms), DE.t)

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

    Returns
    =======

    (delta_a, delta_d, H) : tuple
        ``delta_a`` and ``delta_d`` are Polys in ``DE.t`` such that
        ``delta_a/delta_d`` is the sum of the principal parts of the
        Laurent series of A/D at the real zeros of ``F``.  ``H`` is the
        list of Polys ``[H_1, ..., H_n]`` from Theorem 2.7.1: at each
        zero ``alpha`` of ``F``, ``H_j(alpha)`` is the coefficient of
        ``1/(t - alpha)**(n - j + 1)`` in the principal part of A/D at
        ``alpha``; in particular, ``H[-1]`` evaluated at ``alpha`` is the
        residue of A/D at ``alpha``.  If ``F`` has degree 0 (and so has
        no zeros), ``(Poly(0, DE.t), Poly(1, DE.t), [])`` is returned.

    This is ``LaurentSeries`` from Section 2.7 of Bronstein's book.  The
    book's version returns only the sum of the principal parts (so the
    degree 0 case returns the expression 0); ``H`` is additionally
    returned here for use by ``recognize_derivative()``.
    """
    if F.degree() == 0:
        return (Poly(0, DE.t), Poly(1, DE.t), [])
    Z: list[Symbol] = [*_symbols('z', n)]
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
    for s, n in Sp:
        delta_a, delta_d, H = laurent_series(r, d, s, n, DE)
        if not H[-1].rem(s.as_poly(DE.t)).is_zero:  # Di does not divide H[-1]
            # A conclusive False from a special factor stands regardless
            # of any undecidable normal factors, so check these first.
            flag = False
            break
    else:
        if Np:
            # The Laurent series machinery above is only justified for
            # factors with constant roots (Theorem 2.7.1 is stated for
            # K[x] with the roots algebraic over the constant field K);
            # deciding whether the residues vanish at nonconstant poles
            # of order > 1 is not yet implemented.
            raise NotImplementedError("recognize_derivative() cannot decide "
                "residue vanishing at nonconstant poles of order greater "
                "than 1.")
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

    if any(Poly(s, z).degree() > 0 for s, _ in Np):
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
    z = z or Dummy('z')
    a, d = a.cancel(d, include=True)
    a, d = a.to_field().mul_ground(1/d.LC()), d.to_field().mul_ground(1/d.LC())

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
            s = Poly(s, z).monic()
            h_lc = Poly(h.as_poly(DE.t).LC(), z, field=True)

            h_lc_sqf = h_lc.sqf_list_include(all=True)

            for a, j in h_lc_sqf:
                g = a.gcd(s)
                h = Poly(h, DE.t, field=True).exquo(Poly(g.as_expr()**j, DE.t))

            if invert:
                h_lc = Poly(Poly(h, DE.t).LC(), DE.t, field=True, expand=False)
                try:
                    inv = Poly(h_lc, z, field=True).invert(s)
                except NotInvertible:
                    if DE.transcendental:
                        raise
                    # The leading coefficient shares a factor with this
                    # resultant factor, so this residue term cannot be
                    # normalized.  Give the residues up: the uncaptured
                    # mass stays in the remainder f - Dg, and b == False
                    # sends the level down the give-up path, whose
                    # nonelementary conclusion is degraded to an honest
                    # unevaluated Integral for non-transcendental towers.
                    return (H, False)
                coeffs = [S.One]

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

    Terms whose residues can all be computed explicitly are rewritten as
    real logarithms and arc-tangents using log_to_real() (Rioboo's
    algorithm from Section 2.8 of Bronstein's book, with the roots of
    each s_i as the residues and S_i in place of the Rothstein-Trager
    resultant's remainder); the remaining terms are returned as RootSums.
    log_to_real() is only valid over a real field, so it is not used when
    I appears in a term or in the extension tower.
    """
    # TODO: check what Lambda does with RootOf
    i = Dummy('i')
    s = list(zip(reversed(DE.T), reversed([f(DE.x) for f in DE.Tfuncs])))
    real_tower = not any(f(DE.x).has(I) for f in DE.Tfuncs)

    result = S.Zero
    for a in H:
        real = None
        if real_tower and not (a[0].as_expr().has(I) or a[1].as_expr().has(I)):
            real = log_to_real(a[1], a[0].as_poly(z), DE.t, z)
        if real is not None:
            result += real.subs(s)
        else:
            result += RootSum(a[0].as_poly(z), Lambda(i, i*log(
                a[1].as_expr()).subs({z: i}).subs(s)))

    return result


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


def _nontrans_power_relations(DE):
    """
    The known algebraic relations of a non-transcendental tower.

    Returns (t, rel) pairs, ordered by tower level, where rel == 0 in
    the actual (evaluated) tower and rel is polynomial in t with its
    leading coefficient free of t and of every higher generator.  Only
    generators of the form t == exp(c*tj) with c == p/q a non-integer
    Rational and tj == log(u) are recognized; such a t represents
    u**(p/q) and satisfies t**q == u**p.

    TODO: relations among *several* radicals are not captured -- branch
    choices like t1*t2 == 1 for t1 == exp(t0/2), t2 == exp(-t0/2), or
    radicals with dependent arguments.  The kernel test below is then
    incomplete (it can miss actual zeros), though still sound.

    The branch-ratio constants introduced by the radicand splitting are
    included through their root-of-unity relations s**q == 1 (so that a
    denominator containing a factor of s**q - 1, formally nonzero but
    actually zero, is caught).
    """
    rels = [(s, s**q - 1) for s, _, q in DE.sign_consts or []]
    for i, ext in enumerate(DE.exts, 1):
        if ext != 'exp':
            continue
        coeff, rest = DE.extargs[i - 1].as_coeff_Mul()
        if (coeff.is_Rational and not coeff.is_Integer and rest in DE.T[1:]
                and DE.exts[DE.T.index(rest) - 1] == 'log'):
            u = DE.extargs[DE.T.index(rest) - 1]
            p, q = coeff.p, coeff.q
            if p >= 0:
                rel = DE.T[i]**q - u**p
            else:
                rel = DE.T[i]**q*u**-p - 1
            rels.append((DE.T[i], rel))
    return rels


def _nontrans_branch_corrections(result, DE):
    """
    Substitute the branch-ratio constants back into ``result`` and
    correct the jumps their substitution introduces.

    The candidate was found with each ratio held as a formal constant
    s, so ``result`` with s substituted is an antiderivative on every
    connected region where the ratio really is constant, but it may
    jump by a constant where the ratio changes value (the real roots
    and poles of the radicand factors).  Following Jeffrey, Labahn,
    von Mohrenschildt & Rich (Theorem 5), subtracting
    J*(x - r)/sqrt((x - r)**2) with J == (G(r+) - G(r-))/2 at each
    such breakpoint r restores continuity on the real line, extending
    the answer to Jeffrey's domain of maximum extent.  The correction
    factor equals sign(x - r) for real x but, unlike sign, is locally
    constant on the complex plane away from the vertical line
    Re(x) == r, so the answer stays an antiderivative at complex
    points too (off that line -- the same status as any branch cut).
    The corrections are locally constant, so they never affect the
    derivative; every step here is therefore free to give up
    (unlocatable breakpoints, an unrecognizable ratio value, a
    singular or infinite jump, an uncertifiable sampling offset),
    leaving the plain substitution, which is already correct on each
    region.
    """
    from sympy.core.numbers import pi

    x = DE.x
    sign_subs = {s: R for s, R, _ in DE.sign_consts}
    subbed = result.subs(sign_subs)
    if not any(result.has(s) for s in sign_subs):
        return subbed

    # Breakpoints: real roots of every factor under a ratio (numerator
    # and denominator -- the ratios jump only where some factor crosses
    # zero or infinity).  The set must be *complete*: a missing root
    # would let a side sample cross an unrecognized breakpoint and
    # certify a nonlocal ratio value, so real_roots() (exact isolation)
    # is used and every factor must yield its full real-root set or no
    # correction is attempted at all.
    from sympy.polys.polyerrors import CoercionFailed, DomainError
    pts = set()
    for _, R, _ in DE.sign_consts:
        pieces = [p.base for p in R.atoms(Pow)]
        pieces.extend(lg.args[0] for lg in R.atoms(log))
        for piece in pieces:
            for b in piece.as_numer_denom():
                if not b.has(x):
                    continue
                if b.free_symbols != {x}:
                    return subbed
                try:
                    pts.update(Poly(b, x).real_roots())
                except (PolynomialError, DomainError, CoercionFailed,
                        NotImplementedError):
                    return subbed
    if not pts:
        return subbed
    pts = sorted(pts, key=lambda r: r.evalf(30))

    def ratio_at(R, q, pt):
        # the ratio's value at an exact off-breakpoint point, snapped to
        # an exact q-th root of unity (None if it isn't recognizably one)
        try:
            v = complex(R.subs(x, pt).evalf(30))
        except (TypeError, ValueError):
            return None
        best = None
        for j in range(q):
            w = exp(2*I*pi*Rational(j, q))
            d = abs(v - complex(w.evalf(30)))
            if d < 0.01 and (best is None or d < best[1]):
                best = (w, d)
        return best[0] if best else None

    corrections = []
    for idx, r in enumerate(pts):
        # an exactly-certified sampling offset: r -+ d must not cross
        # the neighboring breakpoints, or the side values are nonlocal
        prev_ = pts[idx - 1] if idx else None
        next_ = pts[idx + 1] if idx + 1 < len(pts) else None
        for d in (Rational(1, 2)**k for k in range(2, 64)):
            if ((prev_ is None or (r - d - prev_).is_positive) and
                    (next_ is None or (next_ - r - d).is_positive)):
                break
        else:
            continue
        left, right = {}, {}
        ok = True
        for s, R, q in DE.sign_consts:
            vl = ratio_at(R, q, r - d)
            vr_ = ratio_at(R, q, r + d)
            if vl is None or vr_ is None:
                ok = False
                break
            left[s], right[s] = vl, vr_
        if not ok or left == right:
            continue
        # cancel before substituting the breakpoint: singularities of
        # the candidate there are often removable, but only in the
        # s-form (the cancellation needs s*w*sqrt(v) still split)
        Gl = cancel(together(result.subs(left))).subs(x, r)
        Gr = cancel(together(result.subs(right))).subs(x, r)
        J = cancel(together((Gr - Gl)/2))
        if J == 0 or J.free_symbols or J.is_real is not True:
            # a complex jump would make the answer complex on whole
            # regions where it is real now; an infinite or undecidable
            # one cannot be corrected by a constant.  Leave the jump:
            # the result is still an antiderivative on each region.
            continue
        corrections.append(J*(x - r)*Pow((x - r)**2, Rational(-1, 2)))
    return subbed - Add(*corrections)


def _nontrans_is_kernel(e, DE):
    """
    Decide whether e, polynomial in the tower generators, evaluates to
    zero in a non-transcendental tower (i.e. lies in the kernel of the
    evaluation map), by reducing modulo the recognized power relations,
    highest generator first.  The relations are triangular (each
    involves only its own and lower generators), so this terminates,
    and e evaluates to zero iff the reduced form is zero -- exactly, if
    the quotient by the recognized relations is an integral domain
    (always true for a single radical), and soundly-but-incompletely
    otherwise (see _nontrans_power_relations()).
    """
    e = cancel(e)
    for t, r in reversed(_nontrans_power_relations(DE)):
        if e.has(t):
            e = rem(e, r, t)
    return cancel(e) == 0


def _nontrans_sign_assignments(DE, exprs):
    """
    Every assignment of the branch-ratio constants occurring in
    ``exprs`` to their possible root-of-unity values, as a list of
    substitution dicts ([{}] when none occur), or None if the
    assignments cannot be enumerated in exact
    radical-free-of-transcendentals form (a root of unity of high
    order whose cos/sin do not evaluate, or too many combinations) --
    the caller must then reject the candidate, since the per-branch
    check cannot be run.  Constants absent from ``exprs`` need no
    assignment: substituting them cannot change the checked
    expressions.
    """
    from itertools import product

    from sympy.core.numbers import pi

    consts = [sc for sc in DE.sign_consts or []
              if any(e.has(sc[0]) for e in exprs)]
    values = []
    for _, _, q in consts:
        vals = []
        for j in range(q):
            v = (cos(2*pi*Rational(j, q)) + I*sin(2*pi*Rational(j, q)))
            if v.has(sin, cos, exp):
                return None
            vals.append(v)
        values.append(vals)
    n = 1
    for vals in values:
        n *= len(vals)
        if n > 64:
            return None
    return [dict(zip([s for s, _, _ in consts], combo))
            for combo in product(*values)]


def _nontrans_vet(result, DE):
    """
    Final check of a sign-carrying result before the ratio constants
    are substituted back, returning False if the result must be
    rejected.

    The acceptance filter runs per tower level, but the base case
    (everything collapsed to a rational function of x and the ratio
    constants) never reaches it, and ratint() over QQ(s) can divide by
    s-polynomials that vanish at attained values.  So the total result
    is vetted under every possible assignment of the ratio constants:
    the denominator must stay out of the relation kernel, the leading
    coefficients of RootSum defining polynomials must not vanish
    (specialization would silently change the root set, which no
    value check can see), and the assigned result must not evaluate
    to nan or an infinity.

    Only constants in positions that can turn singular need
    assigning: denominators at any depth (the top-level one, and
    those inside function arguments and fractional-power bases, which
    as_numer_denom() does not see into), arguments of non-entire
    functions (log(0) is zoo, but no finite exp() argument is
    singular), exponents, and RootSum defining polynomials.  A
    constant occurring purely as a numerator polynomial factor cannot
    produce a singularity.  Nested structures are covered by the
    atoms() walks being recursive.
    """
    den = cancel(result).as_numer_denom()[1]
    risky = [den]
    # entire functions: no finite argument is singular (tan, cot,
    # coth have poles and atan is singular at +-I, so they stay risky)
    entire = (exp, sin, cos, sinh, cosh)
    for fn in result.atoms(Function):
        if isinstance(fn, entire):
            risky.extend(cancel(arg).as_numer_denom()[1]
                         for arg in fn.args)
        else:
            risky.extend(fn.args)
    for p in result.atoms(Pow):
        risky.append(p.exp)
        if not p.exp.is_Integer:
            risky.append(cancel(p.base).as_numer_denom()[1])
    rootsums = list(result.atoms(RootSum))
    risky.extend(rs.poly.as_expr() for rs in rootsums)
    assignments = _nontrans_sign_assignments(DE, risky)
    if assignments is None:
        return False
    for sig in assignments:
        if _nontrans_is_kernel(den.subs(sig), DE):
            return False
        if any(_nontrans_is_kernel(rs.poly.LC().subs(sig), DE)
               for rs in rootsums):
            return False
        if result.subs(sig).has(S.NaN, oo, zoo):
            return False
    return True


def _nontrans_accept(elem, g2, i, a, d, DE, z):
    """
    Acceptance filter for a candidate result of integrate_primitive()
    or integrate_hyperexponential() over a non-transcendental tower:

    1. a/d == D(elem) + D(residue part) + i must hold as a formal
       rational-function identity in the tower generators, decided by
       cancel() -- pure rational arithmetic, no algebraic
       simplification.
    2. The denominators of elem and i must not evaluate to zero under
       the tower's algebraic relations (_nontrans_is_kernel()).

    Why this suffices: the machinery computes in the formal field
    Q(x, t0, ..., tn), where the generators really are transcendental,
    and the actual functions are the image of an evaluation map that is
    defined only on elements whose denominators avoid the kernel of the
    tower's relations (t**2 - x for t representing sqrt(x)).  A formal
    identity between elements on which the evaluation is defined pushes
    forward to an identity of functions, since the evaluation commutes
    with the derivation by construction.  So a candidate passing both
    checks is a genuine antiderivative relation, no matter what invalid
    transcendence assumptions were used to find it.  The checks cannot
    vouch for *negative* conclusions (nonelementary proofs), which is
    why those are separately degraded to unevaluated Integrals.

    Check 1 catches machinery bugs surfacing under the broken
    hypotheses; check 2 catches the subtler failure mode: the formal
    field has genuinely new constants (D(t**2/x) == 0 for t
    representing sqrt(x)), so internal divisions by formally-nonzero
    constants can smuggle kernel factors (t**2/x - 1, which evaluates
    to zero) into denominators while keeping the formal identity
    intact.

    The residue terms are checked through the coefficient denominators
    of their root polynomials and logarithm arguments (a kernel there
    means an infinite residue or an undefined argument; the derivative
    alone cannot show it, since the kernel factors cancel against the
    argument's derivative).  TODO: a logarithm argument that itself
    evaluates to zero is not detected here.

    The branch-ratio constants need more than the generic kernel test:
    s**q - 1 is reducible, so its quotient has zero divisors, and a
    denominator like s - 1 is formally nonzero yet the actual ratio
    equals 1 on whole regions.  So every denominator is additionally
    checked under every assignment of the ratio constants to their
    possible root-of-unity values (conservative: assignments the
    ratios never attain are not excluded, so this can over-reject,
    never under-reject).
    """
    lhs = cancel(a.as_expr()/d.as_expr())
    total = cancel(lhs - derivation(elem, DE, basic=True) -
        residue_reduce_derivation(g2, DE, z) - i)
    if total != 0:
        return False
    checked = [elem, i]
    for qz, sz in g2:
        checked.extend([qz.as_expr(), sz.as_expr()])
    dens = [cancel(e).as_numer_denom()[1] for e in checked]
    assignments = _nontrans_sign_assignments(DE, dens)
    if assignments is None:
        return False
    for den in dens:
        for sig in assignments:
            if _nontrans_is_kernel(den.subs(sig), DE):
                return False
    return True


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


def integrate_primitive(a, d, DE, z=None):
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

    This is ``IntegratePrimitive`` from Section 5.8 of Bronstein's book.
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
        if not DE.transcendental and not _nontrans_accept(
                g1[0].as_expr()/g1[1].as_expr(), g2, i, a, d, DE, z):
            return (S.Zero, Integral(cancel(
                (a.as_expr()/d.as_expr()).subs(s)), DE.x), False)
        if DE.transcendental:
            i = NonElementaryIntegral(cancel(i).subs(s), DE.x)
        else:
            i = Integral(cancel(i).subs(s), DE.x)
        return ((g1[0].as_expr()/g1[1].as_expr()).subs(s) +
            residue_reduce_to_basic(g2, DE, z), i, b)

    # h - Dg2 + r
    p = cancel(h[0].as_expr()/h[1].as_expr() - residue_reduce_derivation(g2,
        DE, z) + r[0].as_expr()/r[1].as_expr())
    p = p.as_poly(DE.t)

    q, i, b = integrate_primitive_polynomial(p, DE)
    if not DE.transcendental and not _nontrans_accept(
            g1[0].as_expr()/g1[1].as_expr() + q.as_expr(), g2, i.as_expr(),
            a, d, DE, z):
        return (S.Zero, Integral(cancel(
            (a.as_expr()/d.as_expr()).subs(s)), DE.x), False)

    ret = ((g1[0].as_expr()/g1[1].as_expr() + q.as_expr()).subs(s) +
        residue_reduce_to_basic(g2, DE, z))
    if not b:
        # i == p - Dq from integrate_primitive_polynomial(), which has been
        # proven to have no elementary integral over k(t), and ret contains
        # the partial q, so f == D(ret) + i holds for the returned values.
        # The nonelementary proof is only valid for transcendental towers.
        if DE.transcendental:
            i = NonElementaryIntegral(cancel(i.as_expr()).subs(s), DE.x)
        else:
            i = Integral(cancel(i.as_expr()).subs(s), DE.x)
    else:
        i = cancel(i.as_expr())

    return (ret, i, b)


def integrate_hyperexponential_polynomial(p, DE, z):
    """
    Integration of hyperexponential polynomials.

    Explanation
    ===========

    Given a hyperexponential monomial t over k and ``p`` in k[t, 1/t], return q in
    k[t, 1/t] and a bool b in {True, False} such that p - Dq in k if b is True,
    or p - Dq does not have an elementary integral over k(t) if b is False.

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
            else:
                # q += v*t**i
                if i > 0:
                    ti = Poly(t1**i, t1)
                else:
                    ti = Poly(z**-i, z)

                qa = qa*vd + va*ti*qd
                qd *= vd

    return (qa, qd, b)


def integrate_hyperexponential(a, d, DE, z=None, conds='piecewise'):
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
        if not DE.transcendental and not _nontrans_accept(
                g1[0].as_expr()/g1[1].as_expr(), g2, i, a, d, DE, z):
            return (S.Zero, Integral(cancel(
                (a.as_expr()/d.as_expr()).subs(s)), DE.x), False)
        if DE.transcendental:
            i = NonElementaryIntegral(cancel(i.subs(s)), DE.x)
        else:
            i = Integral(cancel(i.subs(s)), DE.x)
        return ((g1[0].as_expr()/g1[1].as_expr()).subs(s) +
            residue_reduce_to_basic(g2, DE, z), i, b)

    # p should be a polynomial in t and 1/t, because Sirr == k[t, 1/t]
    # h - Dg2 + r
    p = cancel(h[0].as_expr()/h[1].as_expr() - residue_reduce_derivation(g2,
        DE, z) + r[0].as_expr()/r[1].as_expr())
    pp = as_poly_1t(p, DE.t, z)

    qa, qd, b = integrate_hyperexponential_polynomial(pp, DE, z)

    i = pp.nth(0, 0)
    if b:
        resid = i
    else:
        # the residual that the not-b branch below will return.  The
        # arithmetic is done on expressions: as Polys, qa and qd can
        # disagree about whether the residue dummy z belongs to the
        # ground domain, which fails to unify.
        resid = p - (qd.as_expr()*derivation(qa, DE).as_expr() -
            qa.as_expr()*derivation(qd, DE).as_expr())/qd.as_expr()**2

    if not DE.transcendental:
        if not _nontrans_accept(g1[0].as_expr()/g1[1].as_expr() +
                qa.as_expr()/qd.as_expr(), g2, resid, a, d, DE, z):
            return (S.Zero, Integral(cancel(
                (a.as_expr()/d.as_expr()).subs(s)), DE.x), False)

    ret = ((g1[0].as_expr()/g1[1].as_expr()).subs(s) \
        + residue_reduce_to_basic(g2, DE, z))

    qas = qa.as_expr().subs(s)
    qds = qd.as_expr().subs(s)
    if conds == 'piecewise' and DE.x not in qds.free_symbols:
        # We have to be careful if the exponent is S.Zero!

        # XXX: Does qd = 0 always necessarily correspond to the exponential
        # equaling 1?
        cond = Ne(qds, 0)
        if DE.transcendental:
            fallback = integrate((p - i).subs(DE.t, 1).subs(s), DE.x)
        elif cond is S.true:
            # the denominator provably cannot vanish; there is no
            # degenerate branch
            ret += qas/qds
            cond = None
        else:
            # t == 1 encodes "the exponential degenerates to 1", which
            # has no meaning for an algebraic generator (t representing
            # sqrt(u) does not become 1 when a parameter vanishes).  Nor
            # can the degenerate branch replace only the qas/qds piece:
            # the generic decomposition itself divides by qds, so every
            # piece of it may be invalid there.
            if not b:
                # with a residual returned separately, no branch value
                # can survive the degenerate substitution (both the
                # rest of the level and the residual keep the vanishing
                # denominator); give the level up so the caller retries
                # the original integrand whole
                return (S.Zero, Integral(cancel(
                    (a.as_expr()/d.as_expr()).subs(s)), DE.x), False)
            # everything else this level produced is inside ret, so the
            # branch is the unevaluated rest of the level minus ret,
            # and the assembled degenerate value is exactly
            # Integral(a/d - i) plus the continuing constant's integral
            fallback = Integral(cancel((a.as_expr()/d.as_expr()
                - resid).subs(s)), DE.x) - ret
        if cond is not None:
            ret += Piecewise(
                    (qas/qds, cond),
                    (fallback, True)
                )
    else:
        ret += qas/qds

    if not b:
        if DE.transcendental:
            i = NonElementaryIntegral(cancel(resid).subs(s), DE.x)
        else:
            i = Integral(cancel(resid).subs(s), DE.x)
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
                    conds='piecewise', algebraic=True):
    r"""
    The Risch Integration Algorithm.

    Explanation
    ===========

    Only transcendental functions are supported.  Currently, only exponentials
    and logarithms are supported, but support for trigonometric functions is
    forthcoming.  Radicals are handled through the ``algebraic`` flag
    described below.

    If this function returns a NonElementaryIntegral (a subclass of Integral)
    in the result, it means that it has proven that integral to be
    nonelementary.  A plain unevaluated Integral carries no such proof: it is
    returned for radical integrands (see ``algebraic`` below), where the
    nonelementary conclusions of the transcendental machinery are not valid.
    Any errors will result in raising NotImplementedError.

    handle_first may be either 'exp' or 'log'.  This changes the order in
    which the extension is built, and may result in a different (but
    equivalent) solution (for an example of this, see issue 5109).  It is also
    possible that the integral may be computed with one but not the other,
    because not all cases have been implemented yet.  It defaults to 'log' so
    that the outer extension is exponential when possible, because more of the
    exponential case has been implemented.

    If ``separate_integral`` is ``True``, the result is returned as a tuple (ans, i),
    where the integral is ans + i, ans is elementary, and i is either 0, a
    NonElementaryIntegral, or a plain unevaluated Integral (for radical
    integrands, where nonelementarity is not proven).  This is useful if
    you want to try further integrating the leftover part using other
    algorithms to possibly get a solution in terms of special functions.  It
    is False by default.

    If ``algebraic`` is ``True`` (the default), the transcendental machinery
    is used to solve integrals involving radicals (internally, by representing
    `x**(1/n)` as ``exp(log(x)/n)``). In this case, an unevaluated
    ``Integral`` result is not a proof of nonelementarity. It only means the
    transcendental algorithms aren't able to handle the radical integrand, and
    the full algebraic Risch algorithm may be required.  With
    ``algebraic=False``, radical integrands raise ``NotImplementedError``.

    Examples
    ========

    >>> from sympy.integrals.risch import risch_integrate
    >>> from sympy import exp, log, pprint
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
    with exponentials and logarithms, though note that this can include
    nested exponentials and logarithms, as well as exponentials with bases
    other than E.

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

    if algebraic and f.has(Float):
        # The algebraic towers do exact arithmetic: a Float coefficient
        # becomes a rational with an astronomical denominator
        # (0.333333333333333 == 333333333333333/10**15) and the
        # structure-theorem constant systems grind on it, while nothing
        # exact can be concluded from inexact input anyway.  Leave
        # radical integrands with Floats to the numeric-friendly
        # fallbacks in integrate().
        algebraic = False

    DE = extension or DifferentialExtension(f, x, handle_first=handle_first,
            dummy=True, rewrite_complex=rewrite_complex, algebraic=algebraic)
    fa, fd = DE.fa, DE.fd

    result = S.Zero
    for case in reversed(DE.cases):
        if not fa.expr.has(DE.t) and not fd.expr.has(DE.t) and not case == 'base':
            DE.decrement_level()
            fa, fd = frac_in((fa, fd), DE.t)
            continue

        fa, fd = fa.cancel(fd, include=True)
        if case == 'exp':
            ans, i, b = integrate_hyperexponential(fa, fd, DE, conds=conds)
        elif case == 'primitive':
            ans, i, b = integrate_primitive(fa, fd, DE)
        elif case == 'base':
            # XXX: We can't call ratint() directly here because it doesn't
            # handle polynomials correctly.
            ans = integrate(fa.as_expr()/fd.as_expr(), DE.x, risch=False)
            b = False
            i = S.Zero
        else:
            raise NotImplementedError("Only exponential and logarithmic "
            "extensions are currently supported.")

        result += ans
        if b:
            DE.decrement_level()
            fa, fd = frac_in(i, DE.t)
        else:
            sign_symbols = {s for s, _, _ in DE.sign_consts or []}
            result = result.subs([(o, n) for o, n in DE.backsubs
                                  if o not in sign_symbols])
            if sign_symbols:
                if not _nontrans_vet(result, DE):
                    if not separate_integral:
                        return Integral(f, x)
                    return (S.Zero, Integral(f, x))
                if i == 0:
                    result = _nontrans_branch_corrections(result, DE)
                else:
                    # with a leftover integral the total antiderivative
                    # is unknown, so its jumps cannot be corrected
                    result = result.subs(
                        {s: R for s, R, _ in DE.sign_consts})
            if isinstance(i, Integral):
                if DE.transcendental:
                    i = NonElementaryIntegral(i.function.subs(DE.backsubs), i.limits)
                else:
                    i = Integral(i.function.subs(DE.backsubs), *i.limits)
            leaked = result.free_symbols
            if isinstance(i, Integral):
                leaked = leaked | i.function.free_symbols
            if not DE.transcendental and \
                    leaked - f.free_symbols - {x}:
                # an internal symbol survived the back-substitutions (a
                # tower Dummy leaked into a residue term, say), in the
                # elementary part or in the residual integrand; the
                # result is unusable and must not be returned
                if not separate_integral:
                    return Integral(f, x)
                return (S.Zero, Integral(f, x))
            if not separate_integral:
                result += i
                return result
            else:
                return (result, i)
