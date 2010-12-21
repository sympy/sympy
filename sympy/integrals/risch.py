"""
The Risch Algorithm for transcendental function integration.

The core algorithms for the Risch algorithm are here.  The subproblem
algorithms are in the rde.py and prde.py files for the Risch
Differential Equation solver and the parametric problems solvers,
respectively.  Conventions here is that the base domain is QQ, x == t0, and
each differential extension is t1, t2, ..., t. x is the variable of
integration (Dx == 1), D is a list of the derivatives of x, t1, t2, ..., t,
T is the list of x, t1, t2, ..., t, t is the outer-most variable of the
differential extension, k is the field C(x, t1, ..., tn-1), where t ==
tn.  The numerator of a fraction is denoted by a and the denominator by
d.  If the fraction is named f, fa == numer(f) and fd == denom(f).
Fractions are returned as tuples (fa, fd).  d and t are often used to
represent the topmost derivation and extension variable, respectively.
The docstring of a function signifies whether an argument is in k[t], in
which case it will just return a Poly in t, or in k(t), in which case it
will return the fraction (fa, fd). Other variable names probably come
from the names used in Bronstein's book.
"""
from sympy.core.basic import S
from sympy.core.function import Lambda
from sympy.core.numbers import ilcm
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.symbol import Symbol, Dummy

from sympy.functions import log, exp, sin, cos, tan, asin, acos, atan

from sympy.integrals import Integral, integrate

from sympy.polys import gcd, cancel, PolynomialError, Poly, reduced, RootSum

from sympy.utilities.iterables import numbered_symbols, any, all
#    from pudb import set_trace; set_trace() # Debugging

class DifferentialExtension(object):
    """
    A container for all the information relating to a differential extension.

    The attributes of this object are (see also the docstring of __init__):

    - f: The original (Expr) integrand.
    - x: The variable of integration.
    - T: List of variables in the extension.
    - D: List of derivations in the extension; corresponds to the elements of T.
    - fa: Poly of the numerator of the integrand.
    - fd: Poly of the denominator of the integrand.
    - Tfuncs: Lambda() representations of each element of T, in reverse order.
      For back-substitution after integration.
    - backsubs: A (possibly empty) list of further substitutions to be made on
      the final integral to make it look more like the integrand.
    - E_K: List of the positions of the exponential extensions in T.
    - E_args: The arguments of each of the exponentials in E_K.
    - L_K: List of the positions of the logarithmic extensions in T.
    - L_args: The arguments of each of the logarithms in L_K.
    (See the docstrings of is_deriv_k() and is_log_deriv_k_t_radical() for
    more information on E_K, E_args, L_K, and L_args)
    - cases: List of string representations of the cases of T, in reverse order.
    - t: The top level extension variable, as defined by the current level
      (see level below).
    - d: The top level extension derivation, as defined by the current
      derivation (see level below).
    (Note that self.T and self.D will always contain the complete extension,
    regardless of the level.  Therefore, you should ALWAYS use DE.t and DE.d
    instead of DE.T[-1] and DE.D[-1].)

    The following are also attributes, but will probably not be useful other
    than in internal use:
    - newf: Expr form of fa/fd.
    - level: The number (between -1 and -len(self.T)) such that
      self.T[self.level] == self.t and self.D[self.level] == self.d.
      Use the methods self.increase_level() and self.decrease_level() to change
      the current level.
    """
    __slots__ = ('f', 'x', 'T', 'D', 'fa', 'fd', 'Tfuncs', 'backsubs', 'E_K',
        'E_args', 'L_K', 'L_args', 'cases', 't', 'd', 'newf', 'level', 'ts')

    def __init__(self, f, x, handle_first='log', dummy=True, extension=None):
        """
        Tries to build a transcendental extension tower from f with respect to x.

        If it is successful, creates a DifferentialExtension object with, among
        other attributes, attributes (fa, fd, D, T, Tfuncs, backsubs) such that
        fa and fd are Polys in T[-1] with rational coefficients in T[:-1],
        fa/fd == f, and D[i] is a Poly in T[i] with rational coefficients in
        T[:i] representing the derivative of T[i] for each i from 1 to len(T).
        Tfuncs is a list of Lambda objects for back replacing the functions
        after integrating.  Lambda() is only used (instead of lambda) to make
        them easier to test and debug. Note that Tfuncs corresponds to the
        elements of T (except for T[0] == x) in reverse order, because that is
        the order that they have to be back-substituted in.  backsubs is a
        (possibly empty) back-substitution list that should be applied on the
        completed integral to make it look more like the original integrand.

        If it is unsuccessful, it raises NotImplementedError.

        You can also create an object by manually setting the attributes as a
        dictionary to the extension keyword argument.  You must include at least
        D.  Warning, any attribute that is not given will be set to None. The
        attributes T, t, d, cases, and level are set automatically and do not
        need to be given.  The functions in the Risch Algorithm will NOT check
        to see if an attribute is None before using it.  This also does not
        check to see if the extension is valid (non-algebraic) or even if it is
        self-consistent.  Therefore, this should only be used for
        testing/debugging purposes
        """
        if extension:
            for attr in self.__slots__:
                setattr(self, attr, None)

            if not extension.has_key('D'):
                raise ValueError("At least the key D must be included with " +
                    "the extension flag to DifferentialExtension.")
            for attr in extension:
                setattr(self, attr, extension[attr])

            self._auto_attrs()

            return

        from sympy.integrals.prde import is_deriv_k

        if handle_first not in ['log', 'exp']:
            raise ValueError("handle_first must be 'log' or 'exp', not %s." %
                str(handle_first))

        # f will be the original function, self.f might change if we reset
        self.f = f
        self.x = x

        # Get common cases out of the way:
        if any(i.has(x) for i in self.f.atoms(sin, cos, tan, atan, asin, acos)):
            raise NotImplementedError("Trigonometric extensions are not " +
            "supported (yet!)")
        self.reset(dummy=dummy)
        exp_new_extension, log_new_extension = True, True
        while True:
            restart = False
            if self.newf.is_rational_function(*self.T):
                break

            if not exp_new_extension and not log_new_extension:
                # We couldn't find a new extension on the last pass, so I guess
                # we can't do it.
                raise NotImplementedError("Couldn't find an elementary " +
                    "transcendental extension for %s.  Try using a " % str(f) +
                    "manual extension with the extension flag.")
                # TODO: Actually implement said extension flag

            # Pre-preparsing.
            #################
            # Get all exp arguments, so we can avoid ahead of time doing something
            # like t1 = exp(x), t2 = exp(x/2) == sqrt(t1).

            exps = filter(lambda i: i.exp.is_rational_function(*self.T) and
                i.exp.has_any_symbols(*self.T), self.newf.atoms(exp))
            # 2**x
            pows = filter(lambda i: i.exp.is_rational_function(*self.T) and
                i.exp.has_any_symbols(*self.T), self.newf.atoms(Pow))
            numpows = filter(lambda i: not i.base.has_any_symbols(*self.T),
                pows)
            sympows = filter(lambda i: i.base.is_rational_function(*self.T) and
                not i.exp.is_Integer, list(set(pows) - set(numpows)))

            # The easiest way to deal with non-base E powers is to convert them
            # into base E, integrate, and then convert back.
            for i in pows:
                old = i
                # If exp is ever changed to automatically reduce exp(x*log(2))
                # to 2**x, then this will break.  The solution is to not change
                # exp to do that :)
                if i in sympows:
                    if i.exp.is_Rational:
                        raise NotImplementedError("Algebraic extensions are " +
                            "not supported (%s)." % str(i))
                    # We can add a**b only if log(a) in the extension, because
                    # a**b == exp(b*log(a)).
                    basea, based = frac_in(i.base, self.t)
                    A = is_deriv_k(basea, based, self.L_K, self.E_K,
                        self.L_args, self.E_args, self.D, self.T)
                    if A is None:
                        # Nonelementary monomial (so far)

                        # TODO: Would there ever be any benefit from just adding
                        # log(base) as a new monomial?
                        # ANSWER: Yes, otherwise we can't integrate x**x (or
                        # rather prove that it has no elementary integral)
                        # without first manually rewriting it as exp(x*log(x))
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
                    # http://groups.google.com/group/sympy/browse_thread/thread/a61d48235f16867f

                    self.newf = self.newf.subs(i, newterm)

                elif i not in numpows:
                    continue
                new = exp(i.exp*log(i.base))
                # TODO: Just put it in self.Tfuncs
                self.backsubs.append((new, old))
                self.newf = self.newf.subs(old, newterm)
                exps.append(newterm)

            logs = filter(lambda i: i.args[0].is_rational_function(*self.T) and
                i.args[0].has_any_symbols(*self.T), self.newf.atoms(log))

            if handle_first == 'exp' or not log_new_extension:
                exp_new_extension = self._exp_part(exps)
                if exp_new_extension is None:
                    # reset and restart
                    self.f = self.newf
                    self.reset(dummy=dummy)
                    exp_new_extension = True
                    continue
            if handle_first == 'log' or not exp_new_extension:
                log_new_extension = self._log_part(logs)

        self.fa, self.fd = frac_in(self.newf, self.t)
        self._auto_attrs()

        return

    def _auto_attrs(self):
        """
        Set attributes that are generated automatically.
        """
        if not self.T:
            # i.e., when using the extension flag and T isn't given
            self.T = [i.gen for i in self.D]
        self.cases = [get_case(d, t) for d, t in zip(self.D, self.T)]
        self.cases.reverse()
        self.level = -1
        self.t = self.T[self.level]
        self.d = self.D[self.level]

    def _exp_part(self, exps):
        """
        Try to build an exponential extension.

        Returns True if there was a new extension, False if there was no new
        extension but it was able to rewrite the given exponentials in terms
        of the existing extension, and None if the entire extension building
        process should be restarted.  If the process fails because there is no
        way around an algebraic extension (e.g., exp(log(x)/2)), it will raise
        NotImplementedError.
        """
        from sympy.integrals.prde import is_log_deriv_k_t_radical

        new_extension = False
        restart = False
        expargs = [i.exp for i in exps]
        ip = integer_powers(expargs)
        for arg, others in ip:
            # Minimize potential problems with algebraic substitution
            others.sort(key=lambda i: i[1])

            arga, argd = frac_in(arg, self.t)
            A = is_log_deriv_k_t_radical(arga, argd, self.L_K, self.E_K,
                self.L_args, self.E_args, self.D, self.T)

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
                    self.newf = self.newf.subs(exp(arg), exp(const)*Mul(*[
                        u**power for u, power in ans]))
                    self.newf = self.newf.subs([(exp(p*expargs[i]),
                        exp(const*p)*Mul(*[u**power for u, power in ans]))
                        for i, p in others])
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
                        self.newf = self.newf.subs([(exp(p*expargs[i]),
                            exp(const*p)*rad) for i, p in others])
                        self.newf = self.newf.subs(zip(reversed(self.T),
                            [f(self.x) for f in self.Tfuncs]))
                        restart = True
                        break
                    else:
                        # TODO: give algebraic dependence in error string
                        raise NotImplementedError("Cannot integrate over " +
                        "algebraic extensions.")
            else:
                arga, argd = frac_in(arg, self.t)
                darga = (argd*derivation(Poly(arga, self.t), self.D, self.T) -
                    arga*derivation(Poly(argd, self.t), self.D, self.T))
                dargd = argd**2
                darga, dargd = darga.cancel(dargd, include=True)
                darg = darga.as_basic()/dargd.as_basic()
                self.t = self.ts.next()
                self.T.append(self.t)
                self.E_args.append(arg)
                self.E_K.append(len(self.T) - 1)
                self.D.append(darg.as_poly(self.t, expand=False)*Poly(self.t,
                    self.t, expand=False))
                i = Symbol('i', dummy=True)
                self.Tfuncs = [Lambda(i, exp(arg.subs(self.x, i)))] + self.Tfuncs
                self.newf = self.newf.subs([(exp(expargs[i]), self.t**p) for i,
                    p in others])
                new_extension = True

        if restart:
            return None
        return new_extension

    def _log_part(self, logs):
        """
        Try to build a logarithmic extension.

        Returns True if there was a new extension and False if there was no new
        extension but it was able to rewrite the given logarithms in terms
        of the existing extension.  Unlike with exponential extensions, there
        is no way that a logarithm is not transcendental over and cannot be
        rewritten in terms of an already existing extension in a non-algebraic
        way, so this function does not ever return None or raise
        NotImplementedError.
        """
        from sympy.integrals.prde import is_deriv_k

        new_extension = False
        logargs = [i.args[0] for i in logs]
        for arg in logargs:
            # The log case is easier, because whenever a logarithm is algebraic
            # over the base field, it is of the form a1*t1 + ... an*tn + c,
            # which is a polynomial, so we can just replace it with that.
            # In other words, we don't have to worry about radicals.
            arga, argd = frac_in(arg, self.t)
            A = is_deriv_k(arga, argd, self.L_K, self.E_K, self.L_args,
                self.E_args, self.D, self.T)
            if A is not None:
                ans, u, const = A
                newterm = log(const) + u
                self.newf = self.newf.subs(log(arg), newterm)
                continue

            else:
                arga, argd = frac_in(arg, self.t)
                darga = (argd*derivation(Poly(arga, self.t), self.D, self.T) -
                    arga*derivation(Poly(argd, self.t), self.D, self.T))
                dargd = argd**2
                darg = darga.as_basic()/dargd.as_basic()
                self.t = self.ts.next()
                self.T.append(self.t)
                self.L_args.append(arg)
                self.L_K.append(len(self.T) - 1)
                self.D.append(cancel(darg.as_basic()/arg).as_poly(self.t,
                    expand=False))
                i = Symbol('i', dummy=True)
                self.Tfuncs = [Lambda(i, log(arg.subs(self.x, i)))] + self.Tfuncs
                self.newf = self.newf.subs(log(arg), self.t)
                new_extension = True

        return new_extension

    @property
    def _important_attrs(self):
        """
        Returns some of the more important attributes of self.

        Used for testing and debugging purposes.

        The attributes are (fa, fd, D, T, Tfuncts, backsubs, E_K, E_args,
        L_K, L_args).
        """
        return (self.fa, self.fd, self.D, self.T, self.Tfuncs,
            self.backsubs, self.E_K, self.E_args, self.L_K, self.L_args)

    def __str__(self):
        return str(self._important_attrs)

    def reset(self, dummy=True):
        """
        Reset self to an initial state.  Used by __init__.
        """
        self.t = self.x
        self.T = [self.x]
        self.D = [Poly(1, self.x)]
        self.level = -1
        self.L_K, self.E_K, self.L_args, self.E_args = [], [], [], []
        if dummy:
            self.ts = numbered_symbols('t', cls=Dummy)
        else:
            # For testing
            self.ts = numbered_symbols('t')
        # For various things that we change to make things work that we need to
        # change back when we are done.
        self.backsubs = []
        self.Tfuncs = []
        self.newf = self.f

    def increment_level(self):
        """
        Increase the level of self.

        This makes the working differential extension larger.  self.level is
        given relative to the end of the list (-1, -2, etc.), so we don't need
        do worry about it when building the extension.
        """
        if self.level >= -1:
            raise ValueError("The level of the differential extension cannot " +
                "be increased any further.")

        self.level += 1
        self.t = self.T[self.level]
        self.d = self.D[self.level]
        return None

    def decrement_level(self):
        """
        Decrease the level of self.

        This makes the working differential extension smaller.  self.level is
        given relative to the end of the list (-1, -2, etc.), so we don't need
        do worry about it when building the extension.
        """
        if self.level <= -len(self.T):
            raise ValueError("The level of the differential extension cannot " +
                "be increased any further.")

        self.level -= 1
        self.t = self.T[self.level]
        self.d = self.D[self.level]
        return None

class NonElementaryIntegralException(Exception):
    """
    Exception used by subroutines within the Risch algorithm to indicate to one
    another that the function being integrated does not have an elementary
    integral in the given differential field.
    """
    # TODO: Rewrite algorithms below to use this (?)
    pass

def gcdex_diophantine(a, b, c):
    """
    Extended Euclidean Algorithm, Diophantine version.

    Given a, b in K[x] and c in (a, b), the ideal generated by a and b,
    return (s, t) such that s*a + t*b == c and either s == 0 or s.degree()
    < b.degree().
    """
    # Extended Euclidean Algorithm (Diophantine Version) pg. 13
    # TODO: This should go in densetools.py.
    # XXX: Bettter name?

    s, g = a.half_gcdex(b)
    q = c.quo(g) # Inexact division means c is not in (a, b)
    s = q*s

    if not s.is_zero and b.degree() >= b.degree():
        q, s = s.div(b)

    t = (c - s*a).quo(b)

    return (s, t)

def frac_in(f, t, **kwargs):
    """
    Returns the tuple (fa, fd), where fa and fd are Polys in t.

    This is a common idiom in the Risch Algorithm functions, so we abstract
    it out here.  f should be a basic expresion, a Poly, or a tuple (fa, fd),
    where fa and fd are either basic expressions or Polys, and f == fa/fd.
    **kwargs are applied to Poly.
    """
    cancel = kwargs.pop('cancel', False)
    if type(f) is tuple:
        fa, fd = f
        f = fa.as_basic()/fd.as_basic()
    fa, fd = f.as_basic().as_numer_denom()
    fa, fd = fa.as_poly(t, **kwargs), fd.as_poly(t, **kwargs)
    if cancel:
        fa, fd = fa.cancel(fd, include=True)
    if fa is None or fd is None:
        raise ValueError("Could not turn %s into a fraction in %s." % (f, t))
    return (fa, fd)

def as_poly_1t(p, t, z):
    """
    (Hackish) way to convert an element p of K[t, 1/t] to K[t, z].

    In other words, z == 1/t will be a dummy variable that Poly can handle
    better.

    See issue 2032.

    Doctest
    =======
    >>> from sympy import Symbol, random_poly
    >>> from sympy.integrals.risch import as_poly_1t
    >>> from sympy.abc import x, z

    >>> p1 = random_poly(x, 10, -10, 10)
    >>> p2 = random_poly(x, 10, -10, 10)
    >>> p = p1 + p2.subs(x, 1/x)
    >>> as_poly_1t(p, x, z).as_basic().subs(z, 1/x) == p
    True
    """
    pa, pd = frac_in(p, t, cancel=True)
    if not pd.is_monomial:
        # XXX: Is there a better Poly exception that we could raise here
        # Either way, if you see this (from the Risch Algorithm) it indicates
        # a bug.
        raise PolynomialError("%s is not an element of K[%s, 1/%s]." % (p, t, t))
    d = pd.degree(t)
    one_t_part = pa.slice(0, d + 1) # requires polys11
    r = pd.degree() - pa.degree()
    t_part = pa - one_t_part
    t_part = t_part.to_field().quo(pd)
    # Compute the negative degree parts.  Also requires polys11.
    one_t_part = Poly.from_list(reversed(one_t_part.rep.rep), *one_t_part.gens,
        **{'domain':one_t_part.domain})
    if r > 0:
        one_t_part *= Poly(t**r, t)

    one_t_part = one_t_part.replace(t, z) # z will be 1/t
    if pd.nth(d):
        one_t_part *= Poly(1/pd.nth(d), z, expand=False)
    ans = t_part.as_poly(t, z, expand=False) + one_t_part.as_poly(t, z, expand=False)

    return ans


def derivation(p, D, T, coefficientD=False, basic=False):
    """
    Computes Dp.

    Given the derivation D with D = d/dx and p is a polynomial in t over
    K(x), return Dp.

    If coefficientD is True, it computes the derivation kD
    (kappaD), which is defined as kD(sum(ai*Xi**i, (i, 0, n))) ==
    sum(Dai*Xi**i, (i, 1, n)) (Definition 3.2.2, page 80).  X in this case is
    T[-1], so coefficientD computes the derivative just with respect to T[:-1],
    with T[-1] treated as a constant.

    If basic=True, the returns a Basic expression.  Elements of D can still be
    instances of Poly.
    """
    t = T[-1]

    if basic:
        r = 0
    else:
        r = Poly(0, t)

    if coefficientD:
        D = D[:-1]
        T = T[:-1]

    for d, v in zip(D, T):
        pv = p.as_poly(v)
        if pv is None or basic:
            pv = p.as_basic()

        if basic:
            r += d.as_basic()*pv.diff(v)
        else:
            r += (d*pv.diff(v)).as_poly(t)

    if basic:
        r = cancel(r)
    return r

def get_case(d, t):
    """
    Returns the type of the derivation d.

    Returns one of {'exp', 'tan', 'base', 'primitive', 'other_linear',
    'other_nonlinear'}.
    """
    if not d.has(t):
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

def splitfactor(p, D, T, coefficientD=False, z=None):
    """
    Splitting factorization.

    Given a derivation D on k[t] and p in k[t], return (p_n, p_s) in
    k[t] x k[t] such that p = p_n*p_s, p_s is special, and each square
    factor of p_n is normal.

    Page. 100
    """
    t = T[-1]
    kinv = [1/x for x in T[:-1]]
    if z:
        kinv.append(z)

    One = Poly(1, t, domain=p.get_domain())
    Dp = derivation(p, D, T, coefficientD)
    # XXX: Is this right?
    if p.is_zero:
        return (p, One)

    if not p.has_any_symbols(t):
        s = p.as_poly(*kinv).gcd(Dp.as_poly(*kinv)).as_poly(t)
        n = p.quo(s)
        return (n, s)

    if not Dp.is_zero:
        h = p.gcd(Dp).to_field()
        g = p.gcd(p.diff(t)).to_field()
        s = h.quo(g)

        if s.degree(t) == 0:
            return (p, One)

        q_split = splitfactor(p.quo(s), D, T, coefficientD)

        return (q_split[0], q_split[1]*s)
    else:
        return (p, One)

def splitfactor_sqf(p, D, T, coefficientD=False, z=None):
    """
    Splitting Square-free Factorization

    Given a derivation D on k[t] and p in k[t], returns (N1, ..., Nm)
    and (S1, ..., Sm) in k[t]^m such that p =
    (N1*N2**2*...*Nm**m)*(S1*S2**2*...*Sm**m) is a splitting
    factorization of p and the Ni and Si are square-free and coprime.
    """
    # TODO: This algorithm appears to be faster in every case
    # TODO: Verify this and splitfactor() for multiple extensions
    t = T[-1]
    kkinv = [1/x for x in T[:-1]] + T[:-1]
    if z:
        kkinv = [z]

    S = []
    N = []
    p_sqf = p.sqf_list_include()
    if p.is_zero:
        return (((p, 1),), ())

    for pi, i in p_sqf:
        Si = pi.as_poly(*kkinv).gcd(derivation(pi, D, T,
            coefficientD).as_poly(*kkinv)).as_poly(t)
        pi = Poly(pi, t)
        Si = Poly(Si, t)
        Ni = pi.quo(Si)
        if not Si.is_one:
            S.append((Si, i))
        if not Ni.is_one:
            N.append((Ni, i))

    return (tuple(N), tuple(S))

def canonical_representation(a, d, D, T):
    """
    Canonical Representation.

    Given a derivation D on k[t] and f = a/d in k(t), return (f_p, f_s,
    f_n) in k[t] x k(t) x k(t) such that f = f_p + f_s + f_n is the
    canonical representation of f (f_p is a polynomial, f_s is reduced
    (has a special denominator), and f_n is simple (has a normal
    denominator).
    """
    t = T[-1]

    # Make d monic
    l = Poly(1/d.LC(), t)
    a, d = a.mul(l), d.mul(l)

    q, r = a.div(d)
    dn, ds = splitfactor(d, D, T)

    b, c = gcdex_diophantine(dn.as_poly(t), ds.as_poly(t), r.as_poly(t))
    b, c = b.as_poly(t), c.as_poly(t)

    return (q, (b, ds), (c, dn))

def hermite_reduce(a, d, D, T):
    """
    Hermite Reduction - Quadratic version.

    Given a derivation D on k(t) and f = a/d in k(t), returns g, h, r in
    k(t) such that f = Dg + h + r, h is simple, and r is reduced.
    """
    t = T[-1]

    # TODO: Rewrite this using Mack's linear version
    # Make d monic
    l = Poly(1/d.LC(), t)
    a, d = a.mul(l), d.mul(l)

    fp, fs, fn = canonical_representation(a, d, D, T)

    a, d = fn
    l = Poly(1/d.LC(), t)
    a, d = a.mul(l), d.mul(l)

    d_sqf = d.sqf_list_include()
    ga = Poly(0, t)
    gd = Poly(1, t)

    for v, i in d_sqf:
        if i < 2:
            continue

        u = d.quo(v**i)
        for j in range(i - 1, 0, -1):
            udv = u*derivation(v, D, T)
            b, c = gcdex_diophantine(udv.as_poly(t), v.as_poly(t),
                a.mul(Poly(-S(1)/j, t)).as_poly(t))
            b, c = b.as_poly(t), c.as_poly(t)

            vj = v**j
            ga = ga*vj + b*gd
            gd = gd*vj
            a = c.mul(Poly(-j, t)) - u*derivation(b, D, T)

        d = u*v

    q, r = a.div(d)

    ga, gd = ga.cancel(gd, include=True)
    r, d = r.cancel(d, include=True)

    rra = q*fs[1] + fp*fs[1] + fs[0]
    rrd = fs[1]
    rra, rrd = rra.cancel(rrd, include=True)

    return ((ga, gd), (r, d), (rra, rrd))

def polynomial_reduce(p, D, T):
    """
    Polynomial Reduction.

    Given a derivation D on k(t) and p in k[t] where t is a nonlinear
    monomial over k, return q, r in k[t] such that p = Dq  + r, and
    deg(r) < deg_t(Dt).
    """
    t = T[-1]
    d = D[-1]
    q = Poly(0, t)
    while p.degree(t) >= d.degree(t):
        m = p.degree(t) - d.degree(t) + 1
        q0 = Poly(t**m, t).mul(Poly(p.as_poly(t).LC()/(m*d.as_poly(t).LC()), t))
        q += q0
        p = p - derivation(q0, D, T)

    return (q, p)

def residue_reduce(a, d, D, T, z=None, invert=True):
    """
    Lazard-Rioboo-Rothstein-Trager resultant reduction.

    Given a derivation D on k(t) and f in k(t) simple, return g
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
    """
    # TODO: Use log_to_atan() from rationaltools.py
    # If r = residue_reduce(...), then the logarithmic part is given by:
    # sum([RootSum(a[0].as_poly(z), lambda i: i*log(a[1].as_basic()).subs(z,
    # i)).subs(t, log(x)) for a in r[0]])

    t = T[-1]
    z = z or Symbol('z', dummy=True)
    a, d = a.cancel(d, include=True)
    kkinv = [1/x for x in T[:-1]] + T[:-1]

    if a.is_zero:
        return ([], True)
    p, a = a.div(d)

    pz = Poly(z, t)

    Dd = derivation(d, D, T)
    q = a - pz*Dd

    if Dd.degree(t) <= d.degree(t):
        r, R = d.resultant(q, includePRS=True)
    else:
        r, R = q.resultant(d, includePRS=True)

    R_map, H = {}, []
    for i in R:
        R_map[i.degree()] = i

    r = Poly(r, z)
    Np, Sp = splitfactor_sqf(r, D, T, coefficientD=True, z=z)

    for s, i in Sp:
        if i == d.degree(t):
            s = Poly(s, z).monic()
            H.append((s, d))
        else:
            h = R_map.get(i)
            if h is None:
                continue
            h_lc = Poly(h.as_poly(t).LC(), t, field=True)

            h_lc_sqf = h_lc.sqf_list_include(all=True)

            for a, j in h_lc_sqf:
                h = Poly(h, t, field=True).quo(Poly(gcd(a, s**j, *kkinv), t))

            s = Poly(s, z).monic()

            if invert:
                h_lc = Poly(h.as_poly(t).LC(), t, field=True)
                inv, coeffs = h_lc.as_poly(z, field=True).invert(s), [S(1)]

                for coeff in h.coeffs()[1:]:
                    L = reduced(inv*coeff, [s])[1]
                    coeffs.append(L.as_basic())

                h = Poly(dict(zip(h.monoms(), coeffs)), t)

            H.append((s, h))

    b = all([not cancel(i.as_basic()).has_any_symbols(t, z) for i, _ in Np])

    return (H, b)

def residue_reduce_to_basic(H, T, z, Tfuncs):
    """
    Converts the tuple returned by residue_reduce() into a Basic expression.
    """
    # TODO: check what Lambda does with RootOf
    x = T[0]
    i = Symbol('i', dummy=True)
    s = zip(reversed(T), [f(x) for f in Tfuncs])

    return sum((RootSum(a[0].as_poly(z), Lambda(i, i*log(a[1].as_basic()).subs(
        {z: i}).subs(s))) for a in H))

def residue_reduce_derivation(H, D, T, z):
    """
    Computes the derivation of an expression returned by residue_reduce().

    In general, this is a rational function in t, so this returns an
    as_basic() result.
    """
    # TODO: verify that this is correct for multiple extensions
    i = Symbol('i', dummy=True)
    return S(sum((RootSum(a[0].as_poly(z), Lambda(i, i*derivation(a[1], D,
        T).as_basic().subs(z, i)/a[1].as_basic().subs(z, i))) for a in H)))

def integrate_primitive_polynomial(p, D, T):
    """
    Integration of primitive polynomials.

    Given a primitive monomial t over k, and p in k[t], return q in k[t],
    r in k, and a bool b in {True, False} such that r = p - Dq is in k if b is
    True, or r = p - Dq does not have an elementary integral over k(t) if b is
    False.
    """
    from sympy.integrals.prde import limited_integrate
    t = T[-1]

    if not p.has_any_symbols(t):
        return (Poly(0, t), p, True)

    T1 = T[:-1] # We had better be integrating the lowest extension (x) with
    D1 = D[:-1] # ratint().
    t1 = T1[-1]
    a = p.LC()
    aa, ad = frac_in(a, t1)
    Dt = D[-1]
    Dta, Dtb = frac_in(Dt, t1)

    try:
        (ba, bd), c = limited_integrate(aa, ad, [(Dta, Dtb)], D1, T1)
        assert len(c) == 1
    except NonElementaryIntegralException:
        return (Poly(0, t), p, False)

    m = p.degree(t)
    q0 = c[0].as_poly(t)*Poly(t**(m + 1)/(m + 1), t) + \
        (ba.as_basic()/bd.as_basic()).as_poly(t)*Poly(t**m, t)

    # TODO: Rewrite this non-recursively
    # c.f. risch_integrate(log(x)**1001, x)
    q, r, b = integrate_primitive_polynomial(p - derivation(q0, D, T), D, T)
    return (q + q0, r, b)

def integrate_primitive(a, d, D, T, Tfuncs):
    """
    Integration of primitive functions.

    Given a primitive monomial t over k and f in k(t), return g elementary over
    k(t), i in k(t), and b in {True, False} such that i = f - Dg is in k if b
    is True or i = f - Dg does not have an elementary integral over k(t) if b
    is False.

    This function returns a Basic expression for the first argument.  If b is
    True, the second argument is Basic expression in k to recursively integrate.
    If b is False, the second argument is an unevaluated Integral, which has
    been proven to be nonelementary.
    """
    # XXX: a and d must be canceled, or this might return incorrect results
    t = T[-1]
    x = T[0]
    z = Symbol('z', dummy=True)
    s = zip(reversed(T), [f(x) for f in Tfuncs])

    g1, h, r = hermite_reduce(a, d, D, T)
    g2, b = residue_reduce(h[0], h[1], D, T, z=z)
    if not b:
        i = cancel(a.as_basic()/d.as_basic() - (g1[1]*derivation(g1[0], D, T) -
            g1[0]*derivation(g1[0], D, T)).as_basic()/(g1[1]**2).as_basic() -
            residue_reduce_derivation(g2, D, T, z))
        i = Integral(cancel(i).subs(s), x)
        return ((g1[0].as_basic()/g1[1].as_basic()).subs(s) +
            residue_reduce_to_basic(g2, T, z, Tfuncs), i, b)

    # h - Dg2 + r
    p = cancel(h[0].as_basic()/h[1].as_basic() - residue_reduce_derivation(g2,
        D, T, z) + r[0].as_basic()/r[1].as_basic())
    p = p.as_poly(t)

    q, i, b = integrate_primitive_polynomial(p, D, T)

    ret = ((g1[0].as_basic()/g1[1].as_basic() + q.as_basic()).subs(s) +
        residue_reduce_to_basic(g2, T, z, Tfuncs))
    if not b:
        # TODO: This does not do the right thing when b is False
        i = Integral(cancel(i.as_basic()).subs(s), x)
    else:
        i = cancel(i.as_basic())

    return (ret, i, b)

def integrate_hyperexponential_polynomial(p, D, T, z):
    """
    Integration of hyperexponential polynomials.

    Given a hyperexponential monomial t over k and p in k[t, 1/t], return q in
    k[t, 1/t] and a bool b in {True, False} such that p - Dq in k if b is True,
    or p - Dq does not have an elementary integral over k(t) if b is False.
    """
    from sympy.integrals.rde import rischDE

    t = T[-1]
    d = D[-1]

    D1 = D[:-1]
    T1 = T[:-1]
    t1 = T1[-1]

    dtt = d.quo(Poly(t, t))
    qa = Poly(0, t)
    qd = Poly(1, t)
    b = True
    for i in xrange(-p.degree(z), p.degree(t) + 1):
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
            a = p.as_poly(t, expand=False).nth(i)

        aa, ad = frac_in(a, t1, field=True)
        aa, ad = aa.cancel(ad, include=True)
        iDt = Poly(i, t)*dtt
        iDta, iDtd = frac_in(iDt, t1, field=True)
        try:
            va, vd = rischDE(iDta, iDtd, Poly(aa, t1), Poly(ad, t1), D1, T1)
            va, vd = frac_in((va, vd), t)
        except NonElementaryIntegralException:
            b = False
        else:
            qa = qa*vd + va*Poly(t**i)*qd
            qd *= vd

    return (qa, qd, b)

def integrate_hyperexponential(a, d, D, T, Tfuncs):
    """
    Integration of hyperexponential functions.

    Given a hyperexponential monomial t over k and f in k(t), return g
    elementary over k(t), i in k(t), and a bool b in {True, False} such that
    i = f - Dg is in k if b is True or i = f - Dg does not have an elementary
    integral over k(t) if b is False.

    This function returns a Basic expression for the first argument.  If b is
    True, the second argument is Basic expression in k to recursively integrate.
    If b is False, the second argument is an unevaluated Integral, which has
    been proven to be nonelementary.
    """
    # XXX: a and d must be canceled, or this might return incorrect results
    t = T[-1]
    x = T[0]
    z = Symbol('z', dummy=True)
    s = zip(reversed(T), [f(x) for f in Tfuncs])

    g1, h, r = hermite_reduce(a, d, D, T)
    g2, b = residue_reduce(h[0], h[1], D, T, z=z)
    if not b:
        i = cancel(a.as_basic()/d.as_basic() - (g1[1]*derivation(g1[0], D, T) -
            g1[0]*derivation(g1[0], D, T)).as_basic()/(g1[1]**2).as_basic() -
            residue_reduce_derivation(g2, D, T, z))
        i = Integral(cancel(i.subs(s)), x)
        return ((g1[0].as_basic()/g1[1].as_basic()).subs(s) +
            residue_reduce_to_basic(g2, T, z, Tfuncs), i, b)

    # p should be a polynomial in t and 1/t, because Sirr == k[t, 1/t]
    # h - Dg2 + r
    p = cancel(h[0].as_basic()/h[1].as_basic() - residue_reduce_derivation(g2,
        D, T, z) + r[0].as_basic()/r[1].as_basic())
    pp = as_poly_1t(p, t, z) #

    qa, qd, b = integrate_hyperexponential_polynomial(pp, D, T, z)

    i = pp.nth(0, 0)

    ret = ((g1[0].as_basic()/g1[1].as_basic() + qa.as_basic()/
        qd.as_basic()).subs(s) + residue_reduce_to_basic(g2, T, z, Tfuncs))

    if not b:
        i = p - (qd*derivation(qa, D, T) - qa*derivation(qd, D, T)).as_basic()/\
            (qd**2).as_basic()
        i = Integral(cancel(i).subs(s), x)

    return (ret, i, b)

def integrate_hypertangent_polynomial(p, D, T):
    """
    Integration of hypertangent polynomials.

    Given a differential field k such that sqrt(-1) is not in k, a
    hypertangent monomial t over k, and p in k[t], return q in k[t] and
    c in k such that p - Dq - c*D(t**2 + 1)/(t**1 + 1) is in k and p -
    Dq does not have an elementary integral over k(t) if Dc != 0.
    """
    # XXX: Make sure that sqrt(-1) is not in k.
    t = T[-1]
    q, r = polynomial_reduce(p, D, T)
    a = derivation(t, D, T).quo(Poly(t**2 + 1, t))
    c = Poly(r.nth(1)/(2*a.as_basic()), t)
    return (q, c)

def integrate_nonlinear_no_specials(a, d, D, T, Tfuncs):
    """
    Integration of nonlinear monomials with no specials.

    Given a nonlinear monomial t over k such that Sirr ({p in k[t] | p is
    special, monic, and irreducible}) is empty, and f in k(t), returns g
    elementary over k(t) and a Boolean b in {True, False} such that f - Dg is
    in k if b == True, or f - Dg does not have an elementary integral over k(t)
    if b == False.

    This function is applicable to all nonlinear extensions, but in the case
    where it returns b == False, it will only have proven that the integral of
    f - Dg is nonelementary if Sirr is empty.

    This function returns a Basic expression.
    """
    # TODO: Integral from k?
    # TODO: split out nonelementary integral
    # XXX: a and d must be canceled, or this might not return correct results
    t = T[-1]
    x = T[0]
    z = Symbol('z', dummy=True)
    s = zip(reversed(T), [f(x) for f in Tfuncs])

    g1, h, r = hermite_reduce(a, d, D, T)
    g2, b = residue_reduce(h[0], h[1], D, T, z=z)
    if not b:
        return ((g1[0].as_basic()/g1[1].as_basic()).subs(s) +
            residue_reduce_to_basic(g2, T, z, Tfuncs), b)

    # Because f has no specials, this should be a polynomial in t, or else
    # there is a bug.
    p = cancel(h[0].as_basic()/h[1].as_basic() - residue_reduce_derivation(g2,
        D, T, z).as_basic() + r[0].as_basic()/r[1].as_basic()).as_poly(t)
    q1, q2 = polynomial_reduce(p, D, T)

    if q2.has(t):
        b = False
    else:
        b = True

    ret = (cancel(g1[0].as_basic()/g1[1].as_basic() + q1.as_basic()).subs(s) +
        residue_reduce_to_basic(g2, T, z, Tfuncs))
    return (ret, b)

# TODO: Should this go in the regular namespace?
# If so, index should default to False, I think.
def integer_powers(exprs, index=True):
    """
    Rewrites a list of expressions as integer multiples of each other.

    For example, if you have [x, x/2, x**2 + 1, 2*x/3], then you can rewrite
    this as [(x/6) * 6, (x/6) * 3, (x**2 + 1) * 1, (x/6) * 4].  This is useful
    in the Risch integration algorithm, where we must write exp(x) + exp(x/2) as
    (exp(x/2))**2 + exp(x/2), but we cannot write it as exp(x) + sqrt(exp(x))
    (this is because only the transcendental case is implemented and we
    therefore cannot integrate algebraic extensions).  The integer multiples
    returned by this function for each term are the smallest possible (their
    content equals 1).

    Returns a list of tuples where the first element is the base term and the
    second element is a list of `(index, factor)` terms, where `index` is the
    index of the other terms that can be rewritten in terms of the base term,
    and `factor` is the (rational number) multiplicative factor that must
    multiply the base term to obtain the original item indexed by `index`.

    If index=False, then the original expression will appear, instead of its
    index in the original list.

    The easiest way to understand this is to look at an example:

    >>> from sympy.abc import x
    >>> from sympy.integrals.risch import integer_powers
    >>> integer_powers([x, x/2, x**2 + 1, 2*x/3], index=True)
    [(x/6, [(0, 6), (1, 3), (3, 4)]), (1 + x**2, [(2, 1)])]
    >>> integer_powers([x, x/2, x**2 + 1, 2*x/3], index=False)
    [(x/6, [(x, 6), (x/2, 3), (2*x/3, 4)]), (1 + x**2, [(1 + x**2, 1)])]

    We can see how this relates to the example at the e of the docstring.  It
    chose x/6 as the first base term.  Then, the element 0 (x) can be written
    as (x/2) * 2, so we get (0, 2), and so on.  Now only element 2 (x**2 + 1)
    remains, and there are no other terms that can be written as a rational
    multiple of that, so we get that it can be written as (x**2 + 1) * 1.

    The function only accepts rational number multiples because only those are
    useful for arguments of exponentials, but it could easily be extended to
    support any kind of coefficient.
    """
    # Here is the strategy:

    # First, go through each term and determine if it can be rewritten as a
    # rational multiple of any of the terms gathered so far.  Because we only
    # care about rational number coefficients to rational functions (that is all
    # Risch cares about), cancel(a/b).is_Rational is sufficient for this. If it
    # is a multiple, we add its multiple to the dictionary.

    terms = {}
    for i, term in enumerate(exprs):
        added = False
        for j in terms:
            a = cancel(term/j)
            if a.is_Rational:
                if index:
                    terms[j].append((i, a))
                else:
                    terms[j].append((term, a))
                added = True
                break

        if not added:
            # It wasn't in there
            if index:
                terms[term] = [(i, S(1))]
            else:
                terms[term] = [(term, S(1))]

    # After we have done this, we have all the like terms together, so we just
    # need to find a common denominator so that we can get the base term and
    # integer multiples such that each term can be written as an integer
    # multiple of the base term, and the content of the integers is 1.

    newterms = {}
    for term in terms:
        common_denom = reduce(ilcm, [i.as_numer_denom()[1] for _, i in terms[term]])
        newterm = term/common_denom
        newmults = [(i, j*common_denom) for i, j in terms[term]]
        newterms[newterm] = newmults

    return sorted(list(newterms.iteritems()))

def risch_integrate(f, x, extension=None, handle_first='log'):
    r"""
    Prototype function for the Risch Integration Algorithm.

    Only transcendental functions are supported.  Currently, only exponentials
    and logarithms are supported, but support for trigonometric functions is
    forthcoming.

    If this function returns an unevaluated Integral in the result, it means
    that it has proven that integral to be nonelementary.  Any errors will
    result in raising NotImplementedError.

    Examples
    ========
    >>> from sympy import risch_integrate, exp, log, pprint
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
    log(x + log(x))   log(-x + log(x))    |   1
    --------------- - ---------------- +  | ------ dx
           2                 2            | log(x)
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
    >>> pprint(risch_integrate(x**x*log(x), x))
      /
     |
     |  x
     | x *log(x) dx
     |
    /

    >>> pprint(risch_integrate(-1/(x*log(x)*log(log(x))**2), x))
         1
    -----------
    log(log(x))
    """
    if extension:
       raise NotImplementedError("Manual extensions are not supported yet.")

    DE = DifferentialExtension(f, x, handle_first=handle_first)
    fa, fd, D, T, Tfuncs, backsubs = DE.fa, DE.fd, DE.D, DE.T, DE.Tfuncs, DE.backsubs
    cases = DE.cases
    t = DE.t

    result = 0
    for case in cases:
        if not fa.has(t) and not fd.has(t) and not case == 'base':
            T = T[:-1]
            D = D[:-1]
            Tfuncs = Tfuncs[1:]
            t = T[-1]
            fa, fd = frac_in((fa, fd), t)
            continue
        fa, fd = fa.cancel(fd, include=True)
        if case == 'exp':
            ans, i, b = integrate_hyperexponential(fa, fd, D, T, Tfuncs)
        elif case == 'primitive':
            ans, i, b = integrate_primitive(fa, fd, D, T, Tfuncs)
        elif case == 'base':
            # XXX: We can't call ratint() directly here because it doesn't
            # handle polynomials correctly.
            ans = integrate(fa.as_basic()/fd.as_basic(), x)
            b = False
            i = 0
        else:
            raise NotImplementedError("Only exponential and logarithmic " +
            "extensions are currently supported.")

        result += ans
        if b:
            T = T[:-1]
            D = D[:-1]
            Tfuncs = Tfuncs[1:]
            t = T[-1]
            fa, fd = frac_in(i, t)
        else:
            result += i
            return result.subs(backsubs)
