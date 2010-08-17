from sympy.core import (Basic, Expr, S, C, Symbol, Wild, Add, sympify, diff,
                        oo, Tuple, Dummy)

from sympy.core.symbol import Dummy
from sympy.integrals.trigonometry import trigintegrate
from sympy.integrals.deltafunctions import deltaintegrate
from sympy.integrals.rationaltools import ratint
from sympy.integrals.risch import heurisch
from sympy.utilities import xthreaded, flatten, any, all
from sympy.polys import Poly, PolynomialError
from sympy.solvers import solve
from sympy.functions import Piecewise, sign
from sympy.geometry import Curve
from sympy.functions.elementary.piecewise import piecewise_fold
from sympy.series import limit

class Integral(Expr):
    """Represents unevaluated integral."""

    __slots__ = ['is_commutative']

    def __new__(cls, function, *symbols, **assumptions):
        # Any embedded piecewise functions need to be brought out to the
        # top level so that integration can go into piecewise mode at the
        # earliest possible moment.
        function = piecewise_fold(sympify(function))

        if function.is_Number:
            if function is S.NaN:
                return S.NaN

        if symbols:
            limits = []
            for V in symbols:
                if isinstance(V, Symbol):
                    limits.append(Tuple(V))
                    continue
                elif isinstance(V, (tuple, list, Tuple)):
                    V = sympify(flatten(V))
                    if V[0].is_Symbol:
                        newsymbol = V[0]
                        if len(V) == 3:
                            if V[1] is None and V[2] is not None:
                                nlim = [V[2]]
                            elif V[1] is not None and V[2] is None:
                                function = -function
                                nlim = [V[1]]
                            elif V[1] is None and V[2] is None:
                                nlim = []
                            else:
                                nlim = V[1:]
                            limits.append(Tuple(newsymbol, *nlim ))
                            continue
                        elif len(V) == 1 or (len(V) == 2 and V[1] is None):
                            limits.append(Tuple(newsymbol))
                            continue
                        elif len(V) == 2:
                            limits.append(Tuple(newsymbol, V[1]))
                            continue


                raise ValueError("Invalid integration variable or limits: %s" % str(symbols))
        else:
            # no symbols provided -- let's compute full anti-derivative
            syms = function.atoms(Symbol)
            if not syms:
                raise ValueError('An integration variable is required.')
            limits = [Tuple(symb) for symb in syms]

        obj = Expr.__new__(cls, **assumptions)
        arglist = [function]
        arglist.extend(limits)
        obj._args = tuple(arglist)
        obj.is_commutative = all(s.is_commutative for s in obj.free_symbols)

        return obj

    def __getnewargs__(self):
        return (self.function,) + tuple([tuple(xab) for xab in self.limits])

    @property
    def function(self):
        return self._args[0]

    @property
    def limits(self):
        return self._args[1:]

    @property
    def variables(self):
        """Return a list of the integration variables.

        >>> from sympy import Integral
        >>> from sympy.abc import x, i
        >>> Integral(x**i, (i, 1, 3)).variables
        [i]
        """
        return [l[0] for l in self.limits]

    @property
    def free_symbols(self):
        """
        This method returns the symbols that will exist when the
        integral is evaluated. This is useful if one is trying to
        determine whether an integral is dependent on a certain
        symbol or not.

        >>> from sympy import Integral
        >>> from sympy.abc import x, y
        >>> Integral(x, (x, y, 1)).free_symbols
        set([y])
        """
        # analyze the integral
        # >>> Integral(x*y,(x,1,2),(y,1,3)).args
        # (x*y, Tuple(x, 1, 2), Tuple(y, 1, 3))
        # >>> Integral(x, x, y).args
        # (x, Tuple(x), Tuple(y))
        integrand, limits = self.function, self.limits
        if integrand.is_zero:
            return set()
        isyms = integrand.free_symbols
        for xab in limits:
            if len(xab) == 1:
                isyms.add(xab[0])
                continue
            # take out the target symbol
            if xab[0] in isyms:
                isyms.remove(xab[0])
            if len(xab) == 3 and xab[1] == xab[2]:
                # if two limits are the same the integral is 0
                # and there are no symbols
                return set()
            # add in the new symbols
            for i in xab[1:]:
                isyms.update(i.free_symbols)
        return isyms

    @property
    def is_zero(self):
        """Since Integral doesn't autosimplify it it useful to see if
        it would simplify to zero or not in a trivial manner, i.e. when
        the function is 0 or two limits of a definite integral are the same.

        This is a very naive and quick test, not intended to check for special
        patterns like Integral(sin(m*x)*cos(n*x), (x, 0, 2*pi)) == 0.
        """
        if (self.function.is_zero or
            any(len(xab) == 3 and xab[1] == xab[2] for xab in self.limits)):
            return True
        if not self.free_symbols and self.function.is_number:
            # the integrand is a number and the limits are numerical
            return False

    @property
    def is_number(self):
        """
        Return True if the Integral will result in a number, else False.

        sympy considers anything that will result in a number to have
        is_number == True.

        >>> from sympy import log
        >>> log(2).is_number
        True

        Integrals are a special case since they contain symbols that can
        be replaced with numbers. Whether the integral can be done or not is
        another issue. But answering whether the final result is a number is
        not difficult.

        >>> from sympy import Integral
        >>> from sympy.abc import x, y
        >>> Integral(x).is_number
        False
        >>> Integral(x, y).is_number
        False
        >>> Integral(x, (y, 1, x)).is_number
        False
        >>> Integral(x, (y, 1, 2)).is_number
        False
        >>> Integral(x, (y, 1, 1)).is_number
        True
        >>> Integral(x, (x, 1, 2)).is_number
        True
        >>> Integral(x*y, (x, 1, 2), (y, 1, 3)).is_number
        True
        >>> Integral(1, x, (x, 1, 2)).is_number
        True
        """

        integrand, limits = self.function, self.limits
        isyms = integrand.atoms(Symbol)
        for xab in limits:
            if len(xab) == 1:
                isyms.add(xab[0])
                continue # it may be removed later
            elif len(xab) == 3 and xab[1] == xab[2]: # XXX naive equality test
                return True # integral collapsed
            if xab[0] in isyms:
                # take it out of the symbols since it will be replace
                # with whatever the limits of the integral are
                isyms.remove(xab[0])
            # add in the new symbols
            for i in xab[1:]:
                isyms.update(i.free_symbols)
        # if there are no surviving symbols then the result is a number
        return len(isyms) == 0

    def as_dummy(self):
        """
        Replace instances of the integration variables with their dummy
        counterparts to make clear what are dummy variables and what
        are real-world symbols in an Integral. The "integral at" limit
        that has a length of 1 will be explicated with its length-2
        equivalent.

        >>> from sympy import Integral
        >>> from sympy.abc import x, y
        >>> Integral(x).as_dummy()
        Integral(_x, (_x, x))
        >>> Integral(x, (x, x, y), (y, x, y)).as_dummy()
        Integral(_x, (_x, x, _y), (_y, x, y))

        If there were no dummies in the original expression, then the
        output of this function will show which symbols cannot be
        changed by subs(), those with an underscore prefix.

        """
        reps = {}
        f = self.function
        limits = list(self.limits)
        for i in xrange(-1, -len(limits) - 1, -1):
            xab = list(limits[i])
            if len(xab) == 1:
                xab = xab*2
            x = xab[0]
            xab[0] = x.as_dummy()
            for j in range(1, len(xab)):
                xab[j] = xab[j].subs(reps)
            reps[x] = xab[0]
            limits[i] = xab
        f = f.subs(reps)
        return Integral(f, *limits)

    def transform(self, x, mapping, inverse=False):
        """
        Replace the integration variable x in the integrand with the
        expression given by `mapping`, e.g. 2*x or 1/x. The integrand and
        endpoints are rescaled to preserve the value of the original
        integral.

        In effect, this performs a variable substitution (although the
        symbol remains unchanged; follow up with subs to obtain a
        new symbol.)

        With inverse=True, the inverse transformation is performed.

        The mapping must be uniquely invertible (e.g. a linear or linear
        fractional transformation).
        """
        if x not in self.variables:
            return self
        limits = self.limits
        function = self.function
        y = Dummy('y')
        inverse_mapping = solve(mapping.subs(x, y) - x, y)
        if len(inverse_mapping) != 1 or not inverse_mapping[0].has(x):
            raise ValueError("The mapping must be uniquely invertible")
        inverse_mapping = inverse_mapping[0]
        if inverse:
            mapping, inverse_mapping = inverse_mapping, mapping
        function = function.subs(x, mapping) * mapping.diff(x)

        def calc_limit(a, b):
            """replace x with a, using subs if possible, otherwise limit
            where sign of b is considered"""
            wok = inverse_mapping.subs(x, a)
            if not wok is S.NaN:
                return wok
            return limit(sign(b)*inverse_mapping, x, a)
        newlimits = []
        for xab in limits:
            sym = xab[0]
            if sym == x and len(xab) == 3:
                a, b = xab[1:]
                a, b = calc_limit(a, b), calc_limit(b, a)
                if a == b:
                    raise ValueError("The mapping must transform the "
                        "endpoints into separate points")
                if a > b:
                    a, b = b, a
                    function = -function
                newlimits.append((sym, a, b))
            else:
                newlimits.append(xab)
        return Integral(function, *newlimits)


    def doit(self, **hints):
        if not hints.get('integrals', True):
            return self

        deep = hints.get('deep', True)

        # check for the trivial case of equal upper and lower limits
        if self.is_zero:
            return S.Zero

        # now compute and check the function
        function = self.function
        if deep:
            function = function.doit(**hints)
        if function.is_zero:
            return S.Zero

        # There is no trivial answer, so continue
        for xab in self.limits:
            antideriv = self._eval_integral(function, xab[0])

            if antideriv is None:
                newargs = (function, self.__getnewargs__()[1])
                return self.new(*newargs)
            else:
                if len(xab) == 1:
                    function = antideriv
                else:
                    if len(xab) == 3:
                        x, a, b = xab
                    if len(xab) == 2:
                        x, b = xab
                        a = None

                    if deep:
                        if isinstance(a, Basic):
                            a = a.doit(**hints)
                        if isinstance(b, Basic):
                            b = b.doit(**hints)

                    if antideriv.is_Poly:
                        gens = list(antideriv.gens)
                        gens.remove(x)

                        antideriv = antideriv.as_basic()

                        function = antideriv._eval_interval(x, a, b)
                        function = Poly(function, *gens)
                    else:
                        function = antideriv._eval_interval(x, a, b)

        return function

    def _eval_expand_basic(self, deep=True, **hints):
        from sympy import flatten
        if not deep:
            return self
        else:
            return Integral(self.function.expand(deep=deep, **hints),\
            flatten(*self.limits))

    def _eval_derivative(self, sym):
        """Evaluate the derivative of the current Integral object by
        differentiating under the integral sign [1], using the Fundamental
        Theorem of Calculus [2] when possible.

        Whenever an Integral is encountered that is equivalent to zero or
        has an integrand that is independent of the variable of integration
        those integrals are performed. All others are returned as Integral
        instances which can be resolved with doit() (provided they are integrable).

        References:
           [1] http://en.wikipedia.org/wiki/Differentiation_under_the_integral_sign
           [2] http://en.wikipedia.org/wiki/Fundamental_theorem_of_calculus

        >>> from sympy import Integral
        >>> from sympy.abc import x, y
        >>> i = Integral(x + y, y, (y, 1, x))
        >>> i.diff(x)
        Integral(x + y, (y, x)) + Integral(Integral(1, (y, y)), (y, 1, x))
        >>> i.doit().diff(x) == i.diff(x).doit()
        True
        >>> i.diff(y)
        0

        The previous must be true since there is no y in the evaluated integral:
        >>> i.free_symbols
        set([x])
        >>> i.doit()
        -1/6 - x/2 + 2*x**3/3

        """

        # differentiate under the integral sign; we do not
        # check for regularity conditions (TODO), see issue 1116

        # get limits and the function
        f, limits = self.function, list(self.limits)

        # the order matters if variables of integration appear in the limits
        # so work our way in from the outside to the inside.
        limit = limits.pop(-1)
        if len(limit) == 3:
            x, a, b = limit
        elif len(limit) == 2:
            x, b = limit
            a = None
        else:
            a = b = None
            x = limit[0]

        if limits: # f is the argument to an integral
            f = Integral(f, *tuple(limits))

        # assemble the pieces
        rv = 0
        if b is not None:
            rv += f.subs(x, b)*diff(b, sym)
        if a is not None:
            rv -= f.subs(x, a)*diff(a, sym)
        if len(limit) == 1 and sym == x:
            # the dummy variable *is* also the real-world variable
            arg = f
            rv += arg
        else:
            # the dummy variable might match sym but it's
            # only a dummy and the actual variable is determined
            # by the limits, so mask off the variable of integration
            # while differentiating
            u = Dummy('u')
            arg = f.subs(x, u).diff(sym).subs(u, x)
            rv += Integral(arg, Tuple(x, a, b))
        return rv

    def _eval_integral(self, f, x):
        """Calculate the anti-derivative to the function f(x).

        This is a powerful function that should in theory be able to integrate
        everything that can be integrated. If you find something, that it
        doesn't, it is easy to implement it.

        (1) Simple heuristics (based on pattern matching and integral table):

         - most frequently used functions (e.g. polynomials)
         - functions non-integrable by any of the following algorithms (e.g.
           exp(-x**2))

        (2) Integration of rational functions:

         (a) using apart() - apart() is full partial fraction decomposition
         procedure based on Bronstein-Salvy algorithm. It gives formal
         decomposition with no polynomial factorization at all (so it's fast
         and gives the most general results). However it needs much better
         implementation of RootsOf class (if fact any implementation).
         (b) using Trager's algorithm - possibly faster than (a) but needs
         implementation :)

        (3) Whichever implementation of pmInt (Mateusz, Kirill's or a
        combination of both).

          - this way we can handle efficiently huge class of elementary and
            special functions

        (4) Recursive Risch algorithm as described in Bronstein's integration
        tutorial.

          - this way we can handle those integrable functions for which (3)
            fails

        (5) Powerful heuristics based mostly on user defined rules.

         - handle complicated, rarely used cases
        """

        # if it is a poly(x) then let the polynomial integrate itself (fast)
        #
        # It is important to make this check first, otherwise the other code
        # will return a sympy expression instead of a Polynomial.
        #
        # see Polynomial for details.
        if isinstance(f, Poly):
            return f.integrate(x)

        # Piecewise antiderivatives need to call special integrate.
        if f.func is Piecewise:
            return f._eval_integral(x)

        # let's cut it short if `f` does not depend on `x`
        if not f.has(x):
            return f*x

        # try to convert to poly(x) and then integrate if successful (fast)
        poly = f.as_poly(x)

        if poly is not None:
            return poly.integrate().as_basic()

        # since Integral(f=g1+g2+...) == Integral(g1) + Integral(g2) + ...
        # we are going to handle Add terms separately,
        # if `f` is not Add -- we only have one term
        parts = []
        args = Add.make_args(f)
        for g in args:
            coeff, g = g.as_independent(x)

            # g(x) = const
            if g is S.One:
                parts.append(coeff*x)
                continue

            #               c
            # g(x) = (a*x+b)
            if g.is_Pow and not g.exp.has(x):
                a = Wild('a', exclude=[x])
                b = Wild('b', exclude=[x])

                M = g.base.match(a*x + b)

                if M is not None:
                    if g.exp == -1:
                        h = C.log(g.base)
                    else:
                        h = g.base**(g.exp + 1) / (g.exp + 1)

                    parts.append(coeff * h / M[a])
                    continue

            #        poly(x)
            # g(x) = -------
            #        poly(x)
            if g.is_rational_function(x):
                parts.append(coeff * ratint(g, x))
                continue

            # g(x) = Mul(trig)
            h = trigintegrate(g, x)
            if h is not None:
                parts.append(coeff * h)
                continue

            # g(x) has at least a DiracDelta term
            h = deltaintegrate(g, x)
            if h is not None:
                parts.append(coeff * h)
                continue

            # fall back to the more general algorithm
            try:
                h = heurisch(g, x, hints=[])
            except PolynomialError:
                # XXX: this exception means there is a bug in the
                # implementation of heuristic Risch integration
                # algorithm.
                h = None

            # if we failed maybe it was because we had
            # a product that could have been expanded,
            # so let's try an expansion of the whole
            # thing before giving up; we don't try this
            # out the outset because there are things
            # that cannot be solved unless they are
            # NOT expanded e.g., x**x*(1+log(x)). There
            # should probably be a checker somewhere in this
            # routine to look for such cases and try to do
            # collection on the expressions if they are already
            # in an expanded form
            if not h and len(args) == 1:
                f = f.expand(mul=True, deep=False)
                if f.is_Add:
                    return self._eval_integral(f, x)


            if h is not None:
                parts.append(coeff * h)
            else:
                return None

        return Add(*parts)

    def _eval_lseries(self, x):
        for term in self.function.lseries(x):
            yield integrate(term, *self.limits)

    def _eval_nseries(self, x, n):
        terms, order = self.function.nseries(x, n=n).as_coeff_add(C.Order)
        return integrate(terms, *self.limits) + Add(*order)*x

    def _eval_subs(self, old, new):
        """
        Substitute old with new in the integrand and the limits, but don't
        change anything that is (or corresponds to) a variable of integration.

        The normal substitution semantics -- traversing all arguments looking
        for matching patterns -- should not be applied to the Integrals since
        changing the integration variables should also entail a change in the
        integration limits (which should be done with the transform method). So
        this method just makes changes in the integrand and the limits.

        Not all instances of a given variable are conceptually the same: the
        first argument of the limit tuple and any corresponding variable in
        the integrand are dummy variables while every other symbol is a symbol
        that will be unchanged when the integral is evaluated. For example, in

            Integral(x + a, (a, a, b))

        the dummy variables are shown below with angle-brackets around them and
        will not be changed by this function:

            Integral(x + <a>, (<a>, a, b))

        If you want to change the lower limit to 1 there is no reason to
        prohibit this since it is not conceptually related to the integration
        variable, <a>. Nor is there reason to disallow changing the b to 1.

        If a second limit were added, however, as in:

            Integral(x + a, (a, a, b), (b, 1, 2))

        the dummy variables become:

            Integral(x + <a>, (<a>, a, <b>), (<b>, a, b))

        Note that the `b` of the first limit is now a dummy variable since `b` is a
        dummy variable in the second limit.

        Summary: no variable of the integrand or limit can be the target of
        substitution if it appears as a variable of integration in a limit
        positioned to the right of it.

        >>> from sympy import Integral
        >>> from sympy.abc import a, b, c, x, y

        >>> i = Integral(a + x, (a, a, 3), (b, x, c))
        >>> list(i.free_symbols) # only these can be changed
        [x, a, c]
        >>> i.subs(a, c) # note that the variable of integration is unchanged
        Integral(a + x, (a, c, 3), (b, x, c))
        >>> i.subs(a + x, b) == i # there is no x + a, only x + <a>
        True
        >>> i.subs(x, y - c)
        Integral(a + y - c, (a, a, 3), (b, y - c, c))
        """
        if self == old:
            return new
        integrand, limits = self.function, self.limits
        old_atoms = old.free_symbols
        limits = list(limits)

        # make limits explicit if they are to be targeted by old:
        # Integral(x, x) -> Integral(x, (x, x)) if old = x
        if old.is_Symbol:
            for i, l in enumerate(limits):
                if len(l) == 1 and l[0] == old:
                    limits[i] = Tuple(l[0], l[0])

        dummies = set()
        for i in xrange(-1, -len(limits) - 1, -1):
            xab = limits[i]
            if not dummies.intersection(old_atoms):
                limits[i] = Tuple(xab[0],
                                  *[l.subs(old, new) for l in xab[1:]])
            dummies.add(xab[0])
        if not dummies.intersection(old_atoms):
            integrand = integrand.subs(old, new)
        return Integral(integrand, *limits)

    def as_sum(self, n, method="midpoint"):
        """
        Approximates the integral by a sum.

        method ... one of: left, right, midpoint

        This is basically just the rectangle method [1], the only difference is
        where the function value is taken in each interval.

        [1] http://en.wikipedia.org/wiki/Rectangle_method

        **method = midpoint**:

        Uses the n-order midpoint rule to evaluate the integral.

        Midpoint rule uses rectangles approximation for the given area (e.g.
        definite integral) of the function with heights equal to the point on
        the curve exactly in the middle of each interval (thus midpoint
        method). See [1] for more information.

        Examples:

            >>> from sympy import sqrt
            >>> from sympy.abc import x
            >>> from sympy.integrals import Integral
            >>> e = Integral(sqrt(x**3+1), (x, 2, 10))
            >>> e
            Integral((1 + x**3)**(1/2), (x, 2, 10))
            >>> e.as_sum(4, method="midpoint")
            2*730**(1/2) + 4*7**(1/2) + 4*86**(1/2) + 6*14**(1/2)
            >>> e.as_sum(4, method="midpoint").n()
            124.164447891310
            >>> e.n()
            124.616199194723

        **method=left**:

        Uses the n-order rectangle rule to evaluate the integral, at each
        interval the function value is taken at the left hand side of the
        interval.

        Examples:

            >>> from sympy import sqrt
            >>> from sympy.abc import x
            >>> e = Integral(sqrt(x**3+1), (x, 2, 10))
            >>> e
            Integral((1 + x**3)**(1/2), (x, 2, 10))
            >>> e.as_sum(4, method="left")
            6 + 2*65**(1/2) + 2*217**(1/2) + 6*57**(1/2)
            >>> e.as_sum(4, method="left").n()
            96.8853618335341
            >>> e.n()
            124.616199194723

        """

        limits = self.limits
        if len(limits) > 1:
            raise NotImplementedError("Multidimensional midpoint rule not implemented yet")
        else:
            limit = limits[0]
        if n <= 0:
            raise ValueError("n must be > 0")
        if n == oo:
            raise NotImplementedError("Infinite summation not yet implemented")
        sym, lower_limit, upper_limit = limit
        dx = (upper_limit - lower_limit)/n
        result = 0.
        for i in range(n):
            if method == "midpoint":
                xi = lower_limit + i*dx + dx/2
            elif method == "left":
                xi = lower_limit + i*dx
            elif method == "right":
                xi = lower_limit + i*dx + dx
            else:
                raise NotImplementedError("Unknown method %s" % method)
            result += self.function.subs(sym, xi)
        return result*dx


@xthreaded
def integrate(*args, **kwargs):
    """integrate(f, var, ...)

       Compute definite or indefinite integral of one or more variables
       using Risch-Norman algorithm and table lookup. This procedure is
       able to handle elementary algebraic and transcendental functions
       and also a huge class of special functions, including Airy,
       Bessel, Whittaker and Lambert.

       var can be:

       - a symbol                   -- indefinite integration
       - a tuple (symbol, a, b)     -- definite integration

       Several variables can be specified, in which case the result is multiple
       integration.

       Also, if no var is specified at all, then the full anti-derivative of f is
       returned. This is equivalent to integrating f over all its variables.

       Examples

       >>> from sympy import integrate, log
       >>> from sympy.abc import a, x, y

       >>> integrate(x*y, x)
       y*x**2/2

       >>> integrate(log(x), x)
       -x + x*log(x)

       >>> integrate(log(x), (x, 1, a))
       1 - a + a*log(a)

       >>> integrate(x)
       x**2/2

       >>> integrate(x*y)
       x**2*y**2/4

       See also the doctest of Integral._eval_integral(), which explains
       thoroughly the strategy that SymPy uses for integration.

    """
    integral = Integral(*args, **kwargs)

    if isinstance(integral, Integral):
        return integral.doit(deep = False)
    else:
        return integral


@xthreaded
def line_integrate(field, curve, vars):
    """line_integrate(field, Curve, variables)

       Compute the line integral.

       Examples
       --------
       >>> from sympy import Curve, line_integrate, E, ln
       >>> from sympy.abc import x, y, t
       >>> C = Curve([E**t + 1, E**t - 1], (t, 0, ln(2)))
       >>> line_integrate(x + y, C, [x, y])
       3*2**(1/2)

    """
    F = sympify(field)
    if not F:
        raise ValueError("Expecting function specifying field as first argument.")
    if not isinstance(curve, Curve):
        raise ValueError("Expecting Curve entity as second argument.")
    if not isinstance(vars, (list, tuple)):
        raise ValueError("Expecting list/tuple for variables.")
    if len(curve.functions) != len(vars):
        raise ValueError("Field variable size does not match curve dimension.")

    if curve.parameter in vars:
        raise ValueError("Curve parameter clashes with field parameters.")

    # Calculate derivatives for line parameter functions
    # F(r) -> F(r(t)) and finally F(r(t)*r'(t))
    Ft = F
    dldt = 0
    for i, var in enumerate(vars):
        _f = curve.functions[i]
        _dn = diff(_f, curve.parameter)
        # ...arc length
        dldt = dldt + (_dn * _dn)
        Ft = Ft.subs(var, _f)
    Ft = Ft * dldt**(S(1)/2)

    integral = Integral(Ft, curve.limits).doit(deep = False)
    return integral
