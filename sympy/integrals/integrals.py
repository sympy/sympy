from sympy.core import Basic, Expr, S, C, Symbol, Wild, Add, sympify, diff, oo, Tuple

from sympy.integrals.trigonometry import trigintegrate
from sympy.integrals.deltafunctions import deltaintegrate
from sympy.integrals.rationaltools import ratint
from sympy.integrals.risch import heurisch
from sympy.utilities import threaded, flatten
from sympy.polys import Poly
from sympy.solvers import solve
from sympy.functions import Piecewise, sign
from sympy.geometry import Curve
from sympy.functions.elementary.piecewise import piecewise_fold
from sympy.series import limit

class Integral(Expr):
    """Represents unevaluated integral."""

    def __new__(cls, function, *symbols, **assumptions):
        # Any embedded piecewise functions need to be brought out to the
        # top level so that integration can go into piecewise mode at the
        # earliest possible moment.
        function = piecewise_fold(sympify(function))

        if function.is_Number:
            if function is S.NaN:
                return S.NaN
            elif function is S.Infinity:
                return S.Infinity
            elif function is S.NegativeInfinity:
                return S.NegativeInfinity

        if symbols:
            limits = []

            for V in symbols:
                if isinstance(V, Symbol):
                    limits.append(Tuple(V))
                    continue
                elif isinstance(V, (tuple, list, Tuple)):
                    V = flatten(V)
                    newsymbol = sympify(V[0])
                    if len(V) == 3:
                        if isinstance(newsymbol, Symbol):
                            nlim = map(sympify, V[1:])
                            if V[1] is None and V[2] is not None:
                                nlim = [V[2]]
                            if V[2] is None and V[1] is not None:
                                function = -function
                                nlim = [V[1]]
                            if V[1] is None and V[2] is None:
                                nlim = []
                            limits.append( Tuple(newsymbol, *nlim ))
                            continue
                    elif len(V) == 1 or (len(V) == 2 and V[1] is None):
                        if isinstance(newsymbol, Symbol):
                            limits.append(Tuple(newsymbol))
                            continue
                    elif len(V) == 2:
                        if isinstance(newsymbol, Symbol):
                            limits.append(Tuple(newsymbol,V[1]))
                            continue


                raise ValueError("Invalid integration variable or limits: %s" % str(symbols))
        else:
            # no symbols provided -- let's compute full anti-derivative
            limits = [Tuple(symb) for symb in function.atoms(Symbol)]

            if not limits:
                return function

        obj = Expr.__new__(cls, **assumptions)
        arglist = [function]
        arglist.extend(limits)
        obj._args = tuple(arglist)

        return obj

    def __getnewargs__(self):
        function = self.args[0]
        limits = self.args[1:]
        newlimits = []
        for lim in limits:
            if len(lim) == 1:
                newlimits.append((lim[0]))
            elif len(lim) == 2:
                newlimits.append((lim[0], lim[1]))
            else:
                newlimits.append((lim[0], lim[1], lim[2]))
        return (function,) + tuple(newlimits)

    @property
    def function(self):
        return self._args[0]

    @property
    def limits(self):
        return self._args[1:]

    @property
    def variables(self):
        variables = []

        for xab in self.limits:
            variables.append(xab[0])

        return variables

    @property
    def symbols(self):
        """
        This method returns the symbols that will exist when the
        integral is evaluated. This is useful if one is trying to
        determine whether an integral is dependent on a certain
        symbol or not.

        >>> from sympy import Integral
        >>> from sympy.abc import x, y
        >>> Integral(x, (x, y, 1)).symbols
        set([y])
        """
        # analyze the integral
        # >>> Integral(x*y,(x,1,2),(y,1,3)).args
        # (x*y, Tuple(x, 1, 2), Tuple(y, 1, 3))
        # >>> Integral(x, x, y).args
        # (x, Tuple(x), Tuple(y))
        intgrl = self
        args = intgrl.args
        integrand, limits = args[0], args[1:]
        if integrand.is_zero:
            return set()
        isyms = integrand.atoms(Symbol)
        for ilim in limits:
            if len(ilim) == 1:
                isyms.add(ilim[0])
                continue
            # take out the target symbol
            if ilim[0] in isyms:
                isyms.remove(ilim[0])
            if len(ilim) == 3 and ilim[1] == ilim[2]:
                # if two limits are the same the integral is 0
                # and there are no symbols
                return set()
            # add in the new symbols
            for i in ilim[1:]:
                isyms.update(i.atoms(Symbol))
        return isyms

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
        y = Symbol('y', dummy=True)
        inverse_mapping = solve(mapping.subs(x,y)-x, y)
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
        for lim in limits:
            sym = lim[0]
            if sym == x and len(lim) == 3:
                a, b = lim[1:]
                a, b = calc_limit(a, b), calc_limit(b, a)
                if a == b:
                    raise ValueError("The mapping must transform the "
                        "endpoints into separate points")
                if a > b:
                    a, b = b, a
                    function = -function
                newlimits.append((sym, a, b))
            else:
                newlimits.append((sym, limit))
        return Integral(function, *newlimits)


    def doit(self, **hints):
        if not hints.get('integrals', True):
            return self

        deep = hints.get('deep', True)

        function = self.function
        if deep:
            function = function.doit(**hints)

        for lim in self.limits:
            x = lim[0]
            antideriv = self._eval_integral(function, x)

            if antideriv is None:
                newargs = (function, self.__getnewargs__()[1])
                return self.new(*newargs)
            else:
                if len(lim) == 1:
                    function = antideriv
                else:
                    if len(lim) == 3:
                        a = lim[1]
                        b = lim[2]
                    if len(lim) == 2:
                        a = None
                        b = lim[1]

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
        """Evaluate the derivative of the current Integral object.
        We follow these steps:

        (1) If sym is not part of the function nor the integration limits,
            return 0

        (2) Check for a possible application of the Fundamental Theorem of
            Calculus [1]

        (3) Derive under the integral sign [2]

        References:
           [1] http://en.wikipedia.org/wiki/Fundamental_theorem_of_calculus
           [2] http://en.wikipedia.org/wiki/Differentiation_under_the_integral_sign
        """

        if not self.has(sym):
            return S.Zero

        if Tuple(sym) in self.limits:
            #case undefinite integral
            if len(self.limits) == 1:
                return self.function
            else:
                _limits = list(self.limits)
                _limits.pop( _limits.index(Tuple(sym)) )
                return Integral(self.function, *tuple(_limits))

        #diff under the integral sign
        #we do not check for regularity conditions (TODO), see issue 1116
        if len(self.limits) > 1:
            # TODO:implement the multidimensional case
            raise NotImplementedError
        int_var = self.limits[0][0]
        lower_limit, upper_limit = self.limits[0][1],self.limits[0][2]
        if sym == int_var:
            sym = Symbol(str(int_var), dummy=True)
        return self.function.subs(int_var, upper_limit)*diff(upper_limit, sym) - \
               self.function.subs(int_var, lower_limit)*diff(lower_limit, sym) + \
               integrate(diff(self.function, sym), (int_var, lower_limit, upper_limit))

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
                parts.append(coeff * x)
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
                        h = g.base**(g.exp+1) / (g.exp+1)

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
            h = deltaintegrate(g,x)
            if h is not None:
                parts.append(coeff * h)
                continue

            # fall back to the more general algorithm
            h = heurisch(g, x, hints=[])

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
        terms, order = self.function.nseries(x, n=n
                                             ).as_coeff_factors(C.Order)
        return integrate(terms, *self.limits) + Add(*order)*x

    def _eval_subs(self, old, new):
        if self == old:
            return new
        arg0 = self.args[0].subs(old, new)
        arg1 = []
        for lim in self.args[1:]:
            sym = lim[0]
            if sym == old:
                return self
            if len(lim) == 1:
                arg1.append((sym,))
            elif len(lim) == 2:
                b = lim[1]
                arg1.append((sym, None, b.subs(old, new)))
            else:
                a, b, = lim[1:3]
                arg1.append((sym, a.subs(old, new), b.subs(old, new)))
        return Integral(arg0, *arg1)

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

        if len(self.args[1:]) > 1:
            raise NotImplementedError("Multidimensional midpoint rule not implemented yet")
        if n <= 0:
            raise ValueError("n must be > 0")
        if n == oo:
            raise NotImplementedError("Infinite summation not yet implemented")
        sym, lower_limit,upper_limit = self.args[1]
        dx = (upper_limit-lower_limit)/n
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
            result += self.args[0].subs(sym, xi)
        return result*dx



@threaded(use_add=False)
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


@threaded(use_add=False)
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

