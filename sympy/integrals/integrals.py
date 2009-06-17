
from sympy.core import Basic, S, C, Symbol, Wild, Pow, sympify, diff

from sympy.integrals.trigonometry import trigintegrate
from sympy.integrals.deltafunctions import deltaintegrate
from sympy.integrals.rationaltools import ratint
from sympy.integrals.risch import heurisch
from sympy.utilities import threaded
from sympy.simplify import apart
from sympy.series import limit
from sympy.polys import Poly
from sympy.solvers import solve
from sympy.functions import DiracDelta, Heaviside, Piecewise
from sympy.geometry import Curve

class Integral(Basic):
    """Represents unevaluated integral."""

    def __new__(cls, function, *symbols, **assumptions):
        function = sympify(function)

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
                    limits.append((V,None))
                    continue
                elif isinstance(V, (tuple, list)):
                    if len(V) == 3:
                        if isinstance(V[0], Symbol):
                            nlim = map(sympify, V[1:])
                            limits.append( (V[0], tuple(nlim) ))
                            continue
                    elif len(V) == 1 or (len(V) == 2 and V[1] is None):
                        if isinstance(V[0], Symbol):
                            limits.append((V[0],None))
                            continue

                raise ValueError("Invalid integration variable or limits: %s" % str(symbols))
        else:
            # no symbols provided -- let's compute full antiderivative
            limits = [(symb,None) for symb in function.atoms(Symbol)]

            if not limits:
                return function

        obj = Basic.__new__(cls, **assumptions)
        obj._args = (function, tuple(limits))

        return obj

    def __getnewargs__(self):
        function, limits = self.args
        return (function,) + limits

    @property
    def function(self):
        return self._args[0]

    @property
    def limits(self):
        return self._args[1]

    @property
    def variables(self):
        variables = []

        for x,ab in self.limits:
            variables.append(x)

        return variables

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
        newlimits = []
        for sym, limit in limits:
            if sym == x and limit and len(limit) == 2:
                a, b = limit
                a = inverse_mapping.subs(x, a)
                b = inverse_mapping.subs(x, b)
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

        function = self.function

        for x,ab in self.limits:
            antideriv = self._eval_integral(function, x)

            if antideriv is None:
                return self
            else:
                if ab is None:
                    function = antideriv
                else:
                    a,b = ab
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

        if not sym in self.atoms(Symbol):
            return S.Zero

        if (sym, None) in self.limits:
            #case undefinite integral
            if len(self.limits) == 1:
                return self.function
            else:
                _limits = list(self.limits)
                _limits.pop( _limits.index((sym, None)) )
                return Integral(self.function, *tuple(_limits))

        #diff under the integral sign
        #we do not check for regularity conditions (TODO), see issue 1116
        if len(self.limits) > 1:
            # TODO:implement the multidimensional case
            raise NotImplementedError
        int_var = self.limits[0][0]
        lower_limit, upper_limit = self.limits[0][1]
        if sym == int_var:
            sym = Symbol(str(int_var), dummy=True)
        return self.function.subs(int_var, upper_limit)*diff(upper_limit, sym) - \
               self.function.subs(int_var, lower_limit)*diff(lower_limit, sym) + \
               integrate(diff(self.function, sym), (int_var, lower_limit, upper_limit))

    def _eval_integral(self, f, x):
        """Calculate the antiderivative to the function f(x).

        This is a powerful function that should in theory be able to integrate
        everything that can be integrated. If you find something, that it
        doesn't, it is easy to implement it.

        (1) Simple heuristics (based on pattern matching and integral table):

         - most frequently used functions (eg. polynomials)
         - functions non-integrable by any of the following algorithms (eg.
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
        if isinstance(f,Piecewise):
            return f._eval_integral(x)

        # let's cut it short if `f` does not depend on `x`
        if not f.has(x):
            return f*x

        # try to convert to poly(x) and then integrate if successful (fast)
        poly = f.as_poly(x)

        if poly is not None:
            return poly.integrate(x).as_basic()

        # since Integral(f=g1+g2+...) == Integral(g1) + Integral(g2) + ...
        # we are going to handle Add terms separately,
        # if `f` is not Add -- we only have one term
        if not f.is_Add:
            f = [f]

        parts = []

        if isinstance(f, Basic):
            f = f.args
        for g in f:
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

            if h is not None:
                parts.append(coeff * h)
            else:
                return None

        return C.Add(*parts)

    def _eval_lseries(self, x, x0):
        arg = self.args[0]
        dx, bounds = self.args[1][0]
        for term in arg.lseries(dx, x0):
            if bounds:
                a, b = bounds
                yield integrate(term, (dx, a, b))
            else:
                yield integrate(term, x)

    def _eval_nseries(self, x, x0, n):
        arg = self.args[0]
        x, bounds = self.args[1][0]
        arg = arg.nseries(x, x0, n)
        if bounds:
            a, b = bounds
            return integrate(arg.removeO(), (x, a, b)) + arg.getO()*x
        else:
            return integrate(arg.removeO(), x) + arg.getO()*x

    def _eval_subs(self, old, new):
        arg0 = self.args[0].subs(old, new)
        arg1 = []
        for sym, limits in self.args[1]:
            if sym == old:
                return self
            if limits is not None:
                a, b, = limits
                arg1.append((sym, a.subs(old, new), b.subs(old, new)))
            else:
                arg1.append((sym, limits))
        return Integral(arg0, *arg1)

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

       Also, if no var is specified at all, then full-antiderivative of f is
       returned. This is equivalent of integrating f over all it's variables.

       Examples

       >>> from sympy import *
       >>> x, y = symbols('xy')

       >>> integrate(x*y, x)
       (1/2)*y*x**2

       >>> integrate(log(x), x)
       -x + x*log(x)

       >>> integrate(x)
       (1/2)*x**2

       >>> integrate(x*y)
       (1/4)*x**2*y**2

       See also the doctest of Integral._eval_integral(), which explains
       thoroughly the strategy that SymPy uses for integration.

    """
    new_args = [sympify(arg) for arg in args]
    integral = Integral(*new_args, **kwargs)

    if isinstance(integral, Integral):
        return integral.doit()
    else:
        return integral


@threaded(use_add=False)
def line_integrate(field, curve, vars):
    """line_integrate(field, Curve, variables)

       Compute the line integral.

       Examples
       --------
       >>> from sympy import *
       >>> x, y, t = symbols('xyt')
       >>> C = Curve([E**t + 1, E**t - 1], (t, 0, ln(2)))
       >>> line_integrate(x + y, C, [x, y])
       3*sqrt(2)

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

    integral = Integral(Ft, curve.limits).doit()
    return integral

