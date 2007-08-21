from sympy.core import *
from sympy.core.methods import ArithMeths, RelMeths


class IntegralError(Exception):
    pass

class Integral(Basic, ArithMeths, RelMeths):
    """
    Carries out integration of the given expression with respect to symbols.

    expr must define ._eval_integral(symbol) method that returns integration result or None.

    Integral(Integral(expr, x), y) -> Integral(expr, x, y)

    Defined integral is represented as
      Integral(expr, x==[a,b])
    """

    precedence = Basic.Apply_precedence

    def __new__(cls, expr, *symbols, **assumptions):
        expr = Basic.sympify(expr)
        evaluate = assumptions.get("evaluate", True)
        if not evaluate:
            r = super(Integral, cls).__new__(cls, **assumptions)
            r._args = (expr,)+symbols
            #this is so that diff() works correctly, obviously
            #it will only work for 1D integrals
            r.f = expr
            if isinstance(symbols[0], Basic.Equality):
                r.x = symbols[0].lhs
                r.a = symbols[0].rhs.start
                r.b = symbols[0].rhs.end
            else:
                r.x = symbols[0]
                r.a = r.b = None
            return r
        for s in symbols:
            if isinstance(s, Basic.Equality):
                expr = cls.doit(cls, expr, s.lhs, s.rhs.start, s.rhs.end)
            else:
                expr = cls.doit(cls, expr, s, None, None)
        return expr

    def tostr(self, level=0):
        r = 'Int%s'%repr(tuple(self))
        if self.precedence <= level:
            r = '(%s)' % (r)
        return r

    def diff(self, sym):
        if sym == self.x:
            raise IntegralError("Cannot differentiate the integration variable")
        if self.a is None or self.b is None:
            raise IntegralError("Cannot differentiate because no endpoints specified")
        return (self.b.diff(sym)*self.f.subs(self.x,self.b)-\
            self.a.diff(sym)*self.f.subs(self.x,self.a))

    @staticmethod
    def doit(cls, f, x, a, b):
        """Try to do the integral."""
        F = cls.primitive_function(f, x)
        if isinstance(a, type(None)):
            return F
        else:
            # TODO This works fine as long as no discontinuties in between
            #      the two endpoints. A proper evaluation would involve
            #      finding D = {x_1, x_2, ..., x_n} the [ordered] set of
            #      discontinuities in the interval (a,b) and then doing
            #      something of the form
            #         v = limit(F,x,x_1) - limit(F,x,a)
            #         v += limit(F,x,b) - limit(F,x,x_n)
            #         v += sum( limit(F,x,x_{i+1}) - limit(F,x,x_i), i=1..n-1)
            #         return v
            v1 = F.subs(x, b)
            v2 = F.subs(x, a)
            if isinstance(v1, Basic.NaN):
                from sympy.series.limits import limit
                v1 = limit(F, x, b)
            if isinstance(v2, Basic.NaN):
                from sympy.series.limits import limit
                v2 = limit(F, x, a)
            return (v1 - v2)

    @staticmethod
    def primitive_function(f, x):
        """Try to calculate a primitive function to "f(x)".

        Use heuristics and an integral table.
        """
        from sympy import exp, cos, sin, log, atan, atanh, sqrt, uppergamma

        if not f.has(x): return f*x
        if f == x: return x**2/2

        if isinstance(f, Mul):
            # Pull out coefficient
            coeff, terms = f.as_coeff_terms(x)
            if coeff != 1:
                return coeff*Integral.primitive_function(Mul(*terms), x)
        elif isinstance(f, Add):
            result = 0
            for term in f:
                result += Integral.primitive_function(term, x)
            return result
        elif isinstance(f, Pow):
            if isinstance(f.exp,Number):
                if x == f.base:
                    if f.exp==-1: return log(abs(x))
                    else: return x**(f.exp+1)/(f.exp+1)
                elif isinstance(f.base, Mul) and x in f.base[:]:
                    coeff = 1
                    for b in f.base:
                        if b != x:
                            if x in b:
                                coeff = None
                                break
                            else:
                                coeff *= b
                    if coeff is not None:
                        coeff = coeff ** f.exp
                        if f.exp == -1: return log(abs(x)) / coeff
                        else: return x**(f.exp+1)/(f.exp+1) * coeff

        a = Wild('a', exclude=[x])
        b = Wild('b', exclude=[x])
        c = Wild('c', exclude=[x])
        integral_table = (
            ( x**(c-1)/(a*x**c+b), log(abs(a*x**c+b)) / a / c ),
            ( x/(a*x+b), x/a - b/a**2 * log(abs(a*x+b)) ),
            ( x/(a*x+b)**2, b/(a**2 * (a*x+b)) + log(abs(a*x+b))/a**2 ),
            ( (a*x+b)**c, (a*x+b)**(c+1) / a / (c+1) ),
            ( 1/(x**2 + a), atan(x/sqrt(a)) / sqrt(a) ),
            #( 1/(x**2 - a), -atanh(x/sqrt(a)) / sqrt(a) ),
            ( sin(a*x), -1/a * cos(a*x) ),
            ( cos(a*x), 1/a * sin(a*x) ),
            ( log(a*x+b), (x+b/a)*log(a*x+b)-x ),
            ( x**a * exp(b*x), (-1)*x**(a+1)*(-b*x)**(-a-1)*uppergamma(a+1,-b*x) )
        )

        for k,v in integral_table:
            r = f.match(k)
            if r != None:
                return v.subs_dict(r)

        raise IntegralError("Don't know how to do this integral: " + str(f))


def integrate(f, *args, **kargs):
    """
    Usage
    =====
      Indefinite integrals
      --------------------
      integrate(f, x) -> Returns the indefinite integral S{int} f(x) dx

      integrate(f, x, y) -> Return the indefinite double integral
      S{int} S{int} f(x, y) dy dx

      integrate(f, x, y, z, ...) -> Return the indefinite multiple integral
      (arbitrary number of variables) S{int} S{int} ... S{int} f(x, y, z, ...) dx ... dy dz


      Definite Integrals
      ------------------
      integrate(f, (x, a, b)) -> Returns the definite integral with integration
      limits a, b

      integrate(f, (x, a, b), (y, c, d)) -> Returns the definite double integral

    Notes
    =====
      Currently only very simple integrals are computed.The general algorithm
      for calculating integrals is described U{here<http://sympy.googlecode.com/svn/trunk/doc/issac98.pdf>}
      Someone just needs to implement it. :)

      Has an optional parameter evaluate, which can have value True or False.
      If set to False, the integral will not be evaluated. Default is set to True.

    Further examples
    ================
      >>> from sympy import Symbol
      >>> x, y = Symbol('x'), Symbol('y')
      >>> integrate(y, y)
      (1/2)*y**2
      >>> integrate(y*x, y)
      (1/2)*x*y**2
      >>> integrate(y*x, x, y)
      (1/4)*x**2*y**2

      #python2.4 doctest is missing a really important directive SKIP
      #>>> integrate(2*x*y, (x,0,1), (y,-1,2)) # doctest: +SKIP
      #3/2
      #>>> integrate(x*y**2 , (x,1,2), y) # doctest: +SKIP
      #(1/2)*y**3
      #>>> integrate(x , (x,1,2), evaluate=False) # doctest: +SKIP
      #integrate(x, (x, 1, 2))

    See also
    ========
      - L{limit<sympy.limits.limit>}

      - External links
        - U{Riemman integral<http://planetmath.org/encyclopedia/RiemannIntegral.html>}
    """
    return Integral(f, *args, **kargs)
