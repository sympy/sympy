from sympy.core import Basic, Symbol, Number, Mul, Pow, Add
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

    def __new__(cls, expr, *symbols,  **assumptions):
        expr = Basic.sympify(expr)
        if assumptions.has_key("evaluate"):
            evaluate = assumptions["evaluate"]
            del assumptions["evaluate"]
        else:
            evaluate = True
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
            return (F.subs(x, b)-F.subs(x, a))

    @staticmethod
    def primitive_function(f, x):
        """Try to calculate a primitive function to "f(x)".

        Use heuristics and an integral table.
        """
        from sympy import exp, cos, sin, log
        from sympy.modules.specfun.factorials import upper_gamma

        if isinstance(f,Mul):
            #a,b = f.getab()
            a,b = f[0],Mul(*f[1:])
            if not a.has(x): return a*Integral.primitive_function(b, x)
            if not b.has(x): return b*Integral.primitive_function(a, x)
        if isinstance(f,Add):
            #a,b = f.getab()
            a,b = f[0],Add(*f[1:])
            return Integral.primitive_function(a,x)+Integral.primitive_function(b,x)
        if not f.has(x): return f*x
        if f==x: return x**2/2
        if isinstance(f,Pow):
            if isinstance(f.exp,Number):
                if x == f.base:
                    if f.exp==-1: return log(abs(x))
                    else: return x**(f.exp+1)/(f.exp+1)
                elif x in f.base and isinstance(f.base, Mul):
                    other = 1
                    for b in f.base:
                        if b != x: other *= b
                    other = other ** f.exp

                    if f.exp==-1: return log(abs(x)) * other
                    else: return x**(f.exp+1)/(f.exp+1) * other

        a,b,c = [Symbol(s, dummy = True) for s in ["a","b","c"]]
        integral_table = {
                a/(b*x+c): a/b * log(abs(b*x+c)),
                a*sin(b*x): -a/b * cos(b*x),
                a*cos(b*x): a/b * sin(b*x),
                log(a*x): x*log(a*x)-x,
                # Note: the next two entries are special cases of the
                # third and would be redundant with a more powerful match()
                exp(a*x) : exp(a*x)/a,
                #x * exp(a*x) : exp(a*x) * (a*x-1) / a**2,
                x**a * exp(b*x) : (-1)*x**(a+1)*(-b*x)**(-a-1)*upper_gamma(a+1,-b*x)
                }
        for k in integral_table:
            r = f.match(k, [a,b,c])
            if r != None:
                # Prevent matching nonconstant expressions 
                if [1 for v in r.values() if v.has(x)]:
                    continue
                return integral_table[k].subs_dict(r)

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
      - L{limit<sympy.modules.limits.limit>}

      - External links
        - U{Riemman integral<http://planetmath.org/encyclopedia/RiemannIntegral.html>}
    """
    return Integral(f, *args, **kargs)
