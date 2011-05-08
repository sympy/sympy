from sympy.polys.lpoly import (lgens, LPoly,monomial_tobasic)
from sympy.series.order import O
from sympy.core.singleton import S
from sympy.polys.domains import QQ
from sympy.polys.monomialtools import monomial_lex_key as O_lex
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import (cos,sin,tan,asin,atan,acos,acot)
from sympy.functions.elementary.exponential import (exp,log,LambertW)
from sympy.functions.elementary.hyperbolic import (sinh,cosh,tanh,atanh,asinh,acosh,acoth)
from sympy.core.numbers import (Number,Rational)
from sympy.core import pi
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy import I
from sympy.core.sympify import sympify

# cot(x), coth(x) are Laurent expansions, so they are passed to series

class TaylorEvalError(TypeError):
    pass


def ev_args(te,a):
    if len(a) == 1:
        a = a[0]
        if a == te.var:
            return te.lvar
        if isinstance(a, Number):
            return te.coerce_number(a)
        return te(a)
    else:
        raise NotImplementedError

def taylor(p,var=None,start=0,prec=6,dir="+",pol_pars=[],ov=True,level=0):
    """
    taylor series expansion of p
    kept the same arguments as series, with the addition
    of pol_pars
    var series variable
    start  var=start point of expansion
    prec precision of the series
    dir ...
    pol_pars polynomial parameters
    ov = True return always the series in expanded form
    ov = return a series which must be expanded to be put in canonical
         form; this is faster

    ALGORITHM try first to compute the series in the
    QQ ring; if this fails compute it
    in the symbolic ring SR consisting
    of the sympy expressions which do not depend on
    var and pol_pars; if also this fails compute it
    using the series function

    EXAMPLES

    >>> from sympy.core.symbol import symbols
    >>> from sympy.functions.elementary.trigonometric import (sin,tan,atan)
    >>> from sympy.functions.elementary.exponential import (exp,log)
    >>> from sympy.functions.elementary.miscellaneous import sqrt
    >>> from sympy.polys.ltaylor import taylor
    >>> from sympy import pi
    >>> x,y = symbols('x,y')
    >>> taylor(sin(x*tan(x)),x,0,10)
    x**2 + x**4/3 - x**6/30 - 71*x**8/630 + O(x**10)

    >>> taylor(sqrt(1 + x*sin(pi*x)),x,0,6)
    1 + x**4*(-pi**2/8 - pi**3/12) + pi*x**2/2 + O(x**6)

    >>> taylor(exp(x*log(x)),x,0,3)
    1 + x*log(x) + x**2*log(x)**2/2 + O(x**3*log(x)**3)

    In these examples y is first treated internally
    as a Sympy symbol, then as a polynomial parameter;
    the latter version is faster

    >>> taylor(atan(x*y + x**2),x,0,5)
    x*y + x**2 - x**4*y**2 - x**3*y**3/3 + O(x**5)
    >>> taylor(atan(x*y + x**2),x,0,5,pol_pars=[y])
    x*y + x**2 - x**4*y**2 - x**3*y**3/3 + O(x**5)
    >>> taylor(sum(sin(sin(n*x)) for n in range(1,4)),x,0,10)
    6*x - 12*x**3 + 138*x**5/5 - 6176*x**7/105 + 7293*x**9/70 + O(x**10)
    """

    p = sympify(p)
    # TODO deal with some of these cases within taylor
    if var == None or prec == None or dir != "+" or \
        prec in [S.Infinity, S.NegativeInfinity]:
        return p.series(var,start,prec)

    if var not in p.atoms():
        return p

    gens = [var] + pol_pars
    # taylor(p1+p2,...) = taylor(p1,...) + taylor(p2,...)
    # in the case in which p2=O(x**prec1), if prec1 < prec
    # series gives a value error
    # with p =  1/x + sin(sin(x)) the first term is passed to series
    # the second is computed in the QQ ring
    if p.__class__ == Add:
        p2 = 0
        p2o = 0
        # consider first the Order; in case the order is less than
        # prec series raises ValueError, so it is not necessary
        # to do the rest of the computation
        for q in p.args:
            if q.is_Order:
                p2o += q.series(var,start,prec,dir)
        for q in p.args:
            if q.is_Order:
                continue
            q = sympify(q)
            p2 += taylor(q,var,start,prec,dir,pol_pars,ov=0,level=level)
        if not p2o and ov:
            return p2 + O(var**prec)
        return p2 + p2o
    # case p =  log(2)*sin(sin(x))
    if p.__class__ == Mul and level == 0:
        c = 1
        p1 = 1
        for q in p.args:
            if var not in q.atoms():
                c *= q
            else:
                p1 *= q
        p2 = c*taylor(p1,var,start,prec,dir,pol_pars,ov=0,level=1)
        if ov:
            p2 += O(var**prec)
        return p2

    for ring in [QQ, sympify]:
        te = TaylorEval(gens, ring, prec)
        try:
            p1 = te(p)
            lp = p1.lp
            p2 = 0
            if str(ring) == 'QQ':
                for m1,c1 in p1.iteritems():
                    c1 = Rational(c1.numerator,c1.denominator)
                    #print 'DB10', m1
                    m1 = monomial_tobasic(m1,*gens)
                    p2 += c1*m1
            else:
            # in the symbolic ring case one must expand
            # because the coefficients are often not expanded
                for m1,c1 in p1.iteritems():
                    c1 = c1.expand()
                    m1 = monomial_tobasic(m1,*gens)
                    p2 += c1*m1
            #args = list(p1._args + (O(var**prec),))
            #p1._args = tuple(args)
            # in the symbolic ring case one must expand
            # because the coefficients are often not expanded
            # TODO expand the coefficients
            if ov:
                p2 += O(var**prec)
            return p2
        except TaylorEvalError:
            continue
        except NotImplementedError:
            continue
    p1 = p.series(var,start,prec)
    return p1



class TaylorEval:
    def __init__(self,gens,ring,prec):
        """gens[0] is the series variable

        """
        self.prec = prec
        self.gens = gens
        self.var = gens[0]
        self.ngens = len(gens)
        # try first with ring QQ
        self.ring = ring
        self.lp = LPoly(['X%d' % i for i in range(self.ngens)],ring,O_lex)
        self.lvname = 'X0'
        self.lgens = self.lp.gens()
        self.lvar = self.lgens[0]
        self.dgens = dict(zip(self.gens,self.lgens))

    def coerce_number(self, a):
        ring = self.ring
        if self.lp.SR:
            return a
        if isinstance(a, Rational):
            if ring == QQ:
                return QQ(a.p,a.q)
            return a
        else:
            raise TaylorEvalError


    def __call__(self,f):
        if isinstance(f, Number):
            return self.coerce_number(f)
        if f in self.gens:
            return self.dgens[f]

        head = f.__class__
        if head == Add:
            s = self.lp(0)
            for x in f.args:
                if self.var in x.atoms():
                    x = self(x)
                elif x in self.gens:
                    x = self.dgens[x]
                else:
                    x = self.coerce_number(x)
                s += x
            return s
        if head == Mul:
            s = self.lp(1)
            for x in f.args:
                if self.var in x.atoms():
                    x = self(x)
                    s = s.mul_trunc(x,self.lvname,self.prec)
                elif x in self.gens:
                    x = self.dgens[x]
                    s = s.mul_trunc(x,self.lvname,self.prec)
                else:
                    x = self.coerce_number(x)
                    s = s*x
            return s
        if head == Pow:
            args = f.args
            pw = args[1]
            x = args[0]
            if self.var in x.atoms():
                x = self(x)
                if pw == int(pw):
                    x1 = x.pow_trunc(pw,self.lvname,self.prec)
                else:
                    if isinstance(pw, Rational):
                        num = int(pw.p)
                        den = int(pw.q)
                        x1 = x.pow_trunc(num,self.lvname,self.prec)
                        x1 = x1.nth_root(den,self.lvname,self.prec)
                    else:
                        raise NotImplementedError
                return x1
            x = self.coerce_number(x)
            return x
        if head == cos:
            q = ev_args(self,f.args)
            return q.cos(self.lvname,self.prec)
        if head == sin:
            q = ev_args(self,f.args)
            return q.sin(self.lvname,self.prec)
        if head == exp:
            q = ev_args(self,f.args)
            return q.exp(self.lvname,self.prec)
        if head == log:
            q = ev_args(self,f.args)
            return q.log(self.lvname,self.prec)
        if head == atan:
            q = ev_args(self,f.args)
            return q.atan(self.lvname,self.prec)
        if head == tan:
            q = ev_args(self,f.args)
            return q.tan(self.lvname,self.prec)
        if head == cosh:
            q = ev_args(self,f.args)
            return q.cosh(self.lvname,self.prec)
        if head == sinh:
            q = ev_args(self,f.args)
            return q.sinh(self.lvname,self.prec)
        if head == tanh:
            q = ev_args(self,f.args)
            return q.tanh(self.lvname,self.prec)
        if head == atanh:
            q = ev_args(self,f.args)
            return q.atanh(self.lvname,self.prec)
        if head == asin:
            q = ev_args(self,f.args)
            return q.asin(self.lvname,self.prec)
        if head == asinh:
            q = ev_args(self,f.args)
            return q.asinh(self.lvname,self.prec)
        if head == acos:
            if not self.lp.SR:
                raise NotImplementedError
            q = ev_args(self,f.args)
            return pi/2 - q.asin(self.lvname,self.prec)
        if head == acosh:
            if not self.lp.SR:
                raise NotImplementedError
            q = ev_args(self,f.args)
            return -I*pi/2 + I*q.asin(self.lvname,self.prec)
        if head == acot:
            if not self.lp.SR:
                raise NotImplementedError
            q = ev_args(self,f.args)
            return pi/2 + q.acot1(self.lvname,self.prec)
        if head == acoth:
            # see issue 564
            # sage gives taylor(acoth(x),x,0,9)
            # -1/2*I*pi + 1/9*x^9 + 1/7*x^7 + 1/5*x^5 + 1/3*x^3 + x
            # sage: acoth(0)
            # arccoth(0)
            if not self.lp.SR:
                raise NotImplementedError
            q = ev_args(self,f.args)
            if q == 0:
                return arcoth(0)
            return -I*pi/2 + q.re_acoth(self.lvname,self.prec)
        if head == LambertW:
            q = ev_args(self,f.args)
            return q.lambert(self.lvname,self.prec)
        raise NotImplementedError('case in __call__ not considered f=%s' % f)
