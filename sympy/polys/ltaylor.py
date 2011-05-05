from sympy.polys.lpoly import *
from sympy import *
from sympy.core import C

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

def taylor(p,var=None,start=0,prec=6,dir="+",pol_pars=[],ov=True):
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

    >>> from sympy import *
    >>> from sympy.polys.ltaylor import taylor
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
    """

    if var == None or prec == None or dir != "+" or \
        prec in [S.Infinity, S.NegativeInfinity]:
        return series(p,var,start,prec)
    # case with ov=True; for ov=False taylor is faster than for ov=True
    # in sin(x), etc. series() is faster at high precision
    #              prec
    # sin(x)       100   series 2x faster
    #              200   series 2x faster
    #              1000  series 30x faster
    # similarly for cos
    # tan(x)       100   taylor 50% faster 
    #              200   taylor 20% faster
    #              1000  series 3x faster (chosen series for this)
    # log(1+x)           taylor 50% faster  
    # atan(x)            taylor 20% faster
    # for longer arguments taylor is faster
    # sin(x+x**2)  100   taylor 10x faster
    #              200   taylor 7x faster
    #              500   taylor 6x faster
    if prec > 70 and ov:
        head = p.__class__
        if head in [cos,sin,tan,acos,asin]:
            q = p.args[0]
            if q == var:
                return series(p,var,start,prec)
            if q.__class__ == Mul:
                if var in q.args:
                    b = 1
                    ni = q.args.index(var)
                    for i in range(len(q.args)):
                        if i != ni and isinstance(q.args[i],Number):
                            b = 0
                            break
                    if b:
                          return series(p,var,start,prec)
    if start:
        p0 = p
        p = p.subs(var,var+start)
    gens = [var] + pol_pars
    for ring in [QQ, sympify]:
        te = TaylorEval(gens, ring, prec)
        try:
            p1 = te(p)
            lp = p1.lp
            sb = Subs(lp,lp,{te.lvname:te.lvar-start})
            lvar = sb.subs(p1)
            p1 = p1.tobasic(*gens)
            #args = list(p1._args + (O(var**prec),))
            #p1._args = tuple(args)
            if ov:
                p1 = p1 + O(var**prec)
                if lp.SR:
                    p1 = p1.expand()
            return p1
        except TaylorEvalError:
            continue
        except NotImplementedError:
            continue
    #print 'DB5 used series',p
    if start:
        p = p0
    p1 = series(p,var,start,prec)
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
        raise NotImplementedError('case in __call__ not considered f=%s' % f)
 
