from sympy.polys.lpoly import *
from sympy import *

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

def taylor(p,var,start,prec,pol_pars=[]):
    """
    taylor series expansion of p
    var series variable
    start  var=start point of expansion
    prec precision of the series
    pol_pars polynomial parameters

    ALGORITHM try first to compute the series in the
    QQ ring; if this fails compute it
    in the symbolic ring SR consisting
    of the sympy expressions which do not depend on
    var and pol_pars; if also this fais compute it
    using the series function
    """
    if start:
        raise NotImplementedError
        p0 = p
        p = p.subs(var,var-start)
    gens = [var] + pol_pars
    for ring in [QQ, sympify]:
        te = TaylorEval(gens, ring, prec)
        try:
            p1 = te(p)
            # TODO add subs to Poly
            #p1.subs(X=lp.lvar+te(start))
            p1 = p1.tobasic(*gens)
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
        raise NotImplementedError('case in __call__ not considered f=%s' % f)
 
