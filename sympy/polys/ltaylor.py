from sympy.polys.lpoly import *
from sympy import *

class TaylorEvalError(TypeError):
    pass


def ev_args(te,a):
    #print 'DB20', a
    if len(a) == 1:
        a = a[0]
        #print 'DB21', a,type(a)
        if a == te.var:
            return te.lvar
        if isinstance(a, Number):
            return te.coerce_number(a)
        return te(a)
    else:
        raise NotImplementedError

def taylor(p,var,start,prec):
    #print 'DB9', p, var,prec
    if start:
        raise NotImplementedError
        p0 = p
        p = p.subs(var,var-start)
    atoms = p.atoms()
    #print 'DB9a atoms=',atoms, var
    gens = [_x for _x in atoms if not _x.is_Number and not isinstance(_x, tuple) and _x != var]
    #print 'DB10', gens
    gens = [var] + gens
    for ring in [QQ, sympify]:
        #print 'DB19 ring=', ring
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
        #print 'DB20', self.lp
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
        #print 'DB22', f, self.gens, f in self.gens
        if f in self.gens:
            return self.dgens[f]

        head = f.__class__
        #print 'DB10 head=%s prec=%d' %(head,self.prec)
        if head == Add:
            #print 'DB11', f.args
            s = self.lp(0)
            for x in f.args:
              #print 'DB12', x
                if self.var in x.atoms():
                    x = self(x)
                elif x in self.gens:
                    x = self.dgens[x]
                else:
                    #print 'DB23', x
                    x = self.coerce_number(x)
                s += x
            #s = sum([self(x) for x in f.args])
            return s
        if head == Mul:
            #print 'DB11Mul', f.args
            s = self.lp(1)
            for x in f.args:
                #print 'DB12', x
                if self.var in x.atoms():
                    #print 'DB12 Mul x=',x
                    x = self(x)
                    s = s.mul_trunc(x,self.lvname,self.prec)
                    #print 'DB12a Mul x=',x
                elif x in self.gens:
                    x = self.dgens[x]
                    s = s.mul_trunc(x,self.lvname,self.prec)
                else:
                    x = self.coerce_number(x)
                    s = s*x
                #s = s.mul_trunc(x,self.lvname,self.prec)
            return s
        if head == Pow:
            #print 'DB12'
            args = f.args
            pw = args[1]
            x = args[0]
            #print 'DB13', x
            if self.var in x.atoms():
                x = self(x)
                #print 'DB13a x=%s pw=%s' %(x, pw), self.lvname,self.prec
                if pw == int(pw):
                    #x1 = x.fun('pow_trunc', pw,self.lvname,self.prec)
                    x1 = x.pow_trunc(pw,self.lvname,self.prec)
                    #print 'DB13b',x1
                else:
                    if isinstance(pw, Rational):
                        num = int(pw.p)
                        den = int(pw.q)
                        x1 = x.pow_trunc(num,self.lvname,self.prec)
                        x1 = x1.nth_root(den,self.lvname,self.prec)
                    else:
                        raise NotImplementedError
                #print 'DB13b', x1
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
        raise TypeError('case in __call__ not considered f=%s' % f)
 
