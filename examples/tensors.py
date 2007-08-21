"""
http://www.lncc.br/~portugal/Canon.html
http://www.lncc.br/~portugal/Invar.html
http://www.ginac.de/tutorial/Indexed-objects.html
http://grtensor.phy.queensu.ca/

"""
import sys
sys.path.append(".")
sys.path.append("..")

from sympy import exp, Symbol, sin, Rational, Derivative, dsolve

from sympy.core import Basic, Function
from sympy.matrices import Matrix

class Indexed(Basic):
    def __init__(self, A, idxlist):
        self._args = [A, idxlist]

    def __str__(self):
        r = str(self[0])
        for idx in self[1]:
            r+=str(idx)
        return r

class Idx(Symbol):
    def __init__(self, name, dim = 4, up = True):
        Symbol.__init__(self, name)
        #self._args.extend([dim,up])
        self._name = name
        self._dim = dim
        self._up = up

    def __str__(self):
        if self._up:
            r = "^"
        else:
            r = "_"
        return r+self._name

    @property
    def up(self):
        return Idx(self._name, self._dim, True)

    @property
    def dn(self):
        return Idx(self._name, self._dim, False)

    def values(self):
        return range(self._dim)

t=Symbol("t")
r=Symbol("r")
theta=Symbol(r"\theta")
phi=Symbol(r"\phi")

class nu(Function):
    pass
class lam(Function):
    pass

gdd=Matrix(( 
    (-exp(nu(r)),0,0,0), 
    (0, exp(lam(r)), 0, 0),
    (0, 0, r**2, 0),
    (0, 0, 0, r**2*sin(theta)**2)
    ))

mu = Idx("mu")
nu = Idx("mu")
i = Idx("i")
m = Idx("m")
k = Idx("k")
l = Idx("l")
g = Indexed(Symbol("A"), [mu,nu])
Chr = g[i.up, m.up]/2 * (g[m.dn, k.dn].diff(l.up) + g[m.dn,l.dn].diff(k.up) \
        - g[k.dn, l.dn].diff(m.up))
#G = g.uu(i,m)/2 * (g.dd(m,k).diff(x[l])+g.dd(m,l).diff(x[k]) \
#                    - g.dd(k,l).diff(x[m]))

print g
print Chr
