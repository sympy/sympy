"""
This example calculates the Ricci tensor from the metric and does this
on the example of Schwarzschild solution.
"""
import sys
sys.path.append(".")
sys.path.append("..")

from sympy import exp, Symbol, sin, Rational, Derivative, dsolve

from sympy.core import Basic, Function
from sympy.modules.matrices import Matrix

def grad(f,X):
    a=[]
    for x in X:
        a.append( f.diff(x) )
    return a

def d(m,x):
    return grad(m[0,0],x)

class MT(object):
    def __init__(self,m):
        self.gdd=m
        self.guu=m.inv()

    def __str__(self):
        return "g_dd =\n" + str(self.gdd)

    def dd(self,i,j):
        return self.gdd[i,j]

    def uu(self,i,j):
        return self.guu[i,j]

class G(object):
    def __init__(self,g,x):
        self.g = g
        self.x = x

    def udd(self,i,k,l):
        g=self.g
        x=self.x
        r=0
        for m in [0,1,2,3]:
            r+=g.uu(i,m)/2 * (g.dd(m,k).diff(x[l])+g.dd(m,l).diff(x[k]) \
                    - g.dd(k,l).diff(x[m]))
        return r

class Riemann(object):
    def __init__(self,G,x):
        self.G = G
        self.x = x

    def uddd(self,rho,sigma,mu,nu):
        G=self.G
        x=self.x
        r=G.udd(rho,nu,sigma).diff(x[mu])-G.udd(rho,mu,sigma).diff(x[nu])
        for lam in [0,1,2,3]:
            r+=G.udd(rho,mu,lam)*G.udd(lam,nu,sigma) \
                -G.udd(rho,nu,lam)*G.udd(lam,mu,sigma)
        return r

class Ricci(object):
    def __init__(self,R,x):
        self.R = R
        self.x = x
        self.g = R.G.g

    def dd(self,mu,nu):
        R=self.R
        x=self.x
        r=0
        for lam in [0,1,2,3]:
            r+=R.uddd(lam,mu,lam,nu)
        return r

    def ud(self,mu,nu):
        r=0
        for lam in [0,1,2,3]:
            r+=self.g.uu(mu,lam)*self.dd(lam,nu)
        return r.expand()

def curvature(Rmn):
    return Rmn.ud(0,0)+Rmn.ud(1,1)+Rmn.ud(2,2)+Rmn.ud(3,3)

class nu(Function):
    def getname(self):
        return r"\nu"
        return r"nu"

class lam(Function):
    def getname(self):
        return r"\lambda"
        return r"lambda"

t=Symbol("t")
r=Symbol("r")
theta=Symbol(r"\theta")
phi=Symbol(r"\phi")

#general, spherically symmetric metric
gdd=Matrix(( 
    (-exp(nu(r)),0,0,0), 
    (0, exp(lam(r)), 0, 0),
    (0, 0, r**2, 0),
    (0, 0, 0, r**2*sin(theta)**2)
    ))
#spherical - flat
#gdd=Matrix(( 
#    (-1, 0, 0, 0), 
#    (0, 1, 0, 0),
#    (0, 0, r**2, 0),
#    (0, 0, 0, r**2*sin(theta)**2)
#    ))
#polar - flat
#gdd=Matrix(( 
#    (-1, 0, 0, 0), 
#    (0, 1, 0, 0),
#    (0, 0, 1, 0),
#    (0, 0, 0, r**2)
#    ))
#polar - on the sphere, on the north pole
#gdd=Matrix(( 
#    (-1, 0, 0, 0), 
#    (0, 1, 0, 0),
#    (0, 0, r**2*sin(theta)**2, 0),
#    (0, 0, 0, r**2)
#    ))
g=MT(gdd)
X=(t,r,theta,phi)
Gamma=G(g,X)
Rmn=Ricci(Riemann(Gamma,X),X)

if __name__ == "__main__":
    #print g
    print "-"*40
    print "Christoffel symbols:"
    print Gamma.udd(0,1,0)
    print Gamma.udd(0,0,1)
    print
    print Gamma.udd(1,0,0)
    print Gamma.udd(1,1,1)
    print Gamma.udd(1,2,2)
    print Gamma.udd(1,3,3)
    print
    print Gamma.udd(2,2,1)
    print Gamma.udd(2,1,2)
    print Gamma.udd(2,3,3)
    print
    print Gamma.udd(3,2,3)
    print Gamma.udd(3,3,2)
    print Gamma.udd(3,1,3)
    print Gamma.udd(3,3,1)
    print "-"*40
    print "Ricci tensor:"
    print Rmn.dd(0,0)
    e =  Rmn.dd(1,1)
    print e
    print Rmn.dd(2,2)
    print Rmn.dd(3,3)
    #print
    #print "scalar curvature:"
    #print curvature(Rmn)
    print "-"*40
    print "solve the Einstein's equations:"
    e = e.subs(nu(r), -lam(r))
    l =  dsolve(e, [lam(r)])
    print lam(r)," = ",l
    metric = gdd.subs(lam(r), l).subs(nu(r),-l).combine()
    print "metric:"
    print metric
