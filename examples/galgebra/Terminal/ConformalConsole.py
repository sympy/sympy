
from sympy import Symbol, symbols, sin, cos, Rational, expand, simplify, collect
from sympy import Rational as Rat
from sympy.galgebra.printer import Format, Eprint, Get_Program, Print_Function,xpdf
from sympy.galgebra.ga import Ga, one, zero
from sympy.galgebra.mv import Com, Nga
from math import sqrt

def radius(T):
    '''
    This retrieves the radius from a trivector representing a circle
    '''
    a=(T*T).scalar()
    b=((T^n)*(T^n)).scalar()
    return (-1*a/b)**0.5

def center(T):
    global n
    '''returns the center of a given circle trivector'''
    return T*n*T

def split_bivector(B):
    global ebar
    '''Implements the algorithm described in Doran and Lasenby to recover null vectors wedging to B'''
    print 'B =',B
    print 'B**2 =',B*B
    NB = B.norm()
    print 'NB =',NB
    Bh = B/NB
    ap = ebar - ((ebar^Bh)*Bh)
    a1 = ap + (ap*Bh)
    a2 = ap - (ap*Bh)
    #print '#a1 = ',a1
    #print '#a2 = ',a2
    return [a1,a2]

def norm(X):
    Y=sqrt((X*X).scalar())
    return Y

Get_Program(True)
Eprint()

g='1 0 0 0, \
   0 1 0 0, \
   0 0 0 2, \
   0 0 2 0'

c2d = Ga('e_1 e_2 n \\bar{n}',g=g)
(e1,e2,n,nbar) = c2d.mv()

global n,nbar,I

def F(x):
    global n,nbar
    Fx = ((x*x)*n+2*x-nbar) / 2
    return(Fx)

e = (n+nbar)/2
ebar = n - e
I=e1*e2*e*ebar

def intersect_lines(L1,L2):
    global I
    '''
    Computes the intersection bivector of two conformal lines L1, L2
    '''
    C = I*((I*L1)^(I*L2))
    return C


A=F(Rat(1,2)*e1)
B=F(2*e1)
C=F(Rat(4,5)*e1+Rat(3,5)*e2)
D=F(Rat(4,5)*e1-Rat(3,5)*e2)

print 'A =',A
print 'B =',B
print 'C =',C
print 'D =',D

T=A^B^C
print 'T =',T
U=F(e1)^(F(e2))^F(-1*e1)
print 'U =',U
inter=intersect_lines(U,T)
print 'inter =',inter

x,y = split_bivector(inter)

bases = (e1,e2)
print x.proj(bases)
print y.proj(bases)


print 'One intersection point x = ',x
print 'The other intersection point y = ',y
print 'x**2 = ',x*x
print 'y**2 = ',y*y
print 'T^x = ',T^x
print 'T^y = ',T^y
print 'U^x = ',U^x
print 'U^y = ',U^y
