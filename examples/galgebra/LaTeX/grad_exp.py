# -*- coding: utf-8 -*-

#computer algebra system
from sympy import *
from sympy.matrices import *
from sympy.galgebra.ga import Ga
from sympy.galgebra.printer import Format, xpdf

########################################################################
#ALGEBRA & DEFINITIONS
########################################################################
#Clifford(1,4)
#Flat space, no metric, just signature
#All constants = 1
#Dimensions

########################################################################
#PHYSICS
########################################################################

Format()

print '#Results with all scalar variables declared as real'

vars = t, x, y, z, w = symbols('t x y z w',real=True)
E = symbols('E',real=True)

myBasis='gamma_t gamma_x gamma_y gamma_z gamma_w'
st4d, gt, gx, gy, gz, gw, = Ga.build(myBasis,g=[1,-1,-1,-1,-1],coords=vars)
unit=st4d.i*st4d.i
X = t*gt+x*gx+y*gy+z*gz+w*gw
K = st4d.mv('k','vector')
print 'X =',X
print 'K =',K
print 'K|X =',K|X
print '%I^{2} =',unit
Ixyzw = st4d.i*gx*gy*gz*gw
print r'%I_{xyzw} = I\gamma_{x}\gamma_{y}\gamma_{z}\gamma_{w} =',Ixyzw
print r'%\lp I\gamma_{x}\gamma_{y}\gamma_{z}\gamma_{w}\rp^{2} =',Ixyzw*Ixyzw
grad = st4d.grad
#For symbolic exponent exp() needs hint on whether square of exponent is + or -
the_exponential= (-E*gw*t).exp(hint='-')

#Accepted the 2 following statements
the_exponential.Fmt(fmt=1, title=r'%e^{-E\gamma_{w}t}')
(-E*gw*the_exponential*t).Fmt(fmt=1, title='%-E\gamma_{w}e^{-E\gamma_{w}t}\gamma_{t}')

#why the simplification doesn't take place here ?
(grad*the_exponential).Fmt(fmt=1, title=r'%\bm{\nabla} e^{-E\gamma_{w}t}')

(grad*the_exponential+
E*gw*the_exponential*t).Fmt(fmt=1,
    title=r'%\bm{\nabla}e^{-E\gamma_{w}t}+E\gamma_{w}e^{-E\gamma_{w}t}\gamma_{t}')

EXP = (Ixyzw*(K|X)).exp(hint='-').simplify()

print r'%e^{I_{xyzw}K\cdot X} =', EXP
(grad*EXP).Fmt(2,r'%\bm{\nabla}e^{I_{xyzw}K\cdot X}')


xpdf()
