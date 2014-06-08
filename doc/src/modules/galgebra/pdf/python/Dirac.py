
import sys
from sympy import symbols,sin,cos
from printer import Format,xpdf,Get_Program,Print_Function
from ga import Ga

Format()
coords = symbols('t x y z',real=True)
(st4d,g0,g1,g2,g3) = Ga.build('gamma*t|x|y|z',g=[1,-1,-1,-1],coords=coords)
I = st4d.i

(m,e) = symbols('m e')

psi = st4d.mv('psi','spinor',f=True)
A = st4d.mv('A','vector',f=True)
sig_z = g3*g0

print '\\text{4-Vector Potential\\;\\;}\\bm{A} =',A
print '\\text{8-component real spinor\\;\\;}\\bm{\\psi} =',psi

dirac_eq = (st4d.grad*psi)*I*sig_z-e*A*psi-m*psi*g0
dirac_eq = dirac_eq.simplify()

dirac_eq.Fmt(3,r'%\text{Dirac Equation\;\;}\nabla \bm{\psi}'+\
               r' I \sigma_{z}-e\bm{A}\bm{\psi}-m\bm{\psi}\gamma_{t} = 0')
xpdf(paper='landscape',prog=True)
