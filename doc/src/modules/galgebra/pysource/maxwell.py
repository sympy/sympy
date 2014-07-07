from sympy import symbols
from sympy.galgebra.printer import Format, xpdf
from sympy.galgebra.ga import Ga

Format()
X = symbols('t x y z',real=True)
(st4d,g0,g1,g2,g3) = Ga.build('gamma*t|x|y|z',g=[1,-1,-1,-1],coords=X)

I = st4d.i

B = st4d.mv('B','vector',f=True)
E = st4d.mv('E','vector',f=True)
B.set_coef(1,0,0)
E.set_coef(1,0,0)
B *= g0
E *= g0
J = st4d.mv('J','vector',f=True)
F = E+I*B

print r'\text{Pseudo Scalar\;\;}I =',I
print '\\text{Magnetic Field Bi-Vector\\;\\;} B = \\bm{B\\gamma_{t}} =',B
print '\\text{Electric Field Bi-Vector\\;\\;} E = \\bm{E\\gamma_{t}} =',E
print '\\text{Electromagnetic Field Bi-Vector\\;\\;} F = E+IB =',F
print '%\\text{Four Current Density\\;\\;} J =',J
gradF = st4d.grad*F
print '#Geom Derivative of Electomagnetic Field Bi-Vector'
gradF.Fmt(3,'grad*F')

print '#Maxwell Equations'
print 'grad*F = J'
print '#Div $E$ and Curl $H$ Equations'
(gradF.get_grade(1)-J).Fmt(3,'%\\grade{\\nabla F}_{1} -J = 0')
print '#Curl $E$ and Div $B$ equations'
(gradF.get_grade(3)).Fmt(3,'%\\grade{\\nabla F}_{3} = 0')

xpdf(paper='landscape')
