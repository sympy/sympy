from sympy import symbols, sin, cos
from sympy.galgebra.ga import Ga
from sympy.galgebra.printer import Format, xpdf

Format()

(u, v) = uv = symbols('u,v',real=True)
(g2d, eu, ev) = Ga.build('e_u e_v', coords=uv)

print '#$A$ is a general 2D linear transformation'

A2d = g2d.lt('A')

print 'A =', A2d
print '\\f{\\det}{A} =', A2d.det()
print '\\f{\\Tr}{A} =', A2d.tr()

B2d = g2d.lt('B')

print 'B =', B2d
print 'A + B =', A2d + B2d
print 'AB =', A2d * B2d
print 'A - B =', A2d - B2d

a = g2d.mv('a','vector')
b = g2d.mv('b','vector')

print r'a|\f{\overline{A}}{b}-b|\f{\underline{A}}{a} =',((a|A2d.adj()(b))-(b|A2d(a))).simplify()

m4d = Ga('e_t e_x e_y e_z', g=[1, -1, -1, -1],coords=symbols('t,x,y,z',real=True))

T = m4d.lt('T')

print '#$T$ is a linear transformation in Minkowski space'
print r'\underline{T} =',T
print r'\overline{T} =',T.adj()
print r'\f{\mbox{tr}}{\underline{T}} =',T.tr()

a = m4d.mv('a','vector')
b = m4d.mv('b','vector')

print r'a|\f{\overline{T}}{b}-b|\f{\underline{T}}{a} =',((a|T.adj()(b))-(b|T(a))).simplify()

xpdf(paper='landscape',crop=True)
