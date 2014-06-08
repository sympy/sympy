from sympy import symbols, sin, cos, latex
from ga import Ga
from printer import Format, xpdf

Format()
(x, y, z) = xyz = symbols('x,y,z',real=True)
(o3d, ex, ey, ez) = Ga.build('e_x e_y e_z', g=[1, 1, 1], coords=xyz)
grad = o3d.grad
(u, v) = uv = symbols('u,v',real=True)
(g2d, eu, ev) = Ga.build('e_u e_v', coords=uv)
grad_uv = g2d.grad
A = o3d.lt('A')
print '#3d orthogonal ($A,\\;B$ are linear transformations)'
print 'A =', A
print r'\f{\operatorname{mat}}{A} =', latex(A.matrix())
print '\\f{\\det}{A} =', A.det()
print '\\overline{A} =', A.adj()
print '\\f{\\Tr}{A} =', A.tr()
print '\\f{A}{e_x^e_y} =', A(ex^ey)
print '\\f{A}{e_x}^\\f{A}{e_y} =', A(ex)^A(ey)
B = o3d.lt('B')
print 'A + B =', A + B
print 'AB =', A * B
print 'A - B =', A - B

print '#2d general ($A,\\;B$ are linear transformations)'
A2d = g2d.lt('A')
print 'A =', A2d
print '\\f{\\det}{A} =', A2d.det()
#A2d.adj().Fmt(4,'\\overline{A}')
print '\\f{\\Tr}{A} =', A2d.tr()
print '\\f{A}{e_u^e_v} =', A2d(eu^ev)
print '\\f{A}{e_u}^\\f{A}{e_v} =', A2d(eu)^A2d(ev)
B2d = g2d.lt('B')
print 'B =', B2d
print 'A + B =', A2d + B2d
print 'AB =', A2d * B2d
print 'A - B =', A2d - B2d
a = g2d.mv('a','vector')
b = g2d.mv('b','vector')
print r'a|\f{\overline{A}}{b}-b|\f{\underline{A}}{a} =',((a|A2d.adj()(b))-(b|A2d(a))).simplify()

print '#4d Minkowski spaqce (Space Time)'
m4d = Ga('e_t e_x e_y e_z', g=[1, -1, -1, -1],coords=symbols('t,x,y,z',real=True))
T = m4d.lt('T')
print 'g =', m4d.g
print r'\underline{T} =',T
print r'\overline{T} =',T.adj()
#m4d.mv(T.det()).Fmt(4,r'\f{\det}{\underline{T}}')
print r'\f{\mbox{tr}}{\underline{T}} =',T.tr()
a = m4d.mv('a','vector')
b = m4d.mv('b','vector')
print r'a|\f{\overline{T}}{b}-b|\f{\underline{T}}{a} =',((a|T.adj()(b))-(b|T(a))).simplify()
xpdf(paper='landscape')
