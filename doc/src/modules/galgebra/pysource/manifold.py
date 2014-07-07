import sys
from sympy import symbols, sin, latex, diff, Function, expand
from sympy.galgebra.ga import Ga
from sympy.galgebra.lt import Mlt
from sympy.galgebra.printer import Eprint, Format, xpdf

Format()

#Define spherical coordinate system in 3-d

coords = (r, th, phi) = symbols('r,theta,phi', real=True)

sp3d = Ga('e_r e_th e_ph', g=[1, r**2, r**2*sin(th)**2], coords=coords)
(er, eth, ephi) = sp3d.mv()

#Define coordinates for 2-d (u,v) and 1-d (s) manifolds

u,v,s,alpha = symbols('u v s alpha',real=True)

sub_coords = (u,v)

smap = [1, u, v]  # Coordinate map for sphere of r = 1 in 3-d

print r'(u,v)\rightarrow (r,\theta,\phi) = ',smap

#Define unit sphere manifold

sph2d = sp3d.sm(smap,sub_coords)

print '#Unit Sphere Manifold:'

print 'g =',sph2d.g

(eu,ev) = sph2d.mv()

#Define vector and vector field on unit sphere tangent space

a = sph2d.mv('a','vector')
b = sph2d.mv('b','vector')
c = sph2d.mv('c','vector')
f = sph2d.mv('f','vector',f=True)

print 'a =', a
print 'f =', f
print r'%\nabla =', sph2d.grad

#Define directional derivative in direction a for unit sphere manifold

dd = a|sph2d.grad

print r'%a\cdot\nabla =', dd
print r'%\paren{a\cdot\nabla}\bm{e}_u =', dd * eu
print r'%\paren{a\cdot\nabla}\bm{e}_v =', dd * ev
print r'%\paren{a\cdot\nabla}f =',dd * f

V = Mlt('V',sph2d,nargs=1,fct=True)
T = Mlt('T',sph2d,nargs=2,fct=True)

print '#Tensors on the Unit Sphere'

print 'V =', V
#print 'T =', T
T.Fmt(5,'T')
print '#Tensor Contraction'
print r'T[1,2] =', T.contract(1,2)
print '#Tensor Evaluation'
print r'T(a,b) =', T(a,b)
print r'T(a,b+c) =', T(a,b+c).expand()
print r'T(a,\alpha b) =', T(a,alpha*b)
print '#Geometric Derivative With Respect To Slot'
print r'\nabla_{a_{1}}T =',T.pdiff(1)
print r'\nabla_{a_{2}}T =',T.pdiff(2)
print '#Covariant Derivatives'
print r'\mathcal{D}V =', V.cderiv()
DT = T.cderiv()
#print r'\mathcal{D}T =', DT
DT.Fmt(3,r'\mathcal{D}T')
print r'\mathcal{D}T[1,3](a) =', (DT.contract(1,3))(a)

#Define curve on unit sphere manifold

us = Function('u__s')(s)
vs = Function('v__s')(s)

xpdf(paper='landscape')
