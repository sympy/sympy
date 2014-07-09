from sympy import symbols, sin, pi, latex
from sympy.galgebra.ga import Ga
from sympy.galgebra.printer import Format, xpdf

Format()
coords = (r, th, phi) = symbols('r,theta,phi', real=True)
sp3d = Ga('e_r e_th e_ph', g=[1, r**2, r**2*sin(th)**2], \
       coords=coords, norm=True)

sph_uv = (u, v) = symbols('u,v', real=True)
sph_map = [1, u, v]  # Coordinate map for sphere of r = 1
sph2d = sp3d.sm(sph_map,sph_uv)

print r'(u,v)\rightarrow (r,\theta,\phi) = ',latex(sph_map)
print 'g =',latex(sph2d.g)
F = sph2d.mv('F','vector',f=True) #scalar function
f = sph2d.mv('f','scalar',f=True) #vector function
print r'\nabla f =',sph2d.grad * f
print 'F =',F
print r'\nabla F = ',sph2d.grad * F

cir_s = s = symbols('s',real=True)
cir_map = [pi/8,s]
cir1d = sph2d.sm(cir_map,(cir_s,))

print 'g =',latex(cir1d.g)
h = cir1d.mv('h','scalar',f=True)
H = cir1d.mv('H','vector',f=True)
print r'(s)\rightarrow (u,v) = ',latex(cir_map)
print 'H =', H
print latex(H)
print r'\nabla h =', cir1d.grad * h
print r'\nabla H =', cir1d.grad * H
xpdf(filename='submanifold.tex',paper=(6,5),crop=True)
