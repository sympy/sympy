from sympy import symbols, sin
from sympy.galgebra.printer import Format, xpdf
from sympy.galgebra.ga import Ga

Format()

print r'#\newline Normalized Spherical Coordinates:'
sph_coords = (r, th, phi) = symbols('r theta phi', real=True)
sp3d = Ga('e', g=[1, r**2, r**2 * sin(th)**2], coords=sph_coords, norm=True)
f = sp3d.mv('f', 'scalar', f=True)
F = sp3d.mv('F', 'vector', f=True)
lap = sp3d.grad*sp3d.grad

print 'g =', sp3d.g
print r'%\nabla =', sp3d.grad
print r'%\nabla^{2} =', lap
print r'%\lp\nabla^{2}\rp f =', lap*f
print r'%\nabla\cdot\lp\nabla f\rp =', sp3d.grad | (sp3d.grad * f)
print 'F =', F
print r'%\nabla\cdot F =', sp3d.grad | F


print '#Unnormalized Spherical Coordinates:'

sp3du = Ga('e', g=[1, r**2, r**2 * sin(th)**2], coords=sph_coords)
f = sp3du.mv('f', 'scalar', f=True)
F = sp3du.mv('F', 'vector', f=True)
lap = sp3du.grad*sp3du.grad
print 'g =', sp3du.g
print r'%\nabla =', sp3du.grad
print r'%\nabla^{2} =', lap
print r'%\lp\nabla^{2}\rp f =', (lap*f).trigsimp()
print r'%\nabla\cdot\lp\nabla f\rp =', sp3du.grad | (sp3du.grad * f)
print r'%\nabla\cdot\nabla =', sp3du.grad | sp3du.grad
print 'F =', F
print r'%\nabla\cdot F =', sp3du.grad | F

xpdf(paper='landscape',prog=True)
