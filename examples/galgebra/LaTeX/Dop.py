from sympy import symbols, sin
from sympy.galgebra.printer import Format, xpdf
from sympy.galgebra.ga import Ga

Format()
xyz_coords = (x, y, z) = symbols('x y z', real=True)
(o3d, ex, ey, ez) = Ga.build('e', g=[1, 1, 1], coords=xyz_coords, norm=True)
f = o3d.mv('f', 'scalar', f=True)
lap = o3d.grad*o3d.grad
print r'grad =', o3d.grad
print r'%\nabla^{2} = \nabla\cdot\nabla =', lap
print r'%\lp\nabla^{2}\rp f =', lap*f
print r'%\nabla\cdot\lp\nabla f\rp =', o3d.grad | (o3d.grad * f)

sph_coords = (r, th, phi) = symbols('r theta phi', real=True)
(sp3d, er, eth, ephi) = Ga.build('e', g=[1, r**2, r**2 * sin(th)**2], coords=sph_coords, norm=True)
f = sp3d.mv('f', 'scalar', f=True)
lap = sp3d.grad*sp3d.grad
print r'%\nabla =', sp3d.grad
print r'%\nabla^{2} = \nabla\cdot\nabla =', lap
print r'%\lp\nabla^{2}\rp f =', lap*f
print r'%\nabla\cdot\lp\nabla f\rp =', sp3d.grad | (sp3d.grad * f)
A = o3d.mv('A','vector')
K = o3d.mv('K','vector')
xs = o3d.mv(x)
X = o3d.X()

print o3d.grad*A
print A*o3d.grad
print o3d.grad*xs
print xs*o3d.grad
print o3d.grad*(o3d.grad+xs)
print (o3d.grad+xs)*o3d.grad

print X
print X^o3d.grad
print (X^o3d.grad).components()

xpdf(filename='prod.tex')
