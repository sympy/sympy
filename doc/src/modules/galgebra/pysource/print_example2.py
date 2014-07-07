from sympy.galgebra.printer import Format, xpdf
from sympy.galgebra.ga import Ga
from sympy import symbols

Format()
X = (x,y,z) = symbols('x y z')
o3d = Ga('e_x e_y e_z',g=[1,1,1],coords=X)

f = o3d.mv('f','scalar',f=True)
A = o3d.mv('A','vector',f=True)
B = o3d.mv('B','bivector',f=True)

print r'\bm{A} =',A
print r'\bm{B} =',B

print 'grad*f =',o3d.grad*f
print r'grad|\bm{A} =',o3d.grad|A
(o3d.grad*A).Fmt(2,r'grad*\bm{A}')

print r'-I*(grad^\bm{A}) =',-o3d.mv_I*(o3d.grad^A)
(o3d.grad*B).Fmt(2,r'grad*\bm{B}')
print r'grad^\bm{B} =',o3d.grad^B
print r'grad|\bm{B} =',o3d.grad|B

xpdf(paper='letter')
