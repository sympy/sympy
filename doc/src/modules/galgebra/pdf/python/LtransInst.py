from sympy import symbols, sin, cos, latex, Matrix
from ga import Ga
from printer import Format, xpdf

Format()
(x, y, z) = xyz = symbols('x,y,z',real=True)
(o3d, ex, ey, ez) = Ga.build('e_x e_y e_z', g=[1, 1, 1], coords=xyz)

A = o3d.lt('A')
print r'\mbox{General Instantiation: }A =', A
th = symbols('theta',real=True)
R = cos(th/2)+(ex^ey)*sin(th/2)
B = o3d.lt(R)
print r'\mbox{Rotor: }R =', R
print r'\mbox{Rotor Instantiation: }B =', B
dict1 = {ex:ey+ez,ez:ey+ez,ey:ex+ez}
C = o3d.lt(dict1)
print r'\mbox{Dictionary} =', latex(dict1)
print r'\mbox{Dictionary Instantiation: }C =', C
lst1 = [[1,0,1],[0,1,0],[1,0,1]]
D = o3d.lt(lst1)
print r'\mbox{List} =', latex(lst1)
print r'\mbox{List Instantiation: }D =', D
lst2 = [ey+ez,ex+ez,ex+ey]
E = o3d.lt(lst2)
print r'\mbox{List} =', latex(lst2)
print r'\mbox{List Instantiation: }E =', E
xpdf(paper=(10,12),crop=True)
