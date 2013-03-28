
from sympy import symbols,sin,cos,factor_terms,simplify
from sympy.GA.GAPrint import enhance_print
from sympy.GA.GA import MV

enhance_print()

X = (x,y,z) = symbols('x y z')
(ex,ey,ez,grad) = MV.setup('e_x e_y e_z',metric='[1,1,1]',coords=(x,y,z))

A = x*(ey^ez) + y*(ez^ex) + z*(ex^ey)
print 'A =', A
print 'grad^A =',(grad^A).simplify()
print

f = MV('f','scalar',fct=True)
f = (x**2 + y**2 + z**2)**(-1.5)
print 'f =', f
print 'grad*f =',(grad*f).expand()
print

B = f*A
print 'B =', B
print

Curl_B = grad^B

print 'grad^B =', Curl_B.simplify()

def Symplify(A):
    return(factor_terms(simplify(A)))

print Curl_B.func(Symplify)

