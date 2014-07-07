from sympy import symbols, sin
from sympy.galgebra.printer import Format, xpdf
from sympy.galgebra.ga import Ga

Format()

X = (r,th,phi) = symbols('r theta phi')
s3d = Ga('e_r e_theta e_phi',g=[1,r**2,r**2*sin(th)**2],coords=X,norm=True)
(er,eth,ephi) = s3d.mv()
grad = s3d.grad

f = s3d.mv('f','scalar',f=True)
A = s3d.mv('A','vector',f=True)
B = s3d.mv('B','bivector',f=True)

print 'f =',f
print 'A =',A
print 'B =',B

print 'grad*f =',grad*f
print 'grad|A =',grad|A
print '-I*(grad^A) =',(-s3d.i*(grad^A)).simplify()
print 'grad^B =',grad^B

xpdf(paper='letter')
