from sympy import symbols
from sympy.galgebra.ga import Ga
from sympy.galgebra.printer import Format, xpdf

Format()

(alpha,beta,gamma) = symbols('alpha beta gamma')
(x,t,xp,tp) = symbols("x t x' t'",real=True)
(st2d,g0,g1) = Ga.build('gamma*t|x',g=[1,-1])

from sympy import sinh,cosh

R = cosh(alpha/2)+sinh(alpha/2)*(g0^g1)
X = t*g0+x*g1
Xp = tp*g0+xp*g1
print 'R =',R

print r"#%t\bm{\gamma_{t}}+x\bm{\gamma_{x}} = t'\bm{\gamma'_{t}}+x'\bm{\gamma'_{x}} = R\lp t'\bm{\gamma_{t}}+x'\bm{\gamma_{x}}\rp R^{\dagger}"

Xpp = R*Xp*R.rev()
Xpp = Xpp.collect()
Xpp = Xpp.trigsimp()
print r"%t\bm{\gamma_{t}}+x\bm{\gamma_{x}} =",Xpp
Xpp = Xpp.subs({sinh(alpha):gamma*beta,cosh(alpha):gamma})

print r'%\f{\sinh}{\alpha} = \gamma\beta'
print r'%\f{\cosh}{\alpha} = \gamma'

print r"%t\bm{\gamma_{t}}+x\bm{\gamma_{x}} =",Xpp.collect()

xpdf(paper='letter')
