#!/usr/bin/env python

from __future__ import division, print_function

from sympy import sin, cos, sinh, cosh, symbols, expand, simplify
from sympy.galgebra import xdvi
from sympy.galgebra import MV, Format, Com

def main():
    Format()

    (ex, ey, ez) = MV.setup('e*x|y|z')
    A = MV('A', 'mv')
    print(r'\bm{A} =', A)
    A.Fmt(2, r'\bm{A}')
    A.Fmt(3, r'\bm{A}')

    X = (x, y, z) = symbols('x y z')
    (ex, ey, ez, grad) = MV.setup('e_x e_y e_z', metric='[1,1,1]', coords=X)

    f = MV('f', 'scalar', fct=True)
    A = MV('A', 'vector', fct=True)
    B = MV('B', 'grade2', fct=True)

    print(r'\bm{A} =', A)
    print(r'\bm{B} =', B)

    print('grad*f =', grad*f)
    print(r'grad|\bm{A} =', grad | A)
    print(r'grad*\bm{A} =', grad*A)

    print(r'-I*(grad^\bm{A}) =', -MV.I*(grad ^ A))
    print(r'grad*\bm{B} =', grad*B)
    print(r'grad^\bm{B} =', grad ^ B)
    print(r'grad|\bm{B} =', grad | B)

    (a, b, c, d) = MV.setup('a b c d')

    print('g_{ij} =', MV.metric)

    print('\\bm{a|(b*c)} =', a | (b*c))
    print('\\bm{a|(b^c)} =', a | (b ^ c))
    print('\\bm{a|(b^c^d)} =', a | (b ^ c ^ d))
    print('\\bm{a|(b^c)+c|(a^b)+b|(c^a)} =', (a | (b ^ c)) + (c | (a ^ b)) + (b | (c ^ a)))
    print('\\bm{a*(b^c)-b*(a^c)+c*(a^b)} =', a*(b ^ c) - b*(a ^ c) + c*(a ^ b))
    print('\\bm{a*(b^c^d)-b*(a^c^d)+c*(a^b^d)-d*(a^b^c)} =', a*(b ^ c ^ d) - b*(a ^ c ^ d) + c*(a ^ b ^ d) - d*(a ^ b ^ c))
    print('\\bm{(a^b)|(c^d)} =', (a ^ b) | (c ^ d))
    print('\\bm{((a^b)|c)|d} =', ((a ^ b) | c) | d)
    print('\\bm{(a^b)\\times (c^d)} =', Com(a ^ b, c ^ d))

    metric = '1 # #,' + \
             '# 1 #,' + \
             '# # 1,'

    (e1, e2, e3) = MV.setup('e1 e2 e3', metric)

    E = e1 ^ e2 ^ e3
    Esq = (E*E).scalar()
    print('E =', E)
    print('%E^{2} =', Esq)
    Esq_inv = 1/Esq

    E1 = (e2 ^ e3)*E
    E2 = (-1)*(e1 ^ e3)*E
    E3 = (e1 ^ e2)*E

    print('E1 = (e2^e3)*E =', E1)
    print('E2 =-(e1^e3)*E =', E2)
    print('E3 = (e1^e2)*E =', E3)

    print('E1|e2 =', (E1 | e2).expand())
    print('E1|e3 =', (E1 | e3).expand())
    print('E2|e1 =', (E2 | e1).expand())
    print('E2|e3 =', (E2 | e3).expand())
    print('E3|e1 =', (E3 | e1).expand())
    print('E3|e2 =', (E3 | e2).expand())
    w = ((E1 | e1).expand()).scalar()
    Esq = expand(Esq)
    print('%(E1\\cdot e1)/E^{2} =', simplify(w/Esq))
    w = ((E2 | e2).expand()).scalar()
    print('%(E2\\cdot e2)/E^{2} =', simplify(w/Esq))
    w = ((E3 | e3).expand()).scalar()
    print('%(E3\\cdot e3)/E^{2} =', simplify(w/Esq))

    X = (r, th, phi) = symbols('r theta phi')
    curv = [[r*cos(phi)*sin(th), r*sin(phi)*sin(th), r*cos(th)], [1, r, r*sin(th)]]
    (er, eth, ephi, grad) = MV.setup('e_r e_theta e_phi', metric='[1,1,1]', coords=X, curv=curv)

    f = MV('f', 'scalar', fct=True)
    A = MV('A', 'vector', fct=True)
    B = MV('B', 'grade2', fct=True)

    print('A =', A)
    print('B =', B)

    print('grad*f =', grad*f)
    print('grad|A =', grad | A)
    print('-I*(grad^A) =', -MV.I*(grad ^ A))
    print('grad^B =', grad ^ B)

    vars = symbols('t x y z')
    (g0, g1, g2, g3, grad) = MV.setup('gamma*t|x|y|z', metric='[1,-1,-1,-1]', coords=vars)
    I = MV.I

    B = MV('B', 'vector', fct=True)
    E = MV('E', 'vector', fct=True)
    B.set_coef(1, 0, 0)
    E.set_coef(1, 0, 0)
    B *= g0
    E *= g0
    J = MV('J', 'vector', fct=True)
    F = E + I*B

    print('B = \\bm{B\\gamma_{t}} =', B)
    print('E = \\bm{E\\gamma_{t}} =', E)
    print('F = E+IB =', F)
    print('J =', J)
    gradF = grad*F
    gradF.Fmt(3, 'grad*F')

    print('grad*F = J')
    (gradF.grade(1) - J).Fmt(3, '%\\grade{\\nabla F}_{1} -J = 0')
    (gradF.grade(3)).Fmt(3, '%\\grade{\\nabla F}_{3} = 0')

    (alpha, beta, gamma) = symbols('alpha beta gamma')

    (x, t, xp, tp) = symbols("x t x' t'")
    (g0, g1) = MV.setup('gamma*t|x', metric='[1,-1]')

    R = cosh(alpha/2) + sinh(alpha/2)*(g0 ^ g1)
    X = t*g0 + x*g1
    Xp = tp*g0 + xp*g1
    print('R =', R)

    print(r"#%t\bm{\gamma_{t}}+x\bm{\gamma_{x}} = t'\bm{\gamma'_{t}}+x'\bm{\gamma'_{x}} = R\lp t'\bm{\gamma_{t}}+x'\bm{\gamma_{x}}\rp R^{\dagger}")

    Xpp = R*Xp*R.rev()
    Xpp = Xpp.collect([xp, tp])
    Xpp = Xpp.subs({2*sinh(alpha/2)*cosh(alpha/2): sinh(alpha), sinh(alpha/2)**2 + cosh(alpha/2)**2: cosh(alpha)})
    print(r"%t\bm{\gamma_{t}}+x\bm{\gamma_{x}} =", Xpp)
    Xpp = Xpp.subs({sinh(alpha): gamma*beta, cosh(alpha): gamma})

    print(r'%\f{\sinh}{\alpha} = \gamma\beta')
    print(r'%\f{\cosh}{\alpha} = \gamma')

    print(r"%t\bm{\gamma_{t}}+x\bm{\gamma_{x}} =", Xpp.collect(gamma))

    vars = symbols('t x y z')
    (g0, g1, g2, g3, grad) = MV.setup('gamma*t|x|y|z', metric='[1,-1,-1,-1]', coords=vars)
    I = MV.I
    (m, e) = symbols('m e')

    psi = MV('psi', 'spinor', fct=True)
    A = MV('A', 'vector', fct=True)
    sig_z = g3*g0
    print('\\bm{A} =', A)
    print('\\bm{\\psi} =', psi)

    dirac_eq = (grad*psi)*I*sig_z - e*A*psi - m*psi*g0
    dirac_eq.simplify()

    dirac_eq.Fmt(3, r'\nabla \bm{\psi} I \sigma_{z}-e\bm{A}\bm{\psi}-m\bm{\psi}\gamma_{t} = 0')

    xdvi()
    return

if __name__ == "__main__":
    main()
