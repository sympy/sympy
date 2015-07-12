#!/usr/bin/python

from __future__ import division, print_function

from sympy import Symbol, symbols, sin, cos, Rational, expand, simplify, collect, S
from sympy.galgebra import enhance_print, Get_Program, Print_Function
from sympy.galgebra import MV, Format, Com, Nga
from sympy.galgebra.printing import GA_Printer


def basic_multivector_operations():
    Print_Function()
    (ex, ey, ez) = MV.setup('e*x|y|z')

    A = MV('A', 'mv')

    A.Fmt(1, 'A')
    A.Fmt(2, 'A')
    A.Fmt(3, 'A')

    X = MV('X', 'vector')
    Y = MV('Y', 'vector')

    print('g_{ij} =\n', MV.metric)

    X.Fmt(1, 'X')
    Y.Fmt(1, 'Y')

    (X*Y).Fmt(2, 'X*Y')
    (X ^ Y).Fmt(2, 'X^Y')
    (X | Y).Fmt(2, 'X|Y')

    (ex, ey) = MV.setup('e*x|y')

    print('g_{ij} =\n', MV.metric)

    X = MV('X', 'vector')
    A = MV('A', 'spinor')

    X.Fmt(1, 'X')
    A.Fmt(1, 'A')

    (X | A).Fmt(2, 'X|A')
    (X < A).Fmt(2, 'X<A')
    (A > X).Fmt(2, 'A>X')

    (ex, ey) = MV.setup('e*x|y', metric='[1,1]')

    print('g_{ii} =\n', MV.metric)

    X = MV('X', 'vector')
    A = MV('A', 'spinor')

    X.Fmt(1, 'X')
    A.Fmt(1, 'A')

    (X*A).Fmt(2, 'X*A')
    (X | A).Fmt(2, 'X|A')
    (X < A).Fmt(2, 'X<A')
    (X > A).Fmt(2, 'X>A')

    (A*X).Fmt(2, 'A*X')
    (A | X).Fmt(2, 'A|X')
    (A < X).Fmt(2, 'A<X')
    (A > X).Fmt(2, 'A>X')
    return


def check_generalized_BAC_CAB_formulas():
    Print_Function()

    (a, b, c, d, e) = MV.setup('a b c d e')

    print('g_{ij} =\n', MV.metric)

    print('a|(b*c) =', a | (b*c))
    print('a|(b^c) =', a | (b ^ c))
    print('a|(b^c^d) =', a | (b ^ c ^ d))
    print('a|(b^c)+c|(a^b)+b|(c^a) =', (a | (b ^ c)) +
          (c | (a ^ b)) +
          (b | (c ^ a)))
    print('a*(b^c)-b*(a^c)+c*(a^b) =', a*(b ^ c) - b*(a ^ c) + c*(a ^ b))
    print('a*(b^c^d)-b*(a^c^d)+c*(a^b^d)-d*(a^b^c) =', a *
          (b ^ c ^ d) -
          b *
          (a ^ c ^ d) +
          c *
          (a ^ b ^ d) -
          d *
          (a ^ b ^ c))
    print('(a^b)|(c^d) =', (a ^ b) | (c ^ d))
    print('((a^b)|c)|d =', ((a ^ b) | c) | d)
    print('(a^b)x(c^d) =', Com(a ^ b, c ^ d))
    print('(a|(b^c))|(d^e) =', (a | (b ^ c)) | (d ^ e))

    return


def derivatives_in_rectangular_coordinates():
    Print_Function()

    X = (x, y, z) = symbols('x y z')
    (ex, ey, ez, grad) = MV.setup('e_x e_y e_z', metric='[1,1,1]', coords=X)

    f = MV('f', 'scalar', fct=True)
    A = MV('A', 'vector', fct=True)
    B = MV('B', 'grade2', fct=True)
    C = MV('C', 'mv', fct=True)
    print('f =', f)
    print('A =', A)
    print('B =', B)
    print('C =', C)

    print('grad*f =', grad*f)
    print('grad|A =', grad | A)
    print('grad*A =', grad*A)

    print('-I*(grad^A) =', -MV.I*(grad ^ A))
    print('grad*B =', grad*B)
    print('grad^B =', grad ^ B)
    print('grad|B =', grad | B)

    print('grad<A =', grad < A)
    print('grad>A =', grad > A)
    print('grad<B =', grad < B)
    print('grad>B =', grad > B)
    print('grad<C =', grad < C)
    print('grad>C =', grad > C)

    return


def derivatives_in_spherical_coordinates():
    Print_Function()

    X = (r, th, phi) = symbols('r theta phi')
    curv = [[r *
             cos(phi) *
             sin(th), r *
             sin(phi) *
             sin(th), r *
             cos(th)], [1, r, r *
                        sin(th)]]
    (er,
     eth,
     ephi,
     grad) = MV.setup('e_r e_theta e_phi',
                      metric='[1,1,1]',
                      coords=X,
                      curv=curv)

    f = MV('f', 'scalar', fct=True)
    A = MV('A', 'vector', fct=True)
    B = MV('B', 'grade2', fct=True)

    print('f =', f)
    print('A =', A)
    print('B =', B)

    print('grad*f =', grad*f)
    print('grad|A =', grad | A)
    print('-I*(grad^A) =', -MV.I*(grad ^ A))
    print('grad^B =', grad ^ B)
    return


def rounding_numerical_components():
    Print_Function()

    (ex, ey, ez) = MV.setup('e_x e_y e_z', metric='[1,1,1]')

    X = 1.2*ex + 2.34*ey + 0.555*ez
    Y = 0.333*ex + 4*ey + 5.3*ez

    print('X =', X)
    print('Nga(X,2) =', Nga(X, 2))
    print('X*Y =', X*Y)
    print('Nga(X*Y,2) =', Nga(X*Y, 2))
    return


def noneuclidian_distance_calculation():
    from sympy import solve, sqrt
    Print_Function()

    metric = '0 # #,# 0 #,# # 1'
    (X, Y, e) = MV.setup('X Y e', metric)

    print('g_{ij} =', MV.metric)

    print('(X^Y)**2 =', (X ^ Y)*(X ^ Y))

    L = X ^ Y ^ e
    B = L*e  # D&L 10.152
    print('B =', B)
    Bsq = B*B
    print('B**2 =', Bsq)
    Bsq = Bsq.scalar()
    print('#L = X^Y^e is a non-euclidian line')
    print('B = L*e =', B)

    BeBr = B*e*B.rev()
    print('B*e*B.rev() =', BeBr)
    print('B**2 =', B*B)
    print('L**2 =', L*L)   # D&L 10.153
    (s, c, Binv, M, BigS, BigC, alpha, XdotY, Xdote, Ydote) = symbols(
        's c (1/B) M S C alpha (X.Y) (X.e) (Y.e)')

    Bhat = Binv*B  # D&L 10.154
    R = c + s*Bhat  # Rotor R = exp(alpha*Bhat/2)
    print('s = sinh(alpha/2) and c = cosh(alpha/2)')
    print('exp(alpha*B/(2*|B|)) =', R)

    Z = R*X*R.rev()  # D&L 10.155
    Z.obj = expand(Z.obj)
    Z.obj = Z.obj.collect([Binv, s, c, XdotY])
    Z.Fmt(3, 'R*X*R.rev()')
    W = Z | Y  # Extract scalar part of multivector
    # From this point forward all calculations are with sympy scalars
    print('Objective is to determine value of C = cosh(alpha) such that W = 0')
    W = W.scalar()
    print('Z|Y =', W)
    W = expand(W)
    W = simplify(W)
    W = W.collect([s*Binv])

    M = 1/Bsq
    W = W.subs(Binv**2, M)
    W = simplify(W)
    Bmag = sqrt(XdotY**2 - 2*XdotY*Xdote*Ydote)
    W = W.collect([Binv*c*s, XdotY])

    # Double angle substitutions

    W = W.subs(2*XdotY**2 - 4*XdotY*Xdote*Ydote, 2/(Binv**2))
    W = W.subs(2*c*s, BigS)
    W = W.subs(c**2, (BigC + 1)/2)
    W = W.subs(s**2, (BigC - 1)/2)
    W = simplify(W)
    W = W.subs(1/Binv, Bmag)
    W = expand(W)

    print('S = sinh(alpha) and C = cosh(alpha)')

    print('W =', W)

    Wd = collect(W, [BigC, BigS], exact=True, evaluate=False)

    Wd_1 = Wd[S.One]
    Wd_C = Wd[BigC]
    Wd_S = Wd[BigS]

    print('Scalar Coefficient =', Wd_1)
    print('Cosh Coefficient =', Wd_C)
    print('Sinh Coefficient =', Wd_S)

    print('|B| =', Bmag)
    Wd_1 = Wd_1.subs(Bmag, 1/Binv)
    Wd_C = Wd_C.subs(Bmag, 1/Binv)
    Wd_S = Wd_S.subs(Bmag, 1/Binv)

    lhs = Wd_1 + Wd_C*BigC
    rhs = -Wd_S*BigS
    lhs = lhs**2
    rhs = rhs**2
    W = expand(lhs - rhs)
    W = expand(W.subs(1/Binv**2, Bmag**2))
    W = expand(W.subs(BigS**2, BigC**2 - 1))
    W = W.collect([BigC, BigC**2], evaluate=False)

    a = simplify(W[BigC**2])
    b = simplify(W[BigC])
    c = simplify(W[S.One])

    print('Require a*C**2+b*C+c = 0')

    print('a =', a)
    print('b =', b)
    print('c =', c)

    x = Symbol('x')
    C = solve(a*x**2 + b*x + c, x)[0]
    print('cosh(alpha) = C = -b/(2*a) =', expand(simplify(expand(C))))
    return


def F(x):
    global n, nbar
    Fx = Rational(1, 2)*((x*x)*n + 2*x - nbar)
    return(Fx)


def make_vector(a, n=3):
    if isinstance(a, str):
        sym_str = ''
        for i in range(n):
            sym_str += a + str(i + 1) + ' '
        sym_lst = list(symbols(sym_str))
        sym_lst.append(S.Zero)
        sym_lst.append(S.Zero)
        a = MV(sym_lst, 'vector')
    return(F(a))


def conformal_representations_of_circles_lines_spheres_and_planes():
    global n, nbar
    Print_Function()

    metric = '1 0 0 0 0,0 1 0 0 0,0 0 1 0 0,0 0 0 0 2,0 0 0 2 0'

    (e1, e2, e3, n, nbar) = MV.setup('e_1 e_2 e_3 n nbar', metric)

    print('g_{ij} =\n', MV.metric)

    e = n + nbar
    # conformal representation of points

    A = make_vector(e1)    # point a = (1,0,0)  A = F(a)
    B = make_vector(e2)    # point b = (0,1,0)  B = F(b)
    C = make_vector(-e1)   # point c = (-1,0,0) C = F(c)
    D = make_vector(e3)    # point d = (0,0,1)  D = F(d)
    X = make_vector('x', 3)

    print('F(a) =', A)
    print('F(b) =', B)
    print('F(c) =', C)
    print('F(d) =', D)
    print('F(x) =', X)

    print('a = e1, b = e2, c = -e1, and d = e3')
    print('A = F(a) = 1/2*(a*a*n+2*a-nbar), etc.')
    print('Circle through a, b, and c')
    print('Circle: A^B^C^X = 0 =', (A ^ B ^ C ^ X))
    print('Line through a and b')
    print('Line  : A^B^n^X = 0 =', (A ^ B ^ n ^ X))
    print('Sphere through a, b, c, and d')
    print('Sphere: A^B^C^D^X = 0 =', (((A ^ B) ^ C) ^ D) ^ X)
    print('Plane through a, b, and d')
    print('Plane : A^B^n^D^X = 0 =', (A ^ B ^ n ^ D ^ X))

    L = (A ^ B ^ e) ^ X

    L.Fmt(3, 'Hyperbolic Circle: (A^B^e)^X = 0 =')
    return


def properties_of_geometric_objects():
    Print_Function()

    metric = '# # # 0 0,' + \
             '# # # 0 0,' + \
             '# # # 0 0,' + \
             '0 0 0 0 2,' + \
             '0 0 0 2 0'

    (p1, p2, p3, n, nbar) = MV.setup('p1 p2 p3 n nbar', metric)

    print('g_{ij} =\n', MV.metric)

    P1 = F(p1)
    P2 = F(p2)
    P3 = F(p3)

    print('Extracting direction of line from L = P1^P2^n')

    L = P1 ^ P2 ^ n
    delta = (L | n) | nbar
    print('(L|n)|nbar =', delta)

    print('Extracting plane of circle from C = P1^P2^P3')

    C = P1 ^ P2 ^ P3
    delta = ((C ^ n) | n) | nbar
    print('((C^n)|n)|nbar =', delta)
    print('(p2-p1)^(p3-p1) =', (p2 - p1) ^ (p3 - p1))


def extracting_vectors_from_conformal_2_blade():
    Print_Function()

    metric = ' 0 -1 #,' + \
             '-1  0 #,' + \
             ' #  # #,'

    (P1, P2, a) = MV.setup('P1 P2 a', metric)

    print('g_{ij} =\n', MV.metric)

    B = P1 ^ P2
    Bsq = B*B
    print('B**2 =', Bsq)
    ap = a - (a ^ B)*B
    print("a' = a-(a^B)*B =", ap)

    Ap = ap + ap*B
    Am = ap - ap*B

    print("A+ = a'+a'*B =", Ap)
    print("A- = a'-a'*B =", Am)

    print('(A+)^2 =', Ap*Ap)
    print('(A-)^2 =', Am*Am)

    aB = a | B
    print('a|B =', aB)
    return


def reciprocal_frame_test():
    Print_Function()

    metric = '1 # #,' + \
             '# 1 #,' + \
             '# # 1,'

    (e1, e2, e3) = MV.setup('e1 e2 e3', metric)

    print('g_{ij} =\n', MV.metric)

    E = e1 ^ e2 ^ e3
    Esq = (E*E).scalar()
    print('E =', E)
    print('E**2 =', Esq)
    Esq_inv = 1/Esq

    E1 = (e2 ^ e3)*E
    E2 = (-1)*(e1 ^ e3)*E
    E3 = (e1 ^ e2)*E

    print('E1 = (e2^e3)*E =', E1)
    print('E2 =-(e1^e3)*E =', E2)
    print('E3 = (e1^e2)*E =', E3)

    w = (E1 | e2)
    w = w.expand()
    print('E1|e2 =', w)

    w = (E1 | e3)
    w = w.expand()
    print('E1|e3 =', w)

    w = (E2 | e1)
    w = w.expand()
    print('E2|e1 =', w)

    w = (E2 | e3)
    w = w.expand()
    print('E2|e3 =', w)

    w = (E3 | e1)
    w = w.expand()
    print('E3|e1 =', w)

    w = (E3 | e2)
    w = w.expand()
    print('E3|e2 =', w)

    w = (E1 | e1)
    w = (w.expand()).scalar()
    Esq = expand(Esq)
    print('(E1|e1)/E**2 =', simplify(w/Esq))

    w = (E2 | e2)
    w = (w.expand()).scalar()
    print('(E2|e2)/E**2 =', simplify(w/Esq))

    w = (E3 | e3)
    w = (w.expand()).scalar()
    print('(E3|e3)/E**2 =', simplify(w/Esq))
    return


def dummy():
    return


def main():
    Get_Program(True)
    with GA_Printer():
        enhance_print()
        basic_multivector_operations()
        check_generalized_BAC_CAB_formulas()
        derivatives_in_rectangular_coordinates()
        derivatives_in_spherical_coordinates()
        rounding_numerical_components()
        noneuclidian_distance_calculation()
        conformal_representations_of_circles_lines_spheres_and_planes()
        properties_of_geometric_objects()
        extracting_vectors_from_conformal_2_blade()
        reciprocal_frame_test()
    return

if __name__ == "__main__":
    main()
