from sympy.core import S, Rational, Symbol, Basic
from sympy.combinatorics import Permutation
from sympy.combinatorics.tensor_can import (bsgs_direct_product, riemann_bsgs)
from sympy.tensor.tensor import (TensorIndexType, tensor_indices,
  TensorSymmetry, get_symmetric_group_sgs, TensorType, TensorIndex,
  tensor_mul, canon_bp, TensAdd, riemann_cyclic_replace, riemann_cyclic,
  tensorlist_contract_metric, TensMul, tensorsymmetry, tensorhead,
  TensorManager, TensExpr)
from sympy.utilities.pytest import raises

#################### Tests from tensor_can.py #######################

def test_canonicalize_no_slot_sym():
    # A_d0 * B^d0; T_c = A^d0*B_d0
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a, b, d0, d1 = tensor_indices('a,b,d0,d1', Lorentz)
    sym1 = tensorsymmetry([1])
    S1 = TensorType([Lorentz], sym1)
    A, B = S1('A,B')
    t = A(-d0)*B(d0)
    tc = t.canon_bp()
    assert str(tc) == 'A(L_0)*B(-L_0)'

    # A^a * B^b;  T_c = T
    t = A(a)*B(b)
    tc = t.canon_bp()
    assert tc == t
    # B^b * A^a
    t1 = B(b)*A(a)
    tc = t1.canon_bp()
    assert str(tc) == 'A(a)*B(b)'

    # A symmetric
    # A^{b}_{d0}*A^{d0, a}; T_c = A^{a d0}*A{b}_{d0}
    sym2 = tensorsymmetry([1]*2)
    S2 = TensorType([Lorentz]*2, sym2)
    A = S2('A')
    t = A(b, -d0)*A(d0, a)
    tc = t.canon_bp()
    assert str(tc) == 'A(a, L_0)*A(b, -L_0)'

    # A^{d1}_{d0}*B^d0*C_d1
    # T_c = A^{d0 d1}*B_d0*C_d1
    B, C = S1('B,C')
    t = A(d1, -d0)*B(d0)*C(-d1)
    tc = t.canon_bp()
    assert str(tc) == 'A(L_0, L_1)*B(-L_0)*C(-L_1)'

    # A without symmetry
    # A^{d1}_{d0}*B^d0*C_d1 ord=[d0,-d0,d1,-d1]; g = [2,1,0,3,4,5]
    # T_c = A^{d0 d1}*B_d1*C_d0; can = [0,2,3,1,4,5]
    nsym2 = tensorsymmetry([1],[1])
    NS2 = TensorType([Lorentz]*2, nsym2)
    A = NS2('A')
    B, C = S1('B, C')
    t = A(d1, -d0)*B(d0)*C(-d1)
    tc = t.canon_bp()
    assert str(tc) == 'A(L_0, L_1)*B(-L_1)*C(-L_0)'

    # A, B without symmetry
    # A^{d1}_{d0}*B_{d1}^{d0}
    # T_c = A^{d0 d1}*B_{d0 d1}
    B = NS2('B')
    t = A(d1, -d0)*B(-d1, d0)
    tc = t.canon_bp()
    assert str(tc) == 'A(L_0, L_1)*B(-L_0, -L_1)'
    # A_{d0}^{d1}*B_{d1}^{d0}
    # T_c = A^{d0 d1}*B_{d1 d0}
    t = A(-d0, d1)*B(-d1, d0)
    tc = t.canon_bp()
    assert str(tc) == 'A(L_0, L_1)*B(-L_1, -L_0)'

    # A, B, C without symmetry
    # A^{d1 d0}*B_{a d0}*C_{d1 b}
    # T_c=A^{d0 d1}*B_{a d1}*C_{d0 b}
    C = NS2('C')
    t = A(d1, d0)*B(-a, -d0)*C(-d1, -b)
    tc = t.canon_bp()
    assert str(tc) == 'A(L_0, L_1)*B(-a, -L_1)*C(-L_0, -b)'

    # A symmetric, B and C without symmetry
    # A^{d1 d0}*B_{a d0}*C_{d1 b}
    # T_c = A^{d0 d1}*B_{a d0}*C_{d1 b}
    A = S2('A')
    t = A(d1, d0)*B(-a, -d0)*C(-d1, -b)
    tc = t.canon_bp()
    assert str(tc) == 'A(L_0, L_1)*B(-a, -L_0)*C(-L_1, -b)'

    # A and C symmetric, B without symmetry
    # A^{d1 d0}*B_{a d0}*C_{d1 b} ord=[a,b,d0,-d0,d1,-d1]
    # T_c = A^{d0 d1}*B_{a d0}*C_{b d1}
    C = S2('C')
    t = A(d1, d0)*B(-a, -d0)*C(-d1, -b)
    tc = t.canon_bp()
    assert str(tc) == 'A(L_0, L_1)*B(-a, -L_0)*C(-b, -L_1)'

def test_canonicalize_no_dummies():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a, b, c, d = tensor_indices('a, b, c, d', Lorentz)
    sym1 = tensorsymmetry([1])
    sym2 = tensorsymmetry([1]*2)
    sym2a = tensorsymmetry([2])

    # A commuting
    # A^c A^b A^a
    # T_c = A^a A^b A^c
    S1 = TensorType([Lorentz], sym1)
    A = S1('A')
    t = A(c)*A(b)*A(a)
    tc = t.canon_bp()
    assert str(tc) == 'A(a)*A(b)*A(c)'

    # A anticommuting
    # A^c A^b A^a
    # T_c = -A^a A^b A^c
    A = S1('A', 1)
    t = A(c)*A(b)*A(a)
    tc = t.canon_bp()
    assert str(tc) == '-A(a)*A(b)*A(c)'

    # A commuting and symmetric
    # A^{b,d}*A^{c,a}
    # T_c = A^{a c}*A^{b d}
    S2 = TensorType([Lorentz]*2, sym2)
    A = S2('A')
    t = A(b, d)*A(c, a)
    tc = t.canon_bp()
    assert str(tc) == 'A(a, c)*A(b, d)'

    # A anticommuting and symmetric
    # A^{b,d}*A^{c,a}
    # T_c = -A^{a c}*A^{b d}
    A = S2('A', 1)
    t = A(b, d)*A(c, a)
    tc = t.canon_bp()
    assert str(tc) == '-A(a, c)*A(b, d)'

    # A^{c,a}*A^{b,d}
    # T_c = A^{a c}*A^{b d}
    t = A(c, a)*A(b, d)
    tc = t.canon_bp()
    assert str(tc) == 'A(a, c)*A(b, d)'

def test_no_metric_symmetry():
    # no metric symmetry; A no symmetry
    # A^d1_d0 * A^d0_d1
    # T_c = A^d0_d1 * A^d1_d0
    Lorentz = TensorIndexType('Lorentz', metric=None, dummy_fmt='L')
    d0, d1, d2, d3 = tensor_indices('d:4', Lorentz)
    A = tensorhead('A', [Lorentz]*2, [[1], [1]])
    t = A(d1, -d0)*A(d0, -d1)
    tc = t.canon_bp()
    assert str(tc) == 'A(L_0, -L_1)*A(L_1, -L_0)'

    # A^d1_d2 * A^d0_d3 * A^d2_d1 * A^d3_d0
    # T_c = A^d0_d1 * A^d1_d0 * A^d2_d3 * A^d3_d2
    t = A(d1, -d2)*A(d0, -d3)*A(d2,-d1)*A(d3,-d0)
    tc = t.canon_bp()
    assert str(tc) == 'A(L_0, -L_1)*A(L_1, -L_0)*A(L_2, -L_3)*A(L_3, -L_2)'

    # A^d0_d2 * A^d1_d3 * A^d3_d0 * A^d2_d1
    # T_c = A^d0_d1 * A^d1_d2 * A^d2_d3 * A^d3_d0
    t = A(d0, -d1)*A(d1, -d2)*A(d2, -d3)*A(d3,-d0)
    tc = t.canon_bp()
    assert str(tc) == 'A(L_0, -L_1)*A(L_1, -L_2)*A(L_2, -L_3)*A(L_3, -L_0)'

def test_canonicalize1():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a, a0, a1, a2, a3, b, d0, d1, d2, d3 = \
      tensor_indices('a,a0,a1,a2,a3,b,d0,d1,d2,d3', Lorentz)
    sym1 = tensorsymmetry([1])
    base3, gens3 = get_symmetric_group_sgs(3)
    sym2 = tensorsymmetry([1]*2)
    sym2a = tensorsymmetry([2])
    sym3 = tensorsymmetry([1]*3)
    sym3a = tensorsymmetry([3])

    # A_d0*A^d0; ord = [d0,-d0]
    # T_c = A^d0*A_d0
    S1 = TensorType([Lorentz], sym1)
    A = S1('A')
    t = A(-d0)*A(d0)
    tc = t.canon_bp()
    assert str(tc) == 'A(L_0)*A(-L_0)'

    # A commuting
    # A_d0*A_d1*A_d2*A^d2*A^d1*A^d0
    # T_c = A^d0*A_d0*A^d1*A_d1*A^d2*A_d2
    t = A(-d0)*A(-d1)*A(-d2)*A(d2)*A(d1)*A(d0)
    tc = t.canon_bp()
    assert str(tc) == 'A(L_0)*A(-L_0)*A(L_1)*A(-L_1)*A(L_2)*A(-L_2)'

    # A anticommuting
    # A_d0*A_d1*A_d2*A^d2*A^d1*A^d0
    # T_c 0
    A = S1('A', 1)
    t = A(-d0)*A(-d1)*A(-d2)*A(d2)*A(d1)*A(d0)
    tc = t.canon_bp()
    assert tc == 0

    # A commuting symmetric
    # A^{d0 b}*A^a_d1*A^d1_d0
    # T_c = A^{a d0}*A^{b d1}*A_{d0 d1}
    S2 = TensorType([Lorentz]*2, sym2)
    A = S2('A')
    t = A(d0, b)*A(a, -d1)*A(d1, -d0)
    tc = t.canon_bp()
    assert str(tc) == 'A(a, L_0)*A(b, L_1)*A(-L_0, -L_1)'

    # A, B commuting symmetric
    # A^{d0 b}*A^d1_d0*B^a_d1
    # T_c = A^{b d0}*A_d0^d1*B^a_d1
    B = S2('B')
    t = A(d0, b)*A(d1, -d0)*B(a, -d1)
    tc = t.canon_bp()
    assert str(tc) == 'A(b, L_0)*A(-L_0, L_1)*B(a, -L_1)'

    # A commuting symmetric
    # A^{d1 d0 b}*A^{a}_{d1 d0}; ord=[a,b, d0,-d0,d1,-d1]
    # T_c = A^{a d0 d1}*A^{b}_{d0 d1}
    S3 = TensorType([Lorentz]*3, sym3)
    A = S3('A')
    t = A(d1, d0, b)*A(a, -d1, -d0)
    tc = t.canon_bp()
    assert str(tc) == 'A(a, L_0, L_1)*A(b, -L_0, -L_1)'

    # A^{d3 d0 d2}*A^a0_{d1 d2}*A^d1_d3^a1*A^{a2 a3}_d0
    # T_c = A^{a0 d0 d1}*A^a1_d0^d2*A^{a2 a3 d3}*A_{d1 d2 d3}
    t = A(d3, d0, d2)*A(a0, -d1, -d2)*A(d1, -d3, a1)*A(a2, a3, -d0)
    tc = t.canon_bp()
    assert str(tc) == 'A(a0, L_0, L_1)*A(a1, -L_0, L_2)*A(a2, a3, L_3)*A(-L_1, -L_2, -L_3)'

    # A commuting symmetric, B antisymmetric
    # A^{d0 d1 d2} * A_{d2 d3 d1} * B_d0^d3
    # in this esxample and in the next three,
    # renaming dummy indices and using symmetry of A,
    # T = A^{d0 d1 d2} * A_{d0 d1 d3} * B_d2^d3
    # can = 0
    S2a = TensorType([Lorentz]*2, sym2a)
    A = S3('A')
    B = S2a('B')
    t = A(d0, d1, d2)*A(-d2, -d3, -d1)*B(-d0, d3)
    tc = t.canon_bp()
    assert tc == 0

    # A anticommuting symmetric, B anticommuting
    # A^{d0 d1 d2} * A_{d2 d3 d1} * B_d0^d3
    # T_c = A^{d0 d1 d2} * A_{d0 d1}^d3 * B_{d2 d3}
    A = S3('A', 1)
    B = S2a('B')
    t = A(d0, d1, d2)*A(-d2, -d3, -d1)*B(-d0, d3)
    tc = t.canon_bp()
    assert str(tc) == 'A(L_0, L_1, L_2)*A(-L_0, -L_1, L_3)*B(-L_2, -L_3)'

    # A anticommuting symmetric, B antisymmetric commuting, antisymmetric metric
    # A^{d0 d1 d2} * A_{d2 d3 d1} * B_d0^d3
    # T_c = -A^{d0 d1 d2} * A_{d0 d1}^d3 * B_{d2 d3}
    Spinor = TensorIndexType('Spinor', metric=1, dummy_fmt='S')
    a, a0, a1, a2, a3, b, d0, d1, d2, d3 = \
      tensor_indices('a,a0,a1,a2,a3,b,d0,d1,d2,d3', Spinor)
    S3 = TensorType([Spinor]*3, sym3)
    S2a = TensorType([Spinor]*2, sym2a)
    A = S3('A', 1)
    B = S2a('B')
    t = A(d0, d1, d2)*A(-d2, -d3, -d1)*B(-d0, d3)
    tc = t.canon_bp()
    assert str(tc) == '-A(S_0, S_1, S_2)*A(-S_0, -S_1, S_3)*B(-S_2, -S_3)'

    # A anticommuting symmetric, B antisymmetric anticommuting,
    # no metric symmetry
    # A^{d0 d1 d2} * A_{d2 d3 d1} * B_d0^d3
    # T_c = A^{d0 d1 d2} * A_{d0 d1 d3} * B_d2^d3
    Mat = TensorIndexType('Mat', metric=None, dummy_fmt='M')
    a, a0, a1, a2, a3, b, d0, d1, d2, d3 = \
      tensor_indices('a,a0,a1,a2,a3,b,d0,d1,d2,d3', Mat)
    S3 = TensorType([Mat]*3, sym3)
    S2a = TensorType([Mat]*2, sym2a)
    A = S3('A', 1)
    B = S2a('B')
    t = A(d0, d1, d2)*A(-d2, -d3, -d1)*B(-d0, d3)
    tc = t.canon_bp()
    assert str(tc) == 'A(M_0, M_1, M_2)*A(-M_0, -M_1, -M_3)*B(-M_2, M_3)'

    # Gamma anticommuting
    # Gamma_{mu nu} * gamma^rho * Gamma^{nu mu alpha}
    # T_c = -Gamma^{mu nu} * gamma^rho * Gamma_{alpha mu nu}
    S1 = TensorType([Lorentz], sym1)
    S2a = TensorType([Lorentz]*2, sym2a)
    S3a = TensorType([Lorentz]*3, sym3a)
    alpha, beta, gamma, mu, nu, rho = \
      tensor_indices('alpha,beta,gamma,mu,nu,rho', Lorentz)
    Gamma = S1('Gamma', 2)
    Gamma2 = S2a('Gamma', 2)
    Gamma3 = S3a('Gamma', 2)
    t = Gamma2(-mu,-nu)*Gamma(rho)*Gamma3(nu, mu, alpha)
    tc = t.canon_bp()
    assert str(tc) == '-Gamma(L_0, L_1)*Gamma(rho)*Gamma(alpha, -L_0, -L_1)'

    # Gamma_{mu nu} * Gamma^{gamma beta} * gamma_rho * Gamma^{nu mu alpha}
    # T_c = Gamma^{mu nu} * Gamma^{beta gamma} * gamma_rho * Gamma^alpha_{mu nu}
    t = Gamma2(mu, nu)*Gamma2(beta, gamma)*Gamma(-rho)*Gamma3(alpha, -mu, -nu)
    tc = t.canon_bp()
    assert str(tc) == 'Gamma(L_0, L_1)*Gamma(beta, gamma)*Gamma(-rho)*Gamma(alpha, -L_0, -L_1)'

    # f^a_{b,c} antisymmetric in b,c; A_mu^a no symmetry
    # f^c_{d a} * f_{c e b} * A_mu^d * A_nu^a * A^{nu e} * A^{mu b}
    # g = [8,11,5, 9,13,7, 1,10, 3,4, 2,12, 0,6, 14,15]
    # T_c = -f^{a b c} * f_a^{d e} * A^mu_b * A_{mu d} * A^nu_c * A_{nu e}
    Flavor = TensorIndexType('Flavor', dummy_fmt='F')
    a, b, c, d, e, ff = tensor_indices('a,b,c,d,e,f', Flavor)
    mu, nu = tensor_indices('mu,nu', Lorentz)
    sym_f = tensorsymmetry([1], [2])
    S_f = TensorType([Flavor]*3, sym_f)
    sym_A = tensorsymmetry([1], [1])
    S_A = TensorType([Lorentz, Flavor], sym_A)
    f = S_f('f')
    A = S_A('A')
    t = f(c, -d, -a)*f(-c, -e, -b)*A(-mu, d)*A(-nu, a)*A(nu, e)*A(mu, b)
    tc = t.canon_bp()
    assert str(tc) == '-f(F_0, F_1, F_2)*f(-F_0, F_3, F_4)*A(L_0, -F_1)*A(-L_0, -F_3)*A(L_1, -F_2)*A(-L_1, -F_4)'


def test_riemann_invariants():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11 = \
            tensor_indices('d0:12', Lorentz)
    # R^{d0 d1}_{d1 d0}; ord = [d0,-d0,d1,-d1]
    # T_c = -R^{d0 d1}_{d0 d1}
    R = tensorhead('R', [Lorentz]*4, [[2, 2]])
    t = R(d0, d1, -d1, -d0)
    tc = t.canon_bp()
    assert str(tc) == '-R(L_0, L_1, -L_0, -L_1)'

    # R_d11^d1_d0^d5 * R^{d6 d4 d0}_d5 * R_{d7 d2 d8 d9} *
    # R_{d10 d3 d6 d4} * R^{d2 d7 d11}_d1 * R^{d8 d9 d3 d10}
    # can = [0,2,4,6, 1,3,8,10, 5,7,12,14, 9,11,16,18, 13,15,20,22,
    #        17,19,21<F10,23, 24,25]
    # T_c = R^{d0 d1 d2 d3} * R_{d0 d1}^{d4 d5} * R_{d2 d3}^{d6 d7} *
    # R_{d4 d5}^{d8 d9} * R_{d6 d7}^{d10 d11} * R_{d8 d9 d10 d11}


    t = R(-d11,d1,-d0,d5)*R(d6,d4,d0,-d5)*R(-d7,-d2,-d8,-d9)* \
        R(-d10,-d3,-d6,-d4)*R(d2,d7,d11,-d1)*R(d8,d9,d3,d10)
    tc = t.canon_bp()
    assert str(tc) == 'R(L_0, L_1, L_2, L_3)*R(-L_0, -L_1, L_4, L_5)*R(-L_2, -L_3, L_6, L_7)*R(-L_4, -L_5, L_8, L_9)*R(-L_6, -L_7, L_10, L_11)*R(-L_8, -L_9, -L_10, -L_11)'

def test_riemann_products():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    d0, d1, d2, d3, d4, d5, d6 = tensor_indices('d0:7', Lorentz)
    a0, a1, a2, a3, a4, a5 = tensor_indices('a0:6', Lorentz)
    a, b = tensor_indices('a,b', Lorentz)
    R = tensorhead('R', [Lorentz]*4, [[2, 2]])
    # R^{a b d0}_d0 = 0
    t = R(a, b, d0, -d0)
    tc = t.canon_bp()
    assert tc == 0

    # R^{d0 b a}_d0
    # T_c = -R^{a d0 b}_d0
    t = R(d0, b, a, -d0)
    tc = t.canon_bp()
    assert str(tc) == '-R(a, L_0, b, -L_0)'

    # R^d1_d2^b_d0 * R^{d0 a}_d1^d2; ord=[a,b,d0,-d0,d1,-d1,d2,-d2]
    # T_c = -R^{a d0 d1 d2}* R^b_{d0 d1 d2}
    t = R(d1, -d2, b, -d0)*R(d0, a, -d1, d2)
    tc = t.canon_bp()
    assert str(tc) == '-R(a, L_0, L_1, L_2)*R(b, -L_0, -L_1, -L_2)'

    # A symmetric commuting
    # R^{d6 d5}_d2^d1 * R^{d4 d0 d2 d3} * A_{d6 d0} A_{d3 d1} * A_{d4 d5}
    # g = [12,10,5,2, 8,0,4,6, 13,1, 7,3, 9,11,14,15]
    # T_c = -R^{d0 d1 d2 d3} * R_d0^{d4 d5 d6} * A_{d1 d4}*A_{d2 d5}*A_{d3 d6}
    V = tensorhead('V', [Lorentz]*2, [[1]*2])
    t = R(d6, d5, -d2, d1)*R(d4, d0, d2, d3)*V(-d6, -d0)*V(-d3, -d1)*V(-d4, -d5)
    tc = t.canon_bp()
    assert str(tc) == '-R(L_0, L_1, L_2, L_3)*R(-L_0, L_4, L_5, L_6)*V(-L_1, -L_4)*V(-L_2, -L_5)*V(-L_3, -L_6)'

    # R^{d2 a0 a2 d0} * R^d1_d2^{a1 a3} * R^{a4 a5}_{d0 d1}
    # T_c = R^{a0 d0 a2 d1}*R^{a1 a3}_d0^d2*R^{a4 a5}_{d1 d2}
    t = R(d2, a0, a2, d0)*R(d1, -d2, a1, a3)*R(a4, a5, -d0, -d1)
    tc = t.canon_bp()
    assert str(tc) == 'R(a0, L_0, a2, L_1)*R(a1, a3, -L_0, L_2)*R(a4, a5, -L_1, -L_2)'
######################################################################

def test_canonicalize2():
    D = Symbol('D')
    Eucl = TensorIndexType('Eucl', metric=0, dim=D, dummy_fmt='E')
    i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14 = \
            tensor_indices('i0:15', Eucl)
    A = tensorhead('A', [Eucl]*3, [[3]])

    # two examples from Cvitanovic, Group Theory page 59
    # of identities for antisymmetric tensors of rank 3
    # contracted according to the Kuratowski graph  eq.(6.59)
    t = A(i0,i1,i2)*A(-i1,i3,i4)*A(-i3,i7,i5)*A(-i2,-i5,i6)*A(-i4,-i6,i8)
    t1 = t.canon_bp()
    assert t1 == 0

    # eq.(6.60)
    #t = A(i0,i1,i2)*A(-i1,i3,i4)*A(-i2,i5,i6)*A(-i3,i7,i8)*A(-i6,-i7,i9)*
    #    A(-i8,i10,i13)*A(-i5,-i10,i11)*A(-i4,-i11,i12)*A(-i3,-i12,i14)
    t = A(i0,i1,i2)*A(-i1,i3,i4)*A(-i2,i5,i6)*A(-i3,i7,i8)*A(-i6,-i7,i9)*\
        A(-i8,i10,i13)*A(-i5,-i10,i11)*A(-i4,-i11,i12)*A(-i9,-i12,i14)
    t1 = t.canon_bp()
    assert t1 == 0

def test_canonicalize3():
    D = Symbol('D')
    Spinor = TensorIndexType('Spinor', dim=D, metric=True, dummy_fmt='S')
    a0,a1,a2,a3,a4 = tensor_indices('a0:5', Spinor)
    C = Spinor.metric
    chi, psi = tensorhead('chi,psi', [Spinor], [[1]], 1)

    t = chi(a1)*psi(a0)
    t1 = t.canon_bp()
    assert t1 == t

    t = psi(a1)*chi(a0)
    t1 = t.canon_bp()
    assert t1 == -chi(a0)*psi(a1)

class Metric(Basic):
    def __new__(cls, name, antisym, **kwargs):
        obj = Basic.__new__(cls, name, antisym, **kwargs)
        obj.name = name
        obj.antisym = antisym
        return obj

def test_TensorIndexType():
    D = Symbol('D')
    G = Metric('g', False)
    Lorentz = TensorIndexType('Lorentz', metric=G, dim=D, dummy_fmt='L')
    m0, m1, m2, m3, m4 = tensor_indices('m0:5', Lorentz)
    sym2 = tensorsymmetry([1]*2)
    sym2n = tensorsymmetry(*get_symmetric_group_sgs(2))
    assert sym2 == sym2n
    g = Lorentz.metric
    assert str(g) == 'g(Lorentz,Lorentz)'
    assert Lorentz.eps_dim == Lorentz.dim

    TSpace = TensorIndexType('TSpace')
    i0, i1 = tensor_indices('i0 i1', TSpace)
    g = TSpace.metric
    A = tensorhead('A', [TSpace]*2, [[1]*2])
    assert  str(A(i0,-i0).canon_bp()) == 'A(TSpace_0, -TSpace_0)'

def test_indices():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a, b, c, d = tensor_indices('a,b,c,d', Lorentz)
    assert a.tensortype == Lorentz
    assert a != -a
    A, B = tensorhead('A B', [Lorentz]*2, [[1]*2])
    t = A(a,b)*B(-b,c)
    indices = t.get_indices()
    L_0 = TensorIndex('L_0', Lorentz)
    assert indices == [a, L_0, -L_0, c]
    raises(ValueError, lambda: tensor_indices(3, Lorentz))
    raises(ValueError, lambda: A(a,b,c))

def test_tensorsymmetry():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    sym = tensorsymmetry([1]*2)
    sym1 = TensorSymmetry(get_symmetric_group_sgs(2))
    assert sym == sym1
    sym = tensorsymmetry([2])
    sym1 = TensorSymmetry(get_symmetric_group_sgs(2, 1))
    assert sym == sym1
    sym2 = tensorsymmetry()
    assert sym2.base == [] and sym2.generators == [Permutation(1)]
    raises(NotImplementedError, lambda: tensorsymmetry([2, 1]))

def test_TensorType():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    sym = tensorsymmetry([1]*2)
    A = tensorhead('A', [Lorentz]*2, [[1]*2])
    assert A.typ == TensorType([Lorentz]*2, sym)
    assert A.types == [Lorentz]
    typ =  TensorType([Lorentz]*2, sym)
    assert str(typ) == "TensorType(['Lorentz', 'Lorentz'])"
    raises(ValueError, lambda: typ(2))

def test_TensExpr():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a, b, c, d = tensor_indices('a,b,c,d', Lorentz)
    g = Lorentz.metric
    A, B = tensorhead('A B', [Lorentz]*2, [[1]*2])
    raises(ValueError, lambda: g(c, d)/g(a, b))
    raises(ValueError, lambda: S.One/g(a, b))
    raises(ValueError, lambda: (A(c, d) + g(c, d))/g(a, b))
    raises(ValueError, lambda: S.One/(A(c, d) + g(c, d)))
    raises(ValueError, lambda: A(a, b) + A(a, c))
    t = A(a, b) + B(a, b)
    raises(NotImplementedError, lambda: TensExpr.__mul__(t, 'a'))
    raises(NotImplementedError, lambda: TensExpr.__add__(t, 'a'))
    raises(NotImplementedError, lambda: TensExpr.__radd__(t, 'a'))
    raises(NotImplementedError, lambda: TensExpr.__sub__(t, 'a'))
    raises(NotImplementedError, lambda: TensExpr.__rsub__(t, 'a'))
    raises(NotImplementedError, lambda: TensExpr.__div__(t, 'a'))
    raises(NotImplementedError, lambda: TensExpr.__rdiv__(t, 'a'))
    raises(NotImplementedError, lambda: A(a, b)**2)
    raises(NotImplementedError, lambda: 2**A(a, b))
    raises(NotImplementedError, lambda: abs(A(a, b)))

def test_TensorHead():
    assert TensAdd() == 0
    # simple example of algebraic expression
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a,b = tensor_indices('a,b', Lorentz)
    # A, B symmetric
    A = tensorhead('A', [Lorentz]*2, [[1]*2])
    assert A.rank == 2
    assert A.symmetry == tensorsymmetry([1]*2)

def test_add1():
    assert TensAdd() == 0
    # simple example of algebraic expression
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a,b,d0,d1,i,j,k = tensor_indices('a,b,d0,d1,i,j,k', Lorentz)
    # A, B symmetric
    A, B = tensorhead('A,B', [Lorentz]*2, [[1]*2])
    t1 = A(b,-d0)*B(d0,a)
    assert TensAdd(t1).equals(t1)
    t2a = B(d0,a) + A(d0, a)
    t2 = A(b,-d0)*t2a
    assert str(t2) == 'A(a, L_0)*A(b, -L_0) + A(b, L_0)*B(a, -L_0)'
    t2b = t2 + t1
    assert str(t2b) == '2*A(b, L_0)*B(a, -L_0) + A(a, L_0)*A(b, -L_0)'
    p, q, r = tensorhead('p,q,r', [Lorentz], [[1]])
    t = q(d0)*2
    assert str(t) == '2*q(d0)'
    t = 2*q(d0)
    assert str(t) == '2*q(d0)'
    t1 = p(d0) + 2*q(d0)
    assert str(t1) == '2*q(d0) + p(d0)'
    t2 = p(-d0) + 2*q(-d0)
    assert str(t2) == '2*q(-d0) + p(-d0)'
    t1 = p(d0)
    t3 = t1*t2
    assert str(t3) == '2*p(L_0)*q(-L_0) + p(L_0)*p(-L_0)'
    t3 = t2*t1
    assert str(t3) == '2*p(L_0)*q(-L_0) + p(L_0)*p(-L_0)'
    t1 = p(d0) + 2*q(d0)
    t3 = t1*t2
    assert str(t3) == '4*p(L_0)*q(-L_0) + 4*q(L_0)*q(-L_0) + p(L_0)*p(-L_0)'
    t1 =  p(d0) - 2*q(d0)
    assert str(t1) == '-2*q(d0) + p(d0)'
    t2 = p(-d0) + 2*q(-d0)
    t3 = t1*t2
    assert t3 == p(d0)*p(-d0) - 4*q(d0)*q(-d0)
    t = p(i)*p(j)*(p(k) + q(k)) + p(i)*(p(j) + q(j))*(p(k) - 3*q(k))
    assert t == 2*p(i)*p(j)*p(k) - 2*p(i)*p(j)*q(k) + p(i)*p(k)*q(j) - 3*p(i)*q(j)*q(k)
    t1 = (p(i) + q(i) + 2*r(i))*(p(j) - q(j))
    t2 = (p(j) + q(j) + 2*r(j))*(p(i) - q(i))
    t = t1 + t2
    assert t == 2*p(i)*p(j) + 2*p(i)*r(j) + 2*p(j)*r(i) - 2*q(i)*q(j) - 2*q(i)*r(j) - 2*q(j)*r(i)
    t = p(i)*q(j)/2
    assert 2*t == p(i)*q(j)
    t = (p(i) + q(i))/2
    assert 2*t == p(i) + q(i)

    t = S.One - p(i)*p(-i)
    assert (t + p(-j)*p(j)).equals(1)
    t = S.One + p(i)*p(-i)
    assert (t - p(-j)*p(j)).equals(1)

    t = A(a, b) + B(a, b)
    assert t.rank == 2
    t1 = t - A(a, b) - B(a, b)
    assert t1 == 0
    t = 1 - (A(a, -a) + B(a, -a))
    t1 = 1 + (A(a, -a) + B(a, -a))
    assert (t + t1).equals(2)
    t2 = 1 + A(a, -a)
    assert t1 != t2
    assert t2 != TensMul(0, [],[],[])
    t = p(i) + q(i)
    raises(ValueError, lambda: t(i, j))

def test_add2():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    m, n, p, q = tensor_indices('m,n,p,q', Lorentz)
    R = tensorhead('R', [Lorentz]*4, [[2, 2]])
    A = tensorhead('A', [Lorentz]*3, [[3]])
    t1 = 2*R(m,n,p,q) - R(m,q,n,p) + R(m,p,n,q)
    t2 = t1*A(-n,-p,-q)
    assert t2 == 0
    t1 = S(2)/3*R(m,n,p,q) - S(1)/3*R(m,q,n,p) + S(1)/3*R(m,p,n,q)
    t2 = t1*A(-n,-p,-q)
    assert t2 == 0
    t = A(m, -m, n) + A(n, p, -p)
    assert t == 0

def test_mul():
    from sympy.abc import x
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a, b, c, d = tensor_indices('a,b,c,d', Lorentz)
    sym = tensorsymmetry([1]*2)
    t = TensMul(S.One, [],[],[])
    assert str(t) == '1'
    A, B = tensorhead('A B', [Lorentz]*2, [[1]*2])
    t = (1 + x)*A(a, b)
    assert str(t) == '(x + 1)*A(a, b)'
    assert t.types == [Lorentz]
    assert t.rank == 2
    assert t.dum == []
    assert t.coeff == 1 + x
    assert sorted(t.free) == [(a, 0, 0), (b, 1, 0)]
    assert t.components == [A]

    t = A(-b, a)*B(-a, c)*A(-c, d)
    t1 = tensor_mul(*t.split())
    assert t == t(-b, d)
    assert t == t1
    assert tensor_mul(*[]) == TensMul(S.One, [],[],[])

    t = TensMul(1, [], [], [])
    zsym = tensorsymmetry()
    typ =  TensorType([], zsym)
    C = typ('C')
    assert str(C()) == 'C'
    assert str(t) == '1'
    assert t.split()[0] == t
    raises(ValueError, lambda: TensMul.from_indices(a, a))
    raises(ValueError, lambda: TensMul.from_indices(-a, -a))
    raises(ValueError, lambda: A(a, b)*A(a, c))
    t = A(a, b)*A(-a, c)
    raises(ValueError, lambda: t(a, b, c))

def test_substitute_indices():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    i, j, k, l, m, n, p, q = tensor_indices('i,j,k,l,m,n,p,q', Lorentz)
    A, B = tensorhead('A,B', [Lorentz]*2, [[1]*2])
    t = A(i, k)*B(-k, -j)
    t1 = t.substitute_indices((i,j), (j, k))
    t1a = A(j, l)*B(-l, -k)
    assert t1 == t1a

    p = tensorhead('p', [Lorentz], [[1]])
    t = p(i)
    t1 = t.substitute_indices((j, k))
    assert t1 == t
    t1 = t.substitute_indices((i, j))
    assert t1 == p(j)
    t1 = t.substitute_indices((i, -j))
    assert t1 == p(-j)
    t1 = t.substitute_indices((-i, j))
    assert t1 == p(-j)
    t1 = t.substitute_indices((-i, -j))
    assert t1 == p(j)

def test_riemann_cyclic_replace():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    m0,m1,m2,m3 = tensor_indices('m:4', Lorentz)
    symr = tensorsymmetry([2, 2])
    R = tensorhead('R', [Lorentz]*4, [[2, 2]])
    t = R(m0,m2,m1,m3)
    t1 = riemann_cyclic_replace(t)
    t1a =  -S.One/3*R(m0, m3, m2, m1) + S.One/3*R(m0, m1, m2, m3) + Rational(2,3)*R(m0, m2, m1, m3)
    assert t1 == t1a

def test_riemann_cyclic():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    i, j, k, l, m, n, p, q = tensor_indices('i,j,k,l,m,n,p,q', Lorentz)
    R = tensorhead('R', [Lorentz]*4, [[2, 2]])
    t = R(i,j,k,l) + R(i,l,j,k) + R(i,k,l,j) - \
        R(i,j,l,k) - R(i,l,k,j) - R(i,k,j,l)
    t2 = t*R(-i,-j,-k,-l)
    t3 = riemann_cyclic(t2)
    assert t3 == 0
    t = R(i,j,k,l)*(R(-i,-j,-k,-l) - 2*R(-i,-k,-j,-l))
    t1 = riemann_cyclic(t)
    assert t1 == 0
    t = R(i,j,k,l)
    t1 = riemann_cyclic(t)
    assert t1 == -S(1)/3*R(i, l, j, k) + S(1)/3*R(i, k, j, l) + S(2)/3*R(i, j, k, l)

    t = R(i,j,k,l)*R(-k,-l,m,n)*(R(-m,-n,-i,-j) + 2*R(-m,-j,-n,-i))
    t1 = riemann_cyclic(t)
    assert t1 == 0

def test_div():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    m0,m1,m2,m3 = tensor_indices('m0:4', Lorentz)
    R = tensorhead('R', [Lorentz]*4, [[2, 2]])
    t = R(m0,m1,-m1,m3)
    t1 = t/S(4)
    assert str(t1) == '1/4*R(m0, L_0, -L_0, m3)'
    t = t.canon_bp()
    assert not t1._is_canon_bp
    t1 = t*4
    assert t1._is_canon_bp
    t1 = t1/4
    assert t1._is_canon_bp

def test_contract_metric1():
    D = Symbol('D')
    Lorentz = TensorIndexType('Lorentz', dim=D, dummy_fmt='L')
    a, b, c, d, e = tensor_indices('a,b,c,d,e', Lorentz)
    g = Lorentz.metric
    p = tensorhead('p', [Lorentz], [[1]])
    t = g(a, b)*p(-b)
    t1 = t.contract_metric(g)
    assert t1 == p(a)
    A, B = tensorhead('A,B', [Lorentz]*2, [[1]*2])

    # case with g with all free indices
    t1 = A(a,b)*B(-b,c)*g(d, e)
    t2 = t1.contract_metric(g)
    assert t1 == t2

    # case of g(d, -d)
    t1 = A(a,b)*B(-b,c)*g(-d, d)
    t2 = t1.contract_metric(g)
    assert t2 == D*A(a, d)*B(-d, c)

    # g with one free index
    t1 = A(a,b)*B(-b,-c)*g(c, d)
    t2 = t1.contract_metric(g)
    assert t2 == A(a, c)*B(-c, d)

    # g with both indices contracted with another tensor
    t1 = A(a,b)*B(-b,-c)*g(c, -a)
    t2 = t1.contract_metric(g)
    assert t2 == A(a, b)*B(-b, -a)

    t1 = A(a,b)*B(-b,-c)*g(c, d)*g(-a, -d)
    t2 = t1.contract_metric(g)
    t2 = t2.contract_metric(g)
    assert t2 == A(a,b)*B(-b,-a)

    t1 = A(a,b)*g(-a,-b)
    t2 = t1.contract_metric(g)
    assert t2 == A(a, -a)
    assert not t2._free
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a, b = tensor_indices('a,b', Lorentz)
    g = Lorentz.metric
    raises(ValueError, lambda: g(a, -a).contract_metric(g)) # no dim

def test_contract_metric2():
    D = Symbol('D')
    Lorentz = TensorIndexType('Lorentz', dim=D, dummy_fmt='L')
    a, b, c, d, e, L_0 = tensor_indices('a,b,c,d,e,L_0', Lorentz)
    g = Lorentz.metric
    p, q = tensorhead('p,q', [Lorentz], [[1]])

    t1 = g(a,b)*p(c)*p(-c)
    t2 = 3*g(-a,-b)*q(c)*q(-c)
    t = t1*t2
    t = t.contract_metric(g)
    t = t.contract_metric(g)
    assert t == 3*D*p(a)*p(-a)*q(b)*q(-b)
    t1 = g(a,b)*p(c)*p(-c)
    t2 = 3*q(-a)*q(-b)
    t = t1*t2
    t = t.contract_metric(g)
    t = t.canon_bp()
    assert t == 3*p(a)*p(-a)*q(b)*q(-b)

    t1 = 2*g(a,b)*p(c)*p(-c)
    t2 = - 3*g(-a,-b)*q(c)*q(-c)
    t = t1*t2
    t = t.contract_metric(g)
    t = t.contract_metric(g)
    t = 6*g(a,b)*g(-a,-b)*p(c)*p(-c)*q(d)*q(-d)
    t = t.contract_metric(g)
    t = t.contract_metric(g)

    t1 = 2*g(a,b)*p(c)*p(-c)
    t2 = q(-a)*q(-b) + 3*g(-a,-b)*q(c)*q(-c)
    t = t1*t2
    t = t.contract_metric(g)
    assert t == (2 + 6*D)*p(a)*p(-a)*q(b)*q(-b)

    t1 = p(a)*p(b) + p(a)*q(b) + 2*g(a,b)*p(c)*p(-c)
    t2 = q(-a)*q(-b) - g(-a,-b)*q(c)*q(-c)
    t = t1*t2
    t = t.contract_metric(g)
    t1 = (1 - 2*D)*p(a)*p(-a)*q(b)*q(-b) + p(a)*q(-a)*p(b)*q(-b)
    assert t == t1

    t = g(a,b)*g(c,d)*g(-b,-c)
    t1 = t.contract_metric(g)
    t1 = t1.canon_bp()
    assert t1 == g(a,c)*g(d,-c)

    t2 = t1.contract_metric(g)
    assert t2 == g(a, d)

    t1 = t.contract_metric(g, True)
    assert t1 == g(a, d)

    t1 = g(a,b)*g(c,d) + g(a,c)*g(b,d) + g(a,d)*g(b,c)
    t2 = t1.substitute_indices((a,-a),(b,-b),(c,-c),(d,-d))
    t = t1*t2
    t3 = t.contract_metric(g)
    t3 = t3.contract_metric(g)
    t3 = t3.contract_metric(g)
    assert t3.equals(3*D**2 + 6*D)
    t = t.contract_metric(g, True)
    assert t3.equals(3*D**2 + 6*D)

    t = 2*p(a)*g(b,-b)
    t1 = t.contract_metric(g)
    assert t1.equals(2*D*p(a))

    t = 2*p(a)*g(b,-a)
    t1 = t.contract_metric(g)
    assert t1 == 2*p(b)

    M = Symbol('M')
    t = (p(a)*p(b) + g(a, b)*M**2)*g(-a, -b) - D*M**2
    t1 = t.contract_metric(g)
    assert t1 == p(a)*p(-a)

    v = [p(a), q(b), p(c)]
    v1 = tensorlist_contract_metric(v, g(-a, -b))
    assert v1 == [p(-b), q(b), p(c)]
    v = [p(a), q(b), p(c)]
    v1 = tensorlist_contract_metric(v, g(d, e))
    assert v1 == [p(a), q(b), p(c), g(d, e)]

    A = tensorhead('A', [Lorentz]*2, [[1]*2])
    t = A(a, b)*p(L_0)*g(-a, -b)
    t1 = t.contract_metric(g)
    assert str(t1) == 'A(L_1, -L_1)*p(L_0)' or str(t1) == 'A(-L_1, L_1)*p(L_0)'

def test_metric_contract3():
    D = Symbol('D')
    Spinor = TensorIndexType('Spinor', dim=D, metric=True, dummy_fmt='S')
    a0,a1,a2,a3,a4 = tensor_indices('a0:5', Spinor)
    C = Spinor.metric
    chi, psi = tensorhead('chi,psi', [Spinor], [[1]], 1)
    B = tensorhead('B', [Spinor]*2, [[1],[1]])

    t = C(a0, -a0)
    t1 = t.contract_metric(C)
    assert t1.equals(-D)

    t = C(-a0, a0)
    t1 = t.contract_metric(C)
    assert t1.equals(D)

    t = C(a0,a1)*C(-a0,-a1)
    t1 = t.contract_metric(C)
    assert t1.equals(D)

    t = C(a1,a0)*C(-a0,-a1)
    t1 = t.contract_metric(C)
    assert t1.equals(-D)

    t = C(-a0,a1)*C(a0,-a1)
    t1 = t.contract_metric(C)
    assert t1.equals(-D)

    t = C(a1,-a0)*C(a0,-a1)
    t1 = t.contract_metric(C)
    assert t1.equals(D)

    t = C(a0,a1)*B(-a1,-a0)
    t1 = t.contract_metric(C)
    assert t1 == B(a0, -a0)

    t = C(a1,a0)*B(-a1,-a0)
    t1 = t.contract_metric(C)
    assert t1 == -B(a0, -a0)

    t = C(a0,-a1)*B(a1,-a0)
    t1 = t.contract_metric(C)
    assert t1 == -B(a0, -a0)

    t = C(-a0,a1)*B(-a1,a0)
    t1 = t.contract_metric(C)
    assert t1 == -B(a0, -a0)

    t = C(-a0,-a1)*B(a1,a0)
    t1 = t.contract_metric(C)
    assert t1 == B(a0, -a0)

    t = C(-a1, a0)*B(a1,-a0)
    t1 = t.contract_metric(C)
    assert t1 == B(a0, -a0)

    t = C(a0,a1)*psi(-a1)
    t1 = t.contract_metric(C)
    assert t1 == psi(a0)

    t = C(a1,a0)*psi(-a1)
    t1 = t.contract_metric(C)
    assert t1 == -psi(a0)

    t = C(a0,a1)*chi(-a0)*psi(-a1)
    t1 = t.contract_metric(C)
    assert t1 == -chi(a1)*psi(-a1)

    t = C(a1,a0)*chi(-a0)*psi(-a1)
    t1 = t.contract_metric(C)
    assert t1 == chi(a1)*psi(-a1)

    t = C(-a1,a0)*chi(-a0)*psi(a1)
    t1 = t.contract_metric(C)
    assert t1 == chi(-a1)*psi(a1)

    t = C(a0, -a1)*chi(-a0)*psi(a1)
    t1 = t.contract_metric(C)
    assert t1 == -chi(-a1)*psi(a1)

    t = C(-a0,-a1)*chi(a0)*psi(a1)
    t1 = t.contract_metric(C)
    assert t1 == chi(-a1)*psi(a1)

    t = C(-a1,-a0)*chi(a0)*psi(a1)
    t1 = t.contract_metric(C)
    assert t1 == -chi(-a1)*psi(a1)

    t = C(-a1,-a0)*B(a0,a2)*psi(a1)
    t1 = t.contract_metric(C)
    assert t1 == -B(-a1,a2)*psi(a1)

    t = C(a1,a0)*B(-a2,-a0)*psi(-a1)
    t1 = t.contract_metric(C)
    assert t1 == B(-a2,a1)*psi(-a1)

def test_epsilon():
    Lorentz = TensorIndexType('Lorentz', dim=4, dummy_fmt='L')
    a, b, c, d, e = tensor_indices('a,b,c,d,e', Lorentz)
    g = Lorentz.metric
    epsilon = Lorentz.epsilon
    p, q, r, s = tensorhead('p,q,r,s', [Lorentz], [[1]])

    t = epsilon(b,a,c,d)
    t1 = t.canon_bp()
    assert t1 == -epsilon(a,b,c,d)

    t = epsilon(c,b,d,a)
    t1 = t.canon_bp()
    assert t1 == epsilon(a,b,c,d)

    t = epsilon(c,a,d,b)
    t1 = t.canon_bp()
    assert t1 == -epsilon(a,b,c,d)

    t = epsilon(a,b,c,d)*p(-a)*q(-b)
    t1 = t.canon_bp()
    assert t1 == epsilon(c, d, a, b)*p(-a)*q(-b)

    t = epsilon(c,b,d,a)*p(-a)*q(-b)
    t1 = t.canon_bp()
    assert t1 == epsilon(c, d, a, b)*p(-a)*q(-b)

    t = epsilon(c,a,d,b)*p(-a)*q(-b)
    t1 = t.canon_bp()
    assert t1 == -epsilon(c, d, a, b)*p(-a)*q(-b)

    t = epsilon(c,a,d,b)*p(-a)*p(-b)
    t1 = t.canon_bp()
    assert t1 == 0

    t = epsilon(c,a,d,b)*p(-a)*q(-b) + epsilon(a,b,c,d)*p(-b)*q(-a)
    t1 = t.canon_bp()
    assert t1 == -2*epsilon(c, d, a, b)*p(-a)*q(-b)

def test_contract_delta1():
    # see Group Theory by Cvitanovic page 9
    n = Symbol('n')
    Color = TensorIndexType('Color', metric=None, dim=n, dummy_fmt='C')
    a, b, c, d, e, f = tensor_indices('a,b,c,d,e,f', Color)
    delta = Color.delta

    def idn(a, b, d, c):
        assert a.is_up and d.is_up
        assert not (b.is_up or c.is_up)
        return delta(a, c)*delta(d, b)

    def T(a, b, d, c):
        assert a.is_up and d.is_up
        assert not (b.is_up or c.is_up)
        return delta(a, b)*delta(d, c)

    def P1(a, b, c, d):
        return idn(a,b,c,d) - 1/n*T(a,b,c,d)

    def P2(a, b, c, d):
        return 1/n*T(a,b,c,d)

    t = P1(a, -b, e, -f)*P1(f, -e, d, -c)
    t1 = t.contract_delta(delta)
    assert t1 == P1(a, -b, d, -c)

    t = P2(a, -b, e, -f)*P2(f, -e, d, -c)
    t1 = t.contract_delta(delta)
    assert t1 == P2(a, -b, d, -c)

    t = P1(a, -b, e, -f)*P2(f, -e, d, -c)
    t1 = t.contract_delta(delta)
    assert t1 == 0

    t = P1(a, -b, b, -a)
    t1 = t.contract_delta(delta)
    assert t1.equals(n**2 - 1)

def test_fun():
    D = Symbol('D')
    Lorentz = TensorIndexType('Lorentz', dim=D, dummy_fmt='L')
    a,b,c,d,e = tensor_indices('a,b,c,d,e', Lorentz)
    g = Lorentz.metric

    p, q = tensorhead('p q', [Lorentz], [[1]])
    t = q(c)*p(a)*q(b) + g(a,b)*g(c,d)*q(-d)
    assert t(a,b,c) == t
    assert t - t(b,a,c) == q(c)*p(a)*q(b) - q(c)*p(b)*q(a)
    assert t(b,c,d) == q(d)*p(b)*q(c) + g(b,c)*g(d,e)*q(-e)
    t1 = t.fun_eval((a,b),(b,a))
    assert t1 == q(c)*p(b)*q(a) + g(a,b)*g(c,d)*q(-d)

    # check that g_{a b; c} = 0
    # example taken from  L. Brewin
    # "A brief introduction to Cadabra" arxiv:0903.2085
    # dg_{a b c} = \partial_{a} g_{b c} is symmetric in b, c
    dg = tensorhead('dg', [Lorentz]*3, [[1], [1]*2])
    # gamma^a_{b c} is the Christoffel symbol
    gamma = S.Half*g(a,d)*(dg(-b,-d,-c) + dg(-c,-b,-d) - dg(-d,-b,-c))
    # t = g_{a b; c}
    t = dg(-c,-a,-b) - g(-a,-d)*gamma(d,-b,-c) - g(-b,-d)*gamma(d,-a,-c)
    t = t.contract_metric(g, True)
    assert t == 0
    t = q(c)*p(a)*q(b)
    assert t(b,c,d) == q(d)*p(b)*q(c)

def test_TensorManager():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    LorentzH = TensorIndexType('LorentzH', dummy_fmt='LH')
    i, j = tensor_indices('i,j', Lorentz)
    ih, jh = tensor_indices('ih,jh', LorentzH)
    p, q = tensorhead('p q', [Lorentz], [[1]])
    ph, qh = tensorhead('ph qh', [LorentzH], [[1]])

    Gsymbol = Symbol('Gsymbol')
    GHsymbol = Symbol('GHsymbol')
    TensorManager.set_comm(Gsymbol, GHsymbol, 0)
    G = tensorhead('G', [Lorentz], [[1]], Gsymbol)
    assert TensorManager._comm_i2symbol[G.comm] == Gsymbol
    GH = tensorhead('GH', [LorentzH], [[1]], GHsymbol)
    ps = G(i)*p(-i)
    psh = GH(ih)*ph(-ih)
    t = ps + psh
    t1 = t*t
    assert t1 == ps*ps + 2*ps*psh + psh*psh
    qs = G(i)*q(-i)
    qsh = GH(ih)*qh(-ih)
    assert ps*qsh == qsh*ps
    assert ps*qs != qs*ps
    n = TensorManager.comm_symbols2i(Gsymbol)
    assert TensorManager.comm_i2symbol(n) == Gsymbol

    assert GHsymbol in TensorManager._comm_symbols2i
    raises(ValueError, lambda: TensorManager.set_comm(GHsymbol, 1, 2))
    TensorManager.set_comms((Gsymbol,GHsymbol,0),(Gsymbol,1,1))
    assert TensorManager.get_comm(n, 1) == TensorManager.get_comm(1, n) == 1
    TensorManager.clear()
    assert TensorManager.comm == [{0:0, 1:0, 2:0}, {0:0, 1:1, 2:None}, {0:0, 1:None}]
    assert GHsymbol not in TensorManager._comm_symbols2i
    nh = TensorManager.comm_symbols2i(GHsymbol)
    assert GHsymbol in TensorManager._comm_symbols2i


def test_hash():
    D = Symbol('D')
    Lorentz = TensorIndexType('Lorentz', dim=D, dummy_fmt='L')
    a,b,c,d,e = tensor_indices('a,b,c,d,e', Lorentz)
    g = Lorentz.metric

    p, q = tensorhead('p q', [Lorentz], [[1]])
    t1 = p(a)*q(b)
    t2 = p(a)*p(b)
    assert hash(t1) != hash(t2)
    t3 = p(a)*p(b) + g(a,b)
    t4 = p(a)*p(b) - g(a,b)
    assert hash(t3) != hash(t4)
