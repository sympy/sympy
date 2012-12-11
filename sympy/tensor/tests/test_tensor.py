from sympy.core import S, Rational
from sympy.combinatorics import Permutation
from sympy.combinatorics.tensor_can import (bsgs_direct_product, riemann_bsgs)
from sympy.tensor.tensor import (TensorIndexType, tensor_indices,
  TensorSymmetry, get_symmetric_group_sgs, TensorType, TensorIndex,
  tensor_mul, canon_bp, TensAdd, riemann_cyclic_replaceR, riemann_cyclic)


def test_get_indices():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a, b, c, d = tensor_indices('a,b,c,d', Lorentz)
    sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
    S2 = TensorType([Lorentz]*2, sym2)
    A, B = S2('A,B')
    t = A(a,b)*B(-b,c)
    indices = t.get_indices()
    L_0 = TensorIndex('L_0', Lorentz)
    assert indices == [a, L_0, -L_0, c]
    a = t.split()
    t2 = tensor_mul(*a)
    assert t == t2

def test_canonicalize_no_slot_sym():
    # A_d0 * B^d0; T_c = A^d0*B_d0
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a, b, d0, d1 = tensor_indices('a,b,d0,d1', Lorentz)
    sym1 = TensorSymmetry(get_symmetric_group_sgs(1))
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
    sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
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
    nsym2 = TensorSymmetry(([], [Permutation(range(4))]))
    NS2 = TensorType([Lorentz]*2, nsym2)
    A = NS2('A')
    B, C = S1('B,C')
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
    a, b, c, d = tensor_indices('a,b,c,d', Lorentz)
    sym1 = TensorSymmetry(get_symmetric_group_sgs(1))
    sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
    sym2a = TensorSymmetry(get_symmetric_group_sgs(2, 1))

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
    Lorentz = TensorIndexType('Lorentz', metric_sym=None, dummy_fmt='L')
    d0, d1, d2, d3 = tensor_indices('d0,d1,d2,d3', Lorentz)
    nsym2 = TensorSymmetry(([], [Permutation(range(4))]))
    NS2 = TensorType([Lorentz]*2, nsym2)
    A = NS2('A')
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
    sym1 = TensorSymmetry(get_symmetric_group_sgs(1))
    base3, gens3 = get_symmetric_group_sgs(3)
    sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
    sym2a = TensorSymmetry(get_symmetric_group_sgs(2, 1))
    sym3 = TensorSymmetry(get_symmetric_group_sgs(3))
    sym3a = TensorSymmetry(get_symmetric_group_sgs(3, 1))

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
    Spinor = TensorIndexType('Spinor', metric_sym=1, dummy_fmt='S')
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
    Mat = TensorIndexType('Mat', metric_sym=None, dummy_fmt='M')
    a, a0, a1, a2, a3, b, d0, d1, d2, d3 = \
      tensor_indices('a,a0,a1,a2,a3,b,d0,d1,d2,d3', Spinor)
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
    Gamma = S1('Gamma', None)
    Gamma2 = S2a('Gamma', None)
    Gamma3 = S3a('Gamma', None)
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
    sym_f = TensorSymmetry(bsgs_direct_product(sym1.base, sym1.generators,
            sym2a.base, sym2a.generators))
    S_f = TensorType([Flavor]*3, sym_f)
    sym_A = TensorSymmetry(bsgs_direct_product(sym1.base, sym1.generators, sym1.base, sym1.generators))
    S_A = TensorType([Lorentz, Flavor], sym_A)
    f = S_f('f')
    A = S_A('A')
    t = f(c, -d, -a)*f(-c, -e, -b)*A(-mu, d)*A(-nu, a)*A(nu, e)*A(mu, b)
    tc = t.canon_bp()
    assert str(tc) == '-f(F_0, F_1, F_2)*f(-F_0, F_3, F_4)*A(L_0, -F_1)*A(-L_0, -F_3)*A(L_1, -F_2)*A(-L_1, -F_4)'


def test_riemann_invariants():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11 = \
      tensor_indices(','.join(['d%d' % i for i in range(12)]), Lorentz)
    # R^{d0 d1}_{d1 d0}; ord = [d0,-d0,d1,-d1]
    # T_c = -R^{d0 d1}_{d0 d1}
    symr = TensorSymmetry(riemann_bsgs)
    R4 = TensorType([Lorentz]*4, symr)
    R = R4('R', 0)
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
    d0, d1, d2, d3, d4, d5, d6= \
      tensor_indices(','.join(['d%d' % i for i in range(7)]), Lorentz)
    a0, a1, a2, a3, a4, a5 = \
      tensor_indices(','.join(['a%d' % i for i in range(6)]), Lorentz)
    a, b = tensor_indices('a,b', Lorentz)
    symr = TensorSymmetry(riemann_bsgs)
    R4 = TensorType([Lorentz]*4, symr)
    R = R4('R')
    sym1 = TensorSymmetry(get_symmetric_group_sgs(1))
    sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
    sym2a = TensorSymmetry(get_symmetric_group_sgs(2, 1))
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
    S2 = TensorType([Lorentz]*2, sym2)
    V = S2('V')
    t = R(d6, d5, -d2, d1)*R(d4, d0, d2, d3)*V(-d6, -d0)*V(-d3, -d1)*V(-d4, -d5)
    tc = t.canon_bp()
    assert str(tc) == '-R(L_0, L_1, L_2, L_3)*R(-L_0, L_4, L_5, L_6)*V(-L_1, -L_4)*V(-L_2, -L_5)*V(-L_3, -L_6)'

    # R^{d2 a0 a2 d0} * R^d1_d2^{a1 a3} * R^{a4 a5}_{d0 d1}
    # T_c = R^{a0 d0 a2 d1}*R^{a1 a3}_d0^d2*R^{a4 a5}_{d1 d2}
    t = R(d2, a0, a2, d0)*R(d1, -d2, a1, a3)*R(a4, a5, -d0, -d1)
    tc = t.canon_bp()
    assert str(tc) == 'R(a0, L_0, a2, L_1)*R(a1, a3, -L_0, L_2)*R(a4, a5, -L_1, -L_2)'


def test_add1():
    # simple example of algebraic expression
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    a,b,d0,d1,i,j,k = tensor_indices('a,b,d0,d1,i,j,k', Lorentz)
    # A, B symmetric
    sym1 = TensorSymmetry(get_symmetric_group_sgs(1))
    sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
    S2 = TensorType([Lorentz]*2, sym2)
    A, B = S2('A,B')
    t1 = A(b,-d0)*B(d0,a)
    t2a = B(d0,a) + A(d0, a)
    t2 = A(b,-d0)*t2a
    assert str(t2) == 'A(a, L_0)*A(b, -L_0) + A(b, L_0)*B(a, -L_0)'
    t2b = t2 + t1
    assert str(t2b) == 'A(a, L_0)*A(b, -L_0) + A(b, L_0)*B(a, -L_0) + A(b, L_0)*B(a, -L_0)'
    S1 = TensorType([Lorentz], sym1)
    p, q, r = S1('p,q,r')
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

def test_add2():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    m, n, p, q = tensor_indices('m,n,p,q', Lorentz)
    symr = TensorSymmetry(riemann_bsgs)
    sym3a = TensorSymmetry(get_symmetric_group_sgs(3, 1))
    R4 = TensorType([Lorentz]*4, symr)
    S3a = TensorType([Lorentz]*3, sym3a)
    R = R4('R')
    A = S3a('A')
    t1 = 2*R(m,n,p,q) - R(m,q,n,p) + R(m,p,n,q)
    t2 = t1*A(-n,-p,-q)
    assert t2 == 0
    t1 = S(2)/3*R(m,n,p,q) - S(1)/3*R(m,q,n,p) + S(1)/3*R(m,p,n,q)
    t2 = t1*A(-n,-p,-q)
    assert t2 == 0


def test_substitute_indices():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    i, j, k, l, m, n, p, q = tensor_indices('i,j,k,l,m,n,p,q', Lorentz)
    sym2 = TensorSymmetry(get_symmetric_group_sgs(2))
    S2 = TensorType([Lorentz]*2, sym2)
    A, B = S2('A,B')
    t = A(i, k)*B(-k, -j)
    t1 = t.substitute_indices((i,j), (j, k))
    t1a = A(j, l)*B(-l, -k)
    assert t1 == t1a

def test_riemann_cyclic_replaceR():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    m0,m1,m2,m3 = tensor_indices('m0,m1,m2,m3', Lorentz)
    symr = TensorSymmetry(riemann_bsgs)
    R4 = TensorType([Lorentz]*4, symr)
    R = R4('R')
    t = R(m0,m2,m1,m3)
    t1 = riemann_cyclic_replaceR(t)
    t1a =  -S.One/3*R(m0, m3, m2, m1) + S.One/3*R(m0, m1, m2, m3) + Rational(2,3)*R(m0, m2, m1, m3)
    assert t1 == t1a

def test_riemann_cyclic():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    i, j, k, l, m, n, p, q = tensor_indices('i,j,k,l,m,n,p,q', Lorentz)
    symr = TensorSymmetry(riemann_bsgs)
    R4 = TensorType([Lorentz]*4, symr)
    R = R4('R')
    t = R(i,j,k,l) + R(i,l,j,k) + R(i,k,l,j) - \
        R(i,j,l,k) - R(i,l,k,j) - R(i,k,j,l)
    t2 = t*R(-i,-j,-k,-l)
    t3 = riemann_cyclic(t2)
    assert t3 == 0
    t = R(i,j,k,l)*(R(-i,-j,-k,-l) - 2*R(-i,-k,-j,-l))
    t1 = riemann_cyclic(t)
    assert t1 == 0

def test_div():
    Lorentz = TensorIndexType('Lorentz', dummy_fmt='L')
    m0,m1,m2,m3 = tensor_indices('m0,m1,m2,m3', Lorentz)
    symr = TensorSymmetry(riemann_bsgs)
    R4 = TensorType([Lorentz]*4, symr)
    R = R4('R')
    t = R(m0,m1,-m1,m3)
    t1 = t/S(4)
    assert str(t1) == '1/4*R(m0, L_0, -L_0, m3)'
    t = t.canon_bp()
    assert not t1._is_canon_bp
    t1 = t*4
    assert t1._is_canon_bp
    t1 = t1/4
    assert t1._is_canon_bp
