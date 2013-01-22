from sympy import Symbol, S, I
from sympy.tensor.tensor import (TensorIndexType, tensor_indices, \
TensorSymmetry, get_symmetric_group_sgs, TensorType, TensMul)
from sympy.tensor.dgamma_matr import GammaMatrices
from sympy.utilities.pytest import skip, XFAIL


D = Symbol('D')
Lorentz = TensorIndexType('Lorentz', dim=D, eps_dim=4, dummy_fmt='L')
m0, m1, m2, m3, m4, m5 = tensor_indices('m0,m1,m2,m3,m4,m5', Lorentz)
n0, n1, n2, n3, n4, n5 = tensor_indices('n0,n1,n2,n3,n4,n5', Lorentz)
sym1 = TensorSymmetry(get_symmetric_group_sgs(1))
S1 = TensorType([Lorentz], sym1)

GM = GammaMatrices( Lorentz)
match1_gamma = GM.match1_gamma
G = GM.G
g = GM.g
epsilon = GM.epsilon
rule1_gamma = GM.rule1_gamma
do_rule1_gamma = GM.do_rule1_gamma
match2_gamma = GM.match2_gamma
rule2_gamma = GM.rule2_gamma
do_rule2_gamma = GM.do_rule2_gamma
gamma_trace = GM.gamma_trace
gamma_trace1 = GM.gamma_trace1
gctr = GM.gctr


def test_rule1_gamma():
    t = 2*G(m2)*G(m0)*G(m1)*G(-m0)*G(-m1)
    r = match1_gamma(t, 1)
    t = rule1_gamma(t, 1)
    r = match1_gamma(t, 0)
    t = rule1_gamma(t, 0)
    assert t == (D*(-2*D + 4))*G(m2)

    t = 3*G(m2)*G(m0)*G(m1)*G(-m0)
    t = rule1_gamma(t, 1)

    t = G(m2)*G(m0)*G(m1)*G(-m0)*G(-m2)
    r = match1_gamma(t, 1)
    assert r == (1, 3)
    t = rule1_gamma(t, 1)
    t = rule1_gamma(t, 1)
    assert t == ((-D + 2)**2)*G(m1)

    t = G(m0)*G(m1)*G(m2)*G(m3)*G(-m1)
    r = match1_gamma(t, 2)
    assert r == (1, 4)
    t = rule1_gamma(t, 2)
    assert t == (D - 4)*G(m0)*G(m2)*G(m3) + 4*g(m2, m3)*G(m0)

    t = G(m0)*G(m1)*G(m2)*G(m3)*G(-m1)*G(-m0)
    r = match1_gamma(t, 2)
    t = rule1_gamma(t, 2)

    t = G(m0)*G(m1)*G(m2)*G(m3)*G(-m1)*G(-m0)
    t = do_rule1_gamma(t)
    t = do_rule1_gamma(t)
    assert t == ((D - 4)**2)*G(m2)*G(m3) + (8*D - 16)*g(m2, m3)

    t = G(m2)*G(m0)*G(m1)*G(-m2)*G(-m0)
    t = do_rule1_gamma(t)
    t = do_rule1_gamma(t)
    assert t == ((-D + 2)*(D - 4) + 4)*G(m1)

    t = G(m3)*G(m1)*G(m0)*G(m2)*G(-m3)*G(-m0)*G(-m2)
    t = do_rule1_gamma(t)
    t = do_rule1_gamma(t)
    t = do_rule1_gamma(t)
    assert t == (-4*D + (-D + 2)**2*(D - 4) + 8)*G(m1)

    M = Symbol('M')

    p = S1('p')
    ps = p(m0)*G(-m0)
    t0 = ps + M
    t = G(m0)*(ps + M)*G(-m0)
    t = do_rule1_gamma(t)
    assert t == (2 - D)*ps + D*M

    t = G(m0)*(ps + M)*G(-m0)*(ps + M)
    t = do_rule1_gamma(t)

    t = 2*G(m0)*G(m1)*G(m2)*G(m3)*G(-m0)
    t = do_rule1_gamma(t)
    assert t == (-2*D + 8)*G(m1)*G(m2)*G(m3) - 4*G(m3)*G(m2)*G(m1)

    t = G(m5)*G(m0)*G(m1)*G(m4)*G(m2)*G(-m4)*G(m3)*G(-m0)
    t = do_rule1_gamma(t)
    t = do_rule1_gamma(t)
    assert t == ((-D + 2)*(-D + 4))*G(m5)*G(m1)*G(m2)*G(m3) + (2*D - 4)*G(m5)*G(m3)*G(m2)*G(m1)

    t = -G(m0)*G(m1)*G(m2)*G(m3)*G(-m0)*G(m4)
    t = do_rule1_gamma(t)
    assert t == (D - 4)*G(m1)*G(m2)*G(m3)*G(m4) + 2*G(m3)*G(m2)*G(m1)*G(m4)

    t = G(-m5)*G(m0)*G(m1)*G(m2)*G(m3)*G(m4)*G(-m0)*G(m5)
    t = do_rule1_gamma(t, doall=True)
    assert t == ((-D + 4)**2 + 4)*G(m1)*G(m2)*G(m3)*G(m4) + (4*D - 16)*G(m3)*G(m2)*G(m1)*G(m4) + (4*D - 16)*G(m4)*G(m1)*G(m2)*G(m3) + 4*G(m2)*G(m1)*G(m4)*G(m3) + 4*G(m3)*G(m4)*G(m1)*G(m2) + 4*G(m4)*G(m3)*G(m2)*G(m1)


def test_rule2_gamma():
    p = S1('p')
    M = Symbol('M')
    ps = p(m0)*G(-m0)
    t = ps*ps
    t = t.canon_bp()
    t = rule2_gamma(t, 0)
    assert t == p(m0)*p(-m0)

    t = 2*G(m1)*ps*ps*G(m2)
    t = rule2_gamma(t, 0)
    assert t == 2*G(m1)*G(m2)*p(m0)*p(-m0)

    t0 = ps + M
    t = G(m0)*(ps + M)*G(-m0)
    t = do_rule1_gamma(t)
    assert t == (2 - D)*ps + D*M

    t = G(m0)*(ps + M)*G(-m0)*(ps + M)
    t = do_rule1_gamma(t)
    t = do_rule2_gamma(t)
    assert t == (-D + 2)*p(m0)*p(-m0) + (D*M + M*(-D + 2))*G(m0)*p(-m0) + D*M**2


def test_gamma_trace():
    a = [G(m0), G(m1)]
    t1 = gamma_trace1(*a)
    a = [G(m0), G(m1), G(m2), G(m3)]
    t1 = gamma_trace1(*a)

    t = G(m0)*G(m1)
    t1 = gamma_trace(t)
    assert t1 == 4*g(m0, m1)
    t = G(m0)*G(m1)*G(m2)*G(m3)
    t1 = gamma_trace(t)
    t2 = -4*g(m0, m2)*g(m1, m3) + 4*g(m0, m1)*g(m2, m3) + 4*g(m0, m3)*g(m1, m2)
    st2 = str(t2)
    assert t1 == t2

    t = G(m0)*G(m1)*G(m2)*G(m3)*G(m4)*G(-m0)
    t1 = gamma_trace(t)
    assert t1 == (-4*D)*g(m1, m3)*g(m2, m4) + (4*D)*g(m1, m2)*g(m3, m4) + \
                 (4*D)*g(m1, m4)*g(m2, m3)

    t = G(-m5)*G(m0)*G(m1)*G(m2)*G(m3)*G(m4)*G(-m0)*G(m5)
    t1 = gamma_trace(t)
    assert t1 == (32*D + 4*(-D + 4)**2 - 64)*(g(m1, m2)*g(m3, m4) - \
            g(m1, m3)*g(m2, m4) + g(m1, m4)*g(m2, m3))

    t = G(m0)*G(m1)*G(-m0)*G(m3)
    t1 = gamma_trace(t)
    assert t1 == (-4*D + 8)*g(m1, m3)

    t = G(m0)*G(m1)*G(m2)*G(m3)*G(m4)*G(m5)
    t1 = gamma_trace(t)
    t2 = t1*g(-m0, -m5)
    t2 = t2.contract_metric(g)
    assert t2 == D*gamma_trace(G(m1)*G(m2)*G(m3)*G(m4))

    p, q = S1('p,q')
    ps = p(m0)*G(-m0)
    qs = q(m0)*G(-m0)
    t = ps*qs*ps*qs
    t1 = gamma_trace(t)
    assert t1 == 8*p(m0)*q(-m0)*p(m1)*q(-m1) - 4*p(m0)*p(-m0)*q(m1)*q(-m1)

def test_gamma5():
    G5 = GM.G5
    epsilon = GM.epsilon
    t1 = GM.G5_to_right(2*G5)
    assert t1 == 2*G5
    assert gamma_trace(t1) == 0

    t = G5*G(m0)
    t1 = GM.G5_to_right(t)
    assert t1 == - G(m0)*G5

    t = G5*G(m0)*G5*G(m1)
    t1 = GM.G5_to_right(t)
    assert t1 == -G(m0)*G(m1)
    t = G5*G(m0)*G5*G(m1)*G5
    t1 = GM.G5_to_right(t)
    assert t1 == -G(m0)*G(m1)*G5

    t = (1 + G5)*(1 - G5)
    t1 = GM.G5_to_right(t)
    assert t1 == 0

    t = G5*G(m0)
    t1 = gamma_trace(t)
    assert t1 == 0

    t = G(m0)*G5*G(m1)
    t1 = gamma_trace(t)
    assert t1 == 0

    t= G5*G(m0)*G(m1)*G(m2)*G(m3)
    t1 = gamma_trace(t)
    assert t1 == 4*I*epsilon(m0, m1, m2, m3)

    t= G5*G(m0)*G(-m0)*G(m1)*G(m2)*G(m3)
    t1 = gamma_trace(t)
    assert t1 == 0

    t= G5*G(m0)*G(-m0)*G(m1)*G(m2)*G(m3)*G(m4)
    t1 = gamma_trace(t)
    assert t1 == 4*D*I*epsilon(m1, m2, m3, m4)

    t= G5*G(m0)*G(m1)*G(-m0)*G(m2)*G(m3)*G(m4)
    t1 = gamma_trace(t)
    assert t1 == (I*(-4*D + 8))*epsilon(m1, m2, m3, m4)

    t = G(m0)*(1 - G5)
    t1 = gamma_trace(t)
    assert t1 == 0

    t = G(m0)*(1 - G5)*G(m1)*G(m2)*G(m3)
    t1 = gamma_trace(t)
    assert t1 == 4*I*epsilon(m0, m1, m2, m3) + 4*(g(m0,m1)*g(m2,m3) - \
            g(m0,m2)*g(m1,m3) + g(m0,m3)*g(m1,m2))

    t = G(m0)*(1 - G5)*G(m1)*(1 + G5)*G(m2)*(1 - G5)*G(m3)
    t1 = gamma_trace(t)
    assert t1 == 16*I*epsilon(m0, m1, m2, m3) + 16*(g(m0,m1)*g(m2,m3) - \
            g(m0,m2)*g(m1,m3) + g(m0,m3)*g(m1,m2))

    M = Symbol('M')
    p, q = S1('p,q')
    ps = p(m0)*G(-m0)
    qs = q(m0)*G(-m0)
    p2 = p(m0)*p(-m0)
    q2 = q(m0)*q(-m0)
    pq = p(m0)*q(-m0)
    t = (ps + M)*G5*(qs + M)*G5
    t1 = gamma_trace(t)
    assert t1 == -4*p(m0)*q(-m0) + 4*M**2

    t = G(m0)*(ps + M)*G(m2)*(ps + qs + M)*G5
    t1 = gamma_trace(t)
    assert t1 == -4*I*epsilon(m0,m2,m1,m3)*p(-m1)*q(-m3)

    t = G(m0)*(ps + M)*G(m1)*(ps + qs)*G(m2)*(1 - G5)
    t1 = gamma_trace(t)
    assert t1 == 4*M*(-g(m1, m2)*p(m0) -g(m1, m2)*q(m0) + \
      I*epsilon(m0, m1, m2, m3)*p(-m3) + I*epsilon(m0, m1, m2, m3)*q(-m3) + \
        g(m0, m1)*p(m2) + g(m0, m1)*q(m2) + g(m0, m2)*p(m1) + g(m0, m2)*q(m1))

    t = G(m0)*G(m1)*G(m2)*G(m3)*G(m4)*G(m5)*G5
    t1 = gamma_trace(t)
    t2 = t1*g(-m0, -m2)
    t2 = t2.contract_metric(g, contract_all=True)
    assert t2 == (-4*I*D + 8*I)*epsilon(m1, m3, m4, m5)

def test_trace1():
    M = Symbol('M')
    p, q = S1('p,q')
    ps = p(m0)*G(-m0)
    qs = q(m0)*G(-m0)
    p2 = p(m0)*p(-m0)
    q2 = q(m0)*q(-m0)
    pq = p(m0)*q(-m0)

    t = ps*qs*ps*qs*ps*qs
    t1 = gamma_trace(t)
    assert t1 == -12*p2*pq*q2 + 16*pq*pq*pq

    t = ps*qs*ps*qs*ps*qs*ps*qs
    t1 = gamma_trace(t)
    assert t1 == -32*pq*pq*p2*q2 + 32*pq*pq*pq*pq + 4*p2*p2*q2*q2

def test_epsilon():
    t = G(m0)*G(m1)*G(m2)*G(m3)*G(n0)*G(n1)*G(n2)*G(n3)*epsilon(-n0,-n1,-n2,-n3)

    t1 = gamma_trace(t)
    assert t1 == 96*epsilon(m0, m1, m2, m3)
