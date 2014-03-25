from sympy.tensor.tensor import tensor_indices, TensorIndexType, tensorhead, TensorManager, TensMul, TensAdd, get_lines, TensExpr
from sympy import simplify, trace
from sympy.physics.hep.gamma_matrices import GammaMatrix as G, GammaMatrixHead, DiracSpinorIndex
from sympy.utilities.pytest import XFAIL, raises

def _is_tensor_eq(arg1, arg2):
    if isinstance(arg1, TensExpr):
        return arg1.equals(arg2)
    elif isinstance(arg2, TensExpr):
        return arg2.equals(arg1)
    return arg1 == arg2

def execute_gamma_simplify_tests_for_function(tfunc, D):
    """
    Perform tests to check if sfunc is able to simplify gamma matrix expressions.

    Parameters
    ==========

    `sfunc`     a function to simplify a `TIDS`, shall return the simplified `TIDS`.
    `D`         the number of dimension (in most cases `D=4`).

    """

    mu, nu, rho, sigma = tensor_indices("mu, nu, rho, sigma", G.LorentzIndex)
    a1, a2, a3, a4, a5, a6 = tensor_indices("a1:7", G.LorentzIndex)
    mu11, mu12, mu21, mu31, mu32, mu41, mu51, mu52 = tensor_indices("mu11, mu12, mu21, mu31, mu32, mu41, mu51, mu52", G.LorentzIndex)
    mu61, mu71, mu72 = tensor_indices("mu61, mu71, mu72", G.LorentzIndex)
    m0, m1, m2, m3, m4, m5, m6 = tensor_indices("m0:7", G.LorentzIndex)

    def g(xx, yy):
        return (G(xx)*G(yy) + G(yy)*G(xx))/2

    # Some examples taken from Kahane's paper, 4 dim only:
    if D == 4:
        t = (G(a1)*G(mu11)*G(a2)*G(mu21)*G(-a1)*G(mu31)*G(-a2))
        assert _is_tensor_eq(tfunc(t), -4*G(mu11)*G(mu31)*G(mu21) - 4*G(mu31)*G(mu11)*G(mu21))

        t = (G(a1)*G(mu11)*G(mu12)*\
                              G(a2)*G(mu21)*\
                              G(a3)*G(mu31)*G(mu32)*\
                              G(a4)*G(mu41)*\
                              G(-a2)*G(mu51)*G(mu52)*\
                              G(-a1)*G(mu61)*\
                              G(-a3)*G(mu71)*G(mu72)*\
                              G(-a4))
        assert _is_tensor_eq(tfunc(t), \
            16*G(mu31)*G(mu32)*G(mu72)*G(mu71)*G(mu11)*G(mu52)*G(mu51)*G(mu12)*G(mu61)*G(mu21)*G(mu41) + 16*G(mu31)*G(mu32)*G(mu72)*G(mu71)*G(mu12)*G(mu51)*G(mu52)*G(mu11)*G(mu61)*G(mu21)*G(mu41) + 16*G(mu71)*G(mu72)*G(mu32)*G(mu31)*G(mu11)*G(mu52)*G(mu51)*G(mu12)*G(mu61)*G(mu21)*G(mu41) + 16*G(mu71)*G(mu72)*G(mu32)*G(mu31)*G(mu12)*G(mu51)*G(mu52)*G(mu11)*G(mu61)*G(mu21)*G(mu41))

    # Fully G.Lorentz-contracted expressions, these return scalars:

    def add_delta(ne):
        return ne * DiracSpinorIndex.delta(DiracSpinorIndex.auto_left, -DiracSpinorIndex.auto_right)

    t = (G(mu)*G(-mu))
    ts = add_delta(D)
    assert _is_tensor_eq(tfunc(t), ts)

    t = (G(mu)*G(nu)*G(-mu)*G(-nu))
    ts = add_delta(2*D - D**2)  # -8
    assert _is_tensor_eq(tfunc(t), ts)

    t = (G(mu)*G(nu)*G(-nu)*G(-mu))
    ts = add_delta(D**2)  # 16
    assert _is_tensor_eq(tfunc(t), ts)

    t = (G(mu)*G(nu)*G(-rho)*G(-nu)*G(-mu)*G(rho))
    ts = add_delta(4*D - 4*D**2 + D**3)  # 16
    assert _is_tensor_eq(tfunc(t), ts)

    t = (G(mu)*G(nu)*G(rho)*G(-rho)*G(-nu)*G(-mu))
    ts = add_delta(D**3)  # 64
    assert _is_tensor_eq(tfunc(t), ts)

    t = (G(a1)*G(a2)*G(a3)*G(a4)*G(-a3)*G(-a1)*G(-a2)*G(-a4))
    ts = add_delta(-8*D + 16*D**2 - 8*D**3 + D**4)  # -32
    assert _is_tensor_eq(tfunc(t), ts)

    t = (G(-mu)*G(-nu)*G(-rho)*G(-sigma)*G(nu)*G(mu)*G(sigma)*G(rho))
    ts = add_delta(-16*D + 24*D**2 - 8*D**3 + D**4)  # 64
    assert _is_tensor_eq(tfunc(t), ts)

    t = (G(-mu)*G(nu)*G(-rho)*G(sigma)*G(rho)*G(-nu)*G(mu)*G(-sigma))
    ts = add_delta(8*D - 12*D**2 + 6*D**3 - D**4)  # -32
    assert _is_tensor_eq(tfunc(t), ts)

    t = (G(a1)*G(a2)*G(a3)*G(a4)*G(a5)*G(-a3)*G(-a2)*G(-a1)*G(-a5)*G(-a4))
    ts = add_delta(64*D - 112*D**2 + 60*D**3 - 12*D**4 + D**5)  # 256
    assert _is_tensor_eq(tfunc(t), ts)

    t = (G(a1)*G(a2)*G(a3)*G(a4)*G(a5)*G(-a3)*G(-a1)*G(-a2)*G(-a4)*G(-a5))
    ts = add_delta(64*D - 120*D**2 + 72*D**3 - 16*D**4 + D**5)  # -128
    assert _is_tensor_eq(tfunc(t), ts)

    t = (G(a1)*G(a2)*G(a3)*G(a4)*G(a5)*G(a6)*G(-a3)*G(-a2)*G(-a1)*G(-a6)*G(-a5)*G(-a4))
    ts = add_delta(416*D - 816*D**2 + 528*D**3 - 144*D**4 + 18*D**5 - D**6)  # -128
    assert _is_tensor_eq(tfunc(t), ts)

    t = (G(a1)*G(a2)*G(a3)*G(a4)*G(a5)*G(a6)*G(-a2)*G(-a3)*G(-a1)*G(-a6)*G(-a4)*G(-a5))
    ts = add_delta(416*D - 848*D**2 + 584*D**3 - 172*D**4 + 22*D**5 - D**6)  # -128
    assert _is_tensor_eq(tfunc(t), ts)

    # Expressions with free indices:

    t = (G(mu)*G(nu)*G(rho)*G(sigma)*G(-mu))
    assert _is_tensor_eq(tfunc(t), (-2*G(sigma)*G(rho)*G(nu) + (4-D)*G(nu)*G(rho)*G(sigma)))

    t = (G(mu)*G(nu)*G(-mu))
    assert _is_tensor_eq(tfunc(t), (2-D)*G(nu))

    t = (G(mu)*G(nu)*G(rho)*G(-mu))
    assert _is_tensor_eq(tfunc(t), 2*G(nu)*G(rho) + 2*G(rho)*G(nu) - (4-D)*G(nu)*G(rho))

    t = 2*G(m2)*G(m0)*G(m1)*G(-m0)*G(-m1)
    st = tfunc(t)
    assert _is_tensor_eq(st, (D*(-2*D + 4))*G(m2))

    t = G(m2)*G(m0)*G(m1)*G(-m0)*G(-m2)
    st = tfunc(t)
    assert _is_tensor_eq(st, ((-D + 2)**2)*G(m1))

    t = G(m0)*G(m1)*G(m2)*G(m3)*G(-m1)
    st = tfunc(t)
    assert _is_tensor_eq(st, (D - 4)*G(m0)*G(m2)*G(m3) + 4*G(m0)*g(m2, m3))

    t = G(m0)*G(m1)*G(m2)*G(m3)*G(-m1)*G(-m0)
    st = tfunc(t)
    assert _is_tensor_eq(st, ((D - 4)**2)*G(m2)*G(m3) + (8*D - 16)*g(m2, m3))

    t = G(m2)*G(m0)*G(m1)*G(-m2)*G(-m0)
    st = tfunc(t)
    assert _is_tensor_eq(st, ((-D + 2)*(D - 4) + 4)*G(m1))

    t = G(m3)*G(m1)*G(m0)*G(m2)*G(-m3)*G(-m0)*G(-m2)
    st = tfunc(t)
    assert _is_tensor_eq(st, (-4*D + (-D + 2)**2*(D - 4) + 8)*G(m1))

    t = 2*G(m0)*G(m1)*G(m2)*G(m3)*G(-m0)
    st = tfunc(t)
    assert _is_tensor_eq(st, ((-2*D + 8)*G(m1)*G(m2)*G(m3) - 4*G(m3)*G(m2)*G(m1)))

    t = G(m5)*G(m0)*G(m1)*G(m4)*G(m2)*G(-m4)*G(m3)*G(-m0)
    st = tfunc(t)
    assert _is_tensor_eq(st, (((-D + 2)*(-D + 4))*G(m5)*G(m1)*G(m2)*G(m3) + (2*D - 4)*G(m5)*G(m3)*G(m2)*G(m1)))

    t = -G(m0)*G(m1)*G(m2)*G(m3)*G(-m0)*G(m4)
    st = tfunc(t)
    assert _is_tensor_eq(st, ((D - 4)*G(m1)*G(m2)*G(m3)*G(m4) + 2*G(m3)*G(m2)*G(m1)*G(m4)))

    t = G(-m5)*G(m0)*G(m1)*G(m2)*G(m3)*G(m4)*G(-m0)*G(m5)
    st = tfunc(t)

    result1 = ((-D + 4)**2 + 4)*G(m1)*G(m2)*G(m3)*G(m4) +\
        (4*D - 16)*G(m3)*G(m2)*G(m1)*G(m4) + (4*D - 16)*G(m4)*G(m1)*G(m2)*G(m3)\
        + 4*G(m2)*G(m1)*G(m4)*G(m3) + 4*G(m3)*G(m4)*G(m1)*G(m2) +\
        4*G(m4)*G(m3)*G(m2)*G(m1)

    # Kahane's algorithm yields this result, which is equivalent to `result1`
    # in four dimensions, but is not automatically recognized as equal:
    result2 = 8*G(m1)*G(m2)*G(m3)*G(m4) + 8*G(m4)*G(m3)*G(m2)*G(m1)

    if D == 4:
        assert _is_tensor_eq(st, (result1)) or _is_tensor_eq(st, (result2))
    else:
        assert _is_tensor_eq(st, (result1))

    # and a few very simple cases, with no contracted indices:

    t = G(m0)
    st = tfunc(t)
    assert _is_tensor_eq(st, t)

    t = -7*G(m0)
    st = tfunc(t)
    assert _is_tensor_eq(st, t)

    t = 224*G(m0)*G(m1)*G(-m2)*G(m3)
    st = tfunc(t)
    assert _is_tensor_eq(st, t)


def test_kahane_algorithm():
    # Wrap this function to convert to and from TIDS:

    def tfunc(e):
        return GammaMatrixHead._simplify_single_line(e)

    execute_gamma_simplify_tests_for_function(tfunc, D=4)

def test_kahane_simplify1():
    i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15 = tensor_indices('i0:16', G.LorentzIndex)
    mu, nu, rho, sigma = tensor_indices("mu, nu, rho, sigma", G.LorentzIndex)
    KD = DiracSpinorIndex.delta
    sl = DiracSpinorIndex.auto_left
    sr = DiracSpinorIndex.auto_right
    D = 4
    s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16 = \
                   tensor_indices('s0:17', DiracSpinorIndex)
    t = DiracSpinorIndex.delta(s0,s1)
    t = G(i0)*G(i1)
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals(t)

    t = G(i0)*G(i1)*G(-i0)
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals(-2*G(i1))
    t = G(i0,s0,-s1)*G(i1,s1,-s2)*G(-i0,s2,-s3)
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals(-2*G(i1, s0, -s3))

    t = G(i0, s0, -s1)*G(i1, s1, -s2)
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals(t)
    t = G(i0, s0, -s1)*G(i1, s1, -s0)
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals(t)
    t = G(i0)*G(-i0)
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals(4*KD(sl, -sr))
    t = G(i0,s0,-s1)*G(-i0,s1,-s2)
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals(4*KD(s0, -s2))
    t = G(i0,s0,-s1)*G(-i0,s1,-s0)
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals(16)
    t = G(i0)*G(i1)*G(-i0)
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals(-2*G(i1))
    t = G(i0)*G(i1)*G(-i0)*G(-i1)
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals((2*D - D**2)*KD(sl, -sr))
    t = G(i0,s0,-s1)*G(i1,s1,-s2)*G(-i0,s2,-s3)*G(-i1,s3,-s0)
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals(4*(2*D - D**2))
    t = G(i0,s0,-s1)*G(-i0,s2,-s3)*G(i1,s1,-s2)*G(-i1,s3,-s0)
    raises(ValueError, lambda: G._kahane_simplify(t.coeff, t._tids))
    t = (G(mu)*G(nu)*G(-nu)*G(-mu))
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals(D**2*KD(sl, -sr))
    t = (G(mu,s0,-s1)*G(nu,s1,-s2)*G(-nu,s2,-s3)*G(-mu,s3,-s4))
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals(D**2*KD(s0, -s4))
    t = (G(mu,s0,-s1)*G(nu,s1,-s2)*G(-nu,s2,-s3)*G(-mu,s3,-s0))
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals(4*D**2)
    t = (G(mu)*G(nu)*G(-rho)*G(-nu)*G(-mu)*G(rho))
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals((4*D - 4*D**2 + D**3)*KD(sl, -sr))
    t = (G(-mu)*G(-nu)*G(-rho)*G(-sigma)*G(nu)*G(mu)*G(sigma)*G(rho))
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals((-16*D + 24*D**2 - 8*D**3 + D**4)*KD(sl, -sr))
    t = (G(-mu)*G(nu)*G(-rho)*G(sigma)*G(rho)*G(-nu)*G(mu)*G(-sigma))
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals((8*D - 12*D**2 + 6*D**3 - D**4)*KD(sl, -sr))

    # Expressions with free indices:
    t = (G(mu)*G(nu)*G(rho)*G(sigma)*G(-mu))
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals(-2*G(sigma)*G(rho)*G(nu))
    t = (G(mu,s0,-s1)*G(nu,s1,-s2)*G(rho,s2,-s3)*G(sigma,s3,-s4)*G(-mu,s4,-s5))
    r = G._kahane_simplify(t.coeff, t._tids)
    assert r.equals(-2*G(sigma,s0,-s1)*G(rho,s1,-s2)*G(nu,s2,-s5))


def test_gamma_matrix_class():
    i, j, k = tensor_indices('i,j,k', G.LorentzIndex)

    # define another type of TensorHead to see if exprs are correctly handled:
    A = tensorhead('A', [G.LorentzIndex], [[1]])

    t = A(k)*G(i)*G(-i)
    ts = simplify(t)
    assert _is_tensor_eq(ts, 4*A(k)*DiracSpinorIndex.delta(DiracSpinorIndex.auto_left, -DiracSpinorIndex.auto_right))

    t = G(i)*A(k)*G(j)
    ts = simplify(t)
    assert _is_tensor_eq(ts, A(k)*G(i)*G(j))

    execute_gamma_simplify_tests_for_function(simplify, D=4)

def test_gamma_matrix_trace():
    gamma_trace = G.gamma_trace
    g = G.LorentzIndex.metric

    m0, m1, m2, m3, m4, m5, m6 = tensor_indices('m0:7', G.LorentzIndex)
    n0, n1, n2, n3, n4, n5 = tensor_indices('n0:6', G.LorentzIndex)

    # working in D=4 dimensions
    D = 4

    # traces of odd number of gamma matrices are zero:
    t = G(m0)
    t1 = gamma_trace(t)
    assert t1.equals(0)

    t = G(m0)*G(m1)*G(m2)
    t1 = gamma_trace(t)
    assert t1.equals(0)

    t = G(m0)*G(m1)*G(-m0)
    t1 = gamma_trace(t)
    assert t1.equals(0)

    t = G(m0)*G(m1)*G(m2)*G(m3)*G(m4)
    t1 = gamma_trace(t)
    assert t1.equals(0)

    # traces without internal contractions:
    t = G(m0)*G(m1)
    t1 = gamma_trace(t)
    assert _is_tensor_eq(t1, 4*g(m0, m1))

    t = G(m0)*G(m1)*G(m2)*G(m3)
    t1 = gamma_trace(t)
    t2 = -4*g(m0, m2)*g(m1, m3) + 4*g(m0, m1)*g(m2, m3) + 4*g(m0, m3)*g(m1, m2)
    st2 = str(t2)
    assert _is_tensor_eq(t1, t2)

    t = G(m0)*G(m1)*G(m2)*G(m3)*G(m4)*G(m5)
    t1 = gamma_trace(t)
    t2 = t1*g(-m0, -m5)
    t2 = t2.contract_metric(g)
    assert _is_tensor_eq(t2, D*gamma_trace(G(m1)*G(m2)*G(m3)*G(m4)))

    # traces of expressions with internal contractions:
    t = G(m0)*G(-m0)
    t1 = gamma_trace(t)
    assert t1.equals(4*D)

    t = G(m0)*G(m1)*G(-m0)*G(-m1)
    t1 = gamma_trace(t)
    assert t1.equals(8*D - 4*D**2)

    t = G(m0)*G(m1)*G(m2)*G(m3)*G(m4)*G(-m0)
    t1 = gamma_trace(t)
    t2 = (-4*D)*g(m1, m3)*g(m2, m4) + (4*D)*g(m1, m2)*g(m3, m4) + \
                 (4*D)*g(m1, m4)*g(m2, m3)
    assert t1.equals(t2)

    t = G(-m5)*G(m0)*G(m1)*G(m2)*G(m3)*G(m4)*G(-m0)*G(m5)
    t1 = gamma_trace(t)
    t2 = (32*D + 4*(-D + 4)**2 - 64)*(g(m1, m2)*g(m3, m4) - \
            g(m1, m3)*g(m2, m4) + g(m1, m4)*g(m2, m3))
    assert t1.equals(t2)

    t = G(m0)*G(m1)*G(-m0)*G(m3)
    t1 = gamma_trace(t)
    assert t1.equals((-4*D + 8)*g(m1, m3))

#    p, q = S1('p,q')
#    ps = p(m0)*G(-m0)
#    qs = q(m0)*G(-m0)
#    t = ps*qs*ps*qs
#    t1 = gamma_trace(t)
#    assert t1 == 8*p(m0)*q(-m0)*p(m1)*q(-m1) - 4*p(m0)*p(-m0)*q(m1)*q(-m1)

    t = G(m0)*G(m1)*G(m2)*G(m3)*G(m4)*G(m5)*G(-m0)*G(-m1)*G(-m2)*G(-m3)*G(-m4)*G(-m5)
    t1 = gamma_trace(t)
    assert t1.equals(-4*D**6 + 120*D**5 - 1040*D**4 + 3360*D**3 - 4480*D**2 + 2048*D)

    t = G(m0)*G(m1)*G(n1)*G(m2)*G(n2)*G(m3)*G(m4)*G(-n2)*G(-n1)*G(-m0)*G(-m1)*G(-m2)*G(-m3)*G(-m4)
    t1 = gamma_trace(t)
    tresu = -7168*D + 16768*D**2 - 14400*D**3 + 5920*D**4 - 1232*D**5 + 120*D**6 - 4*D**7
    assert t1.equals(tresu)

    # checked with Mathematica
    # In[1]:= <<Tracer.m
    # In[2]:= Spur[l];
    # In[3]:= GammaTrace[l, {m0},{m1},{n1},{m2},{n2},{m3},{m4},{n3},{n4},{m0},{m1},{m2},{m3},{m4}]
    t = G(m0)*G(m1)*G(n1)*G(m2)*G(n2)*G(m3)*G(m4)*G(n3)*G(n4)*G(-m0)*G(-m1)*G(-m2)*G(-m3)*G(-m4)
    t1 = gamma_trace(t)
#    t1 = t1.expand_coeff()
    c1 = -4*D**5 + 120*D**4 - 1200*D**3 + 5280*D**2 - 10560*D + 7808
    c2 = -4*D**5 + 88*D**4 - 560*D**3 + 1440*D**2 - 1600*D + 640
    assert _is_tensor_eq(t1, c1*g(n1, n4)*g(n2, n3) + c2*g(n1, n2)*g(n3, n4) + \
            (-c1)*g(n1, n3)*g(n2, n4))

    p, q = tensorhead('p,q', [G.LorentzIndex], [[1]])
    ps = p(m0)*G(-m0)
    qs = q(m0)*G(-m0)
    p2 = p(m0)*p(-m0)
    q2 = q(m0)*q(-m0)
    pq = p(m0)*q(-m0)
    t = ps*qs*ps*qs
    r = gamma_trace(t)
    assert _is_tensor_eq(r, 8*pq*pq - 4*p2*q2)
    t = ps*qs*ps*qs*ps*qs
    r = gamma_trace(t)
    assert r.equals(-12*p2*pq*q2 + 16*pq*pq*pq)
    t = ps*qs*ps*qs*ps*qs*ps*qs
    r = gamma_trace(t)
    assert r.equals(-32*pq*pq*p2*q2 + 32*pq*pq*pq*pq + 4*p2*p2*q2*q2)

    t = 4*p(m1)*p(m0)*p(-m0)*q(-m1)*q(m2)*q(-m2)
    assert _is_tensor_eq(gamma_trace(t), t)
    t = ps*ps*ps*ps*ps*ps*ps*ps
    r = gamma_trace(t)
    assert r.equals(4*p2*p2*p2*p2)

def test_simple_trace_cases_symbolic_dim():
    from sympy import symbols
    D = symbols('D')
    G = GammaMatrixHead(dim=D)

    m0, m1, m2, m3 = tensor_indices('m0:4', G.LorentzIndex)
    g = G.LorentzIndex.metric

    t = G(m0)*G(m1)
    t1 = G._trace_single_line(t)
    assert _is_tensor_eq(t1, 4 * G.LorentzIndex.metric(m0, m1))

    t = G(m0)*G(m1)*G(m2)*G(m3)
    t1 = G._trace_single_line(t)
    t2 = -4*g(m0, m2)*g(m1, m3) + 4*g(m0, m1)*g(m2, m3) + 4*g(m0, m3)*g(m1, m2)
    assert _is_tensor_eq(t1, t2)

def test_get_lines():
    i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13 = \
       tensor_indices('i0:14', G.LorentzIndex)
    s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16 = \
       tensor_indices('s0:17', DiracSpinorIndex)
    t = G(i1,s1,-s2)*G(i2,s3,-s4)*G(i4,s2,-s6)*G(i3,s4,-s3)
    r = get_lines(t, DiracSpinorIndex)
    assert r == ([[0, 2]], [[1, 3]], [])
    t = G(i1,s1,-s2)*G(i2,s2,-s3)*G(i3,s3,-s4)*G(i4,s4,-s5)*\
        G(i5,s6,-s7)*G(i6,s7,-s8)*G(i7,s8,-s9)*G(i8,s9,-s6)
    r = get_lines(t, DiracSpinorIndex)
    assert r == ([[0, 1, 2, 3]], [[4, 5, 6, 7]], [])
    t = G(i1,s1,-s2)*G(i0,s0,-s10)*G(i2,s2,-s3)*G(i3,s3,-s4)*\
    G(i4,s4,-s5)*G(i5,s6,-s7)*G(i6,s7,-s8)*G(i7,s8,-s9)*\
    G(i8,s9,-s6)*G(i9,s10,-s0)
    r = get_lines(t, DiracSpinorIndex)
    assert r == ([[0, 2, 3, 4]], [[5, 6, 7, 8], [1, 9]], [])
    t = G(i1,s1,-s2)*G(i11,s12,-s13)*G(i0,s0,-s10)*G(i2,s2,-s3)*G(i3,s3,-s4)*\
        G(i4,s4,-s5)*G(i5,s6,-s7)*G(i10,s11,-s12)*G(i6,s7,-s8)*G(i7,s8,-s9)*\
        G(i8,s9,-s6)*G(i9,s10,-s0)
    r = get_lines(t, DiracSpinorIndex)
    assert r == ([[0, 3, 4, 5], [7, 1]], [[6, 8, 9, 10], [2, 11]], [])
    t = G(i4,s4,-s5)*G(i5,s6,-s7)*G(i10,s11,-s12)*G(i6,s7,-s8)*G(i7,s8,-s9)*\
        G(i8,s9,-s6)*G(i9,s10,-s0)*\
        G(i1,s1,-s2)*G(i11,s12,-s13)*G(i0,s0,-s10)*G(i2,s2,-s3)*G(i3,s3,-s4)
    r = get_lines(t, DiracSpinorIndex)
    assert r == ([[2, 8], [7, 10, 11, 0]], [[1, 3, 4, 5], [6, 9]], [])
    t = G(i8,s9,-s6)*G(i9,s10,-s0)*G(i4,s4,-s5)*G(i13,s14,-s15)*\
        G(i10,s11,-s12)*G(i1,s1,-s2)*G(i11,s12,-s13)*\
        G(i0,s0,-s10)*G(i6,s7,-s8)*G(i7,s8,-s9)*\
        G(i2,s2,-s3)*G(i12,s13,-s14)*G(i3,s3,-s4)*G(i5,s6,-s7)
    r = get_lines(t, DiracSpinorIndex)
    assert r == ([[4, 6, 11, 3], [5, 10, 12, 2]], [[1, 7], [0, 13, 8, 9]], [])

def test_simplify_lines():
    i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12 = tensor_indices('i0:13', G.LorentzIndex)
    s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16 = \
       tensor_indices('s0:17', DiracSpinorIndex)

    g = G.LorentzIndex.metric
    Sdelta = DiracSpinorIndex.delta
    t = G(i1,s1,-s2)*G(i2,s2,-s1)*G(i3,s4,-s5)*G(i4,s5,-s6)*G(i5,s7,-s8)
    r = G.simplify_lines(t)
    assert r.equals(4*G(i5, s7, -s8)*G(i3, s4, -s0)*G(i4, s0, -s6)*g(i1, i2))

    t = G(i1,s1,-s2)*G(i2,s2,-s1)*G(i3,s4,-s5)*G(-i3,s5,-s6)*G(i5,s7,-s8)
    r = G.simplify_lines(t)
    assert r.equals(16*G(i5, s7, -s8)*Sdelta(s4, -s6)*g(i1, i2))
    t = G(i1,s1,-s2)*G(i2,s2,-s1)*G(i3,s4,-s5)*G(i4,s5,-s6)*G(i5,s7,-s8)
    r = G.simplify_lines(t)
    assert r.equals(4*G(i5, s7, -s8)*G(i3, s4, s0)*G(i4, -s0, -s6)*g(i1, i2))
    t = G(i5,s7,-s8)*G(i6,s9,-s10)*G(i1,s1,-s2)*G(i3,s4,-s5)*G(i2,s2,-s1)*G(i4,s5,-s6)*G(-i6,s10,-s9)
    r = G.simplify_lines(t)
    assert r.equals(64*G(i5, s7, -s8)*G(i3, s4, s0)*G(i4, -s0, -s6)*g(i1, i2))
    t = G(i5,s7,-s8)*G(i6,s9,-s10)*G(i1,s1,-s2)*G(i7,s12,-s11)*G(i3,s4,-s5)*\
        G(i2,s2,-s1)*G(i4,s5,-s6)*G(-i6,s10,-s9)*G(-i7,s11,-s13)
    r = G.simplify_lines(t)
    assert r.equals(256*G(i5, s7, -s8)*G(i3, s4, s0)*G(i4, -s0, -s6)*\
           g(i1, i2)*Sdelta(s12,-s13))
