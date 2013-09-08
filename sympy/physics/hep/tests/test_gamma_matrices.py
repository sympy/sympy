from sympy.tensor.tensor import tensor_indices, TensorIndexType, tensorhead, TensorManager
from sympy.physics.hep.gamma_matrices import GammaMatrix4D as G
from sympy.physics.hep.gamma_matrices import kahane_simplify, Lorentz


def test_kahane_algorithm():
    mu, nu, rho, sigma = tensor_indices("mu, nu, rho, sigma", Lorentz)
    a1, a2, a3, a4, a5, a6 = tensor_indices("a1:7", Lorentz)
    mu11, mu12, mu21, mu31, mu32, mu41, mu51, mu52 = tensor_indices("mu11, mu12, mu21, mu31, mu32, mu41, mu51, mu52", Lorentz)
    mu61, mu71, mu72 = tensor_indices("mu61, mu71, mu72", Lorentz)
    m0, m1, m2, m3, m4, m5, m6 = tensor_indices("m0:7", Lorentz)

    def g(xx, yy):
        return (G(xx)*G(yy) + G(yy)*G(xx))/2

    # kahane_simplify only works in four dimensions:
    D = 4

    # Some examples taken from Kahane's paper, 4 dim only:
    if D == 4:
        t = (G(a1)*G(mu11)*G(a2)*G(mu21)*G(-a1)*G(mu31)*G(-a2))
        assert kahane_simplify(t) == -4*G(mu11)*G(mu31)*G(mu21) - 4*G(mu31)*G(mu11)*G(mu21)

        t = (G(a1)*G(mu11)*G(mu12)*\
                              G(a2)*G(mu21)*\
                              G(a3)*G(mu31)*G(mu32)*\
                              G(a4)*G(mu41)*\
                              G(-a2)*G(mu51)*G(mu52)*\
                              G(-a1)*G(mu61)*\
                              G(-a3)*G(mu71)*G(mu72)*\
                              G(-a4))
        assert kahane_simplify(t) == \
            16*G(mu31)*G(mu32)*G(mu72)*G(mu71)*G(mu11)*G(mu52)*G(mu51)*G(mu12)*G(mu61)*G(mu21)*G(mu41) + 16*G(mu31)*G(mu32)*G(mu72)*G(mu71)*G(mu12)*G(mu51)*G(mu52)*G(mu11)*G(mu61)*G(mu21)*G(mu41) + 16*G(mu71)*G(mu72)*G(mu32)*G(mu31)*G(mu11)*G(mu52)*G(mu51)*G(mu12)*G(mu61)*G(mu21)*G(mu41) + 16*G(mu71)*G(mu72)*G(mu32)*G(mu31)*G(mu12)*G(mu51)*G(mu52)*G(mu11)*G(mu61)*G(mu21)*G(mu41)

    # Fully Lorentz-contracted expressions, these return scalars:

    t = (G(mu)*G(-mu))
    assert kahane_simplify(t) == D

    t = (G(mu)*G(nu)*G(-mu)*G(-nu))
    assert kahane_simplify(t) == 2*D - D**2  # -8

    t = (G(mu)*G(nu)*G(-nu)*G(-mu))
    assert kahane_simplify(t) == D**2  # 16

    t = (G(mu)*G(nu)*G(-rho)*G(-nu)*G(-mu)*G(rho))
    assert kahane_simplify(t) == 4*D - 4*D**2 + D**3  # 16

    t = (G(mu)*G(nu)*G(rho)*G(-rho)*G(-nu)*G(-mu))
    assert kahane_simplify(t) == D**3  # 64

    t = (G(a1)*G(a2)*G(a3)*G(a4)*G(-a3)*G(-a1)*G(-a2)*G(-a4))
    assert kahane_simplify(t) == -8*D + 16*D**2 - 8*D**3 + D**4  # -32

    t = (G(-mu)*G(-nu)*G(-rho)*G(-sigma)*G(nu)*G(mu)*G(sigma)*G(rho))
    assert kahane_simplify(t) == -16*D + 24*D**2 - 8*D**3 + D**4  # 64

    t = (G(-mu)*G(nu)*G(-rho)*G(sigma)*G(rho)*G(-nu)*G(mu)*G(-sigma))
    assert kahane_simplify(t) == 8*D - 12*D**2 + 6*D**3 - D**4  # -32

    t = (G(a1)*G(a2)*G(a3)*G(a4)*G(a5)*G(-a3)*G(-a2)*G(-a1)*G(-a5)*G(-a4))
    assert kahane_simplify(t) == 64*D - 112*D**2 + 60*D**3 - 12*D**4 + D**5  # 256

    t = (G(a1)*G(a2)*G(a3)*G(a4)*G(a5)*G(-a3)*G(-a1)*G(-a2)*G(-a4)*G(-a5))
    assert kahane_simplify(t) == 64*D - 120*D**2 + 72*D**3 - 16*D**4 + D**5  # -128

    t = (G(a1)*G(a2)*G(a3)*G(a4)*G(a5)*G(a6)*G(-a3)*G(-a2)*G(-a1)*G(-a6)*G(-a5)*G(-a4))
    assert kahane_simplify(t) == 416*D - 816*D**2 + 528*D**3 - 144*D**4 + 18*D**5 - D**6  # -128

    t = (G(a1)*G(a2)*G(a3)*G(a4)*G(a5)*G(a6)*G(-a2)*G(-a3)*G(-a1)*G(-a6)*G(-a4)*G(-a5))
    assert kahane_simplify(t) == 416*D - 848*D**2 + 584*D**3 - 172*D**4 + 22*D**5 - D**6  # -128

    # Expressions with free indices:

    t = (G(mu)*G(nu)*G(rho)*G(sigma)*G(-mu))
    assert kahane_simplify(t).equals(-2*G(sigma)*G(rho)*G(nu) + (4-D)*G(nu)*G(rho)*G(sigma))

    t = (G(mu)*G(nu)*G(-mu))
    assert kahane_simplify (t) == (2-D)*G(nu)

    t = (G(mu)*G(nu)*G(rho)*G(-mu))
    assert kahane_simplify(t) == 2*G(nu)*G(rho) + 2*G(rho)*G(nu) - (4-D)*G(nu)*G(rho)

    t = 2*G(m2)*G(m0)*G(m1)*G(-m0)*G(-m1)
    st = kahane_simplify(t)
    assert st == (D*(-2*D + 4))*G(m2)

    t = G(m2)*G(m0)*G(m1)*G(-m0)*G(-m2)
    st = kahane_simplify(t)
    assert st == ((-D + 2)**2)*G(m1)

    t = G(m0)*G(m1)*G(m2)*G(m3)*G(-m1)
    st = kahane_simplify(t)
    assert st == (D - 4)*G(m0)*G(m2)*G(m3) + 4*G(m0)*g(m2, m3)

    t = G(m0)*G(m1)*G(m2)*G(m3)*G(-m1)*G(-m0)
    st = kahane_simplify(t)
    assert st == ((D - 4)**2)*G(m2)*G(m3) + (8*D - 16)*g(m2, m3)

    t = G(m2)*G(m0)*G(m1)*G(-m2)*G(-m0)
    st = kahane_simplify(t)
    assert st == ((-D + 2)*(D - 4) + 4)*G(m1)

    t = G(m3)*G(m1)*G(m0)*G(m2)*G(-m3)*G(-m0)*G(-m2)
    st = kahane_simplify(t)
    assert st == (-4*D + (-D + 2)**2*(D - 4) + 8)*G(m1)

    t = 2*G(m0)*G(m1)*G(m2)*G(m3)*G(-m0)
    st = kahane_simplify(t)
    assert st.equals((-2*D + 8)*G(m1)*G(m2)*G(m3) - 4*G(m3)*G(m2)*G(m1))

    t = G(m5)*G(m0)*G(m1)*G(m4)*G(m2)*G(-m4)*G(m3)*G(-m0)
    st = kahane_simplify(t)
    assert st.equals(((-D + 2)*(-D + 4))*G(m5)*G(m1)*G(m2)*G(m3) + (2*D - 4)*G(m5)*G(m3)*G(m2)*G(m1))

    t = -G(m0)*G(m1)*G(m2)*G(m3)*G(-m0)*G(m4)
    st = kahane_simplify(t)
    assert st.equals((D - 4)*G(m1)*G(m2)*G(m3)*G(m4) + 2*G(m3)*G(m2)*G(m1)*G(m4))

    t = G(-m5)*G(m0)*G(m1)*G(m2)*G(m3)*G(m4)*G(-m0)*G(m5)
    st = kahane_simplify(t)

    result1 = ((-D + 4)**2 + 4)*G(m1)*G(m2)*G(m3)*G(m4) +\
        (4*D - 16)*G(m3)*G(m2)*G(m1)*G(m4) + (4*D - 16)*G(m4)*G(m1)*G(m2)*G(m3)\
        + 4*G(m2)*G(m1)*G(m4)*G(m3) + 4*G(m3)*G(m4)*G(m1)*G(m2) +\
        4*G(m4)*G(m3)*G(m2)*G(m1)

    # Kahane's algorithm yields this result, which is equivalent to `result1`
    # in four dimensions, but is not automatically recognized as equal:
    result2 = 8*G(m1)*G(m2)*G(m3)*G(m4) + 8*G(m4)*G(m3)*G(m2)*G(m1)

    if D == 4:
        assert st.equals(result1) or st.equals(result2)
    else:
        assert st.equals(result1)

    # and a few very simple cases, with no contracted indices:

    t = G(m0)
    st = kahane_simplify(t)
    assert st == t

    t = -7*G(m0)
    st = kahane_simplify(t)
    assert st == t

    t = 224*G(m0)*G(m1)*G(-m2)*G(m3)
    st = kahane_simplify(t)
    assert st == t
