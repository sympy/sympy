from sympy import var, sqrt, exp, simplify, S, integrate, oo
from sympy.physics.hydrogen import R_nl

var("r Z")

def test_wavefunction():
    a = 1/Z
    R = {
            (1, 0): 2*sqrt(1/a**3) * exp(-r/a),
            (2, 0): sqrt(1/(2*a**3)) * exp(-r/(2*a)) * (1-r/(2*a)),
            (2, 1): S(1)/2 * sqrt(1/(6*a**3)) * exp(-r/(2*a)) * r/a,
            (3, 0): S(2)/3 * sqrt(1/(3*a**3)) * exp(-r/(3*a)) * \
                    (1-2*r/(3*a) + S(2)/27 * (r/a)**2),
            (3, 1): S(4)/27 * sqrt(2/(3*a**3)) * exp(-r/(3*a)) * \
                    (1-r/(6*a)) * r/a,
            (3, 2): S(2)/81 * sqrt(2/(15*a**3)) * exp(-r/(3*a)) * (r/a)**2,
            (4, 0): S(1)/4 * sqrt(1/a**3) * exp(-r/(4*a)) * \
                    (1-3*r/(4*a)+S(1)/8 * (r/a)**2-S(1)/192 * (r/a)**3),
            (4, 1): S(1)/16 * sqrt(5/(3*a**3)) * exp(-r/(4*a)) * \
                    (1-r/(4*a)+S(1)/80 * (r/a)**2) * (r/a),
            (4, 2): S(1)/64 * sqrt(1/(5*a**3)) * exp(-r/(4*a)) * \
                    (1-r/(12*a)) * (r/a)**2,
            (4, 3): S(1)/768 * sqrt(1/(35*a**3)) * exp(-r/(4*a)) * (r/a)**3,
            }
    for n, l in R:
        assert simplify(R_nl(n, l, r, Z) - R[(n, l)]) == 0

def test_norm():
    # Maximum "n" which is tested:
    n_max = 2
    # you can test any n and it works, but it's slow, so it's commented out:
    #n_max = 4
    for n in range(n_max+1):
        for l in range(n):
            assert integrate(R_nl(n, l, r)**2 * r**2, (r, 0, oo)) == 1
