"""Greatest common divisor for the Polynomial class"""

from sympy.modules.polynomials.base import *

def uv(f, g):
    """Euclidean algorithm for univariate polynomials.

    Coefficients assumed to be in a field.
    """
    from sympy.modules.polynomials import div_

    while True:
        if g.sympy_expr is not S.Zero:
            lc, g = g.as_monic()
            f, g = g, div_.mv(f, g)[-1]
        else:
            break
    return f

def uv_int(f, g):
    """For integer coefficients.
    """
    from sympy.modules.polynomials import div_

    cf, f = f.as_primitive()
    cg, g = g.as_primitive()
    c = Integer(numbers.gcd(int(cf.p), int(cg.p)))

    while True:
        if g.sympy_expr is not S.Zero:
            f, g = g, div_.mv(f, g)[-1]
        else:
            break
    return Polynomial(coeffs=tuple(map(lambda t:(t[0]*c) + t[1:], f.coeffs)),
                      var=f.var, order=f.order)
    
def mv(f, g):
    """Computes the gcd of 2 polynomials by dividing their product by their lcm.

    It is assumed that f and g are instances of the Polynomial class with
    matching variables and orders.
    """
    from sympy.modules.polynomials import lcm_
    from sympy.modules.polynomials import div_

    lcm = lcm_.mv(f, g)
    q, r = div_.mv(f*g, lcm)
    assert r.sympy_expr is S.Zero
    q = q[0] # q is a list!
    return q.as_monic()[1]
