
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.symbol import Symbol
from sympy.core.numbers import Integer
from sympy.core.sympify import sympify
from sympy.core.basic import Basic, S, C, Atom

from sympy.polys.polynomial import Poly, PolynomialError
from sympy.polys.monomial import monomial_div

def poly_div(f, g, *symbols):
    """Generalized polynomial division with remainder.

       Given polynomial f and a set of polynomials g = (g_1, ..., g_n)
       compute a set of quotients q = (q_1, ..., q_n) and remainder r
       such that f = q_1*f_1 + ... + q_n*f_n + r, where r = 0 or r is
       a completely reduced polynomial with respect to g.

       In particular g can be a tuple, list or a singleton. All g_i
       and f can be given as Poly class instances or as expressions.

       For more information on the implemented algorithm refer to:

       [1] D. Cox, J. Little, D. O'Shea, Ideals, Varieties and
           Algorithms, Springer, Second Edition, 1997, pp. 62

       [2] I.A. Ajwa, Z. Liu, P.S. Wang, Groebner Bases Algorithm,
           http://citeseer.ist.psu.edu/ajwa95grbner.html, 1995

    """
    from sympy.simplify import cancel

    f = sympify(f)

    if not isinstance(f, Poly):
        if not symbols:
            symbols = f.atoms(Symbol)

        f = Poly(f, *symbols)
    elif symbols:
        f = Poly(f, *symbols, **f.flags)

    symbols, flags = f.symbols, f.flags

    r = Poly((), *symbols, **flags)

    if isinstance(g, (tuple, list)):
        q = [r] * len(g)
    else:
        g, q = [g], [r]

    g = [ Poly(h, *symbols, **flags) for h in g ]

    while not f.is_zero:
        for i, h in enumerate(g):
            M = monomial_div(f.LM, h.LM)

            if M is not None:
                # cancel() should be enough fast and general solution
                T = Poly((cancel(f.LC / h.LC), M), *symbols, **flags)

                q[i] = q[i] + T
                f = f - T*h
                break
        else:
            r = r.add_term(*f.LT)
            f = f.kill_lead_term()

    if len(q) != 1:
        return tuple(q), r
    else:
        return q[0], r

def poly_quo(f, g, *symbols):
    return poly_div(f, g, *symbols)[0]

def poly_rem(f, g, *symbols):
    return poly_div(f, g, *symbols)[1]
