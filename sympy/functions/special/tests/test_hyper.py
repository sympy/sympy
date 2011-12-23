from sympy import (hyper, meijerg, S, Tuple, pi, I, exp, log,
                   cos, sqrt, symbols, oo, Derivative)
from sympy.abc import x, z, k
from sympy.utilities.pytest import raises
from sympy.utilities.randtest import (
        random_complex_number as randcplx,
        test_numerically as tn,
        test_derivative_numerically as td)

def test_TupleParametersBase():
    # test that our implementation of the chain rule works
    p = hyper((), (), z**2)
    assert p.diff(z) == p*2*z

def test_hyper():
    raises(TypeError, 'hyper(1, 2, z)')

    assert hyper((1, 2),(1,), z) == hyper(Tuple(1, 2), Tuple(1), z)

    h = hyper((1, 2), (3, 4, 5), z)
    assert h.ap == Tuple(1, 2)
    assert h.bq == Tuple(3, 4, 5)
    assert h.argument == z
    assert h.is_commutative is True

    # just a few checks to make sure that all arguments go where they should
    assert tn(hyper(Tuple(), Tuple(), z), exp(z), z)
    assert tn(z*hyper((1, 1), Tuple(2), -z), log(1 + z), z)

    # differentiation
    h = hyper((randcplx(), randcplx(), randcplx()), (randcplx(), randcplx()), z)
    assert td(h, z)

    a1, a2, b1, b2, b3 = symbols('a1:3, b1:4')
    assert hyper((a1, a2), (b1, b2, b3), z).diff(z) == \
             a1*a2/(b1*b2*b3) * hyper((a1+1, a2+1), (b1+1, b2+1, b3+1), z)

    # differentiation wrt parameters is not supported
    assert hyper([z], [], z).diff(z) == Derivative(hyper([z], [], z), z)

    # hyper is unbranched wrt parameters
    from sympy import polar_lift
    assert hyper([polar_lift(z)], [polar_lift(k)], polar_lift(x)) == \
           hyper([z], [k], polar_lift(x))

def test_expand_func():
    # evaluation at 1 of Gauss' hypergeometric function:
    from sympy.abc import a, b, c
    from sympy import gamma, expand_func
    a1, b1, c1 = randcplx(), randcplx(), randcplx() + 5
    assert expand_func(hyper([a, b], [c], 1)) == \
           gamma(c)*gamma(-a - b + c)/(gamma(-a + c)*gamma(-b + c))
    assert abs(expand_func(hyper([a1, b1], [c1], 1)).n()
               - hyper([a1, b1], [c1], 1).n()) < 1e-10

    # hyperexpand wrapper for hyper:
    assert expand_func(hyper([], [], z)) == exp(z)
    assert expand_func(hyper([1, 2, 3], [], z)) == hyper([1, 2, 3], [], z)
    assert expand_func(meijerg([[1,1],[]], [[1],[0]], z)) == log(z + 1)
    assert expand_func(meijerg([[1,1],[]], [[],[]], z)) \
           == meijerg([[1,1],[]], [[],[]], z)

def test_radius_of_convergence():
    assert hyper((1, 2), [3], z).radius_of_convergence == 1
    assert hyper((1, 2), [3, 4], z).radius_of_convergence == oo
    assert hyper((1, 2, 3), [4], z).radius_of_convergence == 0
    assert hyper((0, 1, 2), [4], z).radius_of_convergence == oo
    assert hyper((-1, 1, 2), [-4], z).radius_of_convergence == 0
    assert hyper((-1, -2, 2), [-1], z).radius_of_convergence == oo
    assert hyper((-1, 2), [-1, -2], z).radius_of_convergence == 0
    assert hyper([-1, 1, 3], [-2, 2], z).radius_of_convergence == 1
    assert hyper([-1, 1], [-2, 2], z).radius_of_convergence == oo
    assert hyper([-1, 1, 3], [-2], z).radius_of_convergence == 0
    assert hyper((-1, 2, 3, 4), [], z).radius_of_convergence == oo

    assert hyper([1, 1], [3], 1).convergence_statement is True
    assert hyper([1, 1], [2], 1).convergence_statement is False
    assert hyper([1, 1], [2], -1).convergence_statement is True
    assert hyper([1, 1], [1], -1).convergence_statement is False


def test_meijer():
    raises(TypeError, 'meijerg(1, z)')
    raises(TypeError, 'meijerg(((1,), (2,)), (3,), (4,), z)')

    assert meijerg(((1, 2), (3,)), ((4,), (5,)), z) == \
           meijerg(Tuple(1, 2), Tuple(3), Tuple(4), Tuple(5), z)

    g = meijerg((1, 2), (3, 4, 5), (6, 7, 8, 9), (10, 11, 12, 13, 14), z)
    assert g.an == Tuple(1, 2)
    assert g.ap == Tuple(1, 2, 3, 4, 5)
    assert g.aother == Tuple(3, 4, 5)
    assert g.bm == Tuple(6, 7, 8, 9)
    assert g.bq == Tuple(6, 7, 8, 9, 10, 11, 12, 13, 14)
    assert g.bother == Tuple(10, 11, 12, 13, 14)
    assert g.argument == z
    assert g.nu == 75
    assert g.delta == -1
    assert g.is_commutative is True

    assert meijerg([1, 2], [3], [4], [5], z).delta == S(1)/2

    # just a few checks to make sure that all arguments go where they should
    assert tn(meijerg(Tuple(), Tuple(), Tuple(0), Tuple(), -z), exp(z), z)
    assert tn(sqrt(pi)*meijerg(Tuple(), Tuple(),
                               Tuple(0), Tuple(S(1)/2), z**2/4), cos(z), z)
    assert tn(meijerg(Tuple(1, 1),Tuple(), Tuple(1), Tuple(0), z),
              log(1 + z), z)

    # differentiation
    g = meijerg((randcplx(),), (randcplx() + 2*I,), Tuple(),
                (randcplx(), randcplx()), z)
    assert td(g, z)

    g = meijerg(Tuple(), (randcplx(),), Tuple(),
                (randcplx(), randcplx()), z)
    assert td(g, z)

    g = meijerg(Tuple(), Tuple(), Tuple(randcplx()),
                Tuple(randcplx(), randcplx()), z)
    assert td(g, z)

    a1, a2, b1, b2, c1, c2, d1, d2 = symbols('a1:3, b1:3, c1:3, d1:3')
    assert meijerg((a1, a2), (b1, b2), (c1, c2), (d1, d2), z).diff(z) == \
        (meijerg((a1-1, a2), (b1, b2), (c1, c2), (d1, d2), z) \
         + (a1 - 1)*meijerg((a1, a2), (b1, b2), (c1, c2), (d1, d2), z))/z

    assert meijerg([z, z], [], [], [], z).diff(z) == \
           Derivative(meijerg([z, z], [], [], [], z), z)

    # meijerg is unbranched wrt parameters
    from sympy import polar_lift as pl
    assert meijerg([pl(a1)], [pl(a2)], [pl(b1)], [pl(b2)], pl(z)) == \
           meijerg([a1], [a2], [b1], [b2], pl(z))

def test_meijerg_derivative():
    assert meijerg([], [1, 1], [0, 0, x], [], z).diff(x) == \
           log(z)*meijerg([], [1, 1], [0, 0, x], [], z) \
           + 2*meijerg([], [1, 1, 1], [0, 0, x, 0], [], z)

    y = randcplx()
    a = 5 # mpmath chokes with non-real numbers, and Mod1 with floats
    assert td(meijerg([x], [], [], [], y), x)
    assert td(meijerg([x**2], [], [], [], y), x)
    assert td(meijerg([], [x], [], [], y), x)
    assert td(meijerg([], [], [x], [], y), x)
    assert td(meijerg([], [], [], [x], y), x)
    assert td(meijerg([x], [a], [a + 1], [], y), x)
    assert td(meijerg([x], [a + 1], [a], [], y), x)
    assert td(meijerg([x, a], [], [], [a + 1], y), x)
    assert td(meijerg([x, a + 1], [], [], [a], y), x)
    b = S(3)/2
    assert td(meijerg([a + 2], [b], [b - 3, x], [a], y), x)

def test_meijerg_period():
    assert meijerg([], [1], [0], [], x).get_period() == 2*pi
    assert meijerg([1], [], [], [0], x).get_period() == 2*pi
    assert meijerg([], [], [0], [], x).get_period() == 2*pi # exp(x)
    assert meijerg([], [], [0], [S(1)/2], x).get_period() == 2*pi # cos(sqrt(x))
    assert meijerg([], [], [S(1)/2], [0], x).get_period() == 4*pi # sin(sqrt(x))
    assert meijerg([1, 1], [], [1], [0], x).get_period() == oo # log(1 + x)
