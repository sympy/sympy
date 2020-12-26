from sympy import var, integrate, simplify, expand, factor, log, Integral, \
    diff, FiniteSet, Matrix, exp
from ..relational import Equality
from ..equation import Equation, Eqn
from pytest import raises, warns


def test_define_equation():
    a, b, c = var('a b c')
    raises(TypeError, lambda: Eqn(FiniteSet(a), FiniteSet(b, c)))
    warns(UserWarning, lambda: Eqn(1,0))
    assert Eqn(1, 0, check=False) == Equation(1, 0, check=False)
    tsteqn = Eqn(a, b / c)
    assert tsteqn == Equation(a, b / c)
    assert tsteqn.lhs == Equation(a, b / c).args[0]
    assert tsteqn.rhs == Equation(a, b / c).args[1]
    assert tsteqn.free_symbols == {a, b, c}


def test_convert_equation():
    a, b, c = var('a b c')
    tsteqn = Eqn(a, b / c)
    assert tsteqn.as_Boolean() == Equality(a, b / c)
    assert tsteqn.reversed == Equation(b / c, a)


def test_binary_op():
    a, b, c = var('a b c')
    tsteqn = Eqn(a, b / c)
    assert tsteqn + c == Equation(a + c, b / c + c)
    assert c + tsteqn == Equation(c + a, c + b / c)
    assert tsteqn * c == Equation(a * c, b)
    assert c * tsteqn == Equation(c * a, b)
    assert tsteqn - c == Equation(a - c, b / c - c)
    assert c - tsteqn == Equation(c - a, c - b / c)
    assert tsteqn / c == Equation(a / c, b / c ** 2)
    assert c / tsteqn == Equation(c / a, c ** 2 / b)
    assert tsteqn % c == Equation(a % c, (b / c) % c)
    assert c % tsteqn == Equation(c % a, c % (b / c))
    assert tsteqn ** c == Equation(a ** c, (b / c) ** c)
    assert c ** tsteqn == Equation(c ** a, c ** (b / c))
    assert tsteqn + tsteqn == Equation(2 * a, 2*b / c)
    assert tsteqn * tsteqn == Equation(a ** 2, b ** 2 / c ** 2)
    assert tsteqn - tsteqn == Equation(0, 0)
    assert tsteqn / tsteqn == Equation(1, 1)
    assert tsteqn % tsteqn == Equation(0, 0)
    assert tsteqn ** tsteqn == Equation(a ** a, (b / c) ** (b / c))


def test_outputs():
    a, b, c = var('a b c')
    tsteqn = Eqn(a, b / c)
    assert tsteqn.__repr__() == 'a=b/c'
    assert tsteqn.__str__() == 'a=b/c'
    assert tsteqn._latex(tsteqn) == 'a=\\frac{b}{c}'


def test_helper_functions():
    a, b, c = var('a b c')
    tsteqn = Eqn(a, b / c)
    raises(ValueError, lambda: integrate(tsteqn, c))
    raises(AttributeError, lambda: integrate(tsteqn, c, side='right'))
    raises(ValueError, lambda: tsteqn._applyfunc(log))
    assert tsteqn.applyfunc(log) == Equation(log(a), log(b / c))
    assert tsteqn.applylhs(log) == Equation(log(a), b / c)
    assert tsteqn.applyrhs(log) == Equation(a, log(b / c))
    assert tsteqn.evalf(4, {b: 2.0, c: 4}) == Eqn(a, 0.5000)
    assert diff(tsteqn, c) == Equation(diff(a, c, evaluate=False), -b / c ** 2)
    tsteqn = Eqn(a * c, b / c)
    assert diff(tsteqn, c) == Equation(a, -b / c ** 2)
    assert integrate(tsteqn, c, side='rhs') == integrate(tsteqn.rhs, c)
    assert integrate(tsteqn, c, side='lhs') == integrate(tsteqn.lhs, c)
    assert (tsteqn.integ(c)).lhs == Integral(a * c, c)
    assert tsteqn.integ(c).rhs.simplify() == integrate(tsteqn.rhs, c)

    def adsq(expr):
        # Arbitrary python function
        return expr + expr ** 2

    assert (adsq(tsteqn) == Eqn(a ** 2 * c ** 2 + a * c, b ** 2 / c ** 2 +
                                b / c))
    tsteqn = Eqn((a - 1) * (a + 1) * (a + 3), (a + 3) * (2 * b + c) ** 2)
    tsteqn2 = Eqn(a ** 3 + 3 * a ** 2 - a - 3,
                  4 * a * b ** 2 + 4 * a * b * c + a * c ** 2 + 12 * b ** 2 +
                  12 * b * c + 3 * c ** 2)
    tsteqn3 = Eqn(a ** 3 + 3 * a ** 2 - a - 3,
                  a * c ** 2 + b ** 2 * (4 * a + 12) + b * (
                          4 * a * c + 12 * c) + 3 * c ** 2)
    tsteqn4 = Eqn(a ** 3 + 3 * a ** 2 - a - 3,
                  a * c ** 2 + 4 * b ** 2 * (a + 3) + 4 * b * c * (
                          a + 3) + 3 * c ** 2)
    assert tsteqn.expand() == tsteqn2
    assert expand(tsteqn) == tsteqn2
    assert tsteqn2.collect(b) == tsteqn3
    # collect(Eqn,...) not working yet
    #assert collect(tsteqn2,b) == tsteqn3
    assert tsteqn3.simplify() == tsteqn4
    assert simplify(tsteqn3) == tsteqn4
    assert tsteqn4.factor() == tsteqn
    assert factor(tsteqn4) == tsteqn
    # Check matrix exponentiation is not overridden.
    tsteqn5 = Eqn(a, Matrix([[1, 1], [1, 1]]))
    result5 = 'Eqn(exp(a), Matrix([\n[0.5 + exp(2) / 2, -0.5 + exp(2) / 2],\n['
    result5 +='-0.5 + exp(2) / 2, 0.5 + exp(2) / 2]]))'
    assert (exp(tsteqn5) == result5Eqn))
