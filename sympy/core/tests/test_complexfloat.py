from sympy import (ComplexFloat, Number, Rational, Symbol, Float, I, S, Pow)
from sympy.utilities.pytest import XFAIL, raises
import mpmath

t = Symbol('t', real=False)

def same_and_same_prec(a, b):
    # stricter matching for Floats
    return a == b and a._prec == b._prec


def test_ComplexFloat():
    z = ComplexFloat(1 + 3j)
    assert isinstance(z, ComplexFloat)
    z = ComplexFloat(z)
    assert isinstance(z, ComplexFloat)
    z = ComplexFloat(1, 0)
    assert isinstance(z, ComplexFloat)
    z = ComplexFloat(0, 1)
    assert isinstance(z, ComplexFloat)
    z = S(1.0 + 3.0j)
    assert isinstance(z, ComplexFloat)
    z = Number(1.0 + 3.0j)
    assert isinstance(z, ComplexFloat)


def test_evalf_makes_ComplexFloat():
    from sympy import erf
    z = erf(S(1)+S(8)*I).n()
    assert isinstance(z, ComplexFloat)
    z = ComplexFloat(S(1)/3, S(2)/3)
    a = z.n(3)
    assert isinstance(a, ComplexFloat)
    #assert str(a) == '0.333 + 0.667*I'
    assert str(a) == '0.333 + 0.667j'


def test_ComplexFloat_from_mpc():
    a = Float(1.2)
    b = Float(3.4)
    zmpc = mpmath.mpc(a._as_mpf_val(a._prec), b._as_mpf_val(b._prec))
    z = ComplexFloat(zmpc)
    assert isinstance(z, ComplexFloat)
    z = S(zmpc)
    assert isinstance(z, ComplexFloat)


def test_ComplexFloat_re_im():
    from sympy import re, im
    z = S(2.0 + 3.0j)
    assert re(z) == Float(2)
    assert im(z) == Float(3)
    assert complex(z) == 2.0 + 3.0j


def test_ComplexFloat_assumptions():
    assert ComplexFloat(1, 2).is_complex
    assert ComplexFloat(1, 0).is_complex
    assert ComplexFloat(1, 0).is_real
    assert not ComplexFloat(1, 0).is_imaginary
    assert ComplexFloat(0, 1).is_imaginary
    assert not ComplexFloat(0, 1).is_real
    assert ComplexFloat(0, 0).is_real
    assert not ComplexFloat(0, 0).is_imaginary
    assert ComplexFloat(0, 0).is_zero
    assert not ComplexFloat(0, 1).is_zero


def test_ComplexFloat_equality():
    assert ComplexFloat(1, 2) == ComplexFloat(1, 2)
    assert ComplexFloat(1, 2) != ComplexFloat(1, 3)
    assert ComplexFloat(1, 2) != ComplexFloat(2, 2)
    assert ComplexFloat(1, 2) != S(1)
    assert ComplexFloat(1, 2) != S(2)*I
    # TODO: see Float's behaviour
    #assert ComplexFloat(1, 2) != 1 + 2*I
    #assert ComplexFloat(Rational(1, 3), Rational(2, 3)) != \
    #    Rational(1, 3) + I*Rational(2, 3)


def test_ComplexFloat_arithmetic():
    u = ComplexFloat(2.0 + 3.0j)
    assert u*I == ComplexFloat(-3, 2)
    assert u/I == ComplexFloat(3, -2)
    assert u + I == ComplexFloat(2, 4)
    assert u - I == ComplexFloat(2, 2)
    assert -u == ComplexFloat(-2, -3)
    assert u*2 == ComplexFloat(4, 6)
    assert u/3 == ComplexFloat(Rational(2, 3), 1)
    assert u + 3 == ComplexFloat(5, 3)
    assert u - 3 == ComplexFloat(-1, 3)
    assert u + 2*I == ComplexFloat(2, 5)
    assert u - 2*I == ComplexFloat(2, 1)
    assert same_and_same_prec(u*2, ComplexFloat(4 + 6j))
    assert same_and_same_prec(2*u, u*2)
    v = ComplexFloat(5 + 7j)
    assert same_and_same_prec(u*v, ComplexFloat(-11, 29))
    assert same_and_same_prec(u + v, ComplexFloat(7, 10))
    assert same_and_same_prec(u - v, ComplexFloat(-3, -4))
    # TODO: u / v, u**v, same within eps
    r = Rational(2, 3)
    assert isinstance(r*u, ComplexFloat)
    assert isinstance(r + u, ComplexFloat)
    assert isinstance(r - u, ComplexFloat)
    assert isinstance(r/u, ComplexFloat)


def test_ComplexFloat_ops():
    from sympy import sqrt, conjugate, Abs
    u = S(2+3j)
    v = sqrt(u)
    assert isinstance(v, ComplexFloat)
    # TODO: assert v**2 - u smaller than eps
    assert conjugate(u) == S(2 - 3j)
    assert u.adjoint() == conjugate(u)
    assert isinstance(Abs(u), Float)
    assert Abs(S(3 + 4j)) == Float(5)


def test_ComplexFloat_printing():
    from sympy import srepr, pretty, latex
    z = ComplexFloat('1.25', '2.5')
    assert srepr(z) == \
        "ComplexFloat(Float('1.25', prec=15), Float('2.5', prec=15))"
    assert str(z) == "1.25 + 2.5j"
    assert pretty(z) == "1.25 + 2.5ⅈ"
    assert latex(z) == "1.25 + 2.5i"
    z = ComplexFloat('1.25', '-2.5')
    assert str(z) == "1.25 - 2.5j"
    assert pretty(z) == "1.25 - 2.5ⅈ"
    assert latex(z) == "1.25 - 2.5i"
    assert str(Pow(z, 2, evaluate=False)) == '(1.25 - 2.5j)**2'
    assert latex(Pow(z, 2, evaluate=False)) == r'\left(1.25 - 2.5i\right)^{2}'
    x = Symbol('x')
    assert str(z*x) == '(1.25 - 2.5j)*x'


def test_Float_greedy_make_ComplexFloat():
    f = Float(6)
    assert isinstance(f*I, ComplexFloat)
    assert isinstance(f/I, ComplexFloat)
    assert isinstance(f + I, ComplexFloat)
    assert isinstance(f - I, ComplexFloat)
    assert isinstance(I*f, ComplexFloat)
    assert isinstance(I/f, ComplexFloat)
    assert isinstance(f*(I + 2), ComplexFloat)
    assert isinstance((I + 2)*f, ComplexFloat)
    assert isinstance(I + f, ComplexFloat)
    assert isinstance(I - f, ComplexFloat)
    assert isinstance(f + 2*I, ComplexFloat)
    assert isinstance(f - 2*I, ComplexFloat)


def test_greedy_make_ComplexFloat_w_symbols():
    x = Symbol('x')
    a = 4.2 + x + 2*I
    b = 4.2 + (2*I + x)
    c = 4.2 + 2*I + x
    assert a == b == c == (x + ComplexFloat(4.2, 2))
    f = Float(4.2)
    a = f + x + 2*I
    b = f + (2*I + x)
    c = f + 2*I + x
    assert a == b == c == (x + ComplexFloat(4.2, 2))
    f = ComplexFloat(1.1, 2.2)
    a = f + x + 2*I
    b = f + (2*I + x)
    c = f + 2*I + x
    assert a == b == c == (x + ComplexFloat(1.1, 4.2))


def test_I_ops_greedy_make_ComplexFloat():
    assert isinstance(1.0*I, ComplexFloat)
    assert isinstance(I*1.0, ComplexFloat)
    assert isinstance(1.0/I, ComplexFloat)
    assert isinstance(I/1.0, ComplexFloat)
    assert isinstance(4.2*I, ComplexFloat)
    assert isinstance(I*4.2, ComplexFloat)
    assert isinstance(4.2/I, ComplexFloat)
    assert isinstance(I/4.2, ComplexFloat)
    assert isinstance(4.2 + I, ComplexFloat)
    assert isinstance(I + 4.2, ComplexFloat)
    assert isinstance(4.2 - I, ComplexFloat)
    assert isinstance(I - 4.2, ComplexFloat)
    assert isinstance(4.2 + 2*I, ComplexFloat)
    assert isinstance(2*I + 4.2, ComplexFloat)
    assert isinstance(4.2 - 2*I, ComplexFloat)
    assert isinstance(2*I - 4.2, ComplexFloat)


def test_Mul_combines_I_with_Float_and_ComplexFloat():
    from sympy import Mul
    assert Mul(S(2.0), I, S(3)) == S(6.0j)
    assert Mul(S(2), I, S(3.0)) == S(6.0j)
    assert Mul(S(2.0), -I, S(3)) == S(-6.0j)
    assert Mul(S(2), -I, S(3.0)) == S(-6.0j)
    assert Mul(S(2), S(3), I, S(2.0)) == S(12.0j)
    assert Mul(S(2), S(2.0 + 3.0j), I, S(5.0)) == S(-30.0 + 20.0j)


def test_ComplexFloat_conversions():
    raises(TypeError, lambda: float(S(1.0 + 2.0j)))
