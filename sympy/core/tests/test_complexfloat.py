from sympy import (ComplexFloat, Number, Rational, Symbol, Float, I,
                   S, sqrt, Pow)
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


@XFAIL
def test_ComplexFloat_arithmetic_fail():
    u = ComplexFloat(2 + 3j)
    assert u + 2*I == ComplexFloat(2, 5)
    assert u - 2*I == ComplexFloat(2, 1)


def test_ComplexFloat_ops():
    u = S(2+3j)
    v = sqrt(u)
    assert isinstance(v, ComplexFloat)
    # TODO: assert v**2 - u smaller than eps
    # TODO: Abs


def test_ComplexFloat_printing():
    from sympy import srepr, pretty, latex
    z = ComplexFloat('1.25', '2.5')
    assert srepr(z) == \
        "ComplexFloat(Float('1.25', prec=15), Float('2.5', prec=15))"
    assert str(z) == "1.25 + 2.5j"
    assert pretty(z) == "1.25 + 2.5ⅈ"
    assert latex(z) == "1.25 + 2.5i"


@XFAIL
def test_ComplexFloat_printing_fails():
    assert latex(Pow(z, 2, evaluate=False)) == "\left(1.25 + 2.5i\right)^{2}"


@XFAIL
def test_ComplexFloat_printing_coeff():
    x = S('x')
    a = S(1.0 + 3.0j)
    assert str(a*x) == '(1.0 + 3.0j)*x'


@XFAIL
def test_ComplexFloat_printing_power():
    a = S(1.0 + 3.0j)
    assert str(Pow(a,2, evaluate=False)) == '(1.0 + 3.0j)**2'


def test_Float_greedy_make_ComplexFloat():
    f = Float(6)
    assert isinstance(f*I, ComplexFloat)
    assert isinstance(f/I, ComplexFloat)
    assert isinstance(f + I, ComplexFloat)
    assert isinstance(f - I, ComplexFloat)


@XFAIL
def test_Float_greedy_make_ComplexFloat_fails():
    f = Float(6)
    assert isinstance(I*f, ComplexFloat)
    assert isinstance(I/f, ComplexFloat)
    assert isinstance(I + f, ComplexFloat)
    assert isinstance(I - f, ComplexFloat)


@XFAIL
def test_I_ops_greedy_make_ComplexFloat():
    assert isinstance(4.2*I, ComplexFloat)
    assert isinstance(I*4.2, ComplexFloat)
    assert isinstance(4.2/I, ComplexFloat)
    assert isinstance(I/4.2, ComplexFloat)
    assert isinstance(4.2 + I, ComplexFloat)
    assert isinstance(I + 4.2, ComplexFloat)
    assert isinstance(4.2 - I, ComplexFloat)
    assert isinstance(I - 4.2, ComplexFloat)
