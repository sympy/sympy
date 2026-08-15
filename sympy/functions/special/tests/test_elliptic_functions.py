from __future__ import annotations

from sympy.core import S
from sympy.core.numbers import Float, I, Rational, pi
from sympy.core.symbol import Symbol
from sympy.functions.elementary.exponential import exp
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import cos, sin
from sympy.functions.special.elliptic_functions import jtheta
from sympy.series.order import O
from sympy.testing.pytest import raises


n = Symbol('n', integer=True)
q = Symbol('q')
z = Symbol('z')
d = Symbol('d', integer=True, nonnegative=True)


def _assert_close(left, right, tolerance=S(10)**-25):
    assert abs((left - right).evalf(30)) < tolerance


def test_jtheta_definition():
    assert jtheta(n, z, q).args == (n, z, q)
    assert jtheta(n, z, q, d).args == (n, z, q, d)
    assert jtheta(n, z, q, 0) == jtheta(n, z, q)

    raises(ValueError, lambda: jtheta(0, z, q))
    raises(ValueError, lambda: jtheta(5, z, q))
    raises(ValueError, lambda: jtheta(Rational(3, 2), z, q))
    raises(ValueError, lambda: jtheta(Float(1), z, q))
    raises(ValueError, lambda: jtheta(1, z, q, -1))
    raises(ValueError, lambda: jtheta(1, z, q, Rational(1, 2)))
    raises(ValueError, lambda: jtheta(1, z, q, Float(1)))

    # Do not reject an expression that might represent a valid index merely
    # because it is numeric and has not simplified to an Integer.
    unsimplified_three = ((1 + sqrt(2))**2 + (1 - sqrt(2))**2)/2
    assert jtheta(unsimplified_three, z, q).args[0] == unsimplified_three


def test_jtheta_special_values():
    assert jtheta(1, z, 0) == 0
    assert jtheta(2, z, 0) == 0
    assert jtheta(3, z, 0) == 1
    assert jtheta(4, z, 0) == 1
    assert jtheta(3, z, 0, 1) == 0
    assert jtheta(4, z, 0, 2) == 0

    assert jtheta(1, 0, q) == 0
    assert jtheta(1, 0, q, 2) == 0
    assert jtheta(2, 0, q, 1) == 0
    assert jtheta(3, 0, q, 3) == 0
    assert jtheta(4, 0, q, 1) == 0


def test_jtheta_parity():
    assert jtheta(1, -z, q) == -jtheta(1, z, q)
    assert jtheta(1, -z, q, 1) == jtheta(1, z, q, 1)
    assert jtheta(1, -z, q, 2) == -jtheta(1, z, q, 2)

    for index in (2, 3, 4):
        assert jtheta(index, -z, q) == jtheta(index, z, q)
        assert jtheta(index, -z, q, 1) == -jtheta(index, z, q, 1)
        assert jtheta(index, -z, q, 2) == jtheta(index, z, q, 2)


def test_jtheta_diff():
    assert jtheta(n, z, q).diff(z) == jtheta(n, z, q, 1)
    assert jtheta(n, z, q, d).diff(z) == jtheta(n, z, q, d + 1)
    assert jtheta(n, z, q).diff(q) == -jtheta(n, z, q, 2)/(4*q)
    assert jtheta(n, z, q, d).diff(q) == \
        -jtheta(n, z, q, d + 2)/(4*q)


def test_jtheta_nseries():
    assert jtheta(1, z, q).series(q, 0, 3) == (
        2*q**Rational(1, 4)*sin(z)
        - 2*q**Rational(9, 4)*sin(3*z) + O(q**3))
    assert jtheta(2, z, q).series(q, 0, 3) == (
        2*q**Rational(1, 4)*cos(z)
        + 2*q**Rational(9, 4)*cos(3*z) + O(q**3))
    assert jtheta(3, z, q).series(q, 0, 3) == (
        1 + 2*q*cos(2*z) + O(q**3))
    assert jtheta(4, z, q).series(q, 0, 3) == (
        1 - 2*q*cos(2*z) + O(q**3))

    assert jtheta(1, z, q, 1).series(q, 0, 3) == (
        2*q**Rational(1, 4)*cos(z)
        - 6*q**Rational(9, 4)*cos(3*z) + O(q**3))

    x = Symbol('x', positive=True)
    assert jtheta(3, z, x**2).series(x, 0, 7) == (
        1 + 2*x**2*cos(2*z) + O(x**7))


def test_jtheta_limits_at_zero():
    from sympy.series.limits import limit

    assert [limit(jtheta(index, z, q), q, 0)
            for index in (1, 2, 3, 4)] == [0, 0, 1, 1]
    assert limit(jtheta(3, z, q).diff(q), q, 0) == 2*cos(2*z)
    assert limit(jtheta(4, z, q).diff(q), q, 0) == -2*cos(2*z)
    assert limit(jtheta(1, pi/2, q).diff(q), q, 0, dir='+') is S.Infinity
    assert limit(jtheta(2, 0, q).diff(q), q, 0, dir='+') is S.Infinity


def test_jtheta_evalf():
    # Reference values from the mpmath jtheta documentation:
    # https://mpmath.org/doc/current/functions/elliptic.html#jtheta
    value = jtheta(1, Rational(1, 4), Rational(1, 5)).evalf(25)
    expected = Float('0.2945120798627300045053104', 25)
    assert abs(value - expected) < S(10)**-24

    derivative = jtheta(1, 7, Rational(1, 4), 2).evalf(25)
    expected = Float('-0.2598718791650217206533052', 25)
    assert abs(derivative - expected) < S(10)**-24

    value = jtheta(4, 1 + 2*I, (1 + I)/5).evalf(25)
    expected = (Float('7.180331760146805926356634', 25)
                - Float('1.634292858119162417301683', 25)*I)
    assert abs(value - expected) < S(10)**-24

    # These values were independently computed with Mathematica and are also
    # used by mpmath's test_jtheta:
    # https://github.com/mpmath/mpmath/blob/c90e2427418ea9e69674e957394c3ef2bf3ac0bc/mpmath/tests/test_elliptic.py
    references = (
        '0.1069552990104042681962096',
        '1.101385760258855791140606',
        '1.178319743354331061795905',
        '0.8219318954665153577314573',
    )
    for index, reference in enumerate(references, 1):
        value = jtheta(index, Rational(1, 10), Rational(1, 11)).evalf(25)
        assert abs(value - Float(reference, 25)) < S(10)**-24


def test_jtheta_quasiperiodicity():
    # DLMF 20.2.6--20.2.9.
    tau = 2*I/5
    nome = exp(I*pi*tau)
    argument = Rational(1, 3) + I/7
    m, lattice_n = 2, -1
    translated = argument + (m + lattice_n*tau)*pi
    common_factor = (nome**(-lattice_n**2)
                     * exp(-2*I*lattice_n*argument))
    signs = (
        (-1)**(m + lattice_n),
        (-1)**m,
        1,
        (-1)**lattice_n,
    )

    for index, sign in enumerate(signs, 1):
        _assert_close(
            jtheta(index, translated, nome),
            sign*common_factor*jtheta(index, argument, nome),
        )


def test_jtheta_half_period_translations():
    # DLMF 20.2.10--20.2.14.
    tau = 2*I/5
    nome = exp(I*pi*tau)
    argument = Rational(1, 3) + I/7
    half_pi = pi/2
    half_pi_tau = pi*tau/2
    multiplier = exp(I*argument + I*pi*tau/4)

    translations = (
        (jtheta(1, argument, nome), (
            -jtheta(2, argument + half_pi, nome),
            -I*multiplier*jtheta(4, argument + half_pi_tau, nome),
            -I*multiplier*jtheta(
                3, argument + half_pi + half_pi_tau, nome),
        )),
        (jtheta(2, argument, nome), (
            jtheta(1, argument + half_pi, nome),
            multiplier*jtheta(3, argument + half_pi_tau, nome),
            multiplier*jtheta(
                4, argument + half_pi + half_pi_tau, nome),
        )),
        (jtheta(3, argument, nome), (
            jtheta(4, argument + half_pi, nome),
            multiplier*jtheta(2, argument + half_pi_tau, nome),
            multiplier*jtheta(
                1, argument + half_pi + half_pi_tau, nome),
        )),
        (jtheta(4, argument, nome), (
            jtheta(3, argument + half_pi, nome),
            -I*multiplier*jtheta(1, argument + half_pi_tau, nome),
            I*multiplier*jtheta(
                2, argument + half_pi + half_pi_tau, nome),
        )),
    )

    for expected, values in translations:
        for value in values:
            _assert_close(value, expected)


def test_jtheta_sum_of_squares():
    # DLMF 20.7.1--20.7.5.
    nome = Rational(1, 3) + I/7
    argument = Rational(1, 5) + I/8
    theta = {index: jtheta(index, argument, nome)
             for index in (1, 2, 3, 4)}
    theta_zero = {index: jtheta(index, 0, nome)
                  for index in (2, 3, 4)}

    identities = (
        (theta_zero[3]**2*theta[3]**2,
         theta_zero[4]**2*theta[4]**2 + theta_zero[2]**2*theta[2]**2),
        (theta_zero[3]**2*theta[4]**2,
         theta_zero[2]**2*theta[1]**2 + theta_zero[4]**2*theta[3]**2),
        (theta_zero[2]**2*theta[4]**2,
         theta_zero[3]**2*theta[1]**2 + theta_zero[4]**2*theta[2]**2),
        (theta_zero[2]**2*theta[3]**2,
         theta_zero[4]**2*theta[1]**2 + theta_zero[3]**2*theta[2]**2),
        (theta_zero[3]**4, theta_zero[2]**4 + theta_zero[4]**4),
    )

    for left, right in identities:
        _assert_close(left, right)


def test_jtheta_derivative_identities():
    # Identities also used in mpmath's test_djtheta. They exercise derivative
    # orders 1--3.
    nome = Rational(1, 3)

    _assert_close(
        jtheta(1, 0, nome, 1),
        jtheta(2, 0, nome)*jtheta(3, 0, nome)*jtheta(4, 0, nome),
    )
    _assert_close(
        jtheta(1, 0, nome, 3)/jtheta(1, 0, nome, 1),
        sum(jtheta(index, 0, nome, 2)/jtheta(index, 0, nome)
            for index in (2, 3, 4)),
    )
