from sympy import (
    I, N)

from sympy.functions import (
    bernoulli, bell, euler, polylog, zeta)

from sympy.functions.combinatorial.sequences import (
    andre_sequence, bell_sequence, bernoulli_sequence,
    euler_sequence, tangent_sequence)


def test_andre_sequence():
    a = andre_sequence(20)
    for n, b in enumerate(a):
        if n == 0:
            assert b == 1
        else:
            assert b == int(N(I ** (n + 1) * 2 * polylog(-n, -I)))


def test_bell_sequence():
    a = bell_sequence(20)
    for n, b in enumerate(a):
        assert b == bell(n)


def test_bernoulli_sequence():
    a = bernoulli_sequence(20)
    for n, b in enumerate(a):
        assert b == bernoulli(2 * n)


def test_euler_sequence():
    a = euler_sequence(20)
    for n, b in enumerate(a):
        assert b == abs(euler(2 * n))


def test_tangent_sequence():
    a = tangent_sequence(20)
    for n, b in enumerate(a):
        if n == 0:
            assert b == 1
        else:
            f = (-16) ** (n + 1) - (-4) ** (n + 1)
            assert b == f * zeta(- 1 - 2 * n)
