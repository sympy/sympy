from hypothesis import given
from hypothesis import strategies as st
from sympy import divisors
from sympy.functions.combinatorial.numbers import divisor_sigma, totient
from sympy.ntheory.primetest import is_square
from sympy.ntheory import factorint, digits


@given(n=st.integers(1, 10**10))
def test_tau_hypothesis(n):
    div = divisors(n)
    tau_n = len(div)
    assert is_square(n) == (tau_n % 2 == 1)
    sigmas = [divisor_sigma(i) for i in div]
    totients = [totient(n // i) for i in div]
    mul = [a * b for a, b in zip(sigmas, totients)]
    assert n * tau_n == sum(mul)


@given(n=st.integers(1, 10**10))
def test_totient_hypothesis(n):
    assert totient(n) <= n
    div = divisors(n)
    totients = [totient(i) for i in div]
    assert n == sum(totients)


@given(n=st.integers())
def test_factorint(n):
    factors = factorint(n)
    product = 1
    for prime, exp in factors.items():
        product *= prime ** exp
    assert product == n


@given(n=st.integers(min_value=1), b=st.integers(min_value=2))
def test_digits(n, b):
    digits_list = digits(n, b)
    size = len(digits_list)-1
    x = 0
    for i in range(0, size):
        x = x + digits_list[size-i] * b**i
    assert x == n
