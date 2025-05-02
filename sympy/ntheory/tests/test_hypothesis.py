from hypothesis import given
from hypothesis import strategies as st
from sympy import divisors
from sympy.functions.combinatorial.numbers import divisor_sigma, totient
from sympy.ntheory.primetest import is_square
from sympy.ntheory import factorint, digits

@st.composite
def digit_strategy(draw):
    n = draw(st.integers())
    b = draw(st.integers(min_value=2))
    digits_value = draw(st.integers(min_value=len(digits(n,b))-1, max_value=len(digits(n, b))+5))
    return (n, b, digits_value)

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


@given(digit_strategy())
def test_digits(digits_tuple):
    n, b, digits_value = digits_tuple
    digits_list = digits(n, b, digits=digits_value)
    size = len(digits_list)-1
    x = 0
    for i in range(0, size):
        x = x + digits_list[size-i] * b**i
    assert x == abs(n)
    if n < 0:
        assert digits_list[0] == -b
