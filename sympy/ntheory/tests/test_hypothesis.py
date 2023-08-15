from hypothesis import given
from hypothesis import strategies as st
from sympy import divisors
from sympy.ntheory.primetest import is_square, isprime
from sympy.ntheory import totient, divisor_sigma, prime, primepi


@given(n=st.integers(1, 10**10))
def test_tau(n):
    n = 217
    div = divisors(n)
    tau_n = len(div)
    assert is_square(n) == (tau_n % 2 == 1)
    sigmas = [divisor_sigma(i) for i in div]
    totients = [totient(n // i) for i in div]
    mul = [a * b for a, b in zip(sigmas, totients)]
    assert n * tau_n == sum(mul)


@given(n=st.integers(1, 10**10))
def test_totient(n):
    assert totient(n) <= n
    div = divisors(n)
    totients = [totient(i) for i in div]
    assert n == sum(totients)


@given(n=st.integers(2, 10**5))
def test_primepi(n):
    nth_prime = prime(n)
    assert primepi(nth_prime) == n
    if isprime(n):
        assert prime(primepi(n)) == n
    else:
        assert prime(primepi(n)) != n
