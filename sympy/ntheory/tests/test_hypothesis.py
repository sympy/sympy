from __future__ import annotations
from hypothesis import given,settings
from hypothesis import strategies as st
from sympy.functions.combinatorial.numbers import divisor_sigma, totient
from sympy.ntheory.primetest import is_square
from sympy.ntheory.digits import is_palindromic
from sympy.ntheory.factor_ import is_perfect,divisor_count,divisors,proper_divisors
from sympy.ntheory.multinomial import binomial_coefficients
from sympy.functions.combinatorial.numbers import binomial


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


@given(st.integers())
def test_is_palindromic(n):
    assert is_palindromic(n)==is_palindromic(-n)
    if(n%10==0 and n!=0):
        assert is_palindromic(n) == False
    else:
        rev = int(str(abs(n))[::-1])
        assert is_palindromic(n)==is_palindromic(rev)


@given(st.integers(min_value=1))
def test_divisors(n):
    divisors_list = list(divisors(n))
    assert len(divisors_list) == divisor_count(n)
    for i in range(len(divisors_list)//2):
        assert divisors_list[i]*divisors_list[-(i+1)] == n

@given(st.integers(min_value=1))
def test_is_perfect(n):
    divisors_list = list(proper_divisors(n))
    if(is_perfect(n)):
        assert sum(divisors_list)==n
    else:
        assert sum(divisors_list)!=n


@given(st.integers(1,100))
@settings(max_examples=50)
def test_binomial_coefficients(n):
    coefficients = binomial_coefficients(n)
    for (k1,k2),value in coefficients.items():
        assert k1+k2 ==n
        assert (k1,k2) in coefficients
        assert coefficients[k1,k2] == value
        assert value == binomial(n,k1)

    assert sum(coefficients.values()) == 2**n
    assert coefficients[0,n]==1 and coefficients[n,0]==1

    if n>1:
        prev_coefficients = binomial_coefficients(n-1)
        for k in range(1,n):
            assert coefficients[k,n-k] == prev_coefficients[k-1,n-k] + prev_coefficients[k,n-k-1]
    

