from hypothesis import given
from hypothesis import strategies as st
from sympy.ntheory import totient
from sympy import divisors
import random


@given(n=st.integers(1, 10**10))
def test_totient_hypothesis(n):
    assert totient(n) <= n
    div = divisors(n)
    totients = [totient(i) for i in div]
    assert n == sum(totients)
