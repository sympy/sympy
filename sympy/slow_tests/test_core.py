from sympy import integer_nthroot

def test_integer_nthroot():
    assert integer_nthroot(10**(500*500), 500) == (10**500, True)
    assert integer_nthroot(10**1000000, 100000) == (10**10, True)
