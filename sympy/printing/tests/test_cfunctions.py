from sympy.printing.cfunctions import expm1

def test_expm1():
    assert expm1(x).expand(func=True) - sp.exp(x) == -1
    assert expm1(x).rewrite('tractable') - sp.exp(x) == -1
    assert not ((sp.exp(1e-10).evalf() - 1) - 1e-10 - 5e-21) < 1e-22
    assert abs(expm1(1e-10).evalf() - 1e-10 - 5e-21) < 1e-22
    assert expm1(0) == 0
    assert expm1(x).is_real and expm1(x).is_finite
