from sympy import symbols, sin, exp, cos

x, y, z = symbols('xyz')

def test_count_ops_non_symbolic():
    assert x.count_ops(symbolic=True) == 0
    assert y.count_ops(symbolic=True) == 0
    assert (x+1).count_ops(symbolic=False) == 1
    assert (y+x+1).count_ops(symbolic=False) == 2
    assert (z+y+x+1).count_ops(symbolic=False) == 3
    assert (2*z+y+x+1).count_ops(symbolic=False) == 4
    assert (2*z+y**17+x+1).count_ops(symbolic=False) == 5
    assert (2*z+y**17+x+sin(x)).count_ops(symbolic=False) == 6
    assert (2*z+y**17+x+sin(x**2)).count_ops(symbolic=False) == 7
    assert (2*z+y**17+x+sin(x**2)+exp(cos(x))).count_ops(symbolic=False) == 10
