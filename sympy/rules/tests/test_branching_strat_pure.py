from sympy.rules.branching_strat_pure import (exhaust)

def branch5(x):
    if 0 < x < 5:
        yield x-1
    elif 5 < x < 10:
        yield x+1
    elif x == 5:
        yield x+1
        yield x-1
    else:
        yield x
def test_exhaust():
    brl = exhaust(branch5)
    assert set(brl(3)) == {0}
    assert set(brl(7)) == {10}
    assert set(brl(5)) == {0, 10}

