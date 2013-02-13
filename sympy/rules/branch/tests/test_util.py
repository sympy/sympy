from sympy.rules.branch.util import unique, interleave, fnmap

def test_unique():
    assert tuple(unique((1,2,3))) == (1,2,3)
    assert tuple(unique((1,2,1,3))) == (1,2,3)

def test_interleave():
    assert ''.join(interleave(('ABC', '123'))) == 'A1B2C3'
    assert ''.join(interleave(('ABC', '1'))) == 'A1BC'

def test_fnmap():
    inc = lambda x: x + 1
    dec = lambda x: x - 1
    assert tuple(fnmap((inc, dec), 2)) == (3, 1)
