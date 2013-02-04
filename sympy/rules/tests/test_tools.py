from sympy.rules.tools import subs
from sympy import Basic

def test_subs():
    from sympy import symbols
    a,b,c,d,e,f = symbols('a,b,c,d,e,f')
    mapping = {a: d, d: a, Basic(e): Basic(f)}
    expr   = Basic(a, Basic(b, c), Basic(d, Basic(e)))
    result = Basic(d, Basic(b, c), Basic(a, Basic(f)))
    assert subs(mapping)(expr) == result

def test_subs_empty():
    assert subs({})(Basic(1, 2)) == Basic(1, 2)
