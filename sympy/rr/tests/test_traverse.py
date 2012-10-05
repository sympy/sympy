from sympy.rr.rl import rmid
from sympy.rr.traverse import top_down, bottom_up
from sympy import Basic, symbols, Symbol, S

def test_bottom_up():
    _test_global_traversal(bottom_up)

def test_top_down():
    _test_global_traversal(top_down)

def _test_global_traversal(trav):
    zero_symbols = lambda x: S.Zero if isinstance(x, Symbol) else x
    zero_all_symbols = trav(zero_symbols)
    x,y,z = symbols('x,y,z')

    assert zero_all_symbols(Basic(x, y, Basic(x, z))) == \
                            Basic(0, 0, Basic(0, 0))
