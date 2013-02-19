from sympy.rules.traverse import top_down, bottom_up, sall
from sympy import Basic, symbols, Symbol, S

def test_sall():
    zero_symbols = lambda x: S.Zero if isinstance(x, Symbol) else x
    zero_onelevel = sall(zero_symbols)
    x,y,z = symbols('x,y,z')

    assert zero_onelevel(Basic(x, y, Basic(x, z))) == \
                            Basic(0, 0, Basic(x, z))


def test_bottom_up():
    _test_global_traversal(bottom_up)
    _test_stop_on_non_basics(bottom_up)

def test_top_down():
    _test_global_traversal(top_down)
    _test_stop_on_non_basics(top_down)

def _test_global_traversal(trav):
    zero_symbols = lambda x: S.Zero if isinstance(x, Symbol) else x
    zero_all_symbols = trav(zero_symbols)
    x,y,z = symbols('x,y,z')

    assert zero_all_symbols(Basic(x, y, Basic(x, z))) == \
                            Basic(0, 0, Basic(0, 0))

def _test_stop_on_non_basics(trav):
    def add_one_if_can(expr):
        try:    return expr + 1
        except: return expr

    expr     = Basic(1, 'a', Basic(2, 'b'))
    expected = Basic(2, 'a', Basic(3, 'b'))
    rl = trav(add_one_if_can)

    assert rl(expr) == expected
