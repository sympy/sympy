from sympy import Basic, symbols, Symbol, S
from sympy.rules.branch_traverse import top_down

def inc(x):
    if isinstance(x, int):
        yield x + 1

def test_top_down_easy():
    expr     = Basic(1, 2)
    expected = Basic(2, 3)
    brl = top_down(inc)

    assert set(brl(expr)) == {expected}

def test_top_down_big_tree():
    expr     = Basic(1, Basic(2), Basic(3, Basic(4), 5))
    expected = Basic(2, Basic(3), Basic(4, Basic(5), 6))
    brl = top_down(inc)

    assert set(brl(expr)) == {expected}

def test_top_down_harder_function():
    def split5(x):
        if x == 5:
            yield x - 1
            yield x + 1

    expr     = Basic(Basic(5, 6), 1)
    expected = {Basic(Basic(4, 6), 1), Basic(Basic(6, 6), 1)}
    brl = top_down(split5)

    assert set(brl(expr)) == expected
