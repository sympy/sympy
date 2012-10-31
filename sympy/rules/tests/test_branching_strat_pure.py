from sympy.rules.branching_strat_pure import (exhaust, debug, multiplex,
        condition, notempty)


def posdec(x):
    if x > 0:
        yield x-1
    else:
        yield x

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

even = lambda x: x%2 == 0

def ident_if_even(x):
    if even(x):
        yield x

def test_exhaust():
    brl = exhaust(branch5)
    assert set(brl(3)) == {0}
    assert set(brl(7)) == {10}
    assert set(brl(5)) == {0, 10}

def test_debug():
    import StringIO
    file = StringIO.StringIO()
    rl = debug(posdec, file)
    list(rl(5))
    log = file.getvalue()
    file.close()

    assert posdec.func_name in log
    assert '5' in log
    assert '4' in log

def test_multiplex():
    brl = multiplex(posdec, branch5)
    assert set(brl(3)) == {2}
    assert set(brl(7)) == {6, 8}
    assert set(brl(5)) == {4, 6}

def test_condition():
    brl = condition(even, branch5)
    assert set(brl(4)) == set(branch5(4))
    assert set(brl(5)) == set([])

def test_notempty():
    brl = notempty(ident_if_even)
    assert set(brl(4)) == {4}
    assert set(brl(5)) == {5}
