from sympy.rules.branching_strat_pure import (exhaust, debug)


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
