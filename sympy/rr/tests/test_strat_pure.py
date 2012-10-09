from sympy.rr.strat_pure import null_safe

def test_null_safe():
    def rl(expr):
        if expr == 1:
            return 2
    safe_rl = null_safe(rl)
    assert rl(1) == safe_rl(1)

    assert      rl(3) == None
    assert safe_rl(3) == 3
