from sympy.utilities.pytest import raises

def test_raises():
    class My(Exception):
        pass
    class My2(Exception):
        pass
    raises(My, "raise My()")

    try:
        raises(My, "1+1")
        assert False
    except Exception, e:
        assert str(e) == "DID NOT RAISE"

    try:
        raises(My, "raise My2('my text123')")
        assert False
    except My2, e:
        assert str(e) == "my text123"
