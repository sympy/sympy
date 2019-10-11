from sympy.plotting.intervalmath.interval_membership import intervalMembership
from sympy.utilities.pytest import raises


def test_creation():
    assert intervalMembership(True, True)
    raises(TypeError, lambda: intervalMembership(True))
    raises(TypeError, lambda: intervalMembership(True, True, True))


def test_getitem():
    a = intervalMembership(True, False)
    assert a[0] is True
    assert a[1] is False
    raises(IndexError, lambda: a[2])


def test_str():
    a = intervalMembership(True, False)
    assert str(a) == 'intervalMembership(True, False)'
    assert repr(a) == 'intervalMembership(True, False)'


def test_equivalence():
    a = intervalMembership(True, True)
    b = intervalMembership(True, False)
    assert (a == b) is False
    assert (a != b) is True

    a = intervalMembership(True, False)
    b = intervalMembership(True, False)
    assert (a == b) is True
    assert (a != b) is False


def test_boolean():
    # There can be 9*9 test cases in full mapping of the cartesian product.
    # But we only consider 3*3 cases for simplicity.
    s = [
        intervalMembership(False, False),
        intervalMembership(None, None),
        intervalMembership(True, True)
    ]

    # Reduced tests for 'And'
    a1 = [
        intervalMembership(False, False),
        intervalMembership(False, False),
        intervalMembership(False, False),
        intervalMembership(False, False),
        intervalMembership(None, None),
        intervalMembership(None, None),
        intervalMembership(False, False),
        intervalMembership(None, None),
        intervalMembership(True, True)
    ]
    a1_iter = iter(a1)
    for i in range(len(s)):
        for j in range(len(s)):
            assert s[i] & s[j] == next(a1_iter)

    # Reduced tests for 'Or'
    a1 = [
        intervalMembership(False, False),
        intervalMembership(None, None),
        intervalMembership(True, True),
        intervalMembership(None, None),
        intervalMembership(None, None),
        intervalMembership(True, True),
        intervalMembership(True, True),
        intervalMembership(True, True),
        intervalMembership(True, True)
    ]
    a1_iter = iter(a1)
    for i in range(len(s)):
        for j in range(len(s)):
            assert s[i] | s[j] == next(a1_iter)

    # Reduced tests for 'Xor'
    a1 = [
        intervalMembership(False, False),
        intervalMembership(None, None),
        intervalMembership(True, True),
        intervalMembership(None, None),
        intervalMembership(None, None),
        intervalMembership(None, None),
        intervalMembership(True, True),
        intervalMembership(None, None),
        intervalMembership(False, False)
    ]
    a1_iter = iter(a1)
    for i in range(len(s)):
        for j in range(len(s)):
            assert s[i] ^ s[j] == next(a1_iter)

    # Reduced tests for 'Not'
    a1 = [
        intervalMembership(True, True),
        intervalMembership(None, None),
        intervalMembership(False, False)
    ]
    a1_iter = iter(a1)
    for i in range(len(s)):
        assert ~s[i] == next(a1_iter)
