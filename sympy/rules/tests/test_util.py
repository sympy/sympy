from sympy.rules.util import count, RuleDB

def test_count():
    assert count((1, 1, 2, 5, 5, 6)) == {1: 2, 2: 1, 5: 2, 6: 1}
    assert count('a') == {'a': 1}
