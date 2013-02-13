from sympy.rules.util import count, RuleDB

def test_count():
    assert count((1, 1, 2, 5, 5, 6)) == {1: 2, 2: 1, 5: 2, 6: 1}
    assert count('a') == {'a': 1}

def test_RuleDB():
    db = RuleDB()
    db.insert('a', 1)
    db.insert('b', 2)
    db.insert('aa', 3)
    db.insert('aab', 4)
    db.insert('c', 4)
    assert set(db.query('a')) == set([1])
    assert set(db.query('ab')) == set([1, 2])
    assert set(db.query('aab')) == set([1, 2, 3, 4])
