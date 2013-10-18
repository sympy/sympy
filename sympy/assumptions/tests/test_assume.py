from sympy import Predicate

def test_Predicate():
    # Test that nothing unusual happens with the Predicate.__new__ magic
    a = Predicate('a')
    a2 = Predicate('a')
    b = Predicate('b')
    assert a == a
    assert a == a2
    assert hash(a) == hash(a2)
    assert not (a == b)

    assert isinstance(a, Predicate)
    assert type(a) != Predicate
    assert a.__class__.__name__ == 'APredicate'

    class MyPredicate(Predicate):
        pass

    c = MyPredicate('c')
    c2 = MyPredicate('c')
    d = MyPredicate('d')
    assert c == c2
    assert not (c == d)
    assert isinstance(c, Predicate)
    assert isinstance(c, MyPredicate)
    assert type(c) != MyPredicate
    assert c.__class__.__name__ == 'CMyPredicate'

def test_Predicate_doc():
    doc = """
This is a predicate about a
"""

    a = Predicate('a', doc=doc)
    assert a.__class__.__doc__ == doc
