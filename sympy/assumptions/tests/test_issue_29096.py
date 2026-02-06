from sympy import Symbol, ask, Q

def test_issue_29096():
    x = Symbol('x')
    # Regression test for issue 29096
    # Allow LRASolver to handle generic symbols if context implies they are real
    assert ask(x > -1, Q.positive(x)) is True

    # Other examples that rely on this
    assert ask(x < 1, Q.negative(x)) is True

    y = Symbol('y')
    assert ask(x + y > 0, Q.positive(x) & Q.positive(y)) is True

    z = Symbol('z')
    assert ask(z > 5, Q.gt(z, 10)) is True
