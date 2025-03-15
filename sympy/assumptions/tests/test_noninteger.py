from sympy import symbols, Q, ask, Rational, Float, pi, E, I, oo, log

def test_noninteger():
    x = symbols('x')
    y = symbols('y', rational=True)
    z = symbols('z')

    # Integer cases
    # assert ask(Q.noninteger(x)) == False
    assert ask(Q.noninteger(5)) == False
    assert ask(Q.noninteger(-3)) == False
    assert ask(Q.noninteger(0)) == False

    # Non-integer cases
    assert ask(Q.noninteger(5.5)) == True
    assert ask(Q.noninteger(Float(5.5))) == True
    assert ask(Q.noninteger(Rational(11, 2))) == True
    assert ask(Q.noninteger(Rational(4, 2))) == False

    # Special constants
    assert ask(Q.noninteger(pi)) == True
    assert ask(Q.noninteger(E)) == True

    # Complex numbers
    assert ask(Q.noninteger(I)) == True
    assert ask(Q.noninteger(1 + I)) == True

    # Infinity
    assert ask(Q.noninteger(oo)) == True
    assert ask(Q.noninteger(-oo)) == True

    # Expressions
    assert ask(Q.noninteger(2 + 3)) == False
    assert ask(Q.noninteger(2.5 + 3)) == True
    assert ask(Q.noninteger(log(2))) == True

    # Applied predicates
    assert ask(Q.noninteger(x), Q.integer(x)) == False
    assert ask(Q.integer(x), Q.noninteger(x)) == False
    assert ask(Q.integer(x), Q.integer(x)) == True
    assert ask(Q.noninteger(x), Q.noninteger(x)) == True
