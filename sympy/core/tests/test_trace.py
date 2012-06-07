from sympy import S, Symbol, symbols, Matrix
from sympy.core.trace import Tr

def test_trace_new():
    a, b, c, d = symbols('a b c d')
    A, B, C, D = symbols('A B C D', commutative=False)

    #Tr(A+B) = Tr(A) + Tr(B)
    assert Tr(a+b) == Tr(a) + Tr(b)
    # check for mul and adds
    assert Tr((a*b) + ( c*d)) == (a*b) + (c*d)
    # Tr(scalar*A) = scalar*Tr(A)
    assert Tr( a*A ) == a*Tr(A)
    assert Tr(a*b*A*B) == a*b*Tr(A*B)

    #POW
    assert Tr ( pow(a,b) ) == a**b

    M = Matrix([[1,1],[2,2]])
    assert Tr(a*M) == Tr(Matrix ( [[a,a],[2*a,2*a]]))

    #trace indices test
    t = Tr((a+b), (2))
    assert t.args[0].args[1] == (2) and t.args[1].args[1] == (2)

    t = Tr(a*A, (2,3))
    print t.args[1].args[1]
    assert t.args[1].args[1] == (2,3)

def test_trace_doit():
    a, b, c, d = symbols('a b c d')
    A, B, C, D = symbols('A B C D', commutative=False)

    M = Matrix([[1,1],[2,2]])
    assert Tr(M).doit() == 3
