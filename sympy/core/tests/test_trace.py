from sympy import S, Symbol, symbols, Matrix
from sympy.core.trace import Tr

def test_trace_new():
    a, b, c, d, Y = symbols('a b c d Y')
    A, B, C, D = symbols('A B C D', commutative=False)

    assert Tr(a+b) == a + b
    assert Tr(A+B) == Tr(A) + Tr(B)

    # check for mul and adds
    assert Tr((a*b) + ( c*d)) == (a*b) + (c*d)
    # Tr(scalar*A) = scalar*Tr(A)
    assert Tr(a*A) == a*Tr(A)

    #also check if Muls are permuted in canonical form
    assert Tr(a*A*B*b) == a*b*Tr(A*B)
    assert Tr(a*C*D*A*B) == a*Tr(A*B*C*D)
    assert Tr(a*C*A*D*B*2) == 2*a*Tr(A*D*B*C)
    assert Tr(B*A*C*B*A) == Tr(A*C*B*A*B)
    assert Tr(A*C*B*A*B) == Tr(A*B*A*C*B)

    # since A is symbol and not commutative
    assert isinstance(Tr(A), Tr)

    #POW
    assert Tr(pow(a, b)) == a**b
    assert isinstance(Tr(pow(A, a)), Tr)

    #Matrix
    M = Matrix([[1,1], [2,2]])
    assert Tr(M) == 3

    #trace indices test
    t = Tr((A+B), (2))
    assert t.args[0].args[1] == (2) and t.args[1].args[1] == (2)

    t = Tr(a*A, (2,3))
    assert t.args[1].args[1] == (2,3)

def test_trace_doit():
    a, b, c, d = symbols('a b c d')
    A, B, C, D = symbols('A B C D', commutative=False)

    #TODO: needed while testing reduced density operations, etc.
