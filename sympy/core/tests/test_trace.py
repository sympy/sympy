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
    assert Tr(A*A*B*B) == Tr(A**2*B**2)
    assert Tr(A*B*B*A) == Tr(A**2*B**2)

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

#def test_permute():
#    A, B, C, D, E, F, G = symbols('A B C D E F G', commutative=False)
#    t = Tr(A*B*C*D*E*F*G, cycle=False)

#    assert t.permute(0).args[0].args == (A, B, C, D, E, F, G)
#    assert t.permute(2).args[0].args == (F, G, A, B, C, D, E)
#    assert t.permute(4).args[0].args == (D, E, F, G, A, B, C)
#    assert t.permute(6).args[0].args == (B, C, D, E, F, G, A)
#    assert t.permute(8).args[0].args == t.permute(1).args[0].args

#    assert t.permute(-1).args[0].args == (B, C, D, E, F, G, A)
#    assert t.permute(-3).args[0].args == (D, E, F, G, A, B, C)
#    assert t.permute(-5).args[0].args == (F, G, A, B, C, D, E)
#    assert t.permute(-8).args[0].args == t.permute(-1).args[0].args

#    t = Tr((A+B)*(B*B)*C*D,cycle=False)
#    assert t.permute(2).args[0].args == (C, D, (A+B), (B**2))

def test_equals():
    #TODO: Need to test if permute allowed
    pass
