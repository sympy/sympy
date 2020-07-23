from sympy import Set
from sympy.map import (
    BinaryOperator, AssociativeOperator,
    AppliedBinaryOperator, AppliedAssociativeOperator,
)

A = Set('A')
a1, a2 = A.element('a1'), A.element('a2')
B = Set('B', (A,))
b1, b2 = B.element('b1'), B.element('b2')

class F(BinaryOperator):
    def __new__(cls, domain, codomain, **kwargs):
        obj = super().__new__(cls, domain, codomain, **kwargs)
        obj.domain = domain
        obj.codomain = codomain
        return obj
f = F(A**2, A)

class G(F, AssociativeOperator):
    pass
g_A = G(A**2, B)
g_B = G(B**2, B)

def test_BinaryOperator():
    # calling binary operator returns AppliedBinaryOperator
    assert isinstance(f(a1, a2), AppliedBinaryOperator)

def test_AppliedBinaryOperator():
    # check nested application (not flattened)
    expr = f(f(a1, a1), f(a2, a2), evaluate=True)
    assert expr.arguments == (f(a1, a1), f(a2, a2))

def test_AssociativeOperator():
    # calling associative operator returns AppliedAssociativeOperator
    assert isinstance(g_A(a1, a2), AppliedAssociativeOperator)

def test_AppliedAssociativeOperator():
    # check nested application (flattened)
    # also check that domain's subset's element can be argument
    expr = g_A(g_A(a1, b1), g_A(a2, b2), evaluate=True)
    assert expr.arguments == (a1, b1, a2, b2)

    # restricted function in argument is flattened
    expr = g_A(g_A(a1, b1), g_B(a2, b2), evaluate=True)
    assert expr.arguments == (a1, b1, a2, b2)

    # extended function in argument is not flattened
    expr = g_B(g_B(a1, b1), g_A(b2, b2), evaluate=True)
    assert expr.arguments == (a1, b1, g_A(b2, b2))
