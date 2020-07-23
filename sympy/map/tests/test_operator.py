from sympy import Set
from sympy.map import (
    BinaryOperator, AssociativeOperator,
    AppliedBinaryOperator, AppliedAssociativeOperator,
)
from sympy.testing.pytest import raises

A = Set('A')
a1, a2, a3 = [A.element('a%s' % i) for i in range(3)]
B = Set('B', (A,))
b1, b2, b3 = [B.element('b%s' % i) for i in range(3)]

class F(BinaryOperator):
    def __new__(cls, domain, codomain, **kwargs):
        obj = super().__new__(cls, domain, codomain, **kwargs)
        obj.domain = domain
        obj.codomain = codomain
        return obj

    def _eval_left_identity(self, element):
        # assume identity is in A but outside of B
        if element == a3:
            return True
f = F(A**2, A)

class F2(F):
    identity = a3
f2 = F2(A**2, A)

class G(AssociativeOperator):
    def __new__(cls, domain, codomain, **kwargs):
        obj = super().__new__(cls, domain, codomain, **kwargs)
        obj.domain = domain
        obj.codomain = codomain
        return obj

    def _eval_right_identity(self, element):
        # assume identity is in B (thus in A as well)
        if element == b3:
            return True

g_A = G(A**2, B)
g_B = G(B**2, B)

class G2(G):
    identity = a3
g2 = G2(A**2, A)

def test_BinaryOperator():
    # calling binary operator returns AppliedBinaryOperator
    assert isinstance(f(a1, a2), AppliedBinaryOperator)

    # 0 or 1 arguments to binary operator is not supported
    raises(TypeError, lambda: f())
    raises(TypeError, lambda: f(a1))

    # 3 or more arguments to non-associative binary operator is not supported
    raises(TypeError, lambda: f(a1, a2, a2))

    # identity checking
    assert f.left_identity(a3)
    assert not f.right_identity(a3)

    # two side identity is both left and right identity
    assert f2.left_identity(a3)
    assert f2.right_identity(a3)

def test_AppliedBinaryOperator():
    # check nested application (not flattened)
    expr = f(f(a1, a1), f(a2, a2)).doit()
    assert expr.arguments == (f(a1, a1), f(a2, a2))

    # identity removal
    expr = f(f(a3, a2), f(a1, a3)).doit()
    assert expr.arguments == (a2, f(a1, a3))
    expr = f2(f2(a3, a2), f2(a1, a3)).doit()
    assert expr.arguments == (a2, a1)

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

    # identity removal with flattening
    expr = g_A(g_A(b3, b2), g_A(b1, b2)).doit()
    assert expr.arguments == (b3, b2, b1, b2)
    expr = g_A(g_A(b2, b3), g_A(b1, b2)).doit()
    assert expr.arguments == (b2, b1, b2)
    expr = g_A(g_A(b3, b1), g_A(b3, b2)).doit()
    assert expr.arguments == (b3, b1, b2)
    expr = g_A(g_A(b2, b1), g_A(b2, g_A(b2, b3))).doit()
    assert expr.arguments == (b2, b1, b2, b2)

    expr = g2(g2(a3, a1), g2(a3, a2)).doit()
    assert expr.arguments == (a1, a2)
