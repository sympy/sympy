from sympy import Set, assuming, Q
from sympy.map import (
    BinaryOperator, AppliedBinaryOperator,
)
from sympy.testing.pytest import raises

A = Set('A')
a1, a2, a3 = [A.element('a%s' % i) for i in range(1,4)]
B = Set('B', (A,))
b1, b2, b3 = [B.element('b%s' % i) for i in range(1, 4)]

class BaseOp(BinaryOperator):
    def __new__(cls, domain, codomain, **kwargs):
        obj = super().__new__(cls, domain, codomain, **kwargs)
        obj.domain = domain
        obj.codomain = codomain
        return obj

class F(BaseOp):
    def _eval_left_identity(self, element):
        # assume identity is in A but outside of B
        if element == a3:
            return True
f = F(A**2, A)

class F2(F):
    identity = a3
f2 = F2(A**2, A)

class G(BaseOp):
    is_associative=True
    def _eval_right_identity(self, element):
        # assume identity is in B (thus in A as well)
        if element == b3:
            return True

g_A = G(A**2, B)
g_B = G(B**2, B)

class G2(G):
    identity = a3
g2 = G2(A**2, A)

class H(BaseOp):
    pass
h = H(A**2, A)

class H1(BaseOp):
    is_left_divisible = True
    def _eval_left_division(self):
        return H(self.domain, self.codomain)
h1 = H1(A**2, A)

class H2(BaseOp):
    is_right_divisible = True
    def _eval_right_division(self):
        return H(self.domain, self.codomain)
h2 = H2(A**2, A)

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

    # associativity can be assumed.
    with assuming(Q.associative(f)):
        expr = f(f(a1, a1), f(a2, a2)).doit()
        assert expr.arguments == (a1, a1, a2, a2)

    # identity removal
    expr = f(f(a3, a2), f(a1, a3)).doit()
    assert expr.arguments == (a2, f(a1, a3))
    expr = f2(f2(a3, a2), f2(a1, a3)).doit()
    assert expr.arguments == (a2, a1)

    # associative operators
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

def test_LeftDivision():
    # left divisibility can be assumed
    raises(TypeError, lambda: f.left_division())
    with assuming(Q.left_divisible(f)):
        f.left_division()   # not raise error

    # left division operator can be evaluated
    assert h1.left_division(evaluate=True) == h

    h1_ld = h1.left_division()
    # domain and codomain is preserved
    assert h1_ld.domain == h1.domain
    assert h1_ld.codomain == h1.codomain

    # is_left_division
    assert h1_ld.is_left_division(h1)
    assert h.is_left_division(h1)

    # left division cancellation
    assert h1(a1, h1_ld(a1, a2), evaluate=True) == a2
    assert h1_ld(a1, h1(a1, a2), evaluate=True) == a2
    assert h1(a1, h1_ld(a2, a1), evaluate=True) != a2
    assert h1_ld(a1, h1(a2, a1), evaluate=True) != a2
    assert h1(h1_ld(a2, a1), a1, evaluate=True) != a2
    assert h1_ld(h1(a2, a1), a1, evaluate=True) != a2

    assert h1(a1, h(a1, a2), evaluate=True) == a2
    assert h(a1, h1(a1, a2), evaluate=True) == a2
    assert h1(a1, h(a2, a1), evaluate=True) != a2
    assert h(a1, h1(a2, a1), evaluate=True) != a2
    assert h1(h(a2, a1), a1, evaluate=True) != a2
    assert h(h1(a2, a1), a1, evaluate=True) != a2

    assert h1(h1(a1, a1), h1_ld(a1, a2), evaluate=True) != h1(a1, a2)
    with assuming(Q.associative(h1)):
        assert h1(h1(a1, a1), h1_ld(a1, a2), evaluate=True) == h1(a1, a2)

def test_RightDivision():
    # right divisibility can be assumed
    raises(TypeError, lambda: f.right_division())
    with assuming(Q.right_divisible(f)):
        f.right_division()   # not raise error

    # right division operator can be evaluated
    assert h2.right_division(evaluate=True) == h

    h2_rd = h2.right_division()
    # domain and codomain is preserved
    assert h2_rd.domain == h2.domain
    assert h2_rd.codomain == h2.codomain

    # is_right_division
    assert h2_rd.is_right_division(h2)
    assert h.is_right_division(h2)

    # right division cancellation
    assert h2(h2_rd(a1, a2), a2, evaluate=True) == a1
    assert h2_rd(h2(a1, a2), a2, evaluate=True) == a1
    assert h2(h2_rd(a2, a1), a2, evaluate=True) != a1
    assert h2_rd(h2(a2, a1), a2, evaluate=True) != a1
    assert h2(a2, h2_rd(a2, a1), evaluate=True) != a1
    assert h2_rd(a2, h2(a2, a1), evaluate=True) != a1

    assert h2(h(a1, a2), a2, evaluate=True) == a1
    assert h(h2(a1, a2), a2, evaluate=True) == a1
    assert h2(h(a2, a1), a2, evaluate=True) != a1
    assert h(h2(a2, a1), a2, evaluate=True) != a1
    assert h2(a2, h(a2, a1), evaluate=True) != a1
    assert h(a2, h2(a2, a1), evaluate=True) != a1

    assert h2(h2_rd(a1, a2), h2(a2, a2), evaluate=True) != h2(a2, a2)
    with assuming(Q.associative(h2)):
        assert h2(h2_rd(a1, a2), h2(a2, a2), evaluate=True) == h2(a1, a2)
