from sympy import S, Set, assuming, Q
from sympy.map import (
    BinaryOperator,
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
    def _eval_check_left_identity(self, element):
        # assume identity is in A but outside of B
        if element == a3:
            return True
f = F(A**2, A)

class F2(F):
    identity = a3
f2 = F2(A**2, A)

class G(BaseOp):
    associative=True
    def _eval_check_right_identity(self, element):
        # assume identity is in B (thus in A as well)
        if element == b3:
            return True

g_A = G(A**2, A)
g_B = G(B**2, B)

class G2(G):
    identity = a3
g2 = G2(A**2, A)

class H1(BaseOp):
    left_divisible = True
h1 = H1(A**2, A)

class H2(BaseOp):
    right_divisible = True
h2 = H2(A**2, A)

class L(BaseOp):
    identity = a3
l = L(A**2, A)

class M(BaseOp):
    associative = True
    identity = a3
m = M(A**2, A)

class N(BaseOp):
    associative = True
    commutative = True
    identity = a3
n = N(A**2, A)

def test_BinaryOperator():

    # 0 or 1 arguments to binary operator is not supported
    raises(TypeError, lambda: f())
    raises(TypeError, lambda: f(a1))

    # 3 or more arguments to non-associative binary operator is not supported
    raises(TypeError, lambda: f(a1, a2, a2))

    # identity checking
    assert f.check_left_identity(a3)
    assert not f.check_right_identity(a3)

    # two side identity is both left and right identity
    assert f2.check_left_identity(a3)
    assert f2.check_right_identity(a3)

    class T(BaseOp):
        def eval(self, a, b):
            return a/b
    t = T(S.Integers**2, S.Integers)
    # 3 or more arguments is allowed if mapping is associative
    raises(TypeError, lambda: t(1, 2, 3, evaluate=True))
    with assuming(Q.associative(t)):
        t(1, 2, 3) # not raise error
    # result of evaluation must be in codomain
    raises(TypeError, lambda: t(2, 3, evaluate=True))

    # check nested application (not flattened)
    expr = f(f(a1, a1), f(a2, a2)).doit()
    assert expr.arguments == (f(a1, a1), f(a2, a2))

    # associativity can be assumed.
    with assuming(Q.associative(f)):
        expr = f(f(a1, a2), f(a1, a2)).doit()
        assert expr.arguments == (a1, a2, a1, a2)

    # identity removal
    expr = f(f(a3, a2), f(a1, a3)).doit()
    assert expr.arguments == (a2, f(a1, a3))
    expr = f2(f2(a3, a2), f2(a1, a3)).doit()
    assert expr.arguments == (a2, a1)

    # associative operators
    # check nested application (flattened)
    # also check that domain's subset's element can be argument
    expr = g_A(g_A(a1, a2), g_A(b2, b1), evaluate=True)
    assert expr.arguments == (a1, a2, b2, b1)

    # restricted function in argument is flattened
    expr = g_A(g_A(a1, b1), g_B(b2, b1), evaluate=True)
    assert expr.arguments == (a1, b1, b2, b1)

    # identity removal with flattening
    expr = g_A(g_A(b3, b2), g_A(b1, b2)).doit()
    assert expr.arguments == (b3, b2, b1, b2)
    expr = g_A(g_A(b2, b3), g_A(b1, b2)).doit()
    assert expr.arguments == (b2, b1, b2)
    expr = g_A(g_A(b3, b1), g_A(b3, b2)).doit()
    assert expr.arguments == (b3, b1, b2)
    expr = g_A(g_A(b2, b1), g_A(a2, g_A(b2, b3))).doit()
    assert expr.arguments == (b2, b1, a2, b2)

    expr = g2(g2(a3, a1), g2(a3, a2)).doit()
    assert expr.arguments == (a1, a2)

def test_LeftDivision():
    # left divisibility can be assumed
    raises(TypeError, lambda: f.left_division_operator())
    with assuming(Q.left_divisible(f)):
        f.left_division_operator()   # not raise error

    h1_ld = h1.left_division_operator()
    # domain and codomain is preserved
    assert h1_ld.domain == h1.domain
    assert h1_ld.codomain == h1.codomain

    # is_left_division
    assert h1_ld.is_left_division(h1)

    # left division cancellation
    assert h1(a1, h1_ld(a1, a2), evaluate=True) == a2
    assert h1_ld(a1, h1(a1, a2), evaluate=True) == a2
    assert h1(a1, h1_ld(a2, a1), evaluate=True) != a2
    assert h1_ld(a1, h1(a2, a1), evaluate=True) != a2
    assert h1(h1_ld(a2, a1), a1, evaluate=True) != a2
    assert h1_ld(h1(a2, a1), a1, evaluate=True) != a2

    assert h1(h1(a1, a1), h1_ld(a1, a2), evaluate=True) != h1(a1, a2)
    with assuming(Q.associative(h1)):
        assert h1(h1(a1, a1), h1_ld(a1, a2), evaluate=True) == h1(a1, a2)

def test_RightDivision():
    # right divisibility can be assumed
    raises(TypeError, lambda: f.right_division_operator())
    with assuming(Q.right_divisible(f)):
        f.right_division_operator()   # not raise error

    h2_rd = h2.right_division_operator()
    # domain and codomain is preserved
    assert h2_rd.domain == h2.domain
    assert h2_rd.codomain == h2.codomain

    # is_right_division
    assert h2_rd.is_right_division(h2)

    # right division cancellation
    assert h2(h2_rd(a1, a2), a2, evaluate=True) == a1
    assert h2_rd(h2(a1, a2), a2, evaluate=True) == a1
    assert h2(h2_rd(a2, a1), a2, evaluate=True) != a1
    assert h2_rd(h2(a2, a1), a2, evaluate=True) != a1
    assert h2(a2, h2_rd(a2, a1), evaluate=True) != a1
    assert h2_rd(a2, h2(a2, a1), evaluate=True) != a1

    assert h2(h2_rd(a1, a2), h2(a2, a2), evaluate=True) != h2(a2, a2)
    with assuming(Q.associative(h2)):
        assert h2(h2_rd(a1, a2), h2(a2, a2), evaluate=True) == h2(a1, a2)

def test_InverseOperator():
    l_inv = l.inverse_operator()
    # inverse cancellation
    assert l(a1, l_inv(a1), evaluate=True) == a3
    assert l(l_inv(a1), a1, evaluate=True) == a3

    # are_inverse
    assert l.are_inverse(a1, l_inv(a1)) == l.are_inverse(l_inv(a1), a1) == True
    # inverse of identity
    assert l_inv(a3, evaluate=True) == a3
    # inverse operator is idempotent
    assert l_inv(l_inv(a1), evaluate=True) == a1
    assert l_inv(l_inv(l_inv(a1)), evaluate=True) == l_inv(a1)

    # as_base_exp works only when base operator is passed
    assert l_inv(a1).as_base_exp(l) == (a1, -1)
    assert l_inv(a1).as_base_exp(f) == (l_inv(a1), 1)

    # inverse element of associative operator is converted to exponent
    m_expop = m.exponent_operator()
    m_invop = m.inverse_operator()
    assert m_invop(a1) == m_expop(a1, -1)
    assert m_invop(m_expop(a1, 2)) == m_expop(m_expop(a1, 2), -1)

def test_ExponentOperator():
    m_expop = m.exponent_operator()
    m_invop = m.inverse_operator()
    # repetitive arguments are converted to ExponentOperator
    assert m(a1, a1, a1, evaluate=True) == m_expop(a1, 3)
    # n=1 is evaluated to original element
    assert m_expop(a1, 1, evaluate=True) == a1
    # n=0 is evaluated to identity element
    assert m_expop(a1, 0, evaluate=True) == a3

    # exponentiation of exponentiation is combined
    assert m_expop(m_expop(a1, 2), 3, evaluate=True) == m_expop(a1, 6)

    # as_base_exp works only when base operator is passed
    assert m_expop(a1, 2).as_base_exp(m) == (a1, 2)
    assert m_expop(a1, 2).as_base_exp(l) == (m_expop(a1, 2), 1)

def test_assoc_comm_process():
    n_invop = n.inverse_operator()
    n_expop = n.exponent_operator()
    assert n(a1, a3, a2, a1, n_invop(a2), n_invop(a1), evaluate=True) == a1
    assert n(a1, a2, a1, a3, n_expop(a1, -2), evaluate=True) == a2
