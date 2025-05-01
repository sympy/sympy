from sympy.core.mul import Mul
from sympy.core.symbol import symbols
from sympy.testing.pytest import raises

from sympy.physics.quantum.slidingtransform import SlidingTransform

a, b, c, d, e, f = symbols('a b c d e f', commutative=False)


def doubler(x):
    return (x, x)

def drone(x):
    return (a,)


def test_unary():
    st = SlidingTransform(unary=doubler)
    assert st(a*b*c) == Mul._from_args((a, a, b, b, c, c))

    st = SlidingTransform(unary=drone)
    assert st(a*b*c) == Mul._from_args((a, a, a))

    st = SlidingTransform(unary=lambda x: None)
    assert st(a*b*c) == Mul._from_args((a, b, c))


def mapping(lhs, rhs):
    if lhs == a and rhs == b:
        return (c,)
    elif lhs == b and rhs == c:
        return (d,)
    elif lhs == a and rhs == c:
        return None
    elif lhs == c and rhs == d:
        return (e, f)
    

def test_binary():
    st = SlidingTransform(binary=mapping)

    assert st(a*b) == Mul._from_args((c,))
    assert st(b*c) == Mul._from_args((d,))
    assert st(a*c) == Mul._from_args((a,c))
    assert st(c*d) == Mul._from_args((e,f))

    assert st(a*b*c) == Mul._from_args((c, c))
    assert st(b*c*d) == Mul._from_args((d, d))
    assert st(c*d*e) == Mul._from_args((e, f, e))

    assert st(a*b*d) == Mul._from_args((e, f))


def test_reverse():
    st = SlidingTransform(binary=mapping, reverse=True)

    assert st(a*b) == Mul._from_args((c,))
    assert st(b*c) == Mul._from_args((d,))
    assert st(a*c) == Mul._from_args((a,c))
    assert st(c*d) == Mul._from_args((e,f))

    assert st(a*b*c) == Mul._from_args((a, d))
    assert st(d*b*c) == Mul._from_args((d, d))
    assert st(a*c*b*a*b) == Mul._from_args((a, e, f))

