from sympy import Mul, Basic, Q, Expr, And, symbols, Or

from sympy.assumptions.sathandlers import (ClassFactRegistry, AllArgs,
    AnyArgs, ExactlyOneArg, ArgFactHandler,)

from sympy.testing.pytest import raises


x, y, z = symbols('x y z')


def test_class_handler_registry():
    my_handler_registry = ClassFactRegistry()

    # The predicate doesn't matter here, so just pass
    @my_handler_registry.register(Mul)
    def fact1(expr):
        pass
    @my_handler_registry.multiregister(Expr)
    def fact2(expr):
        pass

    assert my_handler_registry[Basic] == (frozenset(), frozenset())
    assert my_handler_registry[Expr] == (frozenset(), frozenset({fact2}))
    assert my_handler_registry[Mul] == (frozenset({fact1}), frozenset({fact2}))


def test_ArgFactHandler():

    class MyArgFactHandler(ArgFactHandler):
        def __call__(self, expr):
            return self.apply((expr,))

    a = MyArgFactHandler((x,), Q.positive(x))
    b = MyArgFactHandler((x,), Q.positive(x) | Q.negative(x))

    assert a(x) == Q.positive(x)
    assert b(x) == Q.positive(x) | Q.negative(x)

    raises(TypeError, lambda: MyArgFactHandler(Q.positive))


def test_AllArgs():
    a = AllArgs(x, Q.zero(x))
    b = AllArgs(x, Q.positive(x) | Q.negative(x))
    assert a(x*y) == And(Q.zero(x), Q.zero(y))
    assert b(x*y) == And(Q.positive(x) | Q.negative(x), Q.positive(y) | Q.negative(y))


def test_AnyArgs():
    a = AnyArgs(x, Q.zero(x))
    b = AnyArgs(x, Q.positive(x) & Q.negative(x))
    assert a(x*y) == Or(Q.zero(x), Q.zero(y))
    assert b(x*y) == Or(Q.positive(x) & Q.negative(x), Q.positive(y) & Q.negative(y))


def test_ExactlyOneArg():
    a = ExactlyOneArg(x, Q.zero(x))
    b = ExactlyOneArg(x, Q.positive(x) | Q.negative(x))
    assert a(x*y) == Or(Q.zero(x) & ~Q.zero(y), Q.zero(y) & ~Q.zero(x))
    assert a(x*y*z) == Or(Q.zero(x) & ~Q.zero(y) & ~Q.zero(z), Q.zero(y)
        & ~Q.zero(x) & ~Q.zero(z), Q.zero(z) & ~Q.zero(x) & ~Q.zero(y))
    assert b(x*y) == Or((Q.positive(x) | Q.negative(x)) &
        ~(Q.positive(y) | Q.negative(y)), (Q.positive(y) | Q.negative(y)) &
        ~(Q.positive(x) | Q.negative(x)))
