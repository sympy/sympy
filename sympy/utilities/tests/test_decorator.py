from __future__ import annotations
from functools import wraps

from sympy.utilities.decorator import _relational_opaque, threaded, xthreaded, memoize_property, deprecated
from sympy.testing.pytest import raises, warns_deprecated_sympy

from sympy.core.basic import Basic
from sympy.core.containers import Tuple
from sympy.core.relational import Eq
from sympy.polys import partfrac, polytools
from sympy.utilities import decorator as decorator_module
from sympy.matrices.dense import Matrix

from sympy.abc import x, y


def test_threaded():
    @threaded
    def function(expr, *args):
        return 2*expr + sum(args)

    assert function(Matrix([[x, y], [1, x]]), 1, 2) == \
        Matrix([[2*x + 3, 2*y + 3], [5, 2*x + 3]])

    assert function(Eq(x, y), 1, 2) == Eq(x, y)
    assert function([x, Eq(x, y)], 1, 2) == [2*x + 3, Eq(x, y)]

    assert function([x, y], 1, 2) == [2*x + 3, 2*y + 3]
    assert function((x, y), 1, 2) == (2*x + 3, 2*y + 3)

    assert function({x, y}, 1, 2) == {2*x + 3, 2*y + 3}

    @threaded
    def function(expr, n):
        return expr**n

    assert function(x + y, 2) == x**2 + y**2
    assert function(x, 2) == x**2


def test_xthreaded():
    @xthreaded
    def function(expr, n):
        return expr**n

    assert function(x + y, 2) == (x + y)**2
    assert function(Eq(x, y), 2) == Eq(x, y)
    assert function((x, Eq(x, y)), 2) == (x**2, Eq(x, y))


def test_threaded_relational_opaque():
    for decorator in (threaded, xthreaded):
        visited = []

        @decorator
        def function(expr, **kwargs):
            visited.append(expr)
            return expr

        rel = Eq(x + y, x*y, evaluate=False)
        assert function(rel, flag=True) is rel
        assert visited == []
        assert function([rel, x], flag=True) == [rel, x]
        assert visited == [x]

        @decorator
        def required(expr, n):
            return expr**n

        raises(TypeError, lambda: required(rel))
        raises(TypeError, lambda: required(rel, 2, 3))
        raises(TypeError, lambda: required(rel, 2, expr=x))
        raises(TypeError, lambda: required(rel, 2, unknown=True))
        assert required(rel, 2) is rel


def test_relational_opaque():
    visited = []

    @_relational_opaque
    def function(value, n=1):
        visited.append(value)
        return value*n

    rel = Eq(x + y, 1, evaluate=False)
    assert function(value=rel, n=2) is rel
    assert visited == []
    assert function(Tuple(rel, x), 2) == Tuple(rel, 2*x)
    assert visited == [x]
    raises(TypeError, lambda: function())
    raises(TypeError, lambda: function(rel, value=x))
    raises(TypeError, lambda: function(rel, 2, 3))
    raises(TypeError, lambda: function(value=rel, unknown=True))

    @_relational_opaque
    def required(value, n):
        return value*n

    raises(TypeError, lambda: required(rel))
    raises(TypeError, lambda: required(value=rel))
    assert required(value=rel, n=2) is rel

    @_relational_opaque(allow_relational=True)
    def reverse(expr):
        return expr.reversed if expr.is_Relational else expr

    assert reverse(Tuple(rel, x)) == Tuple(rel.reversed, x)


def test_relational_opaque_public_exports():
    assert 'apart' in partfrac.__all__
    assert 'cancel' in polytools.__all__
    assert not hasattr(decorator_module, '__all__')


def test_wraps():
    def my_func(x):
        """My function. """

    my_func.is_my_func = True

    new_my_func = threaded(my_func)
    new_my_func = wraps(my_func)(new_my_func)

    assert new_my_func.__name__ == 'my_func'
    assert new_my_func.__doc__ == 'My function. '
    assert hasattr(new_my_func, 'is_my_func')
    assert new_my_func.is_my_func is True


def test_memoize_property():
    class TestMemoize(Basic):
        @memoize_property
        def prop(self):
            return Basic()

    member = TestMemoize()
    obj1 = member.prop
    obj2 = member.prop
    assert obj1 is obj2

def test_deprecated():
    @deprecated('deprecated_function is deprecated',
                deprecated_since_version='1.10',
                # This is the target at the top of the file, which will never
                # go away.
                active_deprecations_target='active-deprecations')
    def deprecated_function(x):
        return x

    with warns_deprecated_sympy():
        assert deprecated_function(1) == 1

    @deprecated('deprecated_class is deprecated',
                deprecated_since_version='1.10',
                active_deprecations_target='active-deprecations')
    class deprecated_class:
        pass

    with warns_deprecated_sympy():
        assert isinstance(deprecated_class(), deprecated_class)

    # Ensure the class decorator works even when the class never returns
    # itself
    @deprecated('deprecated_class_new is deprecated',
                deprecated_since_version='1.10',
                active_deprecations_target='active-deprecations')
    class deprecated_class_new:
        def __new__(cls, arg):
            return arg

    with warns_deprecated_sympy():
        assert deprecated_class_new(1) == 1

    @deprecated('deprecated_class_init is deprecated',
                deprecated_since_version='1.10',
                active_deprecations_target='active-deprecations')
    class deprecated_class_init:
        def __init__(self, arg):
            self.arg = 1

    with warns_deprecated_sympy():
        assert deprecated_class_init(1).arg == 1

    @deprecated('deprecated_class_new_init is deprecated',
                deprecated_since_version='1.10',
                active_deprecations_target='active-deprecations')
    class deprecated_class_new_init:
        def __new__(cls, arg):
            if arg == 0:
                return arg
            return object.__new__(cls)

        def __init__(self, arg):
            self.arg = 1

    with warns_deprecated_sympy():
        assert deprecated_class_new_init(0) == 0

    with warns_deprecated_sympy():
        assert deprecated_class_new_init(1).arg == 1
