# import sys

# from nose import SkipTest

from multipledispatch import dispatch
from multipledispatch.dispatcher import Dispatcher


def test_function_annotation_register():
    f = Dispatcher('f')

    @f.register()
    def inc(x: int):
        return x + 1

    @f.register()
    def inc(x: float):
        return x - 1

    assert f(1) == 2
    assert f(1.0) == 0.0


def test_function_annotation_dispatch():
    @dispatch()
    def inc(x: int):
        return x + 1

    @dispatch()
    def inc(x: float):
        return x - 1

    assert inc(1) == 2
    assert inc(1.0) == 0.0


def test_function_annotation_dispatch_custom_namespace():
    namespace = {}

    @dispatch(namespace=namespace)
    def inc(x: int):
        return x + 2

    @dispatch(namespace=namespace)
    def inc(x: float):
        return x - 2

    assert inc(1) == 3
    assert inc(1.0) == -1.0

    assert namespace['inc'] == inc
    assert set(inc.funcs.keys()) == set([(int,), (float,)])


def test_method_annotations():
    class Foo():
        @dispatch()
        def f(self, x: int):
            return x + 1

        @dispatch()
        def f(self, x: float):
            return x - 1

    foo = Foo()

    assert foo.f(1) == 2
    assert foo.f(1.0) == 0.0


def test_overlaps():
    @dispatch(int)
    def inc(x: int):
        return x + 1

    @dispatch(float)
    def inc(x: float):
        return x - 1

    assert inc(1) == 2
    assert inc(1.0) == 0.0


def test_overlaps_conflict_annotation():
    @dispatch(int)
    def inc(x: str):
        return x + 1

    @dispatch(float)
    def inc(x: int):
        return x - 1

    assert inc(1) == 2
    assert inc(1.0) == 0.0
