from multipledispatch import dispatch
import pytest


@dispatch(int)
def isint(x):
    return True


@dispatch(object)
def isint(x):
    return False


@dispatch(object, object)
def isint(x, y):
    return False


@pytest.mark.parametrize("val", [1, 'a'])
def test_benchmark_call_single_dispatch(benchmark, val):
    benchmark(isint, val)


@pytest.mark.parametrize("val", [(1, 4)])
def test_benchmark_call_multiple_dispatch(benchmark, val):
    benchmark(isint, *val)


def test_benchmark_add_and_use_instance(benchmark):
    namespace = {}

    @benchmark
    def inner():
        @dispatch(int, int, namespace=namespace)
        def mul(x, y):
            return x * y

        @dispatch(str, int, namespace=namespace)
        def mul(x, y):
            return x * y

        @dispatch(int, int, [float], namespace=namespace)
        def mul(x, y, *args):
            return x * y

        @dispatch([int], namespace=namespace)
        def mul(*args):
            return sum(args)

        mul(4, 5)
        mul('x', 5)
        mul(1, 2, 3., 4., 5.)
        mul(1, 2, 3, 4, 5)
