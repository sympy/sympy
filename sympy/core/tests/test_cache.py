from sympy.core.cache import cached_property, cacheit


def test_cacheit_doc():
    @cacheit
    def testfn():
        "test docstring"
        pass

    assert testfn.__doc__ == "test docstring"
    assert testfn.__name__ == "testfn"

def test_cacheit_unhashable():
    @cacheit
    def testit(x):
        return x

    assert testit(1) == 1
    assert testit(1) == 1
    a = {}
    assert testit(a) == {}
    a[1] = 2
    assert testit(a) == {1: 2}


class _ClassWithProperty:
    instances = 0
    calls = 0
    def __init__(self, value):
        self.__class__.instances += 1
        self.stored_value = value
    @cached_property
    def value(self):
        self.__class__.calls += 1
        return self.stored_value

def test_cached_instance_property():
    assert _ClassWithProperty.instances == 0
    assert _ClassWithProperty.calls == 0
    o1 = _ClassWithProperty(1)
    assert _ClassWithProperty.instances == 1
    assert _ClassWithProperty.calls == 0
    o2 = _ClassWithProperty(2)
    assert _ClassWithProperty.instances == 2
    assert _ClassWithProperty.calls == 0
    o1.value
    assert _ClassWithProperty.instances == 2
    assert _ClassWithProperty.calls == 1
    o1.value
    assert _ClassWithProperty.instances == 2
    assert _ClassWithProperty.calls == 1
    o2.value
    assert _ClassWithProperty.instances == 2
    assert _ClassWithProperty.calls == 2
