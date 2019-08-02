from multipledispatch.variadic import isvariadic, Variadic


class A(object): pass
class B(A): pass
class C(object): pass


def test_is_variadic():
    assert isvariadic(Variadic[A])
    assert not isvariadic(A)


def test_equals_simple():
    type_a = Variadic[A]
    type_b = Variadic[B]
    other_type_a = Variadic[A]
    assert type_a == other_type_a
    assert not type_a == type_b


def test_equals_union():
    union_a_b = Variadic[(A, B)]
    union_b_a = Variadic[(B, A)]

    union_a_b_c = Variadic[(A, B, C)]
    union_c_b_a = Variadic[(C, B, A)]

    assert union_a_b == union_b_a
    assert union_a_b_c == union_c_b_a

    assert not union_a_b == union_a_b_c
