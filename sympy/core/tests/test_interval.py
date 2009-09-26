from sympy import Symbol, Union, Interval, oo, S, Inequality, max_, min_

def test_union():
    assert Union(Interval(1, 2), Interval(2, 3)) == Interval(1, 3)
    assert Union(Interval(1, 2), Interval(2, 3, True)) == Interval(1, 3)
    assert Union(Interval(1, 2, False, True), Interval(2, 3, True)) != \
           Interval(1, 3)

def test_infinite_intervals():
    assert Interval(0, oo) == Interval(0, oo, False, True)
    assert Interval(-oo, 0) == Interval(-oo, 0, True, False)

def test_interval_argument_order():
    assert Interval(1, 0) == S.EmptySet

def test_complement():
    assert Interval(0, 1).complement == \
           Union(Interval(-oo, 0, True, True), Interval(1, oo, True, True))
    assert Interval(0, 1, True, False).complement == \
           Union(Interval(-oo, 0, True, False), Interval(1, oo, True, True))
    assert Interval(0, 1, False, True).complement == \
           Union(Interval(-oo, 0, True, True), Interval(1, oo, False, True))
    assert Interval(0, 1, True, True).complement == \
           Union(Interval(-oo, 0, True, False), Interval(1, oo, False, True))

def test_union_complement():
    assert Union(Interval(0, 1), Interval(2, 3)).complement == \
           Union(Interval(-oo, 0, True, True), Interval(1, 2, True, True),
                 Interval(3, oo, True, True))

def test_intersect():
    assert Interval(0, 2).intersect(Interval(1, 2)) == Interval(1, 2)
    assert Interval(0, 2).intersect(Interval(1, 2, True)) == \
           Interval(1, 2, True)
    assert Interval(0, 2, True).intersect(Interval(1, 2)) == \
           Interval(1, 2, False, False)
    assert Interval(0, 2, True, True).intersect(Interval(1, 2)) == \
           Interval(1, 2, False, True)

def test_symbolic_end_points():
    a = Symbol('a', real=True)

    assert Union(Interval(0, a), Interval(0, 3)).sup == max_(a, 3)
    assert Union(Interval(a, 0), Interval(-3, 0)).inf == min_(-3, a)

    assert Interval(0, a).contains(1) == Inequality(1, a)

def test_measure():
    a = Symbol('a', real=True)

    assert Interval(1, 3).measure == 2
    assert Interval(0, a).measure == a
    assert Interval(1, a).measure == a - 1

    assert Union(Interval(1, 2), Interval(3, 4)).measure == 2
