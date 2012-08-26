from sympy import (
    Symbol, Set, Union, Interval, oo, S, sympify, nan,
    GreaterThan, LessThan, Max, Min, And, Or, Eq, Ge, Le, Gt, Lt, Float,
    FiniteSet, Intersection
)
from sympy.mpmath import mpi

from sympy.utilities.pytest import raises
from sympy.utilities.pytest import raises, XFAIL

def test_interval_arguments():
    assert Interval(0, oo) == Interval(0, oo, False, True)
    assert Interval(0, oo).right_open == True
    assert Interval(-oo, 0) == Interval(-oo, 0, True, False)
    assert Interval(-oo, 0).left_open == True

    assert isinstance(Interval(1, 1), FiniteSet)

    assert Interval(1, 0) == S.EmptySet
    assert Interval(1, 1).measure == 0

    assert Interval(1, 1, False, True) == S.EmptySet
    assert Interval(1, 1, True, False) == S.EmptySet
    assert Interval(1, 1, True, True) == S.EmptySet

    raises(ValueError, lambda: Interval(0, S.ImaginaryUnit))
    raises(ValueError, lambda: Interval(0, Symbol('z')))

    assert isinstance(Interval(1, Symbol('a', real=True)), Interval)

def test_interval_symbolic_end_points():
    a = Symbol('a', real=True)

    assert Union(Interval(0, a), Interval(0, 3)).sup == Max(a, 3)
    assert Union(Interval(a, 0), Interval(-3, 0)).inf == Min(-3, a)

    assert Interval(0, a).contains(1) == LessThan(1, a)

def test_union():
    assert Union(Interval(1, 2), Interval(2, 3)) == Interval(1, 3)
    assert Union(Interval(1, 2), Interval(2, 3, True)) == Interval(1, 3)
    assert Union(Interval(1, 3), Interval(2, 4)) == Interval(1, 4)
    assert Union(Interval(1, 2), Interval(1, 3)) == Interval(1, 3)
    assert Union(Interval(1, 3), Interval(1, 2)) == Interval(1, 3)
    assert Union(Interval(1, 3, False, True), Interval(1, 2)) == \
           Interval(1, 3, False, True)
    assert Union(Interval(1, 3), Interval(1, 2, False, True)) == Interval(1, 3)
    assert Union(Interval(1, 2, True), Interval(1, 3)) == Interval(1, 3)
    assert Union(Interval(1, 2, True), Interval(1, 3, True)) == Interval(1, 3, True)
    assert Union(Interval(1, 2, True), Interval(1, 3, True, True)) == \
           Interval(1, 3, True, True)
    assert Union(Interval(1, 2, True, True), Interval(1, 3, True)) == \
           Interval(1, 3, True)
    assert Union(Interval(1, 3), Interval(2, 3)) == Interval(1, 3)
    assert Union(Interval(1, 3, False, True), Interval(2, 3)) == \
           Interval(1, 3)
    assert Union(Interval(1, 2, False, True), Interval(2, 3, True)) != \
           Interval(1, 3)
    assert Union(Interval(1, 2), S.EmptySet) == Interval(1, 2)
    assert Union(S.EmptySet) == S.EmptySet

    assert Union(Interval(0,1), [FiniteSet(1.0/n) for n in range(1,10)]) == \
            Interval(0,1)

    assert Interval(1, 2).union(Interval(2, 3)) == \
           Interval(1, 2) + Interval(2, 3)

    assert Interval(1, 2).union(Interval(2, 3)) == Interval(1, 3)

    assert Union(Set()) == Set()

    assert FiniteSet(1) + FiniteSet(2) + FiniteSet(3) == FiniteSet(1,2,3)
    assert FiniteSet(['ham']) + FiniteSet(['eggs']) == FiniteSet('ham', 'eggs')
    assert FiniteSet(1,2,3) + S.EmptySet == FiniteSet(1,2,3)

    assert FiniteSet(1,2,3) & FiniteSet(2,3,4) == FiniteSet(2,3)
    assert FiniteSet(1,2,3) | FiniteSet(2,3,4) == FiniteSet(1,2,3,4)

    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    assert S.EmptySet | FiniteSet(x, FiniteSet(y, z)) == \
           FiniteSet(x, FiniteSet(y, z))

    # Test that Intervals and FiniteSets play nicely
    assert Interval(1,3) + FiniteSet(2) == Interval(1,3)
    assert Interval(1,3, True,True) + FiniteSet(3) == Interval(1,3, True,False)
    X = Interval(1,3)+FiniteSet(5)
    Y = Interval(1,2)+FiniteSet(3)
    XandY = X.intersect(Y)
    assert 2 in X and 3 in X and 3 in XandY
    assert X.subset(XandY) and Y.subset(XandY)

    raises(TypeError, lambda: Union(1, 2, 3))

    assert X.is_iterable == False

def test_difference():
    assert Interval(1, 3) - Interval(1, 2) == Interval(2, 3, True)
    assert Interval(1, 3) - Interval(2, 3) == Interval(1, 2, False, True)
    assert Interval(1, 3, True) - Interval(2, 3) == Interval(1, 2, True, True)
    assert Interval(1, 3, True) - Interval(2, 3, True) == \
           Interval(1, 2, True, False)
    assert Interval(0,2) - FiniteSet(1) == \
            Union(Interval(0,1,False, True), Interval(1,2, True, False))

    assert FiniteSet(1,2,3) - FiniteSet(2) == FiniteSet(1,3)
    assert FiniteSet('ham', 'eggs') - FiniteSet(['eggs']) == FiniteSet(['ham'])
    assert FiniteSet(1,2,3,4) - Interval(2,10, True, False) == FiniteSet(1,2)
    assert FiniteSet(1,2,3,4) - S.EmptySet == FiniteSet(1,2,3,4)
    assert Union(Interval(0,2), FiniteSet(2,3,4)) - Interval(1,3) == \
            Union(Interval(0,1,False, True), FiniteSet(4))

def test_complement():
    assert Interval(0, 1).complement == \
           Union(Interval(-oo, 0, True, True), Interval(1, oo, True, True))
    assert Interval(0, 1, True, False).complement == \
           Union(Interval(-oo, 0, True, False), Interval(1, oo, True, True))
    assert Interval(0, 1, False, True).complement == \
           Union(Interval(-oo, 0, True, True), Interval(1, oo, False, True))
    assert Interval(0, 1, True, True).complement == \
           Union(Interval(-oo, 0, True, False), Interval(1, oo, False, True))

    assert -S.EmptySet == S.EmptySet.complement
    assert ~S.EmptySet == S.EmptySet.complement

    assert S.EmptySet.complement == S.UniversalSet
    assert S.UniversalSet.complement == S.EmptySet

    assert Union(Interval(0, 1), Interval(2, 3)).complement == \
           Union(Interval(-oo, 0, True, True), Interval(1, 2, True, True),
                 Interval(3, oo, True, True))

    assert FiniteSet(0).complement == Union(Interval(-oo,0, True,True) ,
            Interval(0,oo, True, True))

    assert (FiniteSet(5) + Interval(S.NegativeInfinity, 0)).complement == \
            Interval(0, 5, True, True) + Interval(5, S.Infinity, True,True)

    assert FiniteSet(1,2,3).complement == Interval(S.NegativeInfinity,1, True,True) + Interval(1,2, True,True) + Interval(2,3, True,True) + Interval(3,S.Infinity, True,True)

    X = Interval(1,3)+FiniteSet(5)
    assert X.intersect(X.complement) == S.EmptySet

    square = Interval(0,1) * Interval(0,1)
    notsquare = square.complement

    assert all(pt in square for pt in [(0,0), (.5,.5), (1,0), (1,1)])
    assert not any(pt in notsquare for pt in [(0,0), (.5,.5), (1,0), (1,1)])
    assert not any(pt in square for pt in [(-1,0), (1.5,.5), (10,10)])
    assert all(pt in notsquare for pt in [(-1,0), (1.5,.5), (10,10)])

def test_intersect():
    x = Symbol('x')
    assert Interval(0, 2).intersect(Interval(1, 2)) == Interval(1, 2)
    assert Interval(0, 2).intersect(Interval(1, 2, True)) == \
           Interval(1, 2, True)
    assert Interval(0, 2, True).intersect(Interval(1, 2)) == \
           Interval(1, 2, False, False)
    assert Interval(0, 2, True, True).intersect(Interval(1, 2)) == \
           Interval(1, 2, False, True)
    assert Interval(0, 2).intersect(Union(Interval(0, 1), Interval(2, 3))) == \
           Union(Interval(0, 1), Interval(2, 2))

    assert FiniteSet(1,2,x).intersect(FiniteSet(x)) == FiniteSet(x)
    assert FiniteSet('ham', 'eggs').intersect(FiniteSet(['ham'])) == \
            FiniteSet(['ham'])
    assert FiniteSet(1,2,3,4,5).intersect(S.EmptySet) == S.EmptySet

    assert Interval(0,5).intersect(FiniteSet(1,3)) == FiniteSet(1,3)
    assert Interval(0,1, True, True).intersect(FiniteSet(1)) == S.EmptySet

    assert Union(Interval(0, 1), Interval(2, 3)).intersect(Interval(1, 2)) == \
           Union(Interval(1, 1), Interval(2, 2))
    assert Union(Interval(0, 1), Interval(2, 3)).intersect(Interval(0, 2)) == \
           Union(Interval(0, 1), Interval(2, 2))
    assert Union(Interval(0, 1), Interval(2, 3)).intersect(Interval(1, 2, True, True)) == \
           S.EmptySet
    assert Union(Interval(0, 1), Interval(2, 3)).intersect(S.EmptySet) == \
           S.EmptySet
    assert Union(Interval(0,5), FiniteSet(['Ham'])).intersect(FiniteSet(2,3,4,5,6)) == \
           FiniteSet(2,3,4,5)

def test_intersection():
    # iterable
    i = Intersection(FiniteSet(1,2,3), Interval(2, 5), evaluate=False)
    assert i.is_iterable
    assert set(i) == set([S(2), S(3)])

    # challenging intervals
    x = Symbol('x', real=True)
    i = Intersection(Interval(0, 3), Interval(x, 6))
    assert (5 in i) == False
    raises(TypeError, lambda: 2 in i)

    # Singleton special cases
    assert Intersection(Interval(0, 1), S.EmptySet) == S.EmptySet
    assert Intersection(Interval(0, 1), S.UniversalSet) == Interval(0, 1)

    # Products
    line = Interval(0, 5)
    i = Intersection(line**2, line**3, evaluate=False)
    assert (2,2) not in i
    assert (2,2,2) not in i
    raises(ValueError, lambda: list(i))

def test_interval_subs():
    a = Symbol('a', real=True)

    assert Interval(0, a).subs(a, 2) == Interval(0, 2)
    assert Interval(a, 0).subs(a, 2) == S.EmptySet

def test_interval_to_mpi():
    assert Interval(0, 1).to_mpi() == mpi(0, 1)
    assert Interval(0, 1, True, False).to_mpi() == mpi(0, 1)
    assert type(Interval(0,1).to_mpi()) == type(mpi(0,1))

def test_measure():
    a = Symbol('a', real=True)

    assert Interval(1, 3).measure == 2
    assert Interval(0, a).measure == a
    assert Interval(1, a).measure == a - 1

    assert Union(Interval(1, 2), Interval(3, 4)).measure == 2
    assert Union(Interval(1, 2), Interval(3, 4), FiniteSet(5,6,7)).measure == 2

    assert FiniteSet(1,2,oo,a,-oo,-5).measure == 0

    assert S.EmptySet.measure == 0

    square = Interval(0,10) * Interval(0,10)
    offsetsquare = Interval(5,15) * Interval(5,15)
    band = Interval(-oo,oo) * Interval(2,4)

    assert square.measure == offsetsquare.measure == 100
    assert (square + offsetsquare).measure == 175 # there is some overlap
    assert (square - offsetsquare).measure == 75
    assert (square * FiniteSet(1,2,3)).measure == 0
    assert (square.intersect(band)).measure == 20
    assert (square + band).measure == oo
    assert (band * FiniteSet(1,2,3)).measure == nan

def test_subset():
    assert Interval(0, 2).subset(Interval(0, 1)) == True
    assert Interval(0, 2).subset(Interval(0, 3)) == False

    assert FiniteSet(1,2,3,4).subset(FiniteSet(1,2))
    assert FiniteSet(1,2,3,4).subset(FiniteSet(4,5)) == False
    assert Interval(0,2).subset(FiniteSet(1))
    assert Interval(0,2,True,True).subset(FiniteSet(1,2)) == False
    assert (Interval(0,2,False,True)+FiniteSet(2,3)).subset(
            Interval(1,2)+FiniteSet(3))

    assert Union(Interval(0, 1), Interval(2, 5)).subset(Interval(3, 4)) == True
    assert Union(Interval(0, 1), Interval(2, 5)).subset(Interval(3, 6)) == False

    assert Interval(0,5).subset(FiniteSet(1,2,3,4)) == True
    assert FiniteSet(1,2,3).subset(S.EmptySet) == True

    assert S.EmptySet.subset(Interval(0, 1)) == False
    assert S.EmptySet.subset(S.EmptySet) == True

    raises(ValueError, lambda: S.EmptySet.subset(1))

def test_contains():
    assert Interval(0, 2).contains(1) == True
    assert Interval(0, 2).contains(3) == False
    assert Interval(0, 2, True, False).contains(0) == False
    assert Interval(0, 2, True, False).contains(2) == True
    assert Interval(0, 2, False, True).contains(0) == True
    assert Interval(0, 2, False, True).contains(2) == False
    assert Interval(0, 2, True, True).contains(0) == False
    assert Interval(0, 2, True, True).contains(2) == False

    assert FiniteSet(1,2,3).contains(2)
    assert FiniteSet(1,2,Symbol('x')).contains(Symbol('x'))

    items = [1, 2, S.Infinity, S('ham'), -1.1]
    fset = FiniteSet(*items)
    assert all(item in fset for item in items)
    assert all(fset.contains(item) is True for item in items)

    assert Union(Interval(0, 1), Interval(2, 5)).contains(3) == True
    assert Union(Interval(0, 1), Interval(2, 5)).contains(6) == False
    assert Union(Interval(0, 1), FiniteSet(2, 5)).contains(3) == False

    assert S.EmptySet.contains(1) == False

def test_interval_symbolic():
    x = Symbol('x')
    e = Interval(0, 1)
    assert e.contains(x) == And(0<=x, x<=1)
    raises(TypeError, lambda: x in e)
    e = Interval(0, 1, True, True)
    assert e.contains(x) == And(0<x, x<1)

def test_union_contains():
    x = Symbol('x')
    i1 = Interval(0, 1)
    i2 = Interval(2, 3)
    i3 = Union(i1, i2)
    raises(TypeError, lambda: x in i3)
    e = i3.contains(x)
    assert e == Or(And(0 <= x, x <= 1), And(2 <= x, x <= 3))
    assert e.subs(x, -0.5) is False
    assert e.subs(x, 0.5) is True
    assert e.subs(x, 1.5) is False
    assert e.subs(x, 2.5) is True
    assert e.subs(x, 3.5) is False

    U = Interval(0,2, True,True) + Interval(10,oo) + FiniteSet(-1,2,5,6)
    assert all(el not in U for el in [0,4,-oo])
    assert all(el in U for el in [2,5,10])

def test_is_number():
    assert Interval(0, 1).is_number is False
    assert Set().is_number is False

def test_Interval_is_left_unbounded():
    assert Interval(3, 4).is_left_unbounded == False
    assert Interval(-oo, 3).is_left_unbounded == True
    assert Interval(Float("-inf"), 3).is_left_unbounded == True

def test_Interval_is_right_unbounded():
    assert Interval(3, 4).is_right_unbounded == False
    assert Interval(3, oo).is_right_unbounded == True
    assert Interval(3, Float("+inf")).is_right_unbounded == True

def test_Interval_as_relational():
    x = Symbol('x')

    assert Interval(-1, 2, False, False).as_relational(x) == And(Le(-1, x), Le(x, 2))
    assert Interval(-1, 2, True, False).as_relational(x) == And(Lt(-1, x), Le(x, 2))
    assert Interval(-1, 2, False, True).as_relational(x) == And(Le(-1, x), Lt(x, 2))
    assert Interval(-1, 2, True, True).as_relational(x) == And(Lt(-1, x), Lt(x, 2))

    assert Interval(-oo, 2, right_open=False).as_relational(x) == Le(x, 2)
    assert Interval(-oo, 2, right_open=True).as_relational(x) == Lt(x, 2)

    assert Interval(-2, oo, left_open=False).as_relational(x) == Le(-2, x)
    assert Interval(-2, oo, left_open=True).as_relational(x) == Lt(-2, x)

    assert Interval(-oo, oo).as_relational(x) == True

def test_Finite_as_relational():
    x = Symbol('x')
    y = Symbol('y')

    assert FiniteSet(1,2).as_relational(x) == Or(Eq(x,1), Eq(x,2))
    assert FiniteSet(y,-5).as_relational(x) == Or(Eq(x,y), Eq(x,-5))

def test_Union_as_relational():
    x = Symbol('x')
    assert (Interval(0,1) + FiniteSet(2)).as_relational(x) ==\
            Or(And(Le(0, x), Le(x, 1)), Eq(x, 2))
    assert (Interval(0,1, True, True) + FiniteSet(1)).as_relational(x) ==\
            And(Lt(0, x), Le(x, 1))

def test_Intersection_as_relational():
    x = Symbol('x')
    assert (Intersection(Interval(0,1), FiniteSet(2),
            evaluate=False).as_relational(x)
            == And(And(Le(0, x), Le(x, 1)), Eq(x, 2)))

def test_EmptySet_as_relational():
    assert S.EmptySet.as_relational(Symbol('x')) == False

def test_finite_basic():
    x = Symbol('x')
    A = FiniteSet(1,2,3)
    B = FiniteSet(3,4,5)
    AorB = Union(A,B)
    AandB = A.intersect(B)
    assert AorB.subset(A) and AorB.subset(B)
    assert A.subset(AandB)
    assert AandB == FiniteSet(3)

    assert A.inf == 1 and A.sup == 3
    assert AorB.inf == 1 and AorB.sup ==5
    assert FiniteSet(x, 1, 5).sup == Max(x,5)
    assert FiniteSet(x, 1, 5).inf == Min(x,1)

    # Ensure a variety of types can exist in a FiniteSet
    S = FiniteSet((1,2), Float, A, -5, x, 'eggs', x**2, Interval)

def test_product_basic():
    H,T = 'H', 'T'
    unit_line = Interval(0,1)
    d6 = FiniteSet(1,2,3,4,5,6)
    d4 = FiniteSet(1,2,3,4)
    coin = FiniteSet(H, T)

    square = unit_line * unit_line

    assert (0,0) in square
    assert 0 not in square
    assert (H, T) in coin ** 2
    assert (.5,.5,.5) in square * unit_line
    assert (H, 3, 3) in coin * d6* d6
    HH, TT = sympify(H), sympify(T)
    assert set(coin**2) == set(((HH, HH), (HH, TT), (TT, HH), (TT, TT)))

    assert (d6*d6).subset(d4*d4)

    inf, neginf = S.Infinity, S.NegativeInfinity
    assert square.complement == Union(
       Interval(0,1) * (Interval(neginf,0,True,True)+Interval(1,inf,True,True)),
       (Interval(neginf,0,True,True)+Interval(1,inf,True,True))*Interval(0,1),
       ((Interval(neginf,0,True,True) + Interval(1,inf, True, True))
                * (Interval(neginf,0,True,True) + Interval(1,inf, True,True))))

    assert (Interval(-10,10)**3).subset(Interval(-5,5)**3)
    assert not (Interval(-5,5)**3).subset(Interval(-10,10)**3)
    assert not (Interval(-10,10)**2).subset(Interval(-5,5)**3)

    assert square.subset(Interval(.2,.5)*FiniteSet(.5)) # segment in square

def test_real():
    x = Symbol('x', real=True)

    I = Interval(0, 5)
    J = Interval(10, 20)
    A = FiniteSet(1, 2, 30, x, S.Pi, S.Infinity)
    B = FiniteSet(-4, 0)
    C = FiniteSet(100, S.NegativeInfinity)
    D = FiniteSet('Ham', 'Eggs')

    assert all(s.is_real for s in [I,J,A,B,C])
    assert not D.is_real
    assert all((a+b).is_real for a in [I,J,A,B,C] for b in [I,J,A,B,C])
    assert not any((a+D).is_real for a in [I,J,A,B,C,D])

    assert not (I+A+D).is_real

def test_supinf():
    x = Symbol('x', real=True)
    y = Symbol('y', real=True)

    assert (Interval(0,1) + FiniteSet(2)).sup == 2
    assert (Interval(0,1) + FiniteSet(2)).inf == 0
    assert (Interval(0,1) + FiniteSet(x)).sup == Max(1,x)
    assert (Interval(0,1) + FiniteSet(x)).inf == Min(0,x)
    assert FiniteSet(5,1,x).sup == Max(5,x)
    assert FiniteSet(5,1,x).inf == Min(1,x)
    assert FiniteSet(5,1,x,y).sup == Max(5,x,y)
    assert FiniteSet(5,1,x,y).inf == Min(1,x,y)
    assert FiniteSet(5,1,x,y,S.Infinity, S.NegativeInfinity).sup == S.Infinity
    assert FiniteSet(5,1,x,y,S.Infinity, S.NegativeInfinity).inf == S.NegativeInfinity
    assert FiniteSet('Ham', 'Eggs').sup == Max('Ham', 'Eggs')

def test_universalset():
    U = S.UniversalSet
    x = Symbol('x')
    assert U.as_relational(x) == True
    assert U.union(Interval(2,4)) == U
