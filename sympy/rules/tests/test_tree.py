from sympy.rules.tree import (treeexec, treeexec, greedyexec, allexec,
        exhaustiveexec)
from functools import partial

def test_treeexec():
    tree = ([3, 3], [4, 1], 2)
    assert treeexec(tree, {list: min, tuple: max}) == 3

    add = lambda *args: sum(args)
    mul = lambda *args: reduce(lambda a, b: a*b, args, 1)
    assert treeexec(tree, {list: add, tuple: mul}) == 60

def test_treeexec_leaf():
    assert treeexec(3, {}, leaf=lambda x: x**2) == 9
    tree = ([3, 3], [4, 1], 2)
    treep1 = ([4, 4], [5, 2], 3)
    assert treeexec(tree, {list: min, tuple: max}, leaf=lambda x: x+1) == \
           treeexec(treep1, {list: min, tuple: max})

def test_treeexec_strategies():
    from sympy.rules import chain, minimize
    join = {list: chain, tuple: minimize}
    inc = lambda x: x + 1
    dec = lambda x: x - 1
    double = lambda x: 2*x

    assert treeexec(inc, join) == inc
    assert treeexec((inc, dec), join)(5) == minimize(inc, dec)(5)
    assert treeexec([inc, dec], join)(5) == chain(inc, dec)(5)
    tree = (inc, [dec, double]) # either inc or dec-then-double
    assert treeexec(tree, join)(5) == 6
    assert treeexec(tree, join)(1) == 0

    maximize = partial(minimize, objective=lambda x: -x)
    join = {list: chain, tuple: maximize}
    fn = treeexec(tree, join)
    assert fn(4) == 6  # highest value comes from the dec then double
    assert fn(1) == 2  # highest value comes from the inc

def test_greedyexec():
    inc = lambda x: x + 1
    dec = lambda x: x - 1
    double = lambda x: 2*x
    tree = (inc, [dec, double]) # either inc or dec-then-double

    fn = greedyexec(tree, objective=lambda x: -x)
    assert fn(4) == 6  # highest value comes from the dec then double
    assert fn(1) == 2  # highest value comes from the inc

    tree = (inc, dec, (inc, dec, ([inc, inc], [dec, dec])))
    lowest = greedyexec(tree)
    assert lowest(10) == 8

    highest = greedyexec(tree, objective=lambda x: -x)
    assert highest(10) == 12

def test_allexec():
    inc = lambda x: x+1
    dec = lambda x: x-1
    double = lambda x: x*2
    square = lambda x: x**2

    assert set(allexec(inc)(3)) == set([inc(3)])
    assert set(allexec((inc, dec))(3)) == set([2, 4])
    assert set(allexec([inc, dec])(3)) == set([3])
    assert set(allexec((inc, [dec, double]))(4)) == set([5, 6])

def test_exhaustiveexec():
    inc = lambda x: x+1
    dec = lambda x: x-1
    square = lambda x: x**2
    tree = [(inc, dec), square]
    fn = exhaustiveexec(tree, lambda x: -x)

    assert fn(2) == (2 + 1)**2
    assert fn(-2) == (-2 - 1)**2
