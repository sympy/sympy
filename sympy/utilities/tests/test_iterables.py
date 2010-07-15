from sympy import symbols, Integral, Basic, raises
from sympy.utilities.iterables import (postorder_traversal, preorder_traversal,
    flatten, group, split, subsets, variations, numbered_symbols, take, dict_merge)
from sympy.functions.elementary.piecewise import Piecewise, ExprCondPair

w,x,y,z= symbols('w,x,y,z')

def test_postorder_traversal():
    expr = z+w*(x+y)
    expected1 = [z, w, y, x, x + y, w*(x + y), z + w*(x + y)]
    expected2 = [z, w, x, y, x + y, w*(x + y), z + w*(x + y)]
    expected3 = [w, y, x, x + y, w*(x + y), z, z + w*(x + y)]
    assert list(postorder_traversal(expr)) in [expected1, expected2, expected3]

    expr = Piecewise((x,x<1),(x**2,True))
    assert list(postorder_traversal(expr)) == [
        x, x, 1, x < 1, ExprCondPair(x, x < 1), x, 2, x**2, True,
        ExprCondPair(x**2, True), Piecewise((x, x < 1), (x**2, True))
    ]
    assert list(preorder_traversal(Integral(x**2, (x, 0, 1)))) == [
        Integral(x**2, (x, 0, 1)), x**2, x, 2, ((x, (0, 1)),), (x, (0, 1)),
        x, (0, 1), 0, 1
    ]
    assert list(preorder_traversal(('abc', ('d', 'ef')))) == [
        ('abc', ('d', 'ef')), 'abc', ('d', 'ef'), 'd', 'ef']



def test_preorder_traversal():
    expr = z+w*(x+y)
    expected1 = [z + w*(x + y), z, w*(x + y), w, x + y, y, x]
    expected2 = [z + w*(x + y), z, w*(x + y), w, x + y, x, y]
    expected3 = [z + w*(x + y), w*(x + y), w, x + y, y, x, z]
    assert list(preorder_traversal(expr)) in [expected1, expected2, expected3]

    expr = Piecewise((x,x<1),(x**2,True))
    assert list(preorder_traversal(expr)) == [
        Piecewise((x, x < 1), (x**2, True)), ExprCondPair(x, x < 1), x, x < 1,
        x, 1, ExprCondPair(x**2, True), x**2, x, 2, True
    ]
    assert list(postorder_traversal(Integral(x**2, (x, 0, 1)))) == [
        x, 2, x**2, x, 0, 1, (0, 1), (x, (0, 1)), ((x, (0, 1)),),
        Integral(x**2, (x, 0, 1))
    ]
    assert list(postorder_traversal(('abc', ('d', 'ef')))) == [
        'abc', 'd', 'ef', ('d', 'ef'), ('abc', ('d', 'ef'))]

def test_flatten():
    assert flatten((1, (1,))) == [1, 1]
    assert flatten((x, (x,))) == [x, x]

    ls = [[(-2, -1), (1, 2)], [(0, 0)]]

    assert flatten(ls, levels=0) == ls
    assert flatten(ls, levels=1) == [(-2, -1), (1, 2), (0, 0)]
    assert flatten(ls, levels=2) == [-2, -1, 1, 2, 0, 0]
    assert flatten(ls, levels=3) == [-2, -1, 1, 2, 0, 0]

    raises(ValueError, "flatten(ls, levels=-1)")

    class MyOp(Basic):
        pass

    assert flatten([MyOp(x, y), z]) == [MyOp(x, y), z]
    assert flatten([MyOp(x, y), z], cls=MyOp) == [x, y, z]

def test_group():
    assert group([]) == []
    assert group([], multiple=False) == []

    assert group([1]) == [[1]]
    assert group([1], multiple=False) == [(1, 1)]

    assert group([1,1]) == [[1,1]]
    assert group([1,1], multiple=False) == [(1, 2)]

    assert group([1,1,1]) == [[1,1,1]]
    assert group([1,1,1], multiple=False) == [(1, 3)]

    assert group([1,2,1]) == [[1],[2],[1]]
    assert group([1,2,1], multiple=False) == [(1, 1), (2, 1), (1, 1)]

    assert group([1,1,2,2,2,1,3,3]) == [[1,1], [2,2,2], [1], [3,3]]
    assert group([1,1,2,2,2,1,3,3], multiple=False) == [(1, 2), (2, 3), (1, 1), (3, 2)]

def test_split():
    assert split([], key=lambda a: a % 3) == []

    assert split([16, 8, 3, 1, 2, 5, 7], key=lambda a: a % 3) == [[3], [16, 1, 7], [8, 2, 5]]
    assert split([16, 8, 3, 7, 2, 5, 1], key=lambda a: a % 3) == [[3], [16, 7, 1], [8, 2, 5]]

def test_subsets():
    assert list(subsets([1, 2, 3], 1)) == [[1], [2], [3]]
    assert list(subsets([1, 2, 3], 2)) == [[1, 2], [1,3], [2, 3]]
    assert list(subsets([1, 2, 3], 3)) == [[1, 2, 3]]

def test_variations():
    assert variations([1,2,3], 2) == [[1, 2], [1, 3], [2, 1], [2, 3], [3, 1], [3, 2]]
    assert variations([1,2,3], 2, True) == [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], \
                        [3,1], [3,2], [3,3]]

def test_take():
    X = numbered_symbols()

    assert take(X, 5) == list(symbols('x0:5'))
    assert take(X, 5) == list(symbols('x5:10'))

def test_dict_merge():
    assert dict_merge({}, {1: x, y: z}) == {1: x, y: z}
    assert dict_merge({1: x, y: z}, {}) == {1: x, y: z}

    assert dict_merge({2: z}, {1: x, y: z}) == {1: x, 2: z, y: z}
    assert dict_merge({1: x, y: z}, {2: z}) == {1: x, 2: z, y: z}

    assert dict_merge({1: y, 2: z}, {1: x, y: z}) == {1: x, 2: z, y: z}
    assert dict_merge({1: x, y: z}, {1: y, 2: z}) == {1: y, 2: z, y: z}

