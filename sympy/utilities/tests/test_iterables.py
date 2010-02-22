from sympy import symbols, Integral
from sympy.utilities.iterables import postorder_traversal, \
    preorder_traversal, flatten, subsets, variations
from sympy.functions.elementary.piecewise import Piecewise, ExprCondPair

w,x,y,z= symbols('wxyz')

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
    assert flatten( (1,(1,)) ) == [1,1]
    assert flatten( (x,(x,)) ) == [x,x]

    from sympy.core.basic import Basic
    class MyOp(Basic):
        pass
    assert flatten( [MyOp(x, y), z]) == [MyOp(x, y), z]
    assert flatten( [MyOp(x, y), z], cls=MyOp) == [x, y, z]


def test_subsets():
    assert list(subsets([1, 2, 3], 1)) == [[1], [2], [3]]
    assert list(subsets([1, 2, 3], 2)) == [[1, 2], [1,3], [2, 3]]
    assert list(subsets([1, 2, 3], 3)) == [[1, 2, 3]]

def test_variations():
    assert variations([1,2,3], 2) == [[1, 2], [1, 3], [2, 1], [2, 3], [3, 1], [3, 2]]
    assert variations([1,2,3], 2, True) == [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], \
                        [3,1], [3,2], [3,3]]
