from sympy import symbols
from sympy.utilities.iterables import postorder_traversal, \
    preorder_traversal, flatten, subsets


w,x,y,z= symbols('wxyz')

def test_postorder_traversal():
    expr = z+w*(x+y)
    expected1 = [z, w, y, x, x + y, w*(x + y), z + w*(x + y)]
    expected2 = [z, w, x, y, x + y, w*(x + y), z + w*(x + y)]
    expected3 = [w, y, x, x + y, w*(x + y), z, z + w*(x + y)]
    assert list(postorder_traversal(expr)) in [expected1, expected2, expected3]


def test_preorder_traversal():
    expr = z+w*(x+y)
    expected1 = [z + w*(x + y), z, w*(x + y), w, x + y, y, x]
    expected2 = [z + w*(x + y), z, w*(x + y), w, x + y, x, y]
    expected3 = [z + w*(x + y), w*(x + y), w, x + y, y, x, z]
    assert list(preorder_traversal(expr)) in [expected1, expected2, expected3]


def test_flatten():
    assert flatten( (1,(1,)) ) == [1,1]
    assert flatten( (x,(x,)) ) == [x,x]


def test_subsets():
    assert list(subsets([1, 2, 3], 1)) == [[1], [2], [3]]
    assert list(subsets([1, 2, 3], 2)) == [[1, 2], [1,3], [2, 3]]
    assert list(subsets([1, 2, 3], 3)) == [[1, 2, 3]]

