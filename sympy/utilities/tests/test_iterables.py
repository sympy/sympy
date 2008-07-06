from sympy import symbols
from sympy.utilities.iterables import postorder_traversal, preorder_traversal


w,x,y,z= symbols('wxyz')

def test_postorder_traversal():
    expr = z+w*(x+y)
    expected = [z, w, y, x, x + y, w*(x + y), z + w*(x + y)]
    assert list(postorder_traversal(expr)) == expected


def test_preorder_traversal():
    expr = z+w*(x+y)
    expected = [z + w*(x + y), z, w*(x + y), w, x + y, y, x]
    assert list(preorder_traversal(expr)) == expected


