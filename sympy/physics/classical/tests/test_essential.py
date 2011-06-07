from sympy import symbols

from sympy.physics.classical.essential import Vector, ReferenceFrame
from sympy.physics.classical.functions import dot, cross

phi, x, y, z = symbols('phi x y z')
A = ReferenceFrame('A')

def test_orientnew():
    B = A.orientnew('B', 'Simple', phi, 1)
    assert B.parent is not None
    assert B.parent_orient is not None
    B = A.orientnew('B', 'Simple', phi, 2)
    assert B.parent is not None
    assert B.parent_orient is not None
    B = A.orientnew('B', 'Simple', phi, 3)
    assert B.parent is not None
    assert B.parent_orient is not None

def test_Vector():
    v1 = x*A.x + y*A.y + z*A.z
    v2 = x**2*A.x + y**2*A.y + z**2*A.z
    v3 = v1 + v2
    v4 = v1 - v2

    assert isinstance(v1, Vector)
    assert dot(v1, A.x) == x
    assert dot(v1, A.y) == y
    assert dot(v1, A.z) == z

    assert isinstance(v2, Vector)
    assert dot(v2, A.x) == x**2
    assert dot(v2, A.y) == y**2
    assert dot(v2, A.z) == z**2

    assert isinstance(v3, Vector)
    # We probably shouldn't be using simplify in dot...
    assert dot(v3, A.x) == x + x**2
    assert dot(v3, A.y) == y + y**2
    assert dot(v3, A.z) == z + z**2

    assert isinstance(v4, Vector)
    # We probably shouldn't be using simplify in dot...
    assert dot(v4, A.x) == x - x**2
    assert dot(v4, A.y) == y - y**2
    assert dot(v4, A.z) == z - z**2



    #TODO: Put some tests here, mew.
