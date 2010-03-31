"""Tests for tools for real and complex root isolation and refinement. """

from sympy.polys.rootisolation import (
    dup_count_complex_roots,
)

from sympy.polys.algebratools import ZZ, QQ

a, b = (-QQ(1), -QQ(1)), (QQ(1), QQ(1))
c, d = ( QQ(0),  QQ(0)), (QQ(1), QQ(1))

def test_dup_count_complex_roots_1():
    # z-1
    assert dup_count_complex_roots([1,-1], a, b, ZZ) == 1
    assert dup_count_complex_roots([1,-1], c, d, ZZ) == 1

    # z+1
    assert dup_count_complex_roots([1,1], a, b, ZZ) == 1
    assert dup_count_complex_roots([1,1], c, d, ZZ) == 0

def test_dup_count_complex_roots_2():
    # (z-1)*(z)
    assert dup_count_complex_roots([1,-1,0], a, b, ZZ) == 2
    assert dup_count_complex_roots([1,-1,0], c, d, ZZ) == 2

    # (z-1)*(-z)
    assert dup_count_complex_roots([-1,1,0], a, b, ZZ) == 2
    assert dup_count_complex_roots([-1,1,0], c, d, ZZ) == 2

    # (z+1)*(z)
    assert dup_count_complex_roots([1,1,0], a, b, ZZ) == 2
    assert dup_count_complex_roots([1,1,0], c, d, ZZ) == 1

    # (z+1)*(-z)
    assert dup_count_complex_roots([-1,-1,0], a, b, ZZ) == 2
    assert dup_count_complex_roots([-1,-1,0], c, d, ZZ) == 1

def test_dup_count_complex_roots_3():
    # (z-1)*(z+1)
    assert dup_count_complex_roots([1,0,-1], a, b, ZZ) == 2
    assert dup_count_complex_roots([1,0,-1], c, d, ZZ) == 1

    # (z-1)*(z+1)*(z)
    assert dup_count_complex_roots([1,0,-1,0], a, b, ZZ) == 3
    assert dup_count_complex_roots([1,0,-1,0], c, d, ZZ) == 2

    # (z-1)*(z+1)*(-z)
    assert dup_count_complex_roots([-1,0,1,0], a, b, ZZ) == 3
    assert dup_count_complex_roots([-1,0,1,0], c, d, ZZ) == 2

def test_dup_count_complex_roots_4():
    # (z-I)*(z+I)
    assert dup_count_complex_roots([1,0,1], a, b, ZZ) == 2
    assert dup_count_complex_roots([1,0,1], c, d, ZZ) == 1

    # (z-I)*(z+I)*(z)
    assert dup_count_complex_roots([1,0,1,0], a, b, ZZ) == 3
    assert dup_count_complex_roots([1,0,1,0], c, d, ZZ) == 2

    # (z-I)*(z+I)*(-z)
    assert dup_count_complex_roots([-1,0,-1,0], a, b, ZZ) == 3
    assert dup_count_complex_roots([-1,0,-1,0], c, d, ZZ) == 2

    # (z-I)*(z+I)*(z-1)
    assert dup_count_complex_roots([1,-1,1,-1], a, b, ZZ) == 3
    assert dup_count_complex_roots([1,-1,1,-1], c, d, ZZ) == 2

    # (z-I)*(z+I)*(z-1)*(z)
    assert dup_count_complex_roots([1,-1,1,-1,0], a, b, ZZ) == 4
    assert dup_count_complex_roots([1,-1,1,-1,0], c, d, ZZ) == 3

    # (z-I)*(z+I)*(z-1)*(-z)
    assert dup_count_complex_roots([-1,1,-1,1,0], a, b, ZZ) == 4
    assert dup_count_complex_roots([-1,1,-1,1,0], c, d, ZZ) == 3

    # (z-I)*(z+I)*(z-1)*(z+1)
    assert dup_count_complex_roots([1,0,0,0,-1], a, b, ZZ) == 4
    assert dup_count_complex_roots([1,0,0,0,-1], c, d, ZZ) == 2

    # (z-I)*(z+I)*(z-1)*(z+1)*(z)
    assert dup_count_complex_roots([1,0,0,0,-1,0], a, b, ZZ) == 5
    assert dup_count_complex_roots([1,0,0,0,-1,0], c, d, ZZ) == 3

    # (z-I)*(z+I)*(z-1)*(z+1)*(-z)
    assert dup_count_complex_roots([-1,0,0,0,1,0], a, b, ZZ) == 5
    assert dup_count_complex_roots([-1,0,0,0,1,0], c, d, ZZ) == 3

def test_dup_count_complex_roots_5():
    # (z-I+1)*(z+I+1)
    assert dup_count_complex_roots([1,2,2], a, b, ZZ) == 2
    assert dup_count_complex_roots([1,2,2], c, d, ZZ) == 0

    # (z-I+1)*(z+I+1)*(z-1)
    assert dup_count_complex_roots([1,1,0,-2], a, b, ZZ) == 3
    assert dup_count_complex_roots([1,1,0,-2], c, d, ZZ) == 1

    # (z-I+1)*(z+I+1)*(z-1)*z
    assert dup_count_complex_roots([1,1,0,-2,0], a, b, ZZ) == 4
    assert dup_count_complex_roots([1,1,0,-2,0], c, d, ZZ) == 2

    # (z-I+1)*(z+I+1)*(z+1)
    assert dup_count_complex_roots([1,3,4,2], a, b, ZZ) == 3
    assert dup_count_complex_roots([1,3,4,2], c, d, ZZ) == 0

    # (z-I+1)*(z+I+1)*(z+1)*z
    assert dup_count_complex_roots([1,3,4,2,0], a, b, ZZ) == 4
    assert dup_count_complex_roots([1,3,4,2,0], c, d, ZZ) == 1

    # (z-I+1)*(z+I+1)*(z-1)*(z+1)
    assert dup_count_complex_roots([1,2,1,-2,-2], a, b, ZZ) == 4
    assert dup_count_complex_roots([1,2,1,-2,-2], c, d, ZZ) == 1

    # (z-I+1)*(z+I+1)*(z-1)*(z+1)*z
    assert dup_count_complex_roots([1,2,1,-2,-2,0], a, b, ZZ) == 5
    assert dup_count_complex_roots([1,2,1,-2,-2,0], c, d, ZZ) == 2

def test_dup_count_complex_roots_6():
    # (z-I-1)*(z+I-1)
    assert dup_count_complex_roots([1,-2,2], a, b, ZZ) == 2
    assert dup_count_complex_roots([1,-2,2], c, d, ZZ) == 1

    # (z-I-1)*(z+I-1)*(z-1)
    assert dup_count_complex_roots([1,-3,4,-2], a, b, ZZ) == 3
    assert dup_count_complex_roots([1,-3,4,-2], c, d, ZZ) == 2

    # (z-I-1)*(z+I-1)*(z-1)*z
    assert dup_count_complex_roots([1,-3,4,-2,0], a, b, ZZ) == 4
    assert dup_count_complex_roots([1,-3,4,-2,0], c, d, ZZ) == 3

    # (z-I-1)*(z+I-1)*(z+1)
    assert dup_count_complex_roots([1,-1,0,2], a, b, ZZ) == 3
    assert dup_count_complex_roots([1,-1,0,2], c, d, ZZ) == 1

    # (z-I-1)*(z+I-1)*(z+1)*z
    assert dup_count_complex_roots([1,-1,0,2,0], a, b, ZZ) == 4
    assert dup_count_complex_roots([1,-1,0,2,0], c, d, ZZ) == 2

    # (z-I-1)*(z+I-1)*(z-1)*(z+1)
    assert dup_count_complex_roots([1,-2,1,2,-2], a, b, ZZ) == 4
    assert dup_count_complex_roots([1,-2,1,2,-2], c, d, ZZ) == 2

    # (z-I-1)*(z+I-1)*(z-1)*(z+1)*z
    assert dup_count_complex_roots([1,-2,1,2,-2,0], a, b, ZZ) == 5
    assert dup_count_complex_roots([1,-2,1,2,-2,0], c, d, ZZ) == 3

def test_dup_count_complex_roots_7():
    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)
    assert dup_count_complex_roots([1,0,0,0,4], a, b, ZZ) == 4
    assert dup_count_complex_roots([1,0,0,0,4], c, d, ZZ) == 1

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-2)
    assert dup_count_complex_roots([1,-2,0,0,4,-8], a, b, ZZ) == 4
    assert dup_count_complex_roots([1,-2,0,0,4,-8], c, d, ZZ) == 1

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z**2-2)
    assert dup_count_complex_roots([1,0,-2,0,4,0,-8], a, b, ZZ) == 4
    assert dup_count_complex_roots([1,0,-2,0,4,0,-8], c, d, ZZ) == 1

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)
    assert dup_count_complex_roots([1,-1,0,0,4,-4], a, b, ZZ) == 5
    assert dup_count_complex_roots([1,-1,0,0,4,-4], c, d, ZZ) == 2

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)*z
    assert dup_count_complex_roots([1,-1,0,0,4,-4,0], a, b, ZZ) == 6
    assert dup_count_complex_roots([1,-1,0,0,4,-4,0], c, d, ZZ) == 3

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z+1)
    assert dup_count_complex_roots([1,1,0,0,4,4], a, b, ZZ) == 5
    assert dup_count_complex_roots([1,1,0,0,4,4], c, d, ZZ) == 1

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z+1)*z
    assert dup_count_complex_roots([1,1,0,0,4,4,0], a, b, ZZ) == 6
    assert dup_count_complex_roots([1,1,0,0,4,4,0], c, d, ZZ) == 2

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)*(z+1)
    assert dup_count_complex_roots([1,0,-1,0,4,0,-4], a, b, ZZ) == 6
    assert dup_count_complex_roots([1,0,-1,0,4,0,-4], c, d, ZZ) == 2

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)*(z+1)*z
    assert dup_count_complex_roots([1,0,-1,0,4,0,-4,0], a, b, ZZ) == 7
    assert dup_count_complex_roots([1,0,-1,0,4,0,-4,0], c, d, ZZ) == 3

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)*(z+1)*(z-I)*(z+I)
    assert dup_count_complex_roots([1,0,0,0,3,0,0,0,-4], a, b, ZZ) == 8
    assert dup_count_complex_roots([1,0,0,0,3,0,0,0,-4], c, d, ZZ) == 3

def test_dup_count_complex_roots_8():
    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)*(z+1)*(z-I)*(z+I)*z
    assert dup_count_complex_roots([1,0,0,0,3,0,0,0,-4,0], a, b, ZZ) == 9
    assert dup_count_complex_roots([1,0,0,0,3,0,0,0,-4,0], c, d, ZZ) == 4

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)*(z+1)*(z-I)*(z+I)*(z**2-2)*z
    assert dup_count_complex_roots([1,0,-2,0,3,0,-6,0,-4,0,8,0], a, b, ZZ) == 9
    assert dup_count_complex_roots([1,0,-2,0,3,0,-6,0,-4,0,8,0], c, d, ZZ) == 4

def test_dup_count_complex_roots_exclude():
    a, b = (-QQ(1), QQ(0)), (QQ(1), QQ(1))

    f = [1, 0, 0, 0, -1, 0] # z*(z-1)*(z+1)*(z-I)*(z+I)

    assert dup_count_complex_roots(f, a, b, ZZ) == 4

    assert dup_count_complex_roots(f, a, b, ZZ, exclude=['S']) == 3
    assert dup_count_complex_roots(f, a, b, ZZ, exclude=['N']) == 3

    assert dup_count_complex_roots(f, a, b, ZZ, exclude=['S', 'N']) == 2

    assert dup_count_complex_roots(f, a, b, ZZ, exclude=['E']) == 4
    assert dup_count_complex_roots(f, a, b, ZZ, exclude=['W']) == 4

    assert dup_count_complex_roots(f, a, b, ZZ, exclude=['E', 'W']) == 4

    assert dup_count_complex_roots(f, a, b, ZZ, exclude=['N', 'S', 'E', 'W']) == 2

    assert dup_count_complex_roots(f, a, b, ZZ, exclude=['SW']) == 3
    assert dup_count_complex_roots(f, a, b, ZZ, exclude=['SE']) == 3

    assert dup_count_complex_roots(f, a, b, ZZ, exclude=['SW', 'SE']) == 2
    assert dup_count_complex_roots(f, a, b, ZZ, exclude=['SW', 'SE', 'S']) == 1
    assert dup_count_complex_roots(f, a, b, ZZ, exclude=['SW', 'SE', 'S', 'N']) == 0
