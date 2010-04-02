"""Tests for tools for real and complex root isolation and refinement. """

from sympy.polys.rootisolation import (
    dup_isolate_complex_roots_sqf,
    dup_isolate_all_roots_sqf,
    dup_isolate_all_roots,
    dup_count_complex_roots,
)

from sympy.polys.algebratools import ZZ, QQ

from sympy.utilities.pytest import raises

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

    a, b = (QQ(0), QQ(0)), (QQ(1), QQ(1))

    assert dup_count_complex_roots(f, a, b, ZZ, exclude=True) == 1

def test_dup_isolate_complex_roots_sqf():
    f = [1, -2, 3]

    assert dup_isolate_complex_roots_sqf(f, ZZ) == \
        [((0, 0), (6, 6)), ((0, -6), (6, 0))]

    assert dup_isolate_complex_roots_sqf(f, ZZ, eps=QQ(1,10)) == \
        [((QQ(15,16), QQ(45,32)), (QQ(33,32), QQ(3,2))), ((QQ(15,16), -QQ(3,2)), (QQ(33,32), -QQ(45,32)))]
    assert dup_isolate_complex_roots_sqf(f, ZZ, eps=QQ(1,100)) == \
        [((QQ(255,256), QQ(723,512)), (QQ(513,512), QQ(363,256))), ((QQ(255,256), -QQ(363,256)), (QQ(513,512), -QQ(723,512)))]

    f = [7, -19, 20, 17, 20]

    assert dup_isolate_complex_roots_sqf(f, ZZ) == \
        [((-QQ(40,7), 0), (0, QQ(40,7))), ((-QQ(40,7), -QQ(40,7)), (0, 0)), ((0, 0), (QQ(40,7), QQ(40,7))), ((0, -QQ(40,7)), (QQ(40,7), 0))]

def test_dup_isolate_all_roots_sqf():
    f = [4, -1, 2, 5, 0]

    assert dup_isolate_all_roots_sqf(f, ZZ) == \
        ([(-1, 0), (0, 0)], [((0, 0), (QQ(5,2), QQ(5,2))), ((0, -QQ(5,2)), (QQ(5,2), 0))])

    assert dup_isolate_all_roots_sqf(f, ZZ, eps=QQ(1,10)) == \
        ([(QQ(-7,8), QQ(-6,7)), (0, 0)], [((QQ(35,64), QQ(65,64)), (QQ(5,8), QQ(35,32))), ((QQ(35,64), -QQ(35,32)), (QQ(5,8), -QQ(65,64)))])

def test_dup_isolate_all_roots():
    f = [4, -1, 2, 5, 0]

    assert dup_isolate_all_roots(f, ZZ) == \
        ([((-1, 0), 1), ((0, 0), 1)], [(((0, 0), (QQ(5,2), QQ(5,2))), 1), (((0, -QQ(5,2)), (QQ(5,2), 0)), 1)])

    assert dup_isolate_all_roots(f, ZZ, eps=QQ(1,10)) == \
        ([((QQ(-7,8), QQ(-6,7)), 1), ((0, 0), 1)], [(((QQ(35,64), QQ(65,64)), (QQ(5,8), QQ(35,32))), 1), (((QQ(35,64), -QQ(35,32)), (QQ(5,8), -QQ(65,64))), 1)])

    raises(NotImplementedError, "dup_isolate_all_roots([1, 1, -2, -2, 1, 1], ZZ)")
