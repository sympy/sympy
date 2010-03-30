"""Tests for tools for real and complex root isolation and refinement. """

from sympy.polys.rootisolation import (
    dup_count_complex_roots,
)

from sympy.polys.algebratools import ZZ, QQ

def test_dup_count_complex_roots():
    a, b = (QQ(-1), QQ(-1)), (QQ(1), QQ(1))

    # (z-1)*(z+1)
    assert dup_count_complex_roots([1,0,-1], a, b, ZZ) == 2

    # (z-1)*(z+1)*(z)
    assert dup_count_complex_roots([1,0,-1,0], a, b, ZZ) == 3

    # (z-1)*(z+1)*(-z)
    assert dup_count_complex_roots([-1,0,1,0], a, b, ZZ) == 3

    # (z-I)*(z+I)
    assert dup_count_complex_roots([1,0,1], a, b, ZZ) == 2

    # (z-I)*(z+I)*(z)
    assert dup_count_complex_roots([1,0,1,0], a, b, ZZ) == 3

    # (z-I)*(z+I)*(-z)
    assert dup_count_complex_roots([-1,0,-1,0], a, b, ZZ) == 3

    # (z-I)*(z+I)*(z-1)
    assert dup_count_complex_roots([1,-1,1,-1], a, b, ZZ) == 3

    # (z-I)*(z+I)*(z-1)*(z)
    assert dup_count_complex_roots([1,-1,1,-1,0], a, b, ZZ) == 4

    # (z-I)*(z+I)*(z-1)*(-z)
    assert dup_count_complex_roots([-1,1,-1,1,0], a, b, ZZ) == 4

    # (z-I)*(z+I)*(z-1)*(z+1)
    assert dup_count_complex_roots([1,0,0,0,-1], a, b, ZZ) == 4

    # (z-I)*(z+I)*(z-1)*(z+1)*(z)
    assert dup_count_complex_roots([1,0,0,0,-1,0], a, b, ZZ) == 5

    # (z-I)*(z+I)*(z-1)*(z+1)*(-z)
    assert dup_count_complex_roots([-1,0,0,0,1,0], a, b, ZZ) == 5

    # (z-I+1)*(z+I+1)
    assert dup_count_complex_roots([1,2,2], a, b, ZZ) == 2

    # (z-I+1)*(z+I+1)*(z-1)
    assert dup_count_complex_roots([1,1,0,-2], a, b, ZZ) == 3

    # (z-I+1)*(z+I+1)*(z-1)*z
    assert dup_count_complex_roots([1,1,0,-2,0], a, b, ZZ) == 4

    # (z-I+1)*(z+I+1)*(z+1)
    assert dup_count_complex_roots([1,3,4,2], a, b, ZZ) == 3

    # (z-I+1)*(z+I+1)*(z+1)*z
    assert dup_count_complex_roots([1,3,4,2,0], a, b, ZZ) == 4

    # (z-I+1)*(z+I+1)*(z-1)*(z+1)
    assert dup_count_complex_roots([1,2,1,-2,-2], a, b, ZZ) == 4

    # (z-I+1)*(z+I+1)*(z-1)*(z+1)*z
    assert dup_count_complex_roots([1,2,1,-2,-2,0], a, b, ZZ) == 5

    # (z-I-1)*(z+I-1)
    assert dup_count_complex_roots([1,-2,2], a, b, ZZ) == 2

    # (z-I-1)*(z+I-1)*(z-1)
    assert dup_count_complex_roots([1,-3,4,-2], a, b, ZZ) == 3

    # (z-I-1)*(z+I-1)*(z-1)*z
    assert dup_count_complex_roots([1,-3,4,-2,0], a, b, ZZ) == 4

    # (z-I-1)*(z+I-1)*(z+1)
    assert dup_count_complex_roots([1,-1,0,2], a, b, ZZ) == 3

    # (z-I-1)*(z+I-1)*(z+1)*z
    assert dup_count_complex_roots([1,-1,0,2,0], a, b, ZZ) == 4

    # (z-I-1)*(z+I-1)*(z-1)*(z+1)
    assert dup_count_complex_roots([1,-2,1,2,-2], a, b, ZZ) == 4

    # (z-I-1)*(z+I-1)*(z-1)*(z+1)*z
    assert dup_count_complex_roots([1,-2,1,2,-2,0], a, b, ZZ) == 5

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)
    assert dup_count_complex_roots([1,0,0,0,4], a, b, ZZ) == 4

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-2)
    assert dup_count_complex_roots([1,-2,0,0,4,-8], a, b, ZZ) == 4

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z**2-2)
    assert dup_count_complex_roots([1,0,-2,0,4,0,-8], a, b, ZZ) == 4

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)
    assert dup_count_complex_roots([1,-1,0,0,4,-4], a, b, ZZ) == 5

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)*z
    assert dup_count_complex_roots([1,-1,0,0,4,-4,0], a, b, ZZ) == 6

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z+1)
    assert dup_count_complex_roots([1,1,0,0,4,4], a, b, ZZ) == 5

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z+1)*z
    assert dup_count_complex_roots([1,1,0,0,4,4,0], a, b, ZZ) == 6

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)*(z+1)
    assert dup_count_complex_roots([1,0,-1,0,4,0,-4], a, b, ZZ) == 6

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)*(z+1)*z
    assert dup_count_complex_roots([1,0,-1,0,4,0,-4,0], a, b, ZZ) == 7

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)*(z+1)*(z-I)*(z+I)
    assert dup_count_complex_roots([1,0,0,0,3,0,0,0,-4], a, b, ZZ) == 8

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)*(z+1)*(z-I)*(z+I)*z
    assert dup_count_complex_roots([1,0,0,0,3,0,0,0,-4,0], a, b, ZZ) == 9

    # (z-I-1)*(z+I-1)*(z-I+1)*(z+I+1)*(z-1)*(z+1)*(z-I)*(z+I)*(z**2-2)*z
    assert dup_count_complex_roots([1,0,-2,0,3,0,-6,0,-4,0,8,0], a, b, ZZ) == 9
