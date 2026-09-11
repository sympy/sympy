from __future__ import annotations
from sympy.core.relational import Ne
from sympy.core.symbol import (Symbol, symbols)
from sympy.functions.elementary.complexes import (adjoint, conjugate, transpose)
from sympy.functions.elementary.piecewise import Piecewise
from sympy.functions.special.tensor_functions import (Eijk, KroneckerDelta, LeviCivita)

x, y = symbols('x y')


def test_levicivita():
    assert Eijk(1, 2, 3) == LeviCivita(1, 2, 3)
    assert LeviCivita(1, 2, 3) == 1
    assert LeviCivita(int(1), int(2), int(3)) == 1
    assert LeviCivita(1, 3, 2) == -1
    assert LeviCivita(1, 2, 2) == 0
    i, j, k = symbols('i j k')
    assert LeviCivita(i, j, k) == LeviCivita(i, j, k, evaluate=False)
    assert LeviCivita(i, j, i) == 0
    assert LeviCivita(1, i, i) == 0
    assert LeviCivita(i, j, k).doit() == (j - i)*(k - i)*(k - j)/2
    assert LeviCivita(1, 2, 3, 1) == 0
    assert LeviCivita(4, 5, 1, 2, 3) == 1
    assert LeviCivita(4, 5, 2, 1, 3) == -1

    assert LeviCivita(i, j, k).is_integer is True

    assert adjoint(LeviCivita(i, j, k)) == LeviCivita(i, j, k)
    assert conjugate(LeviCivita(i, j, k)) == LeviCivita(i, j, k)
    assert transpose(LeviCivita(i, j, k)) == LeviCivita(i, j, k)


def test_kronecker_delta():
    i, j = symbols('i j')
    k = Symbol('k', nonzero=True)
    assert KroneckerDelta(1, 1) == 1
    assert KroneckerDelta(1, 2) == 0
    assert KroneckerDelta(k, 0) == 0
    assert KroneckerDelta(x, x) == 1
    assert KroneckerDelta(x**2 - y**2, x**2 - y**2) == 1
    assert KroneckerDelta(i, i) == 1
    assert KroneckerDelta(i, i + 1) == 0
    assert KroneckerDelta(0, 0) == 1
    assert KroneckerDelta(0, 1) == 0
    assert KroneckerDelta(i + k, i) == 0
    assert KroneckerDelta(i + k, i + k) == 1
    assert KroneckerDelta(i + k, i + 1 + k) == 0
    assert KroneckerDelta(i, j).subs({"i": 1, "j": 0}) == 0
    assert KroneckerDelta(i, j).subs({"i": 3, "j": 3}) == 1

    assert KroneckerDelta(i, j)**0 == 1
    for n in range(1, 10):
        assert KroneckerDelta(i, j)**n == KroneckerDelta(i, j)
        assert KroneckerDelta(i, j)**-n == 1/KroneckerDelta(i, j)

    assert KroneckerDelta(i, j).is_integer is True

    assert adjoint(KroneckerDelta(i, j)) == KroneckerDelta(i, j)
    assert conjugate(KroneckerDelta(i, j)) == KroneckerDelta(i, j)
    assert transpose(KroneckerDelta(i, j)) == KroneckerDelta(i, j)
    # to test if canonical
    assert (KroneckerDelta(i, j) == KroneckerDelta(j, i)) == True

    assert KroneckerDelta(i, j).rewrite(Piecewise) == Piecewise((0, Ne(i, j)), (1, True))

    # Tests with range:
    assert KroneckerDelta(i, j, (0, i)).args == (i, j, (0, i))
    assert KroneckerDelta(i, j, (-j, i)).delta_range == (-j, i)

    # If index is out of range, return zero:
    assert KroneckerDelta(i, j, (0, i-1)) == 0
    assert KroneckerDelta(-1, j, (0, i-1)) == 0
    assert KroneckerDelta(j, -1, (0, i-1)) == 0
    assert KroneckerDelta(j, i, (0, i-1)) == 0
