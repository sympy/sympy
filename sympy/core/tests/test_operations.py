from sympy import Integer, S
from sympy.core.operations import LatticeOp
from sympy.utilities.pytest import raises
from sympy.core.sympify import SympifyError
from sympy.core.add import Add

# create the simplest possible Lattice class


class join(LatticeOp):
    zero = Integer(0)
    identity = Integer(1)


def test_lattice_simple():
    assert join(join(2, 3), 4) == join(2, join(3, 4))
    assert join(2, 3) == join(3, 2)
    assert join(0, 2) == 0
    assert join(1, 2) == 2
    assert join(2, 2) == 2

    assert join(join(2, 3), 4) == join(2, 3, 4)
    assert join() == 1
    assert join(4) == 4
    assert join(1, 4, 2, 3, 1, 3, 2) == join(2, 3, 4)


def test_lattice_shortcircuit():
    raises(SympifyError, lambda: join(object))
    assert join(0, object) == 0


def test_lattice_print():
    assert str(join(5, 4, 3, 2)) == 'join(2, 3, 4, 5)'


def test_lattice_make_args():
    assert join.make_args(join(2, 3, 4)) == {S(2), S(3), S(4)}
    assert join.make_args(0) == {0}
    assert list(join.make_args(0))[0] is S.Zero
    assert Add.make_args(0)[0] is S.Zero
