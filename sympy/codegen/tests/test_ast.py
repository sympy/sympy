from sympy import (symbols, MatrixSymbol, Matrix, IndexedBase, Idx, Range,
    Tuple, sin)
from sympy.core.relational import Relational
from sympy.utilities.pytest import raises


from sympy.codegen.ast import (Assignment, aug_assign, CodeBlock, For,
    AddAugmentedAssignment, SubAugmentedAssignment, MulAugmentedAssignment,
    DivAugmentedAssignment, ModAugmentedAssignment)

x, y, z, t, x0 = symbols("x, y, z, t, x0")
n = symbols("n", integer=True)
A = MatrixSymbol('A', 3, 1)
mat = Matrix([1, 2, 3])
B = IndexedBase('B')
i = Idx("i", n)


def test_Assignment():
    x, y = symbols("x, y")
    A = MatrixSymbol('A', 3, 1)
    mat = Matrix([1, 2, 3])
    B = IndexedBase('B')
    n = symbols("n", integer=True)
    i = Idx("i", n)
    # Here we just do things to show they don't error
    Assignment(x, y)
    Assignment(x, 0)
    Assignment(A, mat)
    Assignment(A[1,0], 0)
    Assignment(A[1,0], x)
    Assignment(B[i], x)
    Assignment(B[i], 0)
    a = Assignment(x, y)
    assert a.func(*a.args) == a
    # Here we test things to show that they error
    # Matrix to scalar
    raises(ValueError, lambda: Assignment(B[i], A))
    raises(ValueError, lambda: Assignment(B[i], mat))
    raises(ValueError, lambda: Assignment(x, mat))
    raises(ValueError, lambda: Assignment(x, A))
    raises(ValueError, lambda: Assignment(A[1,0], mat))
    # Scalar to matrix
    raises(ValueError, lambda: Assignment(A, x))
    raises(ValueError, lambda: Assignment(A, 0))
    # Non-atomic lhs
    raises(TypeError, lambda: Assignment(mat, A))
    raises(TypeError, lambda: Assignment(0, x))
    raises(TypeError, lambda: Assignment(x*x, 1))
    raises(TypeError, lambda: Assignment(A + A, mat))
    raises(TypeError, lambda: Assignment(B, 0))

    assert Relational(x, y, ':=') == Assignment(x, y)

def test_AugAssign():
    # Here we just do things to show they don't error
    aug_assign(x, '+', y)
    aug_assign(x, '+', 0)
    aug_assign(A, '+', mat)
    aug_assign(A[1, 0], '+', 0)
    aug_assign(A[1, 0], '+', x)
    aug_assign(B[i], '+', x)
    aug_assign(B[i], '+', 0)

    a = aug_assign(x, '+', y)
    b = AddAugmentedAssignment(x, y)
    assert a.func(*a.args) == a == b

    a = aug_assign(x, '-', y)
    b = SubAugmentedAssignment(x, y)
    assert a.func(*a.args) == a == b

    a = aug_assign(x, '*', y)
    b = MulAugmentedAssignment(x, y)
    assert a.func(*a.args) == a == b

    a = aug_assign(x, '/', y)
    b = DivAugmentedAssignment(x, y)
    assert a.func(*a.args) == a == b

    a = aug_assign(x, '%', y)
    b = ModAugmentedAssignment(x, y)
    assert a.func(*a.args) == a == b

    # Here we test things to show that they error
    # Matrix to scalar
    raises(ValueError, lambda: aug_assign(B[i], '+', A))
    raises(ValueError, lambda: aug_assign(B[i], '+', mat))
    raises(ValueError, lambda: aug_assign(x, '+', mat))
    raises(ValueError, lambda: aug_assign(x, '+', A))
    raises(ValueError, lambda: aug_assign(A[1, 0], '+', mat))
    # Scalar to matrix
    raises(ValueError, lambda: aug_assign(A, '+', x))
    raises(ValueError, lambda: aug_assign(A, '+', 0))
    # Non-atomic lhs
    raises(TypeError, lambda: aug_assign(mat, '+', A))
    raises(TypeError, lambda: aug_assign(0, '+', x))
    raises(TypeError, lambda: aug_assign(x * x, '+', 1))
    raises(TypeError, lambda: aug_assign(A + A, '+', mat))
    raises(TypeError, lambda: aug_assign(B, '+', 0))

def test_CodeBlock():
    c = CodeBlock(Assignment(x, 1), Assignment(y, x + 1))
    assert c.func(*c.args) == c

    assert c.left_hand_sides == Tuple(x, y)
    assert c.right_hand_sides == Tuple(1, x + 1)

def test_CodeBlock_topological_sort():
    assignments = [
        Assignment(x, y + z),
        Assignment(z, 1),
        Assignment(t, x),
        Assignment(y, 2),
        ]

    ordered_assignments = [
        # Note that the unrelated z=1 and y=2 are kept in that order
        Assignment(z, 1),
        Assignment(y, 2),
        Assignment(x, y + z),
        Assignment(t, x),
        ]
    c = CodeBlock.topological_sort(assignments)
    assert c == CodeBlock(*ordered_assignments)

    # Cycle
    invalid_assignments = [
        Assignment(x, y + z),
        Assignment(z, 1),
        Assignment(y, x),
        Assignment(y, 2),
        ]

    raises(ValueError, lambda: CodeBlock.topological_sort(invalid_assignments))

    # Undefined variable
    invalid_assignments = [
        Assignment(x, y)
        ]

    raises(ValueError, lambda: CodeBlock.topological_sort(invalid_assignments))

def test_CodeBlock_cse():
    c = CodeBlock(
        Assignment(y, 1),
        Assignment(x, sin(y)),
        Assignment(z, sin(y)),
        Assignment(t, x*z),
        )
    assert c.cse() == CodeBlock(
        Assignment(y, 1),
        Assignment(x0, sin(y)),
        Assignment(x, x0),
        Assignment(z, x0),
        Assignment(t, x*z),
        )

    raises(NotImplementedError, lambda: CodeBlock(Assignment(x, 1),
        Assignment(y, 1), Assignment(y, 2)).cse())

def test_For():
    f = For(n, Range(0, 3), (Assignment(A[n, 0], x + n), aug_assign(x, '+', y)))
    f = For(n, (1, 2, 3, 4, 5), (Assignment(A[n, 0], x + n),))
    assert f.func(*f.args) == f
    raises(TypeError, lambda: For(n, x, (x + y,)))
