from sympy import symbols, MatrixSymbol, Matrix, IndexedBase, Idx, Range
from sympy.core.relational import Relational
from sympy.utilities.pytest import raises


from sympy.codegen.ast import (Assignment, aug_assign, datatype, Bool, Int, Float,
        Double, Void, For, InArgument, OutArgument, InOutArgument, Variable,
        Result, Return, FunctionDef, AddAugmentedAssignment,
        SubAugmentedAssignment, MulAugmentedAssignment,
        DivAugmentedAssignment, ModAugmentedAssignment)

x, y = symbols("x, y")
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


def test_datatype():
    assert Bool == datatype('bool')
    assert Int == datatype('int')
    assert Float == datatype('float')
    assert Double == datatype('double')
    assert Void == datatype('void')
    # Check inferred types
    assert datatype(x) == Double
    assert datatype(n) == Int
    # This should work (I think), but doesn't due to how SymPy handles bools.
    # assert datatype(b) == Bool
    assert datatype(A) == Double
    assert datatype(mat) == Int
    d = datatype('int')
    assert d.func(*d.args) == d


def test_For():
    f = For(n, Range(0, 3), (Assignment(A[n, 0], x + n), aug_assign(x, '+', y)))
    f = For(n, (1, 2, 3, 4, 5), (Assignment(A[n, 0], x + n),))
    assert f.func(*f.args) == f
    raises(TypeError, lambda: For(n, x, (x + y,)))


def test_Variable():
    v = Variable('int', x)
    assert v.func(*v.args) == v
    Variable('double', A)
    raises(TypeError, lambda: Variable('int', x + y))


def test_Result():
    r = Result('double')
    assert r.func(*r.args) == r
    r = Result('double', 'var_name')
    assert r.func(*r.args) == r
    raises(TypeError, lambda: Result('int', x + y))


def test_Arguments():
    a = InArgument('int', x)
    b = OutArgument('double', x)
    c = InOutArgument('float', A)
    assert a.func(*a.args) == a
    assert b.func(*b.args) == b
    assert c.func(*c.args) == c


def test_FunctionDef():
    ax = InArgument('double', x)
    ay = InArgument('double', y)
    f = FunctionDef('test', (ax, ay), (Return(x + y),), (Result('double'),))
    assert f.func(*f.args) == f
    raises(TypeError, lambda: FunctionDef('test', (x, y), (Return(x + y),), (Result('double'),)))
    raises(TypeError, lambda: FunctionDef('test', (ax, ay), (Return(x + y),), (x + y,)))


def test_Return():
    r = Return(x + y)
    assert r.func(*r.args) == r
