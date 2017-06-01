from sympy import (
    Float, Idx, IndexedBase, Integer, Matrix, MatrixSymbol, Range, sin, symbols, Tuple
)
from sympy.core.relational import Relational
from sympy.utilities.pytest import raises


from sympy.codegen.ast import (
    Assignment, aug_assign, CodeBlock, For, Type, Variable, Pointer, Declaration,
    AddAugmentedAssignment, SubAugmentedAssignment, MulAugmentedAssignment,
    DivAugmentedAssignment, ModAugmentedAssignment
)

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


def test_Type():
    t = Type('float64')
    assert t.name == 'float64'
    assert t.from_expr(7) == Type('integer')


def test_Type__from_expr():
    assert Type.from_expr(None, i) == Type('integer')
    assert Type.from_expr(None, x) == Type('real')
    assert Type.from_expr(x, n) == Type('integer')
    assert Type.from_expr(3, x) == Type('integer')
    assert Type.from_expr(3.0, x) == Type('real')
    assert Type.from_expr(3) == Type('integer')
    assert Type.from_expr(3+1j) == Type('complex')


def test_Type__cast_check__integers():
    # Rounding
    integer = Type('integer')
    raises(ValueError, lambda: integer.cast_check(3.5))

    # Range
    int8 = Type('int8')
    assert int8.cast_check(127.0) is 127
    raises(ValueError, lambda: int8.cast_check(128))
    assert int8.cast_check(-128) is -128
    raises(ValueError, lambda: int8.cast_check(-129))

    uint8 = Type('uint8')
    assert uint8.cast_check(0) is 0
    assert uint8.cast_check(128) is 128
    raises(ValueError, lambda: uint8.cast_check(256.0))
    raises(ValueError, lambda: uint8.cast_check(-1))

def test_Type__cast_check__floating_point():
    f32 = Type('float32')
    assert abs(0.1234567890499 - f32.cast_check(0.1234567890499) - 4.99e-11) < 1e-16

    f64 = Type('float64')
    dcm21 = Float('0.123456789012345604999')  # 21 decimals
    assert abs(dcm21 - f64.cast_check(dcm21)) < 4999e-21


def test_Variable():
    v = Variable(x)
    assert v.symbol == x
    assert v.type == Type('real')
    assert v.const == False
    w = Variable(y, Type('float32'), True)
    assert w.symbol == y
    assert w.type == Type('float32')
    assert w.const
    v_n = Variable(n)
    assert v_n.type == Type('integer')
    v_i = Variable(i)
    assert v_i.type == Type('integer')


def test_Pointer():
    p = Pointer(x)
    assert p.symbol == x
    assert p.type == Type('real')
    assert not p.value_const
    assert not p.pointer_const
    assert not p.restrict
    pB = Pointer(B)
    assert pB.symbol is B
    assert pB.type == Type('real')


def test_Declaration():
    assert Declaration(x).variable.type.name == 'real'
    assert Declaration(n).variable.type.name == 'integer'

    decl = Declaration(x, 3.0, True)
    assert decl.variable == Variable(x, const=True)
    assert type(decl.value) == Float
    assert decl.value == 3.0

    decl2 = Declaration(y, 3)
    assert decl2.variable == Variable(y, Type('integer'))
    assert decl2.value is Integer(3)

    decl3 = Declaration(i, 3.0)
    assert decl3.variable.type == Type('integer')
    assert decl3.value == 3.0

    decl4 = raises(ValueError, lambda: Declaration(n, 3.5))

    var = Variable(x)
    raises(ValueError, lambda: Declaration(var, const=True))
    assert Declaration(z, 3).variable.type == Type('integer')
    assert Declaration(z, 3.0).variable.type == Type('real')
    assert Declaration(z, 3.0+1j).variable.type == Type('complex')
    assert Declaration(B).variable == Pointer(B)
