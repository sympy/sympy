from sympy import (
    Float, Idx, IndexedBase, Integer, Matrix, MatrixSymbol, Range, sin, symbols, Tuple
)
from sympy.core.relational import Relational
from sympy.utilities.pytest import raises


from sympy.codegen.ast import (
    Assignment, aug_assign, CodeBlock, For, Type, Variable, Pointer, Declaration,
    AddAugmentedAssignment, SubAugmentedAssignment, MulAugmentedAssignment,
    DivAugmentedAssignment, ModAugmentedAssignment, value_const, pointer_const,
    integer, real, complex_, int8, uint8, float16 as f16, float32 as f32,
    float64 as f64, float80 as f80, float128 as f128, complex64 as c64, complex128 as c128
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
    t = Type('MyType')
    assert t.name == 'MyType'
    assert str(t) == 'MyType'
    assert repr(t) == "Type(name=MyType)"


def test_Type_eq():
    t1 = Type('t1')
    t2 = Type('t2')
    assert t1 != t2
    assert t1 == t1 and t2 == t2
    t1b = Type('t1')
    assert t1 == t1b
    assert t2 != t1b


def test_Type__from_expr():
    assert Type.from_expr(i) == integer
    u = symbols('u', real=True)
    assert Type.from_expr(u) == real
    assert Type.from_expr(n) == integer
    assert Type.from_expr(3) == integer
    assert Type.from_expr(3.0) == real
    assert Type.from_expr(3+1j) == complex_
    raises(ValueError, lambda: Type.from_expr(sum))


def test_Type__cast_check__integers():
    # Rounding
    raises(ValueError, lambda: integer.cast_check(3.5))
    assert integer.cast_check('3') == 3
    assert integer.cast_check(Float('3.0000000000000000000')) == 3
    assert integer.cast_check(Float('3.0000000000000000001')) == 3  # unintuitive maybe?

    # Range
    assert int8.cast_check(127.0) == 127
    raises(ValueError, lambda: int8.cast_check(128))
    assert int8.cast_check(-128) == -128
    raises(ValueError, lambda: int8.cast_check(-129))

    assert uint8.cast_check(0) == 0
    assert uint8.cast_check(128) == 128
    raises(ValueError, lambda: uint8.cast_check(256.0))
    raises(ValueError, lambda: uint8.cast_check(-1))


def test_Variable():
    v = Variable(x, type_=Type('real'))
    assert v.symbol == x
    assert v.type == real
    assert v.value_const == False
    w = Variable(y, {value_const}, f32)
    assert w.symbol == y
    assert w.type == f32
    assert w.value_const
    v_n = Variable(n, type_=Type.from_expr(n))
    assert v_n.type == integer
    v_i = Variable(i, type_=Type.from_expr(n))
    assert v_i.type == integer

    a_i = Variable.deduced(i)
    assert a_i.type == integer


def test_Variable__deduced():
    v_i = Variable.deduced(i)
    assert v_i.type == integer


def test_Pointer():
    p = Pointer(x)
    assert p.symbol == x
    assert p.type == None
    assert not p.value_const
    assert not p.pointer_const
    u = symbols('u', real=True)
    py = Pointer(u, {value_const, pointer_const}, Type.from_expr(u))
    assert py.symbol is u
    assert py.type == real
    assert py.value_const
    assert py.pointer_const


def test_Declaration():
    u = symbols('u', real=True)
    vu = Variable(u, type_=Type.from_expr(u))
    assert Declaration(vu).variable.type == real
    vn = Variable(n, type_=Type.from_expr(n))
    assert Declaration(vn).variable.type == integer

    vuc = Variable(u, {value_const}, Type.from_expr(u))
    decl = Declaration(vuc, 3.0)
    assert decl.variable == vuc
    assert isinstance(decl.value, Float)
    assert decl.value == 3.0

    vy = Variable(y, type_=integer)
    decl2 = Declaration(vy, 3)
    assert decl2.variable == vy
    assert decl2.value is Integer(3)

    vi = Variable(i, type_=Type.from_expr(i))
    decl3 = Declaration(vi, 3.0)
    assert decl3.variable.type == integer
    assert decl3.value == 3.0

    decl4 = raises(ValueError, lambda: Declaration.deduced(n, 3.5, cast=True))



def test_Declaration__deduced():
    assert Declaration.deduced(n).variable.type == integer
    assert Declaration.deduced(z, 3).variable.type == integer
    assert Declaration.deduced(z, 3.0).variable.type == real
    assert Declaration.deduced(z, 3.0+1j).variable.type == complex_


def test_FloatType():
    assert f16.dig == 3
    assert f32.dig == 6
    assert f64.dig == 15
    assert f80.dig == 18
    assert f128.dig == 33

    assert f16.decimal_dig == 5
    assert f32.decimal_dig == 9
    assert f64.decimal_dig == 17
    assert f80.decimal_dig == 21
    assert f128.decimal_dig == 36

    assert f16.max_exponent == 16
    assert f32.max_exponent == 128
    assert f64.max_exponent == 1024
    assert f80.max_exponent == 16384
    assert f128.max_exponent == 16384

    assert f16.min_exponent == -13
    assert f32.min_exponent == -125
    assert f64.min_exponent == -1021
    assert f80.min_exponent == -16381
    assert f128.min_exponent == -16381

    assert abs(f16.eps / Float('0.00097656', precision=16) - 1) < 0.1*10**-f16.dig
    assert abs(f32.eps / Float('1.1920929e-07', precision=32) - 1) < 0.1*10**-f32.dig
    assert abs(f64.eps / Float('2.2204460492503131e-16', precision=64) - 1) < 0.1*10**-f64.dig
    assert abs(f80.eps / Float('1.08420217248550443401e-19', precision=80) - 1) < 0.1*10**-f80.dig
    assert abs(f128.eps / Float(' 1.92592994438723585305597794258492732e-34', precision=128) - 1) < 0.1*10**-f128.dig

    assert abs(f16.max / Float('65504', precision=16) - 1) < .1*10**-f16.dig
    assert abs(f32.max / Float('3.40282347e+38', precision=32) - 1) < 0.1*10**-f32.dig
    assert abs(f64.max / Float('1.79769313486231571e+308', precision=64) - 1) < 0.1*10**-f64.dig  # cf. np.finfo(np.float64).max
    assert abs(f80.max / Float('1.18973149535723176502e+4932', precision=80) - 1) < 0.1*10**-f80.dig
    assert abs(f128.max / Float('1.18973149535723176508575932662800702e+4932', precision=128) - 1) < 0.1*10**-f128.dig

    # cf. np.finfo(np.float32).tiny
    assert abs(f16.tiny / Float('6.1035e-05', precision=16) - 1) < 0.1*10**-f16.dig
    assert abs(f32.tiny / Float('1.17549435e-38', precision=32) - 1) < 0.1*10**-f32.dig
    assert abs(f64.tiny / Float('2.22507385850720138e-308', precision=64) - 1) < 0.1*10**-f64.dig
    assert abs(f80.tiny / Float('3.36210314311209350626e-4932', precision=80) - 1) < 0.1*10**-f80.dig
    assert abs(f128.tiny / Float('3.3621031431120935062626778173217526e-4932', precision=128) - 1) < 0.1*10**-f128.dig


def test_Type__cast_check__floating_point():
    raises(ValueError, lambda: f32.cast_check(123.45678949))
    raises(ValueError, lambda: f32.cast_check(12.345678949))
    raises(ValueError, lambda: f32.cast_check(1.2345678949))
    raises(ValueError, lambda: f32.cast_check(.12345678949))
    assert abs(123.456789049 - f32.cast_check(123.456789049) - 4.9e-8) < 1e-8
    assert abs(0.12345678904 - f32.cast_check(0.12345678904) - 4e-11) < 1e-11

    dcm21 = Float('0.123456789012345670499')  # 21 decimals
    assert abs(dcm21 - f64.cast_check(dcm21) - 4.99e-19) < 1e-19

    f80.cast_check(Float('0.12345678901234567890103', precision=88))
    raises(ValueError, lambda: f80.cast_check(Float('0.12345678901234567890149', precision=88)))

    v10 = 12345.67894
    raises(ValueError, lambda: f32.cast_check(v10))
    assert abs(Float(str(v10), precision=64+8) - f64.cast_check(v10)) < v10*1e-16

    assert abs(f32.cast_check(2147483647) - 2147483650) < 1


def test_Type__cast_check__complex_floating_point():
    val9_11 = 123.456789049 + 0.123456789049j
    raises(ValueError, lambda: c64.cast_check(.12345678949 + .12345678949j))
    assert abs(val9_11 - c64.cast_check(val9_11) - 4.9e-8) < 1e-8

    dcm21 = Float('0.123456789012345670499') + 1e-20j  # 21 decimals
    assert abs(dcm21 - c128.cast_check(dcm21) - 4.99e-19) < 1e-19
    v19 = Float('0.1234567890123456749') + 1j*Float('0.1234567890123456749')
    raises(ValueError, lambda: c128.cast_check(v19))
