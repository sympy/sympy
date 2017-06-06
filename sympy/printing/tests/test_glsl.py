from sympy.core import pi, oo, symbols, Rational, Integer, GoldenRatio, EulerGamma, Catalan, Lambda, Dummy
from sympy.functions import Piecewise, sin, cos, Abs, exp, ceiling, sqrt
from sympy.utilities.pytest import raises
from sympy.printing.glsl import GLSLPrinter
from sympy.utilities.lambdify import implemented_function
from sympy.tensor import IndexedBase, Idx
from sympy.matrices import Matrix, MatrixSymbol

from sympy import glsl_code

x, y, z = symbols('x,y,z')


def test_printmethod():
    assert glsl_code(Abs(x)) == "abs(x)"

def test_print_without_operators():
    var('x,y,z')
    assert glsl_code(x*y,use_operators = False) == 'mul(x, y)'
    assert glsl_code(x**y+z,use_operators = False) == 'add(pow(x, y), z)'
    assert glsl_code(x*(y+z),use_operators = False) == 'mul(x, add(y, z))'
    assert glsl_code(x*(y+z),use_operators = False) == 'mul(x, add(y, z))'
    assert glsl_code(x*(y+z**y**0.5),use_operators = False) == 'mul(x, add(y, pow(z, sqrt(y))))'
def test_glsl_code_sqrt():
    assert glsl_code(sqrt(x)) == "sqrt(x)"
    assert glsl_code(x**0.5) == "sqrt(x)"
    assert glsl_code(sqrt(x)) == "sqrt(x)"


def test_glsl_code_Pow():
    g = implemented_function('g', Lambda(x, 2*x))
    assert glsl_code(x**3) == "pow(x, 3.0)"
    assert glsl_code(x**(y**3)) == "pow(x, pow(y, 3.0))"
    assert glsl_code(1/(g(x)*3.5)**(x - y**x)/(x**2 + y)) == \
        "pow(3.5*2*x, -x + pow(y, x))/(pow(x, 2.0) + y)"
    assert glsl_code(x**-1.0) == '1.0/x'


def test_glsl_code_constants_mathh():
    assert glsl_code(exp(1)) == "float E = 2.7182818284590452353602874713527;\nE"
    assert glsl_code(pi) == "float pi = 3.1415926535897932384626433832795;\npi"
    # assert glsl_code(oo) == "Number.POSITIVE_INFINITY"
    # assert glsl_code(-oo) == "Number.NEGATIVE_INFINITY"


def test_glsl_code_constants_other():
    assert glsl_code(2*GoldenRatio) == "float GoldenRatio = 1.6180339887498948482045868343656;\n2*GoldenRatio"
    assert glsl_code(2*Catalan) == "float Catalan = 0.91596559417721901505460351493238;\n2*Catalan"
    assert glsl_code(2*EulerGamma) == "float EulerGamma = 0.5772156649015328606065120900824;\n2*EulerGamma"


def test_glsl_code_Rational():
    assert glsl_code(Rational(3, 7)) == "3/7"
    assert glsl_code(Rational(18, 9)) == "2"
    assert glsl_code(Rational(3, -7)) == "-3/7"
    assert glsl_code(Rational(-3, -7)) == "3/7"


def test_glsl_code_Integer():
    assert glsl_code(Integer(67)) == "67"
    assert glsl_code(Integer(-1)) == "-1"


def test_glsl_code_functions():
    assert glsl_code(sin(x) ** cos(x)) == "pow(sin(x), cos(x))"


def test_glsl_code_inline_function():
    x = symbols('x')
    g = implemented_function('g', Lambda(x, 2*x))
    assert glsl_code(g(x)) == "2*x"
    g = implemented_function('g', Lambda(x, 2*x/Catalan))
    assert glsl_code(g(x)) == "float Catalan = 0.91596559417721901505460351493238;\n2*x/Catalan"
    A = IndexedBase('A')
    i = Idx('i', symbols('n', integer=True))
    g = implemented_function('g', Lambda(x, x*(1 + x)*(2 + x)))
    assert glsl_code(g(A[i]), assign_to=A[i]) == (
        "for (int i=0; i<n; i++){\n"
        "   A[i] = (A[i] + 1)*(A[i] + 2)*A[i];\n"
        "}"
    )


def test_glsl_code_exceptions():
    assert glsl_code(ceiling(x)) == "ceil(x)"
    assert glsl_code(Abs(x)) == "abs(x)"


def test_glsl_code_boolean():
    assert glsl_code(x & y) == "x && y"
    assert glsl_code(x | y) == "x || y"
    assert glsl_code(~x) == "!x"
    assert glsl_code(x & y & z) == "x && y && z"
    assert glsl_code(x | y | z) == "x || y || z"
    assert glsl_code((x & y) | z) == "z || x && y"
    assert glsl_code((x | y) & z) == "z && (x || y)"


def test_glsl_code_Piecewise():
    expr = Piecewise((x, x < 1), (x**2, True))
    p = glsl_code(expr)
    s = \
"""\
((x < 1) ? (
   x
)
: (
   pow(x, 2.0)
))\
"""
    assert p == s
    assert glsl_code(expr, assign_to="c") == (
    "if (x < 1) {\n"
    "   c = x;\n"
    "}\n"
    "else {\n"
    "   c = pow(x, 2.0);\n"
    "}")
    # Check that Piecewise without a True (default) condition error
    expr = Piecewise((x, x < 1), (x**2, x > 1), (sin(x), x > 0))
    raises(ValueError, lambda: glsl_code(expr))


def test_glsl_code_Piecewise_deep():
    p = glsl_code(2*Piecewise((x, x < 1), (x**2, True)))
    s = \
"""\
2*((x < 1) ? (
   x
)
: (
   pow(x, 2.0)
))\
"""
    assert p == s


def test_glsl_code_settings():
    raises(TypeError, lambda: glsl_code(sin(x), method="garbage"))


def test_glsl_code_Indexed():
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols
    n, m, o = symbols('n m o', integer=True)
    i, j, k = Idx('i', n), Idx('j', m), Idx('k', o)
    p = GLSLPrinter()
    p._not_c = set()

    x = IndexedBase('x')[j]
    assert p._print_Indexed(x) == 'x[j]'
    A = IndexedBase('A')[i, j]
    assert p._print_Indexed(A) == 'A[%s]' % (m*i+j)
    B = IndexedBase('B')[i, j, k]
    assert p._print_Indexed(B) == 'B[%s]' % (i*o*m+j*o+k)

    assert p._not_c == set()


def test_glsl_code_loops_matrix_vector():
    n, m = symbols('n m', integer=True)
    A = IndexedBase('A')
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx('i', m)
    j = Idx('j', n)

    s = (
        'for (int i=0; i<m; i++){\n'
        '   y[i] = 0.0;\n'
        '}\n'
        'for (int i=0; i<m; i++){\n'
        '   for (int j=0; j<n; j++){\n'
        '      y[i] = A[n*i + j]*x[j] + y[i];\n'
        '   }\n'
        '}'
    )

    c = glsl_code(A[i, j]*x[j], assign_to=y[i])
    assert c == s


def test_dummy_loops():
    i, m = symbols('i m', integer=True, cls=Dummy)
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx(i, m)

    expected = (
        'for (int i_%(icount)i=0; i_%(icount)i<m_%(mcount)i; i_%(icount)i++){\n'
        '   y[i_%(icount)i] = x[i_%(icount)i];\n'
        '}'
    ) % {'icount': i.label.dummy_index, 'mcount': m.dummy_index}
    code = glsl_code(x[i], assign_to=y[i])
    assert code == expected


def test_glsl_code_loops_add():
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols
    n, m = symbols('n m', integer=True)
    A = IndexedBase('A')
    x = IndexedBase('x')
    y = IndexedBase('y')
    z = IndexedBase('z')
    i = Idx('i', m)
    j = Idx('j', n)

    s = (
        'for (int i=0; i<m; i++){\n'
        '   y[i] = x[i] + z[i];\n'
        '}\n'
        'for (int i=0; i<m; i++){\n'
        '   for (int j=0; j<n; j++){\n'
        '      y[i] = A[n*i + j]*x[j] + y[i];\n'
        '   }\n'
        '}'
    )
    c = glsl_code(A[i, j]*x[j] + x[i] + z[i], assign_to=y[i])
    assert c == s


def test_glsl_code_loops_multiple_contractions():
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols
    n, m, o, p = symbols('n m o p', integer=True)
    a = IndexedBase('a')
    b = IndexedBase('b')
    y = IndexedBase('y')
    i = Idx('i', m)
    j = Idx('j', n)
    k = Idx('k', o)
    l = Idx('l', p)

    s = (
        'for (int i=0; i<m; i++){\n'
        '   y[i] = 0.0;\n'
        '}\n'
        'for (int i=0; i<m; i++){\n'
        '   for (int j=0; j<n; j++){\n'
        '      for (int k=0; k<o; k++){\n'
        '         for (int l=0; l<p; l++){\n'
        '            y[i] = a[%s]*b[%s] + y[i];\n' % (i*n*o*p + j*o*p + k*p + l, j*o*p + k*p + l) +\
        '         }\n'
        '      }\n'
        '   }\n'
        '}'
    )
    c = glsl_code(b[j, k, l]*a[i, j, k, l], assign_to=y[i])
    assert c == s


def test_glsl_code_loops_addfactor():
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols
    n, m, o, p = symbols('n m o p', integer=True)
    a = IndexedBase('a')
    b = IndexedBase('b')
    c = IndexedBase('c')
    y = IndexedBase('y')
    i = Idx('i', m)
    j = Idx('j', n)
    k = Idx('k', o)
    l = Idx('l', p)

    s = (
        'for (int i=0; i<m; i++){\n'
        '   y[i] = 0.0;\n'
        '}\n'
        'for (int i=0; i<m; i++){\n'
        '   for (int j=0; j<n; j++){\n'
        '      for (int k=0; k<o; k++){\n'
        '         for (int l=0; l<p; l++){\n'
        '            y[i] = (a[%s] + b[%s])*c[%s] + y[i];\n' % (i*n*o*p + j*o*p + k*p + l, i*n*o*p + j*o*p + k*p + l, j*o*p + k*p + l) +\
        '         }\n'
        '      }\n'
        '   }\n'
        '}'
    )
    c = glsl_code((a[i, j, k, l] + b[i, j, k, l])*c[j, k, l], assign_to=y[i])
    assert c == s


def test_glsl_code_loops_multiple_terms():
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols
    n, m, o, p = symbols('n m o p', integer=True)
    a = IndexedBase('a')
    b = IndexedBase('b')
    c = IndexedBase('c')
    y = IndexedBase('y')
    i = Idx('i', m)
    j = Idx('j', n)
    k = Idx('k', o)

    s0 = (
        'for (int i=0; i<m; i++){\n'
        '   y[i] = 0.0;\n'
        '}\n'
    )
    s1 = (
        'for (int i=0; i<m; i++){\n'
        '   for (int j=0; j<n; j++){\n'
        '      for (int k=0; k<o; k++){\n'
        '         y[i] = b[j]*b[k]*c[%s] + y[i];\n' % (i*n*o + j*o + k) +\
        '      }\n'
        '   }\n'
        '}\n'
    )
    s2 = (
        'for (int i=0; i<m; i++){\n'
        '   for (int k=0; k<o; k++){\n'
        '      y[i] = a[%s]*b[k] + y[i];\n' % (i*o + k) +\
        '   }\n'
        '}\n'
    )
    s3 = (
        'for (int i=0; i<m; i++){\n'
        '   for (int j=0; j<n; j++){\n'
        '      y[i] = a[%s]*b[j] + y[i];\n' % (i*n + j) +\
        '   }\n'
        '}\n'
    )
    c = glsl_code(
        b[j]*a[i, j] + b[k]*a[i, k] + b[j]*b[k]*c[i, j, k], assign_to=y[i])
    assert (c == s0 + s1 + s2 + s3[:-1] or
            c == s0 + s1 + s3 + s2[:-1] or
            c == s0 + s2 + s1 + s3[:-1] or
            c == s0 + s2 + s3 + s1[:-1] or
            c == s0 + s3 + s1 + s2[:-1] or
            c == s0 + s3 + s2 + s1[:-1])


def test_Matrix_printing():
    # Test returning a Matrix
    mat = Matrix([x*y, Piecewise((2 + x, y>0), (y, True)), sin(z)])
    A = MatrixSymbol('A', 3, 1)
    assert glsl_code(mat, assign_to=A) == (
        "A[0] = x*y;\n"
        "if (y > 0) {\n"
        "   A[1] = x + 2;\n"
        "}\n"
        "else {\n"
        "   A[1] = y;\n"
        "}\n"
        "A[2] = sin(z);")
    # Test using MatrixElements in expressions
    expr = Piecewise((2*A[2, 0], x > 0), (A[2, 0], True)) + sin(A[1, 0]) + A[0, 0]
    assert glsl_code(expr) == (
        "((x > 0) ? (\n"
        "   2*A[2]\n"
        ")\n"
        ": (\n"
        "   A[2]\n"
        ")) + sin(A[1]) + A[0]")
    # Test using MatrixElements in a Matrix
    q = MatrixSymbol('q', 5, 1)
    M = MatrixSymbol('M', 3, 3)
    m = Matrix([[sin(q[1,0]), 0, cos(q[2,0])],
        [q[1,0] + q[2,0], q[3, 0], 5],
        [2*q[4, 0]/q[1,0], sqrt(q[0,0]) + 4, 0]])
    assert glsl_code(m,M) == (
        "M[0] = sin(q[1]);\n"
        "M[1] = 0;\n"
        "M[2] = cos(q[2]);\n"
        "M[3] = q[1] + q[2];\n"
        "M[4] = q[3];\n"
        "M[5] = 5;\n"
        "M[6] = 2*q[4]/q[1];\n"
        "M[7] = sqrt(q[0]) + 4;\n"
        "M[8] = 0;"
        )
