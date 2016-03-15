from sympy.core import (pi, oo, symbols, Rational, Integer,
                        GoldenRatio, EulerGamma, Catalan, Lambda, Dummy, Eq)
from sympy.functions import (Piecewise, sin, cos, Abs, exp, ceiling, sqrt,
                             gamma, sign)
from sympy.sets import Range
from sympy.logic import ITE
from sympy.codegen import For, aug_assign, Assignment
from sympy.utilities.pytest import raises
from sympy.printing.ccode import CCodePrinter
from sympy.utilities.lambdify import implemented_function
from sympy.tensor import IndexedBase, Idx
from sympy.matrices import Matrix, MatrixSymbol

from sympy import ccode

x, y, z = symbols('x,y,z')


def test_printmethod():
    class fabs(Abs):
        def _ccode(self, printer):
            return "fabs(%s)" % printer._print(self.args[0])
    assert ccode(fabs(x)) == "fabs(x)"


def test_ccode_sqrt():
    assert ccode(sqrt(x)) == "sqrt(x)"
    assert ccode(x**0.5) == "sqrt(x)"
    assert ccode(sqrt(x)) == "sqrt(x)"


def test_ccode_Pow():
    assert ccode(x**3) == "pow(x, 3)"
    assert ccode(x**(y**3)) == "pow(x, pow(y, 3))"
    g = implemented_function('g', Lambda(x, 2*x))
    assert ccode(1/(g(x)*3.5)**(x - y**x)/(x**2 + y)) == \
        "pow(3.5*2*x, -x + pow(y, x))/(pow(x, 2) + y)"
    assert ccode(x**-1.0) == '1.0/x'
    assert ccode(x**Rational(2, 3)) == 'pow(x, 2.0L/3.0L)'
    _cond_cfunc = [(lambda base, exp: exp.is_integer, "dpowi"),
                   (lambda base, exp: not exp.is_integer, "pow")]
    assert ccode(x**3, user_functions={'Pow': _cond_cfunc}) == 'dpowi(x, 3)'
    assert ccode(x**3.2, user_functions={'Pow': _cond_cfunc}) == 'pow(x, 3.2)'


def test_ccode_constants_mathh():
    assert ccode(exp(1)) == "M_E"
    assert ccode(pi) == "M_PI"
    assert ccode(oo) == "HUGE_VAL"
    assert ccode(-oo) == "-HUGE_VAL"


def test_ccode_constants_other():
    assert ccode(2*GoldenRatio) == "double const GoldenRatio = 1.61803398874989;\n2*GoldenRatio"
    assert ccode(
        2*Catalan) == "double const Catalan = 0.915965594177219;\n2*Catalan"
    assert ccode(2*EulerGamma) == "double const EulerGamma = 0.577215664901533;\n2*EulerGamma"


def test_ccode_Rational():
    assert ccode(Rational(3, 7)) == "3.0L/7.0L"
    assert ccode(Rational(18, 9)) == "2"
    assert ccode(Rational(3, -7)) == "-3.0L/7.0L"
    assert ccode(Rational(-3, -7)) == "3.0L/7.0L"
    assert ccode(x + Rational(3, 7)) == "x + 3.0L/7.0L"
    assert ccode(Rational(3, 7)*x) == "(3.0L/7.0L)*x"


def test_ccode_Integer():
    assert ccode(Integer(67)) == "67"
    assert ccode(Integer(-1)) == "-1"


def test_ccode_functions():
    assert ccode(sin(x) ** cos(x)) == "pow(sin(x), cos(x))"


def test_ccode_inline_function():
    x = symbols('x')
    g = implemented_function('g', Lambda(x, 2*x))
    assert ccode(g(x)) == "2*x"
    g = implemented_function('g', Lambda(x, 2*x/Catalan))
    assert ccode(
        g(x)) == "double const Catalan = %s;\n2*x/Catalan" % Catalan.n()
    A = IndexedBase('A')
    i = Idx('i', symbols('n', integer=True))
    g = implemented_function('g', Lambda(x, x*(1 + x)*(2 + x)))
    assert ccode(g(A[i]), assign_to=A[i]) == (
        "for (int i=0; i<n; i++){\n"
        "   A[i] = (A[i] + 1)*(A[i] + 2)*A[i];\n"
        "}"
    )


def test_ccode_exceptions():
    assert ccode(ceiling(x)) == "ceil(x)"
    assert ccode(Abs(x)) == "fabs(x)"
    assert ccode(gamma(x)) == "tgamma(x)"


def test_ccode_user_functions():
    x = symbols('x', integer=False)
    n = symbols('n', integer=True)
    custom_functions = {
        "ceiling": "ceil",
        "Abs": [(lambda x: not x.is_integer, "fabs"), (lambda x: x.is_integer, "abs")],
    }
    assert ccode(ceiling(x), user_functions=custom_functions) == "ceil(x)"
    assert ccode(Abs(x), user_functions=custom_functions) == "fabs(x)"
    assert ccode(Abs(n), user_functions=custom_functions) == "abs(n)"


def test_ccode_boolean():
    assert ccode(x & y) == "x && y"
    assert ccode(x | y) == "x || y"
    assert ccode(~x) == "!x"
    assert ccode(x & y & z) == "x && y && z"
    assert ccode(x | y | z) == "x || y || z"
    assert ccode((x & y) | z) == "z || x && y"
    assert ccode((x | y) & z) == "z && (x || y)"


def test_ccode_Piecewise():
    expr = Piecewise((x, x < 1), (x**2, True))
    assert ccode(expr) == (
            "((x < 1) ? (\n"
            "   x\n"
            ")\n"
            ": (\n"
            "   pow(x, 2)\n"
            "))")
    assert ccode(expr, assign_to="c") == (
            "if (x < 1) {\n"
            "   c = x;\n"
            "}\n"
            "else {\n"
            "   c = pow(x, 2);\n"
            "}")
    expr = Piecewise((x, x < 1), (x + 1, x < 2), (x**2, True))
    assert ccode(expr) == (
            "((x < 1) ? (\n"
            "   x\n"
            ")\n"
            ": ((x < 2) ? (\n"
            "   x + 1\n"
            ")\n"
            ": (\n"
            "   pow(x, 2)\n"
            ")))")
    assert ccode(expr, assign_to='c') == (
            "if (x < 1) {\n"
            "   c = x;\n"
            "}\n"
            "else if (x < 2) {\n"
            "   c = x + 1;\n"
            "}\n"
            "else {\n"
            "   c = pow(x, 2);\n"
            "}")
    # Check that Piecewise without a True (default) condition error
    expr = Piecewise((x, x < 1), (x**2, x > 1), (sin(x), x > 0))
    raises(ValueError, lambda: ccode(expr))


def test_ccode_Piecewise_deep():
    p = ccode(2*Piecewise((x, x < 1), (x + 1, x < 2), (x**2, True)))
    assert p == (
            "2*((x < 1) ? (\n"
            "   x\n"
            ")\n"
            ": ((x < 2) ? (\n"
            "   x + 1\n"
            ")\n"
            ": (\n"
            "   pow(x, 2)\n"
            ")))")
    expr = x*y*z + x**2 + y**2 + Piecewise((0, x < 0.5), (1, True)) + cos(z) - 1
    assert ccode(expr) == (
            "pow(x, 2) + x*y*z + pow(y, 2) + ((x < 0.5) ? (\n"
            "   0\n"
            ")\n"
            ": (\n"
            "   1\n"
            ")) + cos(z) - 1")
    assert ccode(expr, assign_to='c') == (
            "c = pow(x, 2) + x*y*z + pow(y, 2) + ((x < 0.5) ? (\n"
            "   0\n"
            ")\n"
            ": (\n"
            "   1\n"
            ")) + cos(z) - 1;")


def test_ccode_ITE():
    expr = ITE(x < 1, x, x**2)
    assert ccode(expr) == (
            "((x < 1) ? (\n"
            "   x\n"
            ")\n"
            ": (\n"
            "   pow(x, 2)\n"
            "))")


def test_ccode_settings():
    raises(TypeError, lambda: ccode(sin(x), method="garbage"))


def test_ccode_Indexed():
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols
    n, m, o = symbols('n m o', integer=True)
    i, j, k = Idx('i', n), Idx('j', m), Idx('k', o)
    p = CCodePrinter()
    p._not_c = set()

    x = IndexedBase('x')[j]
    assert p._print_Indexed(x) == 'x[j]'
    A = IndexedBase('A')[i, j]
    assert p._print_Indexed(A) == 'A[%s]' % (m*i+j)
    B = IndexedBase('B')[i, j, k]
    assert p._print_Indexed(B) == 'B[%s]' % (i*o*m+j*o+k)

    assert p._not_c == set()


def test_ccode_Indexed_without_looking_for_contraction():
    len_y = 5
    y = IndexedBase('y', shape=(len_y,))
    x = IndexedBase('x', shape=(len_y,))
    Dy = IndexedBase('Dy', shape=(len_y-1,))
    i = Idx('i', len_y-1)
    e=Eq(Dy[i], (y[i+1]-y[i])/(x[i+1]-x[i]))
    code0 = ccode(e.rhs, assign_to=e.lhs, contract=False)
    assert code0 == 'Dy[i] = (y[%s] - y[i])/(x[%s] - x[i]);' % (i + 1, i + 1)


def test_ccode_loops_matrix_vector():
    n, m = symbols('n m', integer=True)
    A = IndexedBase('A')
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx('i', m)
    j = Idx('j', n)

    s = (
        'for (int i=0; i<m; i++){\n'
        '   y[i] = 0;\n'
        '}\n'
        'for (int i=0; i<m; i++){\n'
        '   for (int j=0; j<n; j++){\n'
        '      y[i] = x[j]*A[%s] + y[i];\n' % (i*n + j) +\
        '   }\n'
        '}'
    )
    c = ccode(A[i, j]*x[j], assign_to=y[i])
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
    code = ccode(x[i], assign_to=y[i])
    assert code == expected


def test_ccode_loops_add():
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
        '      y[i] = x[j]*A[%s] + y[i];\n' % (i*n + j) +\
        '   }\n'
        '}'
    )
    c = ccode(A[i, j]*x[j] + x[i] + z[i], assign_to=y[i])
    assert c == s


def test_ccode_loops_multiple_contractions():
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
        '   y[i] = 0;\n'
        '}\n'
        'for (int i=0; i<m; i++){\n'
        '   for (int j=0; j<n; j++){\n'
        '      for (int k=0; k<o; k++){\n'
        '         for (int l=0; l<p; l++){\n'
        '            y[i] = y[i] + b[%s]*a[%s];\n' % (j*o*p + k*p + l, i*n*o*p + j*o*p + k*p + l) +\
        '         }\n'
        '      }\n'
        '   }\n'
        '}'
    )
    c = ccode(b[j, k, l]*a[i, j, k, l], assign_to=y[i])
    assert c == s


def test_ccode_loops_addfactor():
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
        '   y[i] = 0;\n'
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
    c = ccode((a[i, j, k, l] + b[i, j, k, l])*c[j, k, l], assign_to=y[i])
    assert c == s


def test_ccode_loops_multiple_terms():
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
        '   y[i] = 0;\n'
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
        '      y[i] = b[k]*a[%s] + y[i];\n' % (i*o + k) +\
        '   }\n'
        '}\n'
    )
    s3 = (
        'for (int i=0; i<m; i++){\n'
        '   for (int j=0; j<n; j++){\n'
        '      y[i] = b[j]*a[%s] + y[i];\n' % (i*n + j) +\
        '   }\n'
        '}\n'
    )
    c = ccode(
        b[j]*a[i, j] + b[k]*a[i, k] + b[j]*b[k]*c[i, j, k], assign_to=y[i])
    assert (c == s0 + s1 + s2 + s3[:-1] or
            c == s0 + s1 + s3 + s2[:-1] or
            c == s0 + s2 + s1 + s3[:-1] or
            c == s0 + s2 + s3 + s1[:-1] or
            c == s0 + s3 + s1 + s2[:-1] or
            c == s0 + s3 + s2 + s1[:-1])


def test_dereference_printing():
    expr = x + y + sin(z) + z
    assert ccode(expr, dereference=[z]) == "x + y + (*z) + sin((*z))"


def test_Matrix_printing():
    # Test returning a Matrix
    mat = Matrix([x*y, Piecewise((2 + x, y>0), (y, True)), sin(z)])
    A = MatrixSymbol('A', 3, 1)
    assert ccode(mat, A) == (
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
    assert ccode(expr) == (
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
    assert ccode(m, M) == (
        "M[0] = sin(q[1]);\n"
        "M[1] = 0;\n"
        "M[2] = cos(q[2]);\n"
        "M[3] = q[1] + q[2];\n"
        "M[4] = q[3];\n"
        "M[5] = 5;\n"
        "M[6] = 2*q[4]*1.0/q[1];\n"
        "M[7] = 4 + sqrt(q[0]);\n"
        "M[8] = 0;")


def test_ccode_reserved_words():

    x, y = symbols('x, if')

    assert ccode(y**2) == 'pow(if_, 2)'
    assert ccode(x * y**2, dereference=[y]) == 'pow((*if_), 2)*x'

    expected = 'pow(if_unreserved, 2)'
    assert ccode(y**2, reserved_word_suffix='_unreserved') == expected

    with raises(ValueError):
        ccode(y**2, error_on_reserved=True)


def test_ccode_sign():

    expr = sign(x) * y
    assert ccode(expr) == 'y*(((x) > 0) - ((x) < 0))'
    assert ccode(expr, 'z') == 'z = y*(((x) > 0) - ((x) < 0));'

    assert ccode(sign(2 * x + x**2) * x + x**2) == \
        'pow(x, 2) + x*(((pow(x, 2) + 2*x) > 0) - ((pow(x, 2) + 2*x) < 0))'

    expr = sign(cos(x))
    assert ccode(expr) == '(((cos(x)) > 0) - ((cos(x)) < 0))'

def test_ccode_Assignment():
    assert ccode(Assignment(x, y + z)) == 'x = y + z;'
    assert ccode(aug_assign(x, '+', y + z)) == 'x += y + z;'


def test_ccode_For():
    f = For(x, Range(0, 10, 2), [aug_assign(y, '*', x)])
    sol = ccode(f)
    assert sol == ("for (x = 0; x < 10; x += 2) {\n"
                   "   y *= x;\n"
                   "}")
