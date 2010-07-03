from sympy.core import pi, oo, symbols, Function, Rational, Integer, GoldenRatio, EulerGamma, Catalan, Lambda, Dummy
from sympy.functions import Piecewise, sin, cos, Abs, exp, ceiling, sqrt
from sympy.utilities.pytest import XFAIL, raises
from sympy.printing.ccode import CCodePrinter
from sympy.utilities.lambdify import implemented_function
from sympy.tensor import IndexedBase, Idx

# import test
from sympy import ccode

x, y, z = symbols('x,y,z')
g = Function('g')

def test_printmethod():
    class fabs(Abs):
        def _ccode(self, printer):
            return "fabs(%s)" % printer._print(self.args[0])
    assert ccode(fabs(x)) == "fabs(x)"

def test_ccode_sqrt():
    assert ccode(sqrt(x)) == "sqrt(x)"
    assert ccode(x**0.5) == "sqrt(x)"
    assert ccode(x**Rational(1,2)) == "sqrt(x)"

def test_ccode_Pow():
    assert ccode(x**3) == "pow(x, 3)"
    assert ccode(x**(y**3)) == "pow(x, pow(y, 3))"
    assert ccode(1/(g(x)*3.5)**(x - y**x)/(x**2 + y)) == \
        "pow(3.5*g(x), -x + pow(y, x))/(y + pow(x, 2))"

def test_ccode_constants_mathh():
    assert ccode(exp(1)) == "M_E"
    assert ccode(pi) == "M_PI"
    assert ccode(oo) == "HUGE_VAL"
    assert ccode(-oo) == "-HUGE_VAL"

def test_ccode_constants_other():
    assert ccode(2*GoldenRatio) == "double const GoldenRatio = 1.61803398874989;\n2*GoldenRatio"
    assert ccode(2*Catalan) == "double const Catalan = 0.915965594177219;\n2*Catalan"
    assert ccode(2*EulerGamma) == "double const EulerGamma = 0.577215664901533;\n2*EulerGamma"

def test_ccode_Rational():
    assert ccode(Rational(3,7)) == "3.0/7.0"
    assert ccode(Rational(18,9)) == "2"
    assert ccode(Rational(3,-7)) == "-3.0/7.0"
    assert ccode(Rational(-3,-7)) == "3.0/7.0"

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
    assert ccode(g(x)) == "double const Catalan = %s;\n2*x/Catalan" %Catalan.n()
    A = IndexedBase('A')
    i = Idx('i', symbols('n', integer=True))
    g = implemented_function('g', Lambda(x, x*(1 + x)*(2 + x)))
    assert ccode(g(A[i]), assign_to=A[i]) == (
            "for (int i=0; i<n; i++){\n"
            "   A[i] = (1 + A[i])*(2 + A[i])*A[i];\n"
            "}"
            )

def test_ccode_exceptions():
    assert ccode(ceiling(x)) == "ceil(x)"
    assert ccode(Abs(x)) == "fabs(x)"

def test_ccode_boolean():
    assert ccode(x&y) == "x&&y"
    assert ccode(x|y) == "x||y"
    assert ccode(~x) == "!x"
    assert ccode(x&y&z) == "x&&y&&z"
    assert ccode(x|y|z) == "x||y||z"
    assert ccode((x&y)|z) == "x&&y||z"
    assert ccode((x|y)&z) == "(x||y)&&z"

def test_ccode_Piecewise():
    p = ccode(Piecewise((x,x<1),(x**2,True)))
    s = \
"""\
if (x < 1) {
   x
}
else {
   pow(x, 2)
}\
"""
    assert p == s

def test_ccode_Piecewise_deep():
    p = ccode(2*Piecewise((x,x<1),(x**2,True)))
    s = \
"""\
2*if (x < 1) {
   x
}
else {
   pow(x, 2)
}\
"""
    assert p == s

def test_ccode_settings():
    raises(TypeError, 'ccode(sin(x),method="garbage")')

def test_ccode_Indexed():
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols
    i,j,k,n,m,o = symbols('i j k n m o', integer=True)

    p = CCodePrinter()
    p._not_c = set()

    x = IndexedBase('x')[Idx(j, n)]
    assert p._print_Indexed(x) == 'x[j]'
    A = IndexedBase('A')[Idx(i, m), Idx(j, n)]
    assert p._print_Indexed(A) == 'A[%s]'% str(j + n*i)
    B = IndexedBase('B')[Idx(i, m), Idx(j, n), Idx(k, o)]
    assert p._print_Indexed(B) == 'B[%s]'% str(k + i*n*o + j*o)

    assert p._not_c == set()


def test_ccode_loops_matrix_vector():
    n,m = symbols('n m', integer=True)
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
            '      y[i] = A[j + i*n]*x[j] + y[i];\n'
            '   }\n'
            '}'
            )
    c = ccode(A[i, j]*x[j], assign_to=y[i])
    assert c == s

def test_dummy_loops():
    # the following line could also be
    # [Dummy(s, integer=True) for s in 'im']
    # or [Dummy(integer=True) for s in 'im']
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
            '      y[i] = A[j + i*n]*x[j] + y[i];\n'
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
'            y[i] = b[l + k*p + j*o*p]*a[l + k*p + j*o*p + i*n*o*p] + y[i];\n'
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
'            y[i] = (a[l + k*p + j*o*p + i*n*o*p] + b[l + k*p + j*o*p + i*n*o*p])*c[l + k*p + j*o*p] + y[i];\n'
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
'         y[i] = b[j]*b[k]*c[k + j*o + i*n*o] + y[i];\n'
'      }\n'
'   }\n'
'}\n'
    )
    s2 = (
'for (int i=0; i<m; i++){\n'
'   for (int k=0; k<o; k++){\n'
'      y[i] = b[k]*a[k + i*o] + y[i];\n'
'   }\n'
'}\n'
    )
    s3 = (
'for (int i=0; i<m; i++){\n'
'   for (int j=0; j<n; j++){\n'
'      y[i] = b[j]*a[j + i*n] + y[i];\n'
'   }\n'
'}\n'
            )
    c = ccode(b[j]*a[i, j] + b[k]*a[i, k] + b[j]*b[k]*c[i, j, k], assign_to=y[i])
    assert (c == s0 + s1 + s2 + s3[:-1] or
            c == s0 + s1 + s3 + s2[:-1] or
            c == s0 + s2 + s1 + s3[:-1] or
            c == s0 + s2 + s3 + s1[:-1] or
            c == s0 + s3 + s1 + s2[:-1] or
            c == s0 + s3 + s2 + s1[:-1])
