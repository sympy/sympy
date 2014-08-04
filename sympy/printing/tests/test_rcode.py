from sympy.core import pi, oo, symbols, Function, Rational, Integer, GoldenRatio, EulerGamma, Catalan, Lambda, Dummy, Eq
from sympy.functions import Piecewise, sin, cos, Abs, exp, ceiling, sqrt, gamma
from sympy.utilities.pytest import raises
from sympy.printing.rcode import RCodePrinter
from sympy.utilities.lambdify import implemented_function
from sympy.tensor import IndexedBase, Idx
import difflib #import Differ 
from pprint import pprint

# import test
from sympy import rcode

x, y, z = symbols('x,y,z')
g = Function('g')


def test_printmethod():
    p=rcode(Abs(x))
    assert p  == "abs(x)"


def test_rcode_sqrt():
    assert rcode(sqrt(x)) == "sqrt(x)"
    assert rcode(x**0.5) == "sqrt(x)"
    assert rcode(sqrt(x)) == "sqrt(x)"

def test_rcode_Pow():
    assert rcode(x**3) == "x^3"
    assert rcode(x**(y**3)) == "x^(y^3)"
    assert rcode(1/(g(x)*3.5)**(x - y**x)/(x**2 + y)) == \
        "(3.5*g(x))^(-x + y^x)/(x^2 + y)"
    assert rcode(x**-1.0) == '1.0/x'
    assert rcode(x**Rational(2, 3)) == 'x^(2.0/3.0)'
    _cond_cfunc = [(lambda base, exp: exp.is_integer, "dpowi"),
                   (lambda base, exp: not exp.is_integer, "pow")]
    assert rcode(x**3, user_functions={'Pow': _cond_cfunc}) == 'dpowi(x, 3)'
    assert rcode(x**3.2, user_functions={'Pow': _cond_cfunc}) == 'pow(x, 3.2)'

def test_rcode_constants_mathh():
    p=rcode(exp(1)) 
    assert p == "exp(1)" #nothing changes
    assert rcode(pi) == "pi" #nothing changes
    assert rcode(oo) == "Inf"
    assert rcode(-oo) == "-Inf"

def test_rcode_constants_other():
    assert rcode(2*GoldenRatio) == "GoldenRatio = 1.61803398874989;\n2*GoldenRatio"
    assert rcode(
        2*Catalan) == "Catalan = 0.915965594177219;\n2*Catalan"
    assert rcode(2*EulerGamma) == "EulerGamma = 0.577215664901533;\n2*EulerGamma"


def test_rcode_Rational():
    assert rcode(Rational(3, 7)) == "3.0/7.0"
    assert rcode(Rational(18, 9)) == "2"
    assert rcode(Rational(3, -7)) == "-3.0/7.0"
    assert rcode(Rational(-3, -7)) == "3.0/7.0"
    assert rcode(x + Rational(3, 7)) == "x + 3.0/7.0"
    assert rcode(Rational(3, 7)*x) == "(3.0/7.0)*x"


def test_rcode_Integer():
    assert rcode(Integer(67)) == "67"
    assert rcode(Integer(-1)) == "-1"


def test_rcode_functions():
    assert rcode(sin(x) ** cos(x)) == "sin(x)^cos(x)"


def test_rcode_inline_function():
    x = symbols('x')
    g = implemented_function('g', Lambda(x, 2*x))
    assert rcode(g(x)) == "2*x"
    g = implemented_function('g', Lambda(x, 2*x/Catalan))
    assert rcode( g(x)) == "Catalan = %s;\n2*x/Catalan" % Catalan.n()
    A = IndexedBase('A')
    i = Idx('i', symbols('n', integer=True))
    g = implemented_function('g', Lambda(x, x*(1 + x)*(2 + x)))
    res=rcode(g(A[i]), assign_to=A[i])
    ref=(
        "for (i in 1:n){\n"
        "   A[i] = (A[i] + 1)*(A[i] + 2)*A[i];\n"
        "}"
    ) 
    #d=difflib.Differ()
    #pprint(list(d.compare(res.splitlines(keepends=True),ref.splitlines(keepends=True))))
    assert res == ref


def test_rcode_exceptions():
    assert rcode(Abs(x)) == "abs(x)"


def test_rcode_user_functions():
    x = symbols('x', integer=False)
    n = symbols('n', integer=True)
    custom_functions = {
        "ceiling": "my_ceil",
        "Abs": [(lambda x: not x.is_integer, "my_abs"), (lambda x: x.is_integer, "abs")],
    }
    assert rcode(ceiling(x), user_functions=custom_functions) == "my_ceil(x)"
    assert rcode(Abs(x), user_functions=custom_functions) == "my_abs(x)"
    assert rcode(Abs(n), user_functions=custom_functions) == "abs(n)"


def test_rcode_boolean():
    assert rcode(x & y) == "x && y"
    assert rcode(x | y) == "x || y"
    assert rcode(~x) == "!x"
    assert rcode(x & y & z) == "x && y && z"
    assert rcode(x | y | z) == "x || y || z"
    assert rcode((x & y) | z) == "z || x && y"
    assert rcode((x | y) & z) == "z && (x || y)"


def test_rcode_Piecewise():
    p = rcode(Piecewise((x, x < 1), (x**2, True)))
    s = \
"""\
if (x < 1) {
   x
}
else {
   x^2
}\
"""
    assert p == s

def test_rcode_Piecewise_deep():
    p = rcode(2*Piecewise((x, x < 1), (x**2, True)))
    s = \
"""\
2*ifelse(x < 1,x,x^2)\
"""
    assert p == s
    
    p = rcode(2*Piecewise((x, x < 1), (x**2, x<2), (x**3,True)))
    s = \
"""\
2*ifelse(x < 1,x,ifelse(x < 2,x^2,x^3))\
"""
    assert p == s
    # last condition missing
    p = rcode(2*Piecewise((x, x < 1)))
    s = \
"""\
2*ifelse(x < 1,x,NA)\
"""
    assert p == s
    
    p = rcode(2*Piecewise((x, x < 1), (x**2, x<2)))
    s = \
"""\
2*ifelse(x < 1,x,ifelse(x < 2,x^2,NA))\
"""
    assert p == s

def test_rcode_settings():
    raises(TypeError, lambda: rcode(sin(x), method="garbage"))


def test_rcode_Indexed():
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols
    n, m, o = symbols('n m o', integer=True)
    i, j, k = Idx('i', n), Idx('j', m), Idx('k', o)
    p = RCodePrinter()
    p._not_r = set()

    x = IndexedBase('x')[j]
    assert p._print_Indexed(x) == 'x[j]'
    A = IndexedBase('A')[i, j]
    assert p._print_Indexed(A) == 'A[i, j]' 
    B = IndexedBase('B')[i, j, k]
    assert p._print_Indexed(B) == 'B[i, j, k]' 

    assert p._not_r == set()

def test_rcode_Indexed_without_looking_for_contraction():
    len_y = 5
    y = IndexedBase('y', shape=(len_y,))
    x = IndexedBase('x', shape=(len_y,))
    Dy = IndexedBase('Dy', shape=(len_y-1,))
    i = Idx('i', len_y-1)
    e=Eq(Dy[i], (y[i+1]-y[i])/(x[i+1]-x[i]))
    code0 = rcode(e.rhs, assign_to=e.lhs, contract=False)
    assert code0 == 'Dy[i] = (y[%s] - y[i])/(x[%s] - x[i]);' % (i + 1, i + 1)


def test_rcode_loops_matrix_vector():
    n, m = symbols('n m', integer=True)
    A = IndexedBase('A')
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx('i', m)
    j = Idx('j', n)

    s = (
        'for (i in 1:m){\n'
        '   y[i] = 0;\n'
        '}\n'
        'for (i in 1:m){\n'
        '   for (j in 1:n){\n'
        '      y[i] = x[j]*A[i, j] + y[i];\n' 
        '   }\n'
        '}'
    )
    c = rcode(A[i, j]*x[j], assign_to=y[i])
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
            'for (i_%(icount)i in 1:m_%(mcount)i){\n'
        '   y[i_%(icount)i] = x[i_%(icount)i];\n'
        '}'
    ) % {'icount': i.label.dummy_index, 'mcount': m.dummy_index}
    code = rcode(x[i], assign_to=y[i])
    assert code == expected


def test_rcode_loops_add():
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
        'for (i in 1:m){\n'
        '   y[i] = x[i] + z[i];\n'
        '}\n'
        'for (i in 1:m){\n'
        '   for (j in 1:n){\n'
        '      y[i] = x[j]*A[i, j] + y[i];\n' 
        '   }\n'
        '}'
    )
    c = rcode(A[i, j]*x[j] + x[i] + z[i], assign_to=y[i])
    assert c == s

#
#mm
#

def test_rcode_loops_multiple_contractions():
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
        'for (i in 1:m){\n'
        '   y[i] = 0;\n'
        '}\n'
        'for (i in 1:m){\n'
        '   for (j in 1:n){\n'
        '      for (k in 1:o){\n'
        '         for (l in 1:p){\n'
        '            y[i] = y[i] + b[j, k, l]*a[i, j, k, l];\n' 
        '         }\n'
        '      }\n'
        '   }\n'
        '}'
    )
    c = rcode(b[j, k, l]*a[i, j, k, l], assign_to=y[i])
    assert c == s


def test_rcode_loops_addfactor():
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
        'for (i in 1:m){\n'
        '   y[i] = 0;\n'
        '}\n'
        'for (i in 1:m){\n'
        '   for (j in 1:n){\n'
        '      for (k in 1:o){\n'
        '         for (l in 1:p){\n'
        '            y[i] = (a[i, j, k, l] + b[i, j, k, l])*c[j, k, l] + y[i];\n' 
        '         }\n'
        '      }\n'
        '   }\n'
        '}'
    )
    c = rcode((a[i, j, k, l] + b[i, j, k, l])*c[j, k, l], assign_to=y[i])
    assert c == s


def test_rcode_loops_multiple_terms():
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
        'for (i in 1:m){\n'
        '   y[i] = 0;\n'
        '}\n'
    )
    s1 = (
        'for (i in 1:m){\n'
        '   for (j in 1:n){\n'
        '      for (k in 1:o){\n'
        '         y[i] = b[j]*b[k]*c[i, j, k] + y[i];\n' 
        '      }\n'
        '   }\n'
        '}\n'
    )
    s2 = (
        'for (i in 1:m){\n'
        '   for (k in 1:o){\n'
        '      y[i] = b[k]*a[i, k] + y[i];\n' 
        '   }\n'
        '}\n'
    )
    s3 = (
        'for (i in 1:m){\n'
        '   for (j in 1:n){\n'
        '      y[i] = b[j]*a[i, j] + y[i];\n' 
        '   }\n'
        '}\n'
    )
    c = rcode(
        b[j]*a[i, j] + b[k]*a[i, k] + b[j]*b[k]*c[i, j, k], assign_to=y[i])
    assert (c == s0 + s1 + s2 + s3[:-1] or
            c == s0 + s1 + s3 + s2[:-1] or
            c == s0 + s2 + s1 + s3[:-1] or
            c == s0 + s2 + s3 + s1[:-1] or
            c == s0 + s3 + s1 + s2[:-1] or
            c == s0 + s3 + s2 + s1[:-1])
