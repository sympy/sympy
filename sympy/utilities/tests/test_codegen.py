from sympy.core import symbols, Eq, pi, Catalan, Lambda, Dummy
from sympy.core.compatibility import StringIO
from sympy import erf, Integral
from sympy.utilities.codegen import (CCodeGen, Routine, InputArgument,
    CodeGenError, FCodeGen, codegen, CodeGenArgumentListError, OutputArgument,
    InOutArgument)
from sympy.utilities.pytest import raises
from sympy.utilities.lambdify import implemented_function

# import test:
#FIXME: Fails due to circular import in with core
# from sympy import codegen


def get_string(dump_fn, routines, prefix="file", header=False, empty=False):
    """Wrapper for dump_fn. dump_fn writes its results to a stream object and
       this wrapper returns the contents of that stream as a string. This
       auxiliary function is used by many tests below.

       The header and the empty lines are not generator to facilitate the
       testing of the output.
    """
    output = StringIO()
    dump_fn(routines, output, prefix, header, empty)
    source = output.getvalue()
    output.close()
    return source


def test_Routine_argument_order():
    a, x, y, z = symbols('a x y z')
    expr = (x + y)*z
    raises(CodeGenArgumentListError, lambda: Routine("test", expr,
           argument_sequence=[z, x]))
    raises(CodeGenArgumentListError, lambda: Routine("test", Eq(a,
           expr), argument_sequence=[z, x, y]))
    r = Routine('test', Eq(a, expr), argument_sequence=[z, x, a, y])
    assert [ arg.name for arg in r.arguments ] == [z, x, a, y]
    assert [ type(arg) for arg in r.arguments ] == [
        InputArgument, InputArgument, OutputArgument, InputArgument  ]
    r = Routine('test', Eq(z, expr), argument_sequence=[z, x, y])
    assert [ type(arg) for arg in r.arguments ] == [
        InOutArgument, InputArgument, InputArgument ]

    from sympy.tensor import IndexedBase, Idx
    A, B = map(IndexedBase, ['A', 'B'])
    m = symbols('m', integer=True)
    i = Idx('i', m)
    r = Routine('test', Eq(A[i], B[i]), argument_sequence=[B, A, m])
    assert [ arg.name for arg in r.arguments ] == [B.label, A.label, m]

    expr = Integral(x*y*z, (x, 1, 2), (y, 1, 3))
    r = Routine('test', Eq(a, expr), argument_sequence=[z, x, a, y])
    assert [ arg.name for arg in r.arguments ] == [z, x, a, y]


def test_empty_c_code():
    code_gen = CCodeGen()
    source = get_string(code_gen.dump_c, [])
    assert source == "#include \"file.h\"\n#include <math.h>\n"


def test_empty_c_code_with_comment():
    code_gen = CCodeGen()
    source = get_string(code_gen.dump_c, [], header=True)
    assert source[:82] == (
        "/******************************************************************************\n *"
    )
          #   "                    Code generated with sympy 0.7.2-git                    "
    assert source[158:] == (                                                              "*\n"
            " *                                                                            *\n"
            " *              See http://www.sympy.org/ for more information.               *\n"
            " *                                                                            *\n"
            " *                       This file is part of 'project'                       *\n"
            " ******************************************************************************/\n"
            "#include \"file.h\"\n"
            "#include <math.h>\n"
            )


def test_empty_c_header():
    code_gen = CCodeGen()
    source = get_string(code_gen.dump_h, [])
    assert source == "#ifndef PROJECT__FILE__H\n#define PROJECT__FILE__H\n#endif\n"


def test_simple_c_code():
    x, y, z = symbols('x,y,z')
    expr = (x + y)*z
    routine = Routine("test", expr)
    code_gen = CCodeGen()
    source = get_string(code_gen.dump_c, [routine])
    expected = (
        "#include \"file.h\"\n"
        "#include <math.h>\n"
        "double test(double x, double y, double z) {\n"
        "   return z*(x + y);\n"
        "}\n"
    )
    assert source == expected


def test_numbersymbol_c_code():
    routine = Routine("test", pi**Catalan)
    code_gen = CCodeGen()
    source = get_string(code_gen.dump_c, [routine])
    expected = (
        "#include \"file.h\"\n"
        "#include <math.h>\n"
        "double test() {\n"
        "   double const Catalan = 0.915965594177219;\n"
        "   return pow(M_PI, Catalan);\n"
        "}\n"
    )
    assert source == expected


def test_c_code_argument_order():
    x, y, z = symbols('x,y,z')
    expr = x + y
    routine = Routine("test", expr, argument_sequence=[z, x, y])
    code_gen = CCodeGen()
    source = get_string(code_gen.dump_c, [routine])
    expected = (
        "#include \"file.h\"\n"
        "#include <math.h>\n"
        "double test(double z, double x, double y) {\n"
        "   return x + y;\n"
        "}\n"
    )
    assert source == expected


def test_simple_c_header():
    x, y, z = symbols('x,y,z')
    expr = (x + y)*z
    routine = Routine("test", expr)
    code_gen = CCodeGen()
    source = get_string(code_gen.dump_h, [routine])
    expected = (
        "#ifndef PROJECT__FILE__H\n"
        "#define PROJECT__FILE__H\n"
        "double test(double x, double y, double z);\n"
        "#endif\n"
    )
    assert source == expected


def test_simple_c_codegen():
    x, y, z = symbols('x,y,z')
    expr = (x + y)*z
    result = codegen(("test", expr), "C", "file", header=False, empty=False)
    expected = [
        ("file.c",
        "#include \"file.h\"\n"
        "#include <math.h>\n"
        "double test(double x, double y, double z) {\n"
        "   return z*(x + y);\n"
        "}\n"),
        ("file.h",
        "#ifndef PROJECT__FILE__H\n"
        "#define PROJECT__FILE__H\n"
        "double test(double x, double y, double z);\n"
        "#endif\n")
    ]
    assert result == expected


def test_multiple_results_c():
    x, y, z = symbols('x,y,z')
    expr1 = (x + y)*z
    expr2 = (x - y)*z
    routine = Routine(
        "test",
        [expr1, expr2]
    )
    code_gen = CCodeGen()
    raises(CodeGenError, lambda: get_string(code_gen.dump_h, [routine]))


def test_no_results_c():
    raises(ValueError, lambda: Routine("test", []))


def test_ansi_math1_codegen():
    # not included: log10
    from sympy import (acos, asin, atan, ceiling, cos, cosh, floor, log, ln,
        sin, sinh, sqrt, tan, tanh, Abs)
    x = symbols('x')
    name_expr = [
        ("test_fabs", Abs(x)),
        ("test_acos", acos(x)),
        ("test_asin", asin(x)),
        ("test_atan", atan(x)),
        ("test_ceil", ceiling(x)),
        ("test_cos", cos(x)),
        ("test_cosh", cosh(x)),
        ("test_floor", floor(x)),
        ("test_log", log(x)),
        ("test_ln", ln(x)),
        ("test_sin", sin(x)),
        ("test_sinh", sinh(x)),
        ("test_sqrt", sqrt(x)),
        ("test_tan", tan(x)),
        ("test_tanh", tanh(x)),
    ]
    result = codegen(name_expr, "C", "file", header=False, empty=False)
    assert result[0][0] == "file.c"
    assert result[0][1] == (
        '#include "file.h"\n#include <math.h>\n'
        'double test_fabs(double x) {\n   return fabs(x);\n}\n'
        'double test_acos(double x) {\n   return acos(x);\n}\n'
        'double test_asin(double x) {\n   return asin(x);\n}\n'
        'double test_atan(double x) {\n   return atan(x);\n}\n'
        'double test_ceil(double x) {\n   return ceil(x);\n}\n'
        'double test_cos(double x) {\n   return cos(x);\n}\n'
        'double test_cosh(double x) {\n   return cosh(x);\n}\n'
        'double test_floor(double x) {\n   return floor(x);\n}\n'
        'double test_log(double x) {\n   return log(x);\n}\n'
        'double test_ln(double x) {\n   return log(x);\n}\n'
        'double test_sin(double x) {\n   return sin(x);\n}\n'
        'double test_sinh(double x) {\n   return sinh(x);\n}\n'
        'double test_sqrt(double x) {\n   return sqrt(x);\n}\n'
        'double test_tan(double x) {\n   return tan(x);\n}\n'
        'double test_tanh(double x) {\n   return tanh(x);\n}\n'
    )
    assert result[1][0] == "file.h"
    assert result[1][1] == (
        '#ifndef PROJECT__FILE__H\n#define PROJECT__FILE__H\n'
        'double test_fabs(double x);\ndouble test_acos(double x);\n'
        'double test_asin(double x);\ndouble test_atan(double x);\n'
        'double test_ceil(double x);\ndouble test_cos(double x);\n'
        'double test_cosh(double x);\ndouble test_floor(double x);\n'
        'double test_log(double x);\ndouble test_ln(double x);\n'
        'double test_sin(double x);\ndouble test_sinh(double x);\n'
        'double test_sqrt(double x);\ndouble test_tan(double x);\n'
        'double test_tanh(double x);\n#endif\n'
    )


def test_ansi_math2_codegen():
    # not included: frexp, ldexp, modf, fmod
    from sympy import atan2
    x, y = symbols('x,y')
    name_expr = [
        ("test_atan2", atan2(x, y)),
        ("test_pow", x**y),
    ]
    result = codegen(name_expr, "C", "file", header=False, empty=False)
    assert result[0][0] == "file.c"
    assert result[0][1] == (
        '#include "file.h"\n#include <math.h>\n'
        'double test_atan2(double x, double y) {\n   return atan2(x, y);\n}\n'
        'double test_pow(double x, double y) {\n   return pow(x, y);\n}\n'
    )
    assert result[1][0] == "file.h"
    assert result[1][1] == (
        '#ifndef PROJECT__FILE__H\n#define PROJECT__FILE__H\n'
        'double test_atan2(double x, double y);\n'
        'double test_pow(double x, double y);\n'
        '#endif\n'
    )


def test_complicated_codegen():
    from sympy import sin, cos, tan
    x, y, z = symbols('x,y,z')
    name_expr = [
        ("test1", ((sin(x) + cos(y) + tan(z))**7).expand()),
        ("test2", cos(cos(cos(cos(cos(cos(cos(cos(x + y + z))))))))),
    ]
    result = codegen(name_expr, "C", "file", header=False, empty=False)
    assert result[0][0] == "file.c"
    assert result[0][1] == (
        '#include "file.h"\n#include <math.h>\n'
        'double test1(double x, double y, double z) {\n'
        '   return '
        'pow(sin(x), 7) + '
        '7*pow(sin(x), 6)*cos(y) + '
        '7*pow(sin(x), 6)*tan(z) + '
        '21*pow(sin(x), 5)*pow(cos(y), 2) + '
        '42*pow(sin(x), 5)*cos(y)*tan(z) + '
        '21*pow(sin(x), 5)*pow(tan(z), 2) + '
        '35*pow(sin(x), 4)*pow(cos(y), 3) + '
        '105*pow(sin(x), 4)*pow(cos(y), 2)*tan(z) + '
        '105*pow(sin(x), 4)*cos(y)*pow(tan(z), 2) + '
        '35*pow(sin(x), 4)*pow(tan(z), 3) + '
        '35*pow(sin(x), 3)*pow(cos(y), 4) + '
        '140*pow(sin(x), 3)*pow(cos(y), 3)*tan(z) + '
        '210*pow(sin(x), 3)*pow(cos(y), 2)*pow(tan(z), 2) + '
        '140*pow(sin(x), 3)*cos(y)*pow(tan(z), 3) + '
        '35*pow(sin(x), 3)*pow(tan(z), 4) + '
        '21*pow(sin(x), 2)*pow(cos(y), 5) + '
        '105*pow(sin(x), 2)*pow(cos(y), 4)*tan(z) + '
        '210*pow(sin(x), 2)*pow(cos(y), 3)*pow(tan(z), 2) + '
        '210*pow(sin(x), 2)*pow(cos(y), 2)*pow(tan(z), 3) + '
        '105*pow(sin(x), 2)*cos(y)*pow(tan(z), 4) + '
        '21*pow(sin(x), 2)*pow(tan(z), 5) + '
        '7*sin(x)*pow(cos(y), 6) + '
        '42*sin(x)*pow(cos(y), 5)*tan(z) + '
        '105*sin(x)*pow(cos(y), 4)*pow(tan(z), 2) + '
        '140*sin(x)*pow(cos(y), 3)*pow(tan(z), 3) + '
        '105*sin(x)*pow(cos(y), 2)*pow(tan(z), 4) + '
        '42*sin(x)*cos(y)*pow(tan(z), 5) + '
        '7*sin(x)*pow(tan(z), 6) + '
        'pow(cos(y), 7) + '
        '7*pow(cos(y), 6)*tan(z) + '
        '21*pow(cos(y), 5)*pow(tan(z), 2) + '
        '35*pow(cos(y), 4)*pow(tan(z), 3) + '
        '35*pow(cos(y), 3)*pow(tan(z), 4) + '
        '21*pow(cos(y), 2)*pow(tan(z), 5) + '
        '7*cos(y)*pow(tan(z), 6) + '
        'pow(tan(z), 7);\n'
        '}\n'
        'double test2(double x, double y, double z) {\n'
        '   return cos(cos(cos(cos(cos(cos(cos(cos(x + y + z))))))));\n'
        '}\n'
    )
    assert result[1][0] == "file.h"
    assert result[1][1] == (
        '#ifndef PROJECT__FILE__H\n'
        '#define PROJECT__FILE__H\n'
        'double test1(double x, double y, double z);\n'
        'double test2(double x, double y, double z);\n'
        '#endif\n'
    )


def test_loops_c():
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols
    n, m = symbols('n m', integer=True)
    A = IndexedBase('A')
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx('i', m)
    j = Idx('j', n)

    (f1, code), (f2, interface) = codegen(
        ('matrix_vector', Eq(y[i], A[i, j]*x[j])), "C", "file", header=False, empty=False)

    assert f1 == 'file.c'
    expected = (
        '#include "file.h"\n'
        '#include <math.h>\n'
        'void matrix_vector(double *A, int m, int n, double *x, double *y) {\n'
        '   for (int i=0; i<m; i++){\n'
        '      y[i] = 0;\n'
        '   }\n'
        '   for (int i=0; i<m; i++){\n'
        '      for (int j=0; j<n; j++){\n'
        '         y[i] = %(rhs)s + y[i];\n'
        '      }\n'
        '   }\n'
        '}\n'
    )

    assert (code == expected % {'rhs': 'A[%s]*x[j]' % (i*n + j)} or
            code == expected % {'rhs': 'A[%s]*x[j]' % (j + i*n)} or
            code == expected % {'rhs': 'x[j]*A[%s]' % (i*n + j)} or
            code == expected % {'rhs': 'x[j]*A[%s]' % (j + i*n)})
    assert f2 == 'file.h'
    assert interface == (
        '#ifndef PROJECT__FILE__H\n'
        '#define PROJECT__FILE__H\n'
        'void matrix_vector(double *A, int m, int n, double *x, double *y);\n'
        '#endif\n'
    )


def test_dummy_loops_c():
    from sympy.tensor import IndexedBase, Idx
    # the following line could also be
    # [Dummy(s, integer=True) for s in 'im']
    # or [Dummy(integer=True) for s in 'im']
    i, m = symbols('i m', integer=True, cls=Dummy)
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx(i, m)
    expected = (
        '#include "file.h"\n'
        '#include <math.h>\n'
        'void test_dummies(int m_%(mno)i, double *x, double *y) {\n'
        '   for (int i_%(ino)i=0; i_%(ino)i<m_%(mno)i; i_%(ino)i++){\n'
        '      y[i_%(ino)i] = x[i_%(ino)i];\n'
        '   }\n'
        '}\n'
    ) % {'ino': i.label.dummy_index, 'mno': m.dummy_index}
    r = Routine('test_dummies', Eq(y[i], x[i]))
    c = CCodeGen()
    code = get_string(c.dump_c, [r])
    assert code == expected


def test_partial_loops_c():
    # check that loop boundaries are determined by Idx, and array strides
    # determined by shape of IndexedBase object.
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols
    n, m, o, p = symbols('n m o p', integer=True)
    A = IndexedBase('A', shape=(m, p))
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx('i', (o, m - 5))  # Note: bounds are inclusive
    j = Idx('j', n)          # dimension n corresponds to bounds (0, n - 1)

    (f1, code), (f2, interface) = codegen(
        ('matrix_vector', Eq(y[i], A[i, j]*x[j])), "C", "file", header=False, empty=False)

    assert f1 == 'file.c'
    expected = (
        '#include "file.h"\n'
        '#include <math.h>\n'
        'void matrix_vector(double *A, int m, int n, int o, int p, double *x, double *y) {\n'
        '   for (int i=o; i<%(upperi)s; i++){\n'
        '      y[i] = 0;\n'
        '   }\n'
        '   for (int i=o; i<%(upperi)s; i++){\n'
        '      for (int j=0; j<n; j++){\n'
        '         y[i] = %(rhs)s + y[i];\n'
        '      }\n'
        '   }\n'
        '}\n'
    ) % {'upperi': m - 4, 'rhs': '%(rhs)s'}

    assert (code == expected % {'rhs': 'A[%s]*x[j]' % (i*p + j)} or
            code == expected % {'rhs': 'A[%s]*x[j]' % (j + i*p)} or
            code == expected % {'rhs': 'x[j]*A[%s]' % (i*p + j)} or
            code == expected % {'rhs': 'x[j]*A[%s]' % (j + i*p)})
    assert f2 == 'file.h'
    assert interface == (
        '#ifndef PROJECT__FILE__H\n'
        '#define PROJECT__FILE__H\n'
        'void matrix_vector(double *A, int m, int n, int o, int p, double *x, double *y);\n'
        '#endif\n'
    )


def test_output_arg_c():
    from sympy import sin, cos, Equality
    x, y, z = symbols("x,y,z")
    r = Routine("foo", [Equality(y, sin(x)), cos(x)])
    c = CCodeGen()
    result = c.write([r], "test", header=False, empty=False)
    assert result[0][0] == "test.c"
    expected = (
        '#include "test.h"\n'
        '#include <math.h>\n'
        'double foo(double x, double &y) {\n'
        '   y = sin(x);\n'
        '   return cos(x);\n'
        '}\n'
    )
    assert result[0][1] == expected


def test_empty_f_code():
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_f95, [])
    assert source == ""


def test_empty_f_code_with_header():
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_f95, [], header=True)
    assert source[:82] == (
        "!******************************************************************************\n!*"
    )
          #   "                    Code generated with sympy 0.7.2-git                    "
    assert source[158:] == (                                                              "*\n"
            "!*                                                                            *\n"
            "!*              See http://www.sympy.org/ for more information.               *\n"
            "!*                                                                            *\n"
            "!*                       This file is part of 'project'                       *\n"
            "!******************************************************************************\n"
            )


def test_empty_f_header():
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_h, [])
    assert source == ""


def test_simple_f_code():
    x, y, z = symbols('x,y,z')
    expr = (x + y)*z
    routine = Routine("test", expr)
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_f95, [routine])
    expected = (
        "REAL*8 function test(x, y, z)\n"
        "implicit none\n"
        "REAL*8, intent(in) :: x\n"
        "REAL*8, intent(in) :: y\n"
        "REAL*8, intent(in) :: z\n"
        "test = z*(x + y)\n"
        "end function\n"
    )
    assert source == expected


def test_numbersymbol_f_code():
    routine = Routine("test", pi**Catalan)
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_f95, [routine])
    expected = (
        "REAL*8 function test()\n"
        "implicit none\n"
        "REAL*8, parameter :: Catalan = 0.915965594177219d0\n"
        "REAL*8, parameter :: pi = 3.14159265358979d0\n"
        "test = pi**Catalan\n"
        "end function\n"
    )
    assert source == expected

def test_erf_f_code():
    x = symbols('x')
    routine = Routine("test", erf(x) - erf(-2 * x))
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_f95, [routine])
    expected = (
        "REAL*8 function test(x)\n"
        "implicit none\n"
        "REAL*8, intent(in) :: x\n"
        "test = erf(x) + erf(2.0d0*x)\n"
        "end function\n"
    )
    assert source == expected, source

def test_f_code_argument_order():
    x, y, z = symbols('x,y,z')
    expr = x + y
    routine = Routine("test", expr, argument_sequence=[z, x, y])
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_f95, [routine])
    expected = (
        "REAL*8 function test(z, x, y)\n"
        "implicit none\n"
        "REAL*8, intent(in) :: z\n"
        "REAL*8, intent(in) :: x\n"
        "REAL*8, intent(in) :: y\n"
        "test = x + y\n"
        "end function\n"
    )
    assert source == expected


def test_simple_f_header():
    x, y, z = symbols('x,y,z')
    expr = (x + y)*z
    routine = Routine("test", expr)
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_h, [routine])
    expected = (
        "interface\n"
        "REAL*8 function test(x, y, z)\n"
        "implicit none\n"
        "REAL*8, intent(in) :: x\n"
        "REAL*8, intent(in) :: y\n"
        "REAL*8, intent(in) :: z\n"
        "end function\n"
        "end interface\n"
    )
    assert source == expected


def test_simple_f_codegen():
    x, y, z = symbols('x,y,z')
    expr = (x + y)*z
    result = codegen(
        ("test", expr), "F95", "file", header=False, empty=False)
    expected = [
        ("file.f90",
        "REAL*8 function test(x, y, z)\n"
        "implicit none\n"
        "REAL*8, intent(in) :: x\n"
        "REAL*8, intent(in) :: y\n"
        "REAL*8, intent(in) :: z\n"
        "test = z*(x + y)\n"
        "end function\n"),
        ("file.h",
        "interface\n"
        "REAL*8 function test(x, y, z)\n"
        "implicit none\n"
        "REAL*8, intent(in) :: x\n"
        "REAL*8, intent(in) :: y\n"
        "REAL*8, intent(in) :: z\n"
        "end function\n"
        "end interface\n")
    ]
    assert result == expected


def test_multiple_results_f():
    x, y, z = symbols('x,y,z')
    expr1 = (x + y)*z
    expr2 = (x - y)*z
    routine = Routine(
        "test",
        [expr1, expr2]
    )
    code_gen = FCodeGen()
    raises(CodeGenError, lambda: get_string(code_gen.dump_h, [routine]))


def test_no_results_f():
    raises(ValueError, lambda: Routine("test", []))


def test_intrinsic_math_codegen():
    # not included: log10
    from sympy import (acos, asin, atan, ceiling, cos, cosh, floor, log, ln,
            sin, sinh, sqrt, tan, tanh, Abs)
    x = symbols('x')
    name_expr = [
        ("test_abs", Abs(x)),
        ("test_acos", acos(x)),
        ("test_asin", asin(x)),
        ("test_atan", atan(x)),
        # ("test_ceil", ceiling(x)),
        ("test_cos", cos(x)),
        ("test_cosh", cosh(x)),
        # ("test_floor", floor(x)),
        ("test_log", log(x)),
        ("test_ln", ln(x)),
        ("test_sin", sin(x)),
        ("test_sinh", sinh(x)),
        ("test_sqrt", sqrt(x)),
        ("test_tan", tan(x)),
        ("test_tanh", tanh(x)),
    ]
    result = codegen(name_expr, "F95", "file", header=False, empty=False)
    assert result[0][0] == "file.f90"
    expected = (
        'REAL*8 function test_abs(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_abs = Abs(x)\n'
        'end function\n'
        'REAL*8 function test_acos(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_acos = acos(x)\n'
        'end function\n'
        'REAL*8 function test_asin(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_asin = asin(x)\n'
        'end function\n'
        'REAL*8 function test_atan(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_atan = atan(x)\n'
        'end function\n'
        'REAL*8 function test_cos(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_cos = cos(x)\n'
        'end function\n'
        'REAL*8 function test_cosh(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_cosh = cosh(x)\n'
        'end function\n'
        'REAL*8 function test_log(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_log = log(x)\n'
        'end function\n'
        'REAL*8 function test_ln(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_ln = log(x)\n'
        'end function\n'
        'REAL*8 function test_sin(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_sin = sin(x)\n'
        'end function\n'
        'REAL*8 function test_sinh(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_sinh = sinh(x)\n'
        'end function\n'
        'REAL*8 function test_sqrt(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_sqrt = sqrt(x)\n'
        'end function\n'
        'REAL*8 function test_tan(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_tan = tan(x)\n'
        'end function\n'
        'REAL*8 function test_tanh(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_tanh = tanh(x)\n'
        'end function\n'
    )
    assert result[0][1] == expected

    assert result[1][0] == "file.h"
    expected = (
        'interface\n'
        'REAL*8 function test_abs(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_acos(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_asin(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_atan(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_cos(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_cosh(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_log(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_ln(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_sin(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_sinh(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_sqrt(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_tan(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_tanh(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
    )
    assert result[1][1] == expected


def test_intrinsic_math2_codegen():
    # not included: frexp, ldexp, modf, fmod
    from sympy import atan2
    x, y = symbols('x,y')
    name_expr = [
        ("test_atan2", atan2(x, y)),
        ("test_pow", x**y),
    ]
    result = codegen(name_expr, "F95", "file", header=False, empty=False)
    assert result[0][0] == "file.f90"
    expected = (
        'REAL*8 function test_atan2(x, y)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'test_atan2 = atan2(x, y)\n'
        'end function\n'
        'REAL*8 function test_pow(x, y)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'test_pow = x**y\n'
        'end function\n'
    )
    assert result[0][1] == expected

    assert result[1][0] == "file.h"
    expected = (
        'interface\n'
        'REAL*8 function test_atan2(x, y)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_pow(x, y)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'end function\n'
        'end interface\n'
    )
    assert result[1][1] == expected


def test_complicated_codegen_f95():
    from sympy import sin, cos, tan
    x, y, z = symbols('x,y,z')
    name_expr = [
        ("test1", ((sin(x) + cos(y) + tan(z))**7).expand()),
        ("test2", cos(cos(cos(cos(cos(cos(cos(cos(x + y + z))))))))),
    ]
    result = codegen(name_expr, "F95", "file", header=False, empty=False)
    assert result[0][0] == "file.f90"
    expected = (
        'REAL*8 function test1(x, y, z)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'REAL*8, intent(in) :: z\n'
        'test1 = sin(x)**7 + 7*sin(x)**6*cos(y) + 7*sin(x)**6*tan(z) + 21*sin(x) &\n'
        '      **5*cos(y)**2 + 42*sin(x)**5*cos(y)*tan(z) + 21*sin(x)**5*tan(z) &\n'
        '      **2 + 35*sin(x)**4*cos(y)**3 + 105*sin(x)**4*cos(y)**2*tan(z) + &\n'
        '      105*sin(x)**4*cos(y)*tan(z)**2 + 35*sin(x)**4*tan(z)**3 + 35*sin( &\n'
        '      x)**3*cos(y)**4 + 140*sin(x)**3*cos(y)**3*tan(z) + 210*sin(x)**3* &\n'
        '      cos(y)**2*tan(z)**2 + 140*sin(x)**3*cos(y)*tan(z)**3 + 35*sin(x) &\n'
        '      **3*tan(z)**4 + 21*sin(x)**2*cos(y)**5 + 105*sin(x)**2*cos(y)**4* &\n'
        '      tan(z) + 210*sin(x)**2*cos(y)**3*tan(z)**2 + 210*sin(x)**2*cos(y) &\n'
        '      **2*tan(z)**3 + 105*sin(x)**2*cos(y)*tan(z)**4 + 21*sin(x)**2*tan &\n'
        '      (z)**5 + 7*sin(x)*cos(y)**6 + 42*sin(x)*cos(y)**5*tan(z) + 105* &\n'
        '      sin(x)*cos(y)**4*tan(z)**2 + 140*sin(x)*cos(y)**3*tan(z)**3 + 105 &\n'
        '      *sin(x)*cos(y)**2*tan(z)**4 + 42*sin(x)*cos(y)*tan(z)**5 + 7*sin( &\n'
        '      x)*tan(z)**6 + cos(y)**7 + 7*cos(y)**6*tan(z) + 21*cos(y)**5*tan( &\n'
        '      z)**2 + 35*cos(y)**4*tan(z)**3 + 35*cos(y)**3*tan(z)**4 + 21*cos( &\n'
        '      y)**2*tan(z)**5 + 7*cos(y)*tan(z)**6 + tan(z)**7\n'
        'end function\n'
        'REAL*8 function test2(x, y, z)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'REAL*8, intent(in) :: z\n'
        'test2 = cos(cos(cos(cos(cos(cos(cos(cos(x + y + z))))))))\n'
        'end function\n'
    )
    assert result[0][1] == expected
    assert result[1][0] == "file.h"
    expected = (
        'interface\n'
        'REAL*8 function test1(x, y, z)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'REAL*8, intent(in) :: z\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test2(x, y, z)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'REAL*8, intent(in) :: z\n'
        'end function\n'
        'end interface\n'
    )
    assert result[1][1] == expected


def test_loops():
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols

    n, m = symbols('n,m', integer=True)
    A, x, y = map(IndexedBase, 'Axy')
    i = Idx('i', m)
    j = Idx('j', n)

    (f1, code), (f2, interface) = codegen(
        ('matrix_vector', Eq(y[i], A[i, j]*x[j])), "F95", "file", header=False, empty=False)

    assert f1 == 'file.f90'
    expected = (
        'subroutine matrix_vector(A, m, n, x, y)\n'
        'implicit none\n'
        'INTEGER*4, intent(in) :: m\n'
        'INTEGER*4, intent(in) :: n\n'
        'REAL*8, intent(in), dimension(1:m, 1:n) :: A\n'
        'REAL*8, intent(in), dimension(1:n) :: x\n'
        'REAL*8, intent(out), dimension(1:m) :: y\n'
        'INTEGER*4 :: i\n'
        'INTEGER*4 :: j\n'
        'do i = 1, m\n'
        '   y(i) = 0\n'
        'end do\n'
        'do i = 1, m\n'
        '   do j = 1, n\n'
        '      y(i) = %(rhs)s + y(i)\n'
        '   end do\n'
        'end do\n'
        'end subroutine\n'
    )

    assert code == expected % {'rhs': 'A(i, j)*x(j)'} or\
        code == expected % {'rhs': 'x(j)*A(i, j)'}
    assert f2 == 'file.h'
    assert interface == (
        'interface\n'
        'subroutine matrix_vector(A, m, n, x, y)\n'
        'implicit none\n'
        'INTEGER*4, intent(in) :: m\n'
        'INTEGER*4, intent(in) :: n\n'
        'REAL*8, intent(in), dimension(1:m, 1:n) :: A\n'
        'REAL*8, intent(in), dimension(1:n) :: x\n'
        'REAL*8, intent(out), dimension(1:m) :: y\n'
        'end subroutine\n'
        'end interface\n'
    )


def test_dummy_loops_f95():
    from sympy.tensor import IndexedBase, Idx
    # the following line could also be
    # [Dummy(s, integer=True) for s in 'im']
    # or [Dummy(integer=True) for s in 'im']
    i, m = symbols('i m', integer=True, cls=Dummy)
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx(i, m)
    expected = (
        'subroutine test_dummies(m_%(mcount)i, x, y)\n'
        'implicit none\n'
        'INTEGER*4, intent(in) :: m_%(mcount)i\n'
        'REAL*8, intent(in), dimension(1:m_%(mcount)i) :: x\n'
        'REAL*8, intent(out), dimension(1:m_%(mcount)i) :: y\n'
        'INTEGER*4 :: i_%(icount)i\n'
        'do i_%(icount)i = 1, m_%(mcount)i\n'
        '   y(i_%(icount)i) = x(i_%(icount)i)\n'
        'end do\n'
        'end subroutine\n'
    ) % {'icount': i.label.dummy_index, 'mcount': m.dummy_index}
    r = Routine('test_dummies', Eq(y[i], x[i]))
    c = FCodeGen()
    code = get_string(c.dump_f95, [r])
    assert code == expected


def test_loops_InOut():
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols

    i, j, n, m = symbols('i,j,n,m', integer=True)
    A, x, y = symbols('A,x,y')
    A = IndexedBase(A)[Idx(i, m), Idx(j, n)]
    x = IndexedBase(x)[Idx(j, n)]
    y = IndexedBase(y)[Idx(i, m)]

    (f1, code), (f2, interface) = codegen(
        ('matrix_vector', Eq(y, y + A*x)), "F95", "file", header=False, empty=False)

    assert f1 == 'file.f90'
    expected = (
        'subroutine matrix_vector(A, m, n, x, y)\n'
        'implicit none\n'
        'INTEGER*4, intent(in) :: m\n'
        'INTEGER*4, intent(in) :: n\n'
        'REAL*8, intent(in), dimension(1:m, 1:n) :: A\n'
        'REAL*8, intent(in), dimension(1:n) :: x\n'
        'REAL*8, intent(inout), dimension(1:m) :: y\n'
        'INTEGER*4 :: i\n'
        'INTEGER*4 :: j\n'
        'do i = 1, m\n'
        '   do j = 1, n\n'
        '      y(i) = %(rhs)s + y(i)\n'
        '   end do\n'
        'end do\n'
        'end subroutine\n'
    )

    assert (code == expected % {'rhs': 'A(i, j)*x(j)'} or
            code == expected % {'rhs': 'x(j)*A(i, j)'})
    assert f2 == 'file.h'
    assert interface == (
        'interface\n'
        'subroutine matrix_vector(A, m, n, x, y)\n'
        'implicit none\n'
        'INTEGER*4, intent(in) :: m\n'
        'INTEGER*4, intent(in) :: n\n'
        'REAL*8, intent(in), dimension(1:m, 1:n) :: A\n'
        'REAL*8, intent(in), dimension(1:n) :: x\n'
        'REAL*8, intent(inout), dimension(1:m) :: y\n'
        'end subroutine\n'
        'end interface\n'
    )


def test_partial_loops_f():
    # check that loop boundaries are determined by Idx, and array strides
    # determined by shape of IndexedBase object.
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols
    n, m, o, p = symbols('n m o p', integer=True)
    A = IndexedBase('A', shape=(m, p))
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx('i', (o, m - 5))  # Note: bounds are inclusive
    j = Idx('j', n)          # dimension n corresponds to bounds (0, n - 1)

    (f1, code), (f2, interface) = codegen(
        ('matrix_vector', Eq(y[i], A[i, j]*x[j])), "F95", "file", header=False, empty=False)

    expected = (
        'subroutine matrix_vector(A, m, n, o, p, x, y)\n'
        'implicit none\n'
        'INTEGER*4, intent(in) :: m\n'
        'INTEGER*4, intent(in) :: n\n'
        'INTEGER*4, intent(in) :: o\n'
        'INTEGER*4, intent(in) :: p\n'
        'REAL*8, intent(in), dimension(1:m, 1:p) :: A\n'
        'REAL*8, intent(in), dimension(1:n) :: x\n'
        'REAL*8, intent(out), dimension(1:%(iup-ilow)s) :: y\n'
        'INTEGER*4 :: i\n'
        'INTEGER*4 :: j\n'
        'do i = %(ilow)s, %(iup)s\n'
        '   y(i) = 0\n'
        'end do\n'
        'do i = %(ilow)s, %(iup)s\n'
        '   do j = 1, n\n'
        '      y(i) = %(rhs)s + y(i)\n'
        '   end do\n'
        'end do\n'
        'end subroutine\n'
    ) % {
        'rhs': '%(rhs)s',
        'iup': str(m - 4),
        'ilow': str(1 + o),
        'iup-ilow': str(m - 4 - o)
    }

    assert code == expected % {'rhs': 'A(i, j)*x(j)'} or\
        code == expected % {'rhs': 'x(j)*A(i, j)'}


def test_output_arg_f():
    from sympy import sin, cos, Equality
    x, y, z = symbols("x,y,z")
    r = Routine("foo", [Equality(y, sin(x)), cos(x)])
    c = FCodeGen()
    result = c.write([r], "test", header=False, empty=False)
    assert result[0][0] == "test.f90"
    assert result[0][1] == (
        'REAL*8 function foo(x, y)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(out) :: y\n'
        'y = sin(x)\n'
        'foo = cos(x)\n'
        'end function\n'
    )


def test_inline_function():
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols
    n, m = symbols('n m', integer=True)
    A, x, y = map(IndexedBase, 'Axy')
    i = Idx('i', m)
    p = FCodeGen()
    func = implemented_function('func', Lambda(n, n*(n + 1)))
    routine = Routine('test_inline', Eq(y[i], func(x[i])))
    code = get_string(p.dump_f95, [routine])
    expected = (
        'subroutine test_inline(m, x, y)\n'
        'implicit none\n'
        'INTEGER*4, intent(in) :: m\n'
        'REAL*8, intent(in), dimension(1:m) :: x\n'
        'REAL*8, intent(out), dimension(1:m) :: y\n'
        'INTEGER*4 :: i\n'
        'do i = 1, m\n'
        '   y(i) = %s*%s\n'
        'end do\n'
        'end subroutine\n'
    )
    args = ('x(i)', '(x(i) + 1)')
    assert code == expected % args or\
        code == expected % args[::-1]


def test_check_case():
    x, X = symbols('x,X')
    raises(CodeGenError, lambda: codegen(('test', x*X), 'f95', 'prefix'))


def test_check_case_false_positive():
    # The upper case/lower case exception should not be triggered by SymPy
    # objects that differ only because of assumptions.  (It may be useful to
    # have a check for that as well, but here we only want to test against
    # false positives with respect to case checking.)
    x1 = symbols('x')
    x2 = symbols('x', my_assumption=True)
    try:
        codegen(('test', x1*x2), 'f95', 'prefix')
    except CodeGenError as e:
        if e.args[0].startswith("Fortran ignores case."):
            raise AssertionError("This exception should not be raised!")
