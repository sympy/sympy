from sympy import symbols, raises, Eq
from sympy.utilities.codegen import CCodeGen, Routine, InputArgument, Result, \
    codegen, CodeGenError, FCodeGen
from StringIO import StringIO

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
          #   "                    Code generated with sympy 0.6.7-git                    "
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
    x,y,z = symbols('xyz')
    expr = (x+y)*z
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

def test_simple_c_header():
    x,y,z = symbols('xyz')
    expr = (x+y)*z
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
    x,y,z = symbols('xyz')
    expr = (x+y)*z
    result = codegen(("test", (x+y)*z), "C", "file", header=False, empty=False)
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
    x,y,z = symbols('xyz')
    expr1 = (x+y)*z
    expr2 = (x-y)*z
    routine = Routine(
        "test",
        [expr1,expr2]
    )
    code_gen = CCodeGen()
    raises(CodeGenError, 'get_string(code_gen.dump_h, [routine])')

def test_no_results_c():
    raises(ValueError, 'Routine("test", [])')

def test_ansi_math1_codegen():
    # not included: log10
    from sympy import acos, asin, atan, ceiling, cos, cosh, floor, log, ln, \
        sin, sinh, sqrt, tan, tanh, N
    x = symbols('x')
    name_expr = [
        ("test_fabs", abs(x)),
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
    from sympy import atan2, N
    x, y = symbols('xy')
    name_expr = [
        ("test_atan2", atan2(x,y)),
        ("test_pow", x**y),
    ]
    result = codegen(name_expr, "C", "file", header=False, empty=False)
    assert result[0][0] == "file.c"
    assert result[0][1] == (
        '#include "file.h"\n#include <math.h>\n'
        'double test_atan2(double x, double y) {\n   return atan2(x, y);\n}\n'
        'double test_pow(double x, double y) {\n   return pow(x,y);\n}\n'
    )
    assert result[1][0] == "file.h"
    assert result[1][1] == (
        '#ifndef PROJECT__FILE__H\n#define PROJECT__FILE__H\n'
        'double test_atan2(double x, double y);\n'
        'double test_pow(double x, double y);\n'
        '#endif\n'
    )

def test_complicated_codegen():
    from sympy import sin, cos, tan, N
    x,y,z = symbols('xyz')
    name_expr = [
        ("test1", ((sin(x)+cos(y)+tan(z))**7).expand()),
        ("test2", cos(cos(cos(cos(cos(cos(cos(cos(x+y+z))))))))),
    ]
    result = codegen(name_expr, "C", "file", header=False, empty=False)
    assert result[0][0] == "file.c"
    assert result[0][1] == (
        '#include "file.h"\n#include <math.h>\n'
        'double test1(double x, double y, double z) {\n'
        '   return '
        '7*pow(cos(y),6)*sin(x) + '
        '7*pow(cos(y),6)*tan(z) + '
        '7*pow(sin(x),6)*cos(y) + '
        '7*pow(sin(x),6)*tan(z) + '
        '7*pow(tan(z),6)*cos(y) + '
        '7*pow(tan(z),6)*sin(x) + '
        '42*pow(cos(y),5)*sin(x)*tan(z) + '
        '42*pow(sin(x),5)*cos(y)*tan(z) + '
        '42*pow(tan(z),5)*cos(y)*sin(x) + '
        '105*pow(cos(y),2)*pow(sin(x),4)*tan(z) + '
        '105*pow(cos(y),2)*pow(tan(z),4)*sin(x) + '
        '105*pow(cos(y),4)*pow(sin(x),2)*tan(z) + '
        '105*pow(cos(y),4)*pow(tan(z),2)*sin(x) + '
        '105*pow(sin(x),2)*pow(tan(z),4)*cos(y) + '
        '105*pow(sin(x),4)*pow(tan(z),2)*cos(y) + '
        '140*pow(cos(y),3)*pow(sin(x),3)*tan(z) + '
        '140*pow(cos(y),3)*pow(tan(z),3)*sin(x) + '
        '140*pow(sin(x),3)*pow(tan(z),3)*cos(y) + '
        '21*pow(cos(y),5)*pow(sin(x),2) + '
        '21*pow(cos(y),5)*pow(tan(z),2) + '
        '21*pow(sin(x),5)*pow(tan(z),2) + '
        '210*pow(cos(y),2)*pow(sin(x),3)*pow(tan(z),2) + '
        '210*pow(cos(y),3)*pow(sin(x),2)*pow(tan(z),2) + '
        '35*pow(cos(y),4)*pow(sin(x),3) + '
        '35*pow(cos(y),4)*pow(tan(z),3) + '
        '35*pow(sin(x),4)*pow(tan(z),3) + '
        '210*pow(cos(y),2)*pow(sin(x),2)*pow(tan(z),3) + '
        '35*pow(cos(y),3)*pow(sin(x),4) + '
        '35*pow(cos(y),3)*pow(tan(z),4) + '
        '35*pow(sin(x),3)*pow(tan(z),4) + '
        '21*pow(cos(y),2)*pow(sin(x),5) + '
        '21*pow(cos(y),2)*pow(tan(z),5) + '
        '21*pow(sin(x),2)*pow(tan(z),5) + '
        'pow(cos(y),7) + pow(sin(x),7) + pow(tan(z),7);\n'
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
          #   "                    Code generated with sympy 0.6.7-git                    "
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
    x,y,z = symbols('xyz')
    expr = (x+y)*z
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

def test_simple_f_header():
    x,y,z = symbols('xyz')
    expr = (x+y)*z
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
    x,y,z = symbols('xyz')
    expr = (x+y)*z
    result = codegen(("test", (x+y)*z), "F95", "file", header=False, empty=False)
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
    x,y,z = symbols('xyz')
    expr1 = (x+y)*z
    expr2 = (x-y)*z
    routine = Routine(
        "test",
        [expr1,expr2]
    )
    code_gen = FCodeGen()
    raises(CodeGenError, 'get_string(code_gen.dump_h, [routine])')

def test_no_results_f():
    raises(ValueError, 'Routine("test", [])')

def test_intrinsic_math_codegen():
    # not included: log10
    from sympy import acos, asin, atan, ceiling, cos, cosh, floor, log, ln, \
            sin, sinh, sqrt, tan, tanh, N
    x = symbols('x')
    name_expr = [
            ("test_abs", abs(x)),
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
            'test_abs = abs(x)\n'
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
    expected =  (
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
    from sympy import atan2, N
    x, y = symbols('xy')
    name_expr = [
        ("test_atan2", atan2(x,y)),
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
    from sympy import sin, cos, tan, N
    x,y,z = symbols('xyz')
    name_expr = [
        ("test1", ((sin(x)+cos(y)+tan(z))**7).expand()),
        ("test2", cos(cos(cos(cos(cos(cos(cos(cos(x+y+z))))))))),
    ]
    result = codegen(name_expr, "F95", "file", header=False, empty=False)
    assert result[0][0] == "file.f90"
    expected = (
            'REAL*8 function test1(x, y, z)\n'
            'implicit none\n'
            'REAL*8, intent(in) :: x\n'
            'REAL*8, intent(in) :: y\n'
            'REAL*8, intent(in) :: z\n'
            'test1 = 7*cos(y)**6*sin(x) + 7*cos(y)**6*tan(z) + 7*sin(x)**6*cos(y) + &\n'
            '      7*sin(x)**6*tan(z) + 7*tan(z)**6*cos(y) + 7*tan(z)**6*sin(x) + &\n'
            '      42*cos(y)**5*sin(x)*tan(z) + 42*sin(x)**5*cos(y)*tan(z) + 42*tan( &\n'
            '      z)**5*cos(y)*sin(x) + 105*cos(y)**2*sin(x)**4*tan(z) + 105*cos(y) &\n'
            '      **2*tan(z)**4*sin(x) + 105*cos(y)**4*sin(x)**2*tan(z) + 105*cos(y &\n'
            '      )**4*tan(z)**2*sin(x) + 105*sin(x)**2*tan(z)**4*cos(y) + 105*sin( &\n'
            '      x)**4*tan(z)**2*cos(y) + 140*cos(y)**3*sin(x)**3*tan(z) + 140*cos &\n'
            '      (y)**3*tan(z)**3*sin(x) + 140*sin(x)**3*tan(z)**3*cos(y) + 21*cos &\n'
            '      (y)**5*sin(x)**2 + 21*cos(y)**5*tan(z)**2 + 21*sin(x)**5*tan(z) &\n'
            '      **2 + 210*cos(y)**2*sin(x)**3*tan(z)**2 + 210*cos(y)**3*sin(x) &\n'
            '      **2*tan(z)**2 + 35*cos(y)**4*sin(x)**3 + 35*cos(y)**4*tan(z)**3 + &\n'
            '      35*sin(x)**4*tan(z)**3 + 210*cos(y)**2*sin(x)**2*tan(z)**3 + 35* &\n'
            '      cos(y)**3*sin(x)**4 + 35*cos(y)**3*tan(z)**4 + 35*sin(x)**3*tan(z &\n'
            '      )**4 + 21*cos(y)**2*sin(x)**5 + 21*cos(y)**2*tan(z)**5 + 21*sin(x &\n'
            '      )**2*tan(z)**5 + cos(y)**7 + sin(x)**7 + tan(z)**7\n'
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
    from sympy.tensor import Indexed, Idx
    from sympy import symbols
    i,j,n,m = symbols('i j n m', integer=True)
    A,x,y = symbols('A x y')
    A = Indexed(A)(Idx(i, m), Idx(j, n))
    x = Indexed(x)(Idx(j, n))
    y = Indexed(y)(Idx(i, m))

    (f1, code), (f2, interface) = codegen(
            ('matrix_vector', Eq(y, A*x)), "F95", "file", header=False, empty=False)

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
            'y = 0.d0\n'
            'do i = 1, m\n'
            '   do j = 1, n\n'
            '      y(i) = A(i, j)*x(j) + y(i)\n'
            '   end do\n'
            'end do\n'
            'end subroutine\n'
            )

    assert expected == code
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

def test_loops_InOut():
    from sympy.tensor import Indexed, Idx
    from sympy import symbols
    i,j,n,m = symbols('i j n m', integer=True)
    A,x,y = symbols('A x y')
    A = Indexed(A)(Idx(i, m), Idx(j, n))
    x = Indexed(x)(Idx(j, n))
    y = Indexed(y)(Idx(i, m))

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
            '      y(i) = A(i, j)*x(j) + y(i)\n'
            '   end do\n'
            'end do\n'
            'end subroutine\n'
            )

    assert expected == code
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

def test_output_arg_f():
    from sympy import sin, cos, Equality
    x, y, z = symbols("xyz")
    r = Routine("foo", [Equality(y, sin(x)), cos(x)])
    c = FCodeGen()
    result = c.write([r], "test", header=False, empty=False)
    assert result[0][0] == "test.f90"
    assert result[0][1] == (
        'REAL*8 function foo(x, y)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(out) :: y\n'
        'y = 0.d0\n'
        'foo = 0.d0\n'
        'y = sin(x)\n'
        'foo = cos(x)\n'
        'end function\n'
    )
