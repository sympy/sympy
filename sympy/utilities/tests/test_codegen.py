from sympy import symbols, raises
from sympy.utilities.codegen import CCodeGen, Routine, InputArgument, Result, \
    codegen, CodeGenError
from StringIO import StringIO

def get_string(dump_fn, routines, prefix="file"):
    """Wrapper for dump_fn. dump_fn writes its results to a stream object and
       this wrapper returns the contents of that stream as a string. This
       auxiliary function is used by many tests below.

       The header and the empty lines are not generator to facilitate the
       testing of the output.
    """
    output = StringIO()
    dump_fn(routines, output, prefix, header=False, empty=False)
    source = output.getvalue()
    output.close()
    return source

def test_empty_c_code():
    code_gen = CCodeGen()
    source = get_string(code_gen.dump_c, [])
    assert source == "#include \"file.h\"\n#include <math.h>\n"

def test_empty_c_header():
    code_gen = CCodeGen()
    source = get_string(code_gen.dump_h, [])
    assert source == "#ifndef PROJECT__FILE__H\n#define PROJECT__FILE__H\n#endif\n"

def test_simple_c_code():
    x,y,z = symbols('xyz')
    expr = (x+y)*z
    routine = Routine("test", [InputArgument(symbol) for symbol in x,y,z], [Result(expr)])
    code_gen = CCodeGen()
    source = get_string(code_gen.dump_c, [routine])
    expected = (
        "#include \"file.h\"\n"
        "#include <math.h>\n"
        "double test(double x, double y, double z) {\n"
        "  return z*(x + y);\n"
        "}\n"
    )
    assert source == expected

def test_simple_c_header():
    x,y,z = symbols('xyz')
    expr = (x+y)*z
    routine = Routine("test", [InputArgument(symbol) for symbol in x,y,z], [Result(expr)])
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
        "  return z*(x + y);\n"
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
        [InputArgument(symbol) for symbol in x,y,z],
        [Result(expr1),Result(expr2)]
    )
    code_gen = CCodeGen()
    raises(CodeGenError, 'get_string(code_gen.dump_h, [routine])')

def test_no_results_c():
    x = symbols('x')
    raises(ValueError, 'Routine("test", [InputArgument(x)], [])')

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
        'double test_fabs(double x) {\n  return fabs(x);\n}\n'
        'double test_acos(double x) {\n  return acos(x);\n}\n'
        'double test_asin(double x) {\n  return asin(x);\n}\n'
        'double test_atan(double x) {\n  return atan(x);\n}\n'
        'double test_ceil(double x) {\n  return ceil(x);\n}\n'
        'double test_cos(double x) {\n  return cos(x);\n}\n'
        'double test_cosh(double x) {\n  return cosh(x);\n}\n'
        'double test_floor(double x) {\n  return floor(x);\n}\n'
        'double test_log(double x) {\n  return log(x);\n}\n'
        'double test_ln(double x) {\n  return log(x);\n}\n'
        'double test_sin(double x) {\n  return sin(x);\n}\n'
        'double test_sinh(double x) {\n  return sinh(x);\n}\n'
        'double test_sqrt(double x) {\n  return pow(x,(1.0/2.0));\n}\n'
        'double test_tan(double x) {\n  return tan(x);\n}\n'
        'double test_tanh(double x) {\n  return tanh(x);\n}\n'
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
        'double test_atan2(double x, double y) {\n  return atan2(x, y);\n}\n'
        'double test_pow(double x, double y) {\n  return pow(x,y);\n}\n'
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
        '  return '
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
        '  return cos(cos(cos(cos(cos(cos(cos(cos(x + y + z))))))));\n'
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
