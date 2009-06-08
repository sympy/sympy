
from sympy import symbols
from sympy.utilities.codegen import CCodeGen, Routine, InputArgument, Result
from StringIO import StringIO
import sys

x,y,z = symbols('xyz')


def get_string(dump_fn, routines, prefix="file"):
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


