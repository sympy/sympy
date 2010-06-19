# This tests the compilation and execution of the source code generated with
# utilities.codegen. The compilation takes place in a temporary directory that
# is removed after the test. By default the test directory is always removed,
# but this behavior can be changed by setting the environment variable
# SYMPY_TEST_CLEAN_TEMP to:
#   export SYMPY_TEST_CLEAN_TEMP=always   : the default behavior.
#   export SYMPY_TEST_CLEAN_TEMP=success  : only remove the directories of working tests.
#   export SYMPY_TEST_CLEAN_TEMP=never    : never remove the directories with the test code.
# When a directory is not removed, the necessary information is printed on
# screen to find the files that belong to the (failed) tests. If a test does
# not fail, py.test captures all the output and you will not see the directories
# corresponding to the successful tests. Use the --nocapture option to see all
# the output.

# All tests below have a counterpart in utilities/test/test_codegen.py. In the
# latter file, the resulting code is compared with predefined strings, without
# compilation or execution.

# All the generated Fortran code should conform with the Fortran 95 standard,
# which facilitates the incorporation in various projects. The tests below
# assume that the binary ifort is somewhere in the path and that it can compile
# standard Fortran 95 code.


from sympy import symbols
from sympy.utilities.codegen import codegen, FCodeGen, Routine, InputArgument, \
    Result
import sys, os, tempfile


# templates for the main program that will test the generated code.

main_template = """
program main
  include "codegen.h"
  integer :: result;
  result = 0

  %s

  call exit(result)
end program
"""

numerical_test_template = """
  if (abs(%(call)s)>%(threshold)s) then
    write(6,"('Numerical validation failed:')")
    write(6,"('%(call)s=',e15.5,'threshold=',e15.5)") %(call)s, %(threshold)s
    result = -1;
  end if
"""


def try_run(commands):
    """Run a series of commands and only return True if all ran fine."""
    for command in commands:
        retcode = os.system(command)
        if retcode != 0:
            return False
    return True


def run_f95_test(label, routines, numerical_tests, friendly=True):
    """A driver for the codegen tests.

       This driver assumes that a compiler ifort is present in the PATH and that
       ifort is (at least) a Fortran 90 compiler. The generated code is written in
       a temporary directory, together with a main program that validates the
       generated code. The test passes when the compilation and the validation
       run correctly.
    """
    # Do all the magic to compile, run and validate the test code
    # 1) prepare the temporary working directory, switch to that dir
    work = tempfile.mkdtemp("_sympy_f95_test", "%s_" % label)
    oldwork = os.getcwd()
    os.chdir(work)
    # 2) write the generated code
    if friendly:
        # interpret the routines as a name_expr list and call the friendly
        # function codegen
        codegen(routines, 'F95', "codegen", to_files=True)
    else:
        code_gen = FCodeGen()
        code_gen.write(routines, "codegen", to_files=True)
    # 3) write a simple main program that links to the generated code, and that
    #    includes the numerical tests
    test_strings = []
    for fn_name, args, expected, threshold in numerical_tests:
        call_string = "%s(%s)-(%s)" % (fn_name, ",".join(str(arg) for arg in args), expected)
        call_string = fortranize_double_constants(call_string)
        threshold = fortranize_double_constants(str(threshold))
        test_strings.append(numerical_test_template % {
            "call": call_string,
            "threshold": threshold,
        })
    f = file("main.f90", "w")
    f.write(main_template % "".join(test_strings))
    f.close()
    # 4) Compile and link
    compiled = try_run([
        "ifort -c codegen.f90 -o codegen.o",
        "ifort -c main.f90 -o main.o",
        "ifort main.o codegen.o -o test.exe"
    ])
    # 5) Run if compiled
    if compiled:
        executed = try_run(["./test.exe"])
    else:
        executed = False
    # 6) Clean up stuff
    clean = os.getenv('SYMPY_TEST_CLEAN_TEMP', 'always').lower()
    if clean not in ('always', 'success', 'never'):
        raise ValueError("SYMPY_TEST_CLEAN_TEMP must be one of the following: 'always', 'success' or 'never'.")
    if clean == 'always' or (clean == 'success' and compiled and executed):
        def safe_remove(filename):
            if os.path.isfile(filename):
                os.remove(filename)
        safe_remove("codegen.f90")
        safe_remove("codegen.h")
        safe_remove("codegen.o")
        safe_remove("main.f90")
        safe_remove("main.o")
        safe_remove("test.exe")
        os.chdir(oldwork)
        os.rmdir(work)
    else:
        print >> sys.stderr, "TEST NOT REMOVED: %s" % work
        os.chdir(oldwork)
    # 7) Do the assertions in the end
    assert compiled
    assert executed

def fortranize_double_constants(code_string):
    """
    Replaces every literal float with literal doubles
    """
    import re
    pattern_exp = re.compile('\d+(\.)?\d*[eE]-?\d+')
    pattern_float = re.compile('\d+\.\d*(?!\d*d)')

    def subs_exp(matchobj):
        return re.sub('[eE]','d',matchobj.group(0))

    def subs_float(matchobj):
        return "%sd0" % matchobj.group(0)

    code_string = pattern_exp.sub(subs_exp, code_string)
    code_string = pattern_float.sub(subs_float, code_string)

    return code_string



def is_feasible():
    # This test should always work, otherwise the f95 compiler is not present.
    x,y,z = symbols('xyz')
    expr = (x+y)*z
    routine = Routine("test", expr)
    numerical_tests = [
        ("test", (1.0, 6.0, 3.0), 21.0, 1e-15),
        ("test", (-1.0, 2.0, -2.5), -2.5, 1e-15),
    ]
    run_f95_test("is_feasible", [routine], numerical_tests, friendly=False)


try:
    is_feasible()
except AssertionError:
    disabled = True


def test_basic():
    is_feasible()


def test_basic_codegen():
    x,y,z = symbols('xyz')
    numerical_tests = [
        ("test", (1.0, 6.0, 3.0), 21.0, 1e-15),
        ("test", (-1.0, 2.0, -2.5), -2.5, 1e-15),
    ]
    run_f95_test("basic_codegen", [("test", (x+y)*z)], numerical_tests)

def test_intrinsic_math1_codegen():
    # not included: log10
    from sympy import acos, asin, atan, ceiling, cos, cosh, floor, log, ln, \
        sin, sinh, sqrt, tan, tanh, N
    x = symbols('x')
    name_expr = [
        ("test_fabs", abs(x)),
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
    numerical_tests = []
    for name, expr in name_expr:
        for xval in 0.2, 0.5, 0.8:
            expected = N(expr.subs(x, xval))
            numerical_tests.append((name, (xval,), expected, 1e-14))
    run_f95_test("intrinsic_math1", name_expr, numerical_tests)

def test_instrinsic_math2_codegen():
    # not included: frexp, ldexp, modf, fmod
    from sympy import atan2, N
    x, y = symbols('xy')
    name_expr = [
        ("test_atan2", atan2(x,y)),
        ("test_pow", x**y),
    ]
    numerical_tests = []
    for name, expr in name_expr:
        for xval,yval in (0.2, 1.3), (0.5, -0.2), (0.8, 0.8):
            expected = N(expr.subs(x, xval).subs(y, yval))
            numerical_tests.append((name, (xval,yval), expected, 1e-14))
    run_f95_test("intrinsic_math2", name_expr, numerical_tests)

def test_complicated_codegen():
    from sympy import sin, cos, tan, N
    x,y,z = symbols('xyz')
    name_expr = [
        ("test1", ((sin(x)+cos(y)+tan(z))**7).expand()),
        ("test2", cos(cos(cos(cos(cos(cos(cos(cos(x+y+z))))))))),
    ]
    numerical_tests = []
    for name, expr in name_expr:
        for xval,yval,zval in (0.2, 1.3, -0.3), (0.5, -0.2, 0.0), (0.8, 2.1, 0.8):
            expected = N(expr.subs(x, xval).subs(y, yval).subs(z, zval))
            numerical_tests.append((name, (xval,yval,zval), expected, 1e-12))
    run_f95_test("complicated_codegen", name_expr, numerical_tests)
