# This tests the compilation and execution of the source code generated with
# utilities.codegen. The compilation takes place in a temporary directory that
# is removed after the test. By default the test directory is always removed,
# but this behavior can be changed by setting the environment variable
# SYMPY_TEST_CLEAN_TEMP to:
#   'always': the default behavior
#   'success': only remove the directories of working tests
# When a directory is not removed, the necessary information is printed on
# screen to find the files that belong to the (failed) tests.

# All tests below have a counterpart in utilities/test/test_codegen.py. In the
# latter file, the resulting code is compared with predefined strings, without
# compilation or execution.

# All the generated C code should be ANSI C, which facilitates the incorporation
# in various projects. The tests below assume that the binary cc is somewhere in
# the path and that it can compile ANSI C code.


from sympy import symbols
from sympy.utilities.codegen import codegen, CCodeGen, Routine, InputArgument, \
    Result
import sys, os, tempfile

x,y,z = symbols('xyz')

main_template = """
#include "codegen.h"
#include <stdio.h>
#include <math.h>

int main() {
  int result = 0;
  %s

  return result;
}
"""

numerical_test_template = """
  if (fabs(%(call)s)>%(threshold)s) {
    printf("Numerical validation failed: %(call)s=%%e threshold=%(threshold)s\\n", %(call)s);
    result = -1;
  }
"""


def try_run(commands):
    for command in commands:
        retcode = os.system(command)
        if retcode != 0:
            return False
    return True


def run_cc_test(label, routines, numerical_tests, friendly=True):
    # Do all the magic to compile, run and validate the test code
    # 1) prepare the temporary working directory, swith to that dir
    work = tempfile.mkdtemp("sympy_cc_test", label)
    oldwork = os.getcwd()
    os.chdir(work)
    # 2) write the generated code
    if friendly:
        # interpret the routines as a name_expr list and call the friendly
        # function codegen
        codegen(routines, 'C', "codegen")
    else:
        code_gen = CCodeGen()
        code_gen.write(routines, "codegen")
    # 3) write a simple main program that links to the generated code, and that
    #    includes the numerical tests
    test_strings = []
    for fn_name, args, expected, threshold in numerical_tests:
        call_string = "%s(%s)-(%s)" % (fn_name, ",".join(str(arg) for arg in args), expected)
        test_strings.append(numerical_test_template % {
            "call": call_string,
            "threshold": threshold,
        })
    f = file("main.c", "w")
    f.write(main_template % "".join(test_strings))
    f.close()
    # 4) Compile and link
    compiled = try_run([
        "cc -c codegen.c -o codegen.o",
        "cc -c main.c -o main.o",
        "cc main.o codegen.o -lm -o test"
    ])
    # 5) Run if compiled
    if compiled:
        executed = try_run(["./test"])
    else:
        executed = False
    # 6) Clean up stuff
    clean = os.getenv('SYMPY_TEST_CLEAN_TEMP', 'always').lower()
    if clean not in ('always', 'success'):
        raise ValueError("SYMPY_TEST_CLEAN_TEMP must be one of the following: always or success.")
    if clean == 'always' or (clean == 'success' and compiled and executed):
        def safe_remove(filename):
            if os.path.isfile(filename):
                os.remove(filename)
        safe_remove("codegen.c")
        safe_remove("codegen.h")
        safe_remove("codegen.o")
        safe_remove("main.c")
        safe_remove("main.o")
        safe_remove("test")
        os.chdir(oldwork)
        os.rmdir(work)
    else:
        print >> sys.stderr, "TEST NOT REMOVED: %s" % work
        os.chdir(oldwork)
    # 7) Do the assertions in the end
    assert compiled
    assert executed


def is_feasible():
    # This test should always work, otherwise the cc compiler is not present.
    expr = (x+y)*z
    routine = Routine("test", [InputArgument(symbol) for symbol in x,y,z], [Result(expr)])
    numerical_tests = [
        ("test", (1.0, 6.0, 3.0), 21.0, 1e-15),
        ("test", (-1.0, 2.0, -2.5), -2.5, 1e-15),
    ]
    run_cc_test("is_feasible", [routine], numerical_tests, friendly=False)


try:
    is_feasible()
except AssertionError:
    disabled = True


def test_basic():
    is_feasible()


def test_basic_codegen():
    numerical_tests = [
        ("test", (1.0, 6.0, 3.0), 21.0, 1e-15),
        ("test", (-1.0, 2.0, -2.5), -2.5, 1e-15),
    ]
    run_cc_test("is_feasible", [("test", (x+y)*z)], numerical_tests)


