from os import walk, sep, chdir, pardir
from os.path import split, join, abspath, exists
from glob import glob
import re

# System path separator (usually slash or backslash)
sepd = {"sep": sep}

# Files having at least one of these in their path will be excluded
# Example line:
#    "%(sep)sprinting%(sep)spretty%(sep)s" % sepd,
EXCLUDE = set([
    "%(sep)sthirdparty%(sep)s" % sepd,
])

# Error messages
message_space = "File contains trailing whitespace: %s, line %s."
message_implicit = "File contains an implicit import: %s, line %s."
message_tabs = "File contains tabs instead of spaces: %s, line %s."
message_carriage = "File contains carriage returns at end of line: %s, line %s"

def get_whitespace(s):
    """Returns all whitespace at the beginning of the line."""
    for i in range(len(s)):
        if not s[i].isspace():
            return s[:i]
    return ""

def check_directory_tree(base_path):
    """
    Checks all files in the directory tree (with base_path as starting point)
    for bad coding style.
    """
    for root, dirs, files in walk(base_path):
        for fname in glob(join(root, "*.py")):
            if filter(lambda ex: ex in fname, EXCLUDE):
                continue
            file = open(fname, "r")
            try:
                for idx, line in enumerate(file):
                    #if line.endswith(" \n"):
                    if line.rstrip()+"\n" != line and line.rstrip() != line:
                        assert False, message_space % (fname, idx+1)
                    if "\r\n" in line:
                        assert False, message_carriagereturn % (fname, idx+1)
                    w = get_whitespace(line)
                    if w.expandtabs() != w:
                        assert False, message_tabs % (fname, idx+1)
            finally:
                file.close()


def test_no_trailing_whitespace_and_no_tabs():
    """
    This test tests all files in sympy and makes sure no lines contains a
    trailing whitespace at the end, and also that no line uses tabs instead of
    spaces as the indentation.
    """
    path = split(__file__)[0]
    path = path + sep + pardir + sep + pardir # go to sympy/
    sympy_path = abspath(path)
    examples_path = abspath(join(path + sep + pardir, "examples"))
    paths = (sympy_path, examples_path)
    for p in paths:
        assert exists(p)
        check_directory_tree(p)

def test_implicit_imports():
    """
    Tests that all files except __init__.py use explicit imports.
    """
    path = split(abspath(__file__))[0]
    path = path + sep + pardir + sep + pardir # go to sympy/
    sympy_path = abspath(path)
    examples_path = abspath(join(path + sep + pardir, "examples"))
    paths = (sympy_path, examples_path)
    exclude = set([
        "%(sep)sthirdparty%(sep)s" % sepd,
        "%(sep)s__init__.py" % sepd,
        # these two should be fixed:
        "%(sep)smpmath%(sep)s" % sepd,
        "%(sep)splotting%(sep)s" % sepd,

        # everything below should be fixed soon:
        "sympy/core/tests/test_evalf.py" % sepd,
        "sympy/functions/combinatorial/numbers.py" % sepd,
        "sympy/physics/tests/test_units.py" % sepd,
        "sympy/polynomials/factor_.py" % sepd,
        "sympy/polynomials/groebner_.py" % sepd,
        "sympy/polynomials/wrapper.py" % sepd,
        "sympy/polynomials/base.py" % sepd,
        "sympy/polynomials/div_.py" % sepd,
        "sympy/polynomials/roots_.py" % sepd,
        "sympy/polys/tests/test_polynomial.py" % sepd,
        "sympy/solvers/tests/test_solvers.py" % sepd,
        "sympy/statistics/distributions.py" % sepd,
        "sympy/examples/fem_test.py" % sepd,
        "sympy/examples/pidigits.py" % sepd,
        "sympy/examples/print_gtk.py" % sepd,
        "sympy/examples/print_pretty.py" % sepd,
        "sympy/examples/fem.py" % sepd,
        "sympy/examples/trees.py" % sepd,
        "sympy/galgebra/GAcalc.py" % sepd,
        "sympy/galgebra/GAsympy.py" % sepd,
        "sympy/galgebra/latex_out.py" % sepd,
    ])
    for p in paths:
        assert exists(p)
        for root, dirs, files in walk(p):
            for fname in glob(join(root, "*.py")):
                if filter(lambda ex: ex in fname, exclude):
                    continue
                file = open(fname, "r")
                try:
                    for idx, line in enumerate(file):
                        if re.match("^\s*from.*import.*\*",line):
                            assert False, message_implicit % (fname, idx+1)
                finally:
                    file.close()
