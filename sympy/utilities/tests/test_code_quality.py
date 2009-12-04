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
message_str_raise = "File contains string exception: %s, line %s"
message_gen_raise = "File contains generic exception: %s, line %s"
message_eof = "File does not end with a newline: %s, line %s"

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
    strRaise = re.compile(r'raise(\s+(\'|\")|\s*(\(\s*)+(\'|\"))')
    genRaise = re.compile(r'raise(\s+Exception|\s*(\(\s*)+Exception)')
    for root, dirs, files in walk(base_path):
        for fname in glob(join(root, "*.py")):
            if filter(lambda ex: ex in fname, EXCLUDE):
                continue
            file = open(fname, "r")
            try:
                for idx, line in enumerate(file):
                    if line.endswith(" \n"):
                        assert False, message_space % (fname, idx+1)
                    if line.endswith("\r\n"):
                        assert False, message_carriage % (fname, idx+1)
                    w = get_whitespace(line)
                    if w.expandtabs() != w:
                        assert False, message_tabs % (fname, idx+1)
                    if strRaise.search(line):
                        assert False, message_str_raise % (fname, idx+1)
                    if genRaise.search(line):
                        assert False, message_gen_raise % (fname, idx+1)
            finally:
                # eof newline check
                if not line.endswith('\n'):
                    assert False, message_eof % (fname, idx+1)
                file.close()


def test_no_trailing_whitespace_and_no_tabs():
    """
    This test tests all files in sympy and checks that:
      o no lines contains a trailing whitespace
      o no line uses tabs instead of spaces as the indentation
      o that the file ends with a newline
    """
    path = split(__file__)[0]
    path = path + sep + pardir + sep + pardir # go to sympy/
    sympy_path = abspath(path)
    assert exists(sympy_path)
    check_directory_tree(sympy_path)
    examples_path = abspath(join(path + sep + pardir, "examples"))
    # Remember that this test can be executed when examples are not installed
    # (e.g. after "./setup.py install"), so don't raise an error if the
    # examples are not found:
    if exists(examples_path):
        check_directory_tree(examples_path)

def check_directory_tree_imports(p, exclude):
    for root, dirs, files in walk(p):
        for fname in glob(join(root, "*.py")):
            if filter(lambda ex: ex in fname, exclude):
                continue
            file = open(fname, "r")
            try:
                for idx, line in enumerate(file):
                    if re.match("^\s*(>>>)? from .* import .*\*",line):
                        assert False, message_implicit % (fname, idx+1)
            finally:
                file.close()

def test_implicit_imports():
    """
    Tests that all files except __init__.py use explicit imports,
    even in the docstrings.
    """
    path = split(abspath(__file__))[0]
    path = path + sep + pardir + sep + pardir # go to sympy/
    sympy_path = abspath(path)
    examples_path = abspath(join(path + sep + pardir, "examples"))
    exclude = set([
        "%(sep)sthirdparty%(sep)s" % sepd,
        "%(sep)s__init__.py" % sepd,
        # these two should be fixed:
        "%(sep)smpmath%(sep)s" % sepd,
        "%(sep)splotting%(sep)s" % sepd,
    ])
    assert exists(sympy_path)
    check_directory_tree_imports(sympy_path, exclude)
    # Remember that this test can be executed when examples are not installed
    # (e.g. after "./setup.py install"), so don't raise an error if the
    # examples are not found:
    if exists(examples_path):
        check_directory_tree_imports(examples_path, exclude)
