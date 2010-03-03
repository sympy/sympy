from os import walk, sep, chdir, pardir
from os.path import split, join, abspath, exists
from glob import glob
import re

# System path separator (usually slash or backslash) to be
# used with excluded files, e.g.
#     exclude = set([
#                    "%(sep)sthirdparty%(sep)s" % sepd,
#                   ])
sepd = {"sep": sep}

# path and sympy_path
PATH = reduce(join, [split(__file__)[0], pardir, pardir]) # go to sympy/
SYMPY_PATH = abspath(PATH)
assert exists(SYMPY_PATH)

# Tests can be executed when examples are not installed
# (e.g. after "./setup.py install") so set the examples path
# to null so it will be skipped by the checker if it is not
# there.
EXAMPLES_PATH = abspath(reduce(join, [PATH, pardir, "examples"]))
if not exists(EXAMPLES_PATH):
    EXAMPLES_PATH = ""

# Error messages
message_space = "File contains trailing whitespace: %s, line %s."
message_implicit = "File contains an implicit import: %s, line %s."
message_tabs = "File contains tabs instead of spaces: %s, line %s."
message_carriage = "File contains carriage returns at end of line: %s, line %s"
message_str_raise = "File contains string exception: %s, line %s"
message_gen_raise = "File contains generic exception: %s, line %s"
message_eof = "File does not end with a newline: %s, line %s"

def tab_in_leading(s):
    """Returns True if there are tabs in the leading whitespace of a line,
    including the whitespace of docstring code samples."""
    n = len(s)-len(s.lstrip())
    if not s[n:n+3] in ['...', '>>>']:
        check = s[:n]
    else:
        smore = s[n+3:]
        check = s[:n] + smore[:len(smore)-len(smore.lstrip())]
    return not (check.expandtabs() == check)

def check_directory_tree(base_path, file_check, exclusions=set()):
    """
    Checks all files in the directory tree (with base_path as starting point)
    with the file_check function provided, skipping files that contain
    any of the strings in the set provided by exclusions.
    """
    if not base_path:
        return
    for root, dirs, files in walk(base_path):
        for fname in glob(join(root, "*.py")):
            if filter(lambda ex: ex in fname, exclusions):
                continue
            file_check(fname)

def test_whitespace_and_exceptions():
    """
    This test tests all files in sympy and checks that:
      o no lines contains a trailing whitespace
      o no lines end with \r\n
      o no line uses tabs instead of spaces
      o that the file ends with a newline
      o there are no general or string exceptions
    """
    strRaise = re.compile(r'raise(\s+(\'|\")|\s*(\(\s*)+(\'|\"))')
    genRaise = re.compile(r'raise(\s+Exception|\s*(\(\s*)+Exception)')

    def test(fname):
        file = open(fname, "rb") # without "b" the lines from all systems will appear to be \n terminated
        try:
            line = None # to flag the case where there were no lines in file
            for idx, line in enumerate(file):
                if line.endswith(" \n"):
                    assert False, message_space % (fname, idx+1)
                if line.endswith("\r\n"):
                    assert False, message_carriage % (fname, idx+1)
                if tab_in_leading(line):
                    assert False, message_tabs % (fname, idx+1)
                if strRaise.search(line):
                    assert False, message_str_raise % (fname, idx+1)
                if genRaise.search(line):
                    assert False, message_gen_raise % (fname, idx+1)
        finally:
            if line != None:
                # eof newline check
                if not line.endswith('\n'):
                    assert False, message_eof % (fname, idx+1)
            file.close()

    exclude = set([
        "%(sep)sthirdparty%(sep)s" % sepd,
    ])
    check_directory_tree(SYMPY_PATH, test, exclude)
    check_directory_tree(EXAMPLES_PATH, test, exclude)

def test_implicit_imports():
    """
    Tests that all files except __init__.py use explicit imports,
    even in the docstrings.
    """
    def test(fname):
        file = open(fname, "r")
        try:
            for idx, line in enumerate(file):
                if re.match("^\s*(>>>)? from .* import .*\*",line):
                    assert False, message_implicit % (fname, idx+1)
        finally:
            file.close()

    exclude = set([
        "%(sep)sthirdparty%(sep)s" % sepd,
        "%(sep)s__init__.py" % sepd,
        # these two should be fixed:
        "%(sep)smpmath%(sep)s" % sepd,
        "%(sep)splotting%(sep)s" % sepd,
    ])
    check_directory_tree(SYMPY_PATH, test, exclude)
    check_directory_tree(EXAMPLES_PATH, test, exclude)
