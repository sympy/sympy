from os import walk, sep, pardir
from os.path import split, join, abspath, exists, isfile
from glob import glob
import re
import random

from sympy.core.compatibility import PY3

# System path separator (usually slash or backslash) to be
# used with excluded files, e.g.
#     exclude = set([
#                    "%(sep)smpmath%(sep)s" % sepd,
#                   ])
sepd = {"sep": sep}

# path and sympy_path
SYMPY_PATH = abspath(join(split(__file__)[0], pardir, pardir))  # go to sympy/
assert exists(SYMPY_PATH)

TOP_PATH = abspath(join(SYMPY_PATH, pardir))
BIN_PATH = join(TOP_PATH, "bin")
EXAMPLES_PATH = join(TOP_PATH, "examples")

# Error messages
message_space = "File contains trailing whitespace: %s, line %s."
message_implicit = "File contains an implicit import: %s, line %s."
message_tabs = "File contains tabs instead of spaces: %s, line %s."
message_carriage = "File contains carriage returns at end of line: %s, line %s"
message_str_raise = "File contains string exception: %s, line %s"
message_gen_raise = "File contains generic exception: %s, line %s"
message_old_raise = "File contains old-style raise statement: %s, line %s, \"%s\""
message_eof = "File does not end with a newline: %s, line %s"
message_multi_eof = "File ends with more than 1 newline: %s, line %s"
message_test_suite_def = "Function should start with 'test_' or '_': %s, line %s"

implicit_test_re = re.compile(r'^\s*(>>> )?(\.\.\. )?from .* import .*\*')
str_raise_re = re.compile(
    r'^\s*(>>> )?(\.\.\. )?raise(\s+(\'|\")|\s*(\(\s*)+(\'|\"))')
gen_raise_re = re.compile(
    r'^\s*(>>> )?(\.\.\. )?raise(\s+Exception|\s*(\(\s*)+Exception)')
old_raise_re = re.compile(r'^\s*(>>> )?(\.\.\. )?raise((\s*\(\s*)|\s+)\w+\s*,')
test_suite_def_re = re.compile(r'^def\s+(?!(_|test))[^(]*\(\s*\)\s*:$')
test_file_re = re.compile(r'.*test_.*\.py$')

def tab_in_leading(s):
    """Returns True if there are tabs in the leading whitespace of a line,
    including the whitespace of docstring code samples."""
    n = len(s) - len(s.lstrip())
    if not s[n:n + 3] in ['...', '>>>']:
        check = s[:n]
    else:
        smore = s[n + 3:]
        check = s[:n] + smore[:len(smore) - len(smore.lstrip())]
    return not (check.expandtabs() == check)


def check_directory_tree(base_path, file_check, exclusions=set(), pattern="*.py"):
    """
    Checks all files in the directory tree (with base_path as starting point)
    with the file_check function provided, skipping files that contain
    any of the strings in the set provided by exclusions.
    """
    if not base_path:
        return
    for root, dirs, files in walk(base_path):
        check_files(glob(join(root, pattern)), file_check, exclusions)


def check_files(files, file_check, exclusions=set(), pattern=None):
    """
    Checks all files with the file_check function provided, skipping files
    that contain any of the strings in the set provided by exclusions.
    """
    if not files:
        return
    for fname in files:
        if not exists(fname) or not isfile(fname):
            continue
        if any(ex in fname for ex in exclusions):
            continue
        if pattern is None or re.match(pattern, fname):
            file_check(fname)


def test_files():
    """
    This test tests all files in sympy and checks that:
      o no lines contains a trailing whitespace
      o no lines end with \r\n
      o no line uses tabs instead of spaces
      o that the file ends with a single newline
      o there are no general or string exceptions
      o there are no old style raise statements
      o name of arg-less test suite functions start with _ or test_
    """

    def test(fname):
        if PY3:
            with open(fname, "rt", encoding="utf8") as test_file:
                test_this_file(fname, test_file)
        else:
            with open(fname, "rt") as test_file:
                test_this_file(fname, test_file)

    def test_this_file(fname, test_file):
        line = None  # to flag the case where there were no lines in file
        for idx, line in enumerate(test_file):
            if test_file_re.match(fname) and test_suite_def_re.match(line):
                assert False, message_test_suite_def % (fname, idx + 1)
            if line.endswith(" \n") or line.endswith("\t\n"):
                assert False, message_space % (fname, idx + 1)
            if line.endswith("\r\n"):
                assert False, message_carriage % (fname, idx + 1)
            if tab_in_leading(line):
                assert False, message_tabs % (fname, idx + 1)
            if str_raise_re.search(line):
                assert False, message_str_raise % (fname, idx + 1)
            if gen_raise_re.search(line):
                assert False, message_gen_raise % (fname, idx + 1)
            if (implicit_test_re.search(line) and
                    not filter(lambda ex: ex in fname, import_exclude)):
                assert False, message_implicit % (fname, idx + 1)

            result = old_raise_re.search(line)

            if result is not None:
                assert False, message_old_raise % (
                    fname, idx + 1, result.group(2))

        if line is not None:
            if line == '\n' and idx > 0:
                assert False, message_multi_eof % (fname, idx + 1)
            elif not line.endswith('\n'):
                # eof newline check
                assert False, message_eof % (fname, idx + 1)

    # Files to test at top level
    top_level_files = [join(TOP_PATH, file) for file in [
        "build.py",
        "setup.py",
        "setupegg.py",
    ]]
    # Files to exclude from all tests
    exclude = set([
        "%(sep)smpmath%(sep)s" % sepd,
    ])
    # Files to exclude from the implicit import test
    import_exclude = set([
        # glob imports are allowed in top-level __init__.py:
        "%(sep)ssympy%(sep)s__init__.py" % sepd,
        # these __init__.py should be fixed:
        # XXX: not really, they use useful import pattern (DRY)
        "%(sep)svector%(sep)s__init__.py" % sepd,
        "%(sep)smechanics%(sep)s__init__.py" % sepd,
        "%(sep)squantum%(sep)s__init__.py" % sepd,
        "%(sep)spolys%(sep)s__init__.py" % sepd,
        "%(sep)spolys%(sep)sdomains%(sep)s__init__.py" % sepd,
        # interactive sympy executes ``from sympy import *``:
        "%(sep)sinteractive%(sep)ssession.py" % sepd,
        # isympy executes ``from sympy import *``:
        "%(sep)sbin%(sep)sisympy" % sepd,
        # these two are import timing tests:
        "%(sep)sbin%(sep)ssympy_time.py" % sepd,
        "%(sep)sbin%(sep)ssympy_time_cache.py" % sepd,
        # Taken from Python stdlib:
        "%(sep)sparsing%(sep)ssympy_tokenize.py" % sepd,
        # these two should be fixed:
        "%(sep)splotting%(sep)spygletplot%(sep)s" % sepd,
        "%(sep)splotting%(sep)stextplot.py" % sepd,
    ])
    check_files(top_level_files, test)
    check_directory_tree(BIN_PATH, test, set(["~", ".pyc", ".sh"]), "*")
    check_directory_tree(SYMPY_PATH, test, exclude)
    check_directory_tree(EXAMPLES_PATH, test, exclude)


def _with_space(c):
    # return c with a random amount of leading space
    return random.randint(0, 10)*' ' + c


def test_raise_statement_regular_expression():
    candidates_ok = [
        "some text # raise Exception, 'text'",
        "raise ValueError('text') # raise Exception, 'text'",
        "raise ValueError('text')",
        "raise ValueError",
        "raise ValueError('text')",
        "raise ValueError('text') #,",
        # Talking about an exception in a docstring
        ''''"""This function will raise ValueError, except when it doesn't"""''',
        "raise (ValueError('text')",
    ]
    str_candidates_fail = [
        "raise 'exception'",
        "raise 'Exception'",
        'raise "exception"',
        'raise "Exception"',
        "raise 'ValueError'",
    ]
    gen_candidates_fail = [
        "raise Exception('text') # raise Exception, 'text'",
        "raise Exception('text')",
        "raise Exception",
        "raise Exception('text')",
        "raise Exception('text') #,",
        "raise Exception, 'text'",
        "raise Exception, 'text' # raise Exception('text')",
        "raise Exception, 'text' # raise Exception, 'text'",
        ">>> raise Exception, 'text'",
        ">>> raise Exception, 'text' # raise Exception('text')",
        ">>> raise Exception, 'text' # raise Exception, 'text'",
    ]
    old_candidates_fail = [
        "raise Exception, 'text'",
        "raise Exception, 'text' # raise Exception('text')",
        "raise Exception, 'text' # raise Exception, 'text'",
        ">>> raise Exception, 'text'",
        ">>> raise Exception, 'text' # raise Exception('text')",
        ">>> raise Exception, 'text' # raise Exception, 'text'",
        "raise ValueError, 'text'",
        "raise ValueError, 'text' # raise Exception('text')",
        "raise ValueError, 'text' # raise Exception, 'text'",
        ">>> raise ValueError, 'text'",
        ">>> raise ValueError, 'text' # raise Exception('text')",
        ">>> raise ValueError, 'text' # raise Exception, 'text'",
        "raise(ValueError,",
        "raise (ValueError,",
        "raise( ValueError,",
        "raise ( ValueError,",
        "raise(ValueError ,",
        "raise (ValueError ,",
        "raise( ValueError ,",
        "raise ( ValueError ,",
    ]

    for c in candidates_ok:
        assert str_raise_re.search(_with_space(c)) is None, c
        assert gen_raise_re.search(_with_space(c)) is None, c
        assert old_raise_re.search(_with_space(c)) is None, c
    for c in str_candidates_fail:
        assert str_raise_re.search(_with_space(c)) is not None, c
    for c in gen_candidates_fail:
        assert gen_raise_re.search(_with_space(c)) is not None, c
    for c in old_candidates_fail:
        assert old_raise_re.search(_with_space(c)) is not None, c


def test_implicit_imports_regular_expression():
    candidates_ok = [
        "from sympy import something",
        ">>> from sympy import something",
        "from sympy.somewhere import something",
        ">>> from sympy.somewhere import something",
        "import sympy",
        ">>> import sympy",
        "import sympy.something.something",
        "... import sympy",
        "... import sympy.something.something",
        "... from sympy import something",
        "... from sympy.somewhere import something",
        ">> from sympy import *",  # To allow 'fake' docstrings
        "# from sympy import *",
        "some text # from sympy import *",
    ]
    candidates_fail = [
        "from sympy import *",
        ">>> from sympy import *",
        "from sympy.somewhere import *",
        ">>> from sympy.somewhere import *",
        "... from sympy import *",
        "... from sympy.somewhere import *",
    ]
    for c in candidates_ok:
        assert implicit_test_re.search(_with_space(c)) is None, c
    for c in candidates_fail:
        assert implicit_test_re.search(_with_space(c)) is not None, c


def test_test_suite_defs():
    candidates_ok = [
        "    def foo():\n",
        "def foo(arg):\n",
        "def _foo():\n",
        "def test_foo():\n",
    ]
    candidates_fail = [
        "def foo():\n",
        "def foo() :\n",
        "def foo( ):\n",
        "def  foo():\n",
    ]
    for c in candidates_ok:
        assert test_suite_def_re.search(c) is None, c
    for c in candidates_fail:
        assert test_suite_def_re.search(c) is not None, c
