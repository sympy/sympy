from os import walk, sep, chdir, pardir
from os.path import split, join, abspath
from glob import glob

# System path separator (usually slash or backslash)
sepd = {"sep": sep}

# Files having at least one of these in their path will be excluded
# Example line:
#    "%(sep)sprinting%(sep)spretty%(sep)s" % sepd,
EXCLUDE = set([
    "%(sep)sthirdparty%(sep)s" % sepd,
])

def get_whitespace(s):
    """Returns all whitespace at the beginning of the line."""
    for i in range(len(s)):
        if not s[i].isspace():
            return s[:i]
    return ""


def test_no_trailing_whitespace_and_no_tabs():
    """
    This test tests all files in sympy and makes sure no lines contains a
    trailing whitespace at the end, and also that no line uses tabs instead of
    spaces as the indentation.
    """
    message_space = "File contains trailing whitespace: %s, line %s."
    message_tabs = "File contains tabs instead of spaces: %s, line %s."
    base_path = split(__file__)[0]
    base_path += sep + pardir + sep + pardir # go to sympy/
    base_path = abspath(base_path)
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
                    w = get_whitespace(line)
                    if w.expandtabs() != w:
                        assert False, message_tabs % (fname, idx+1)
            finally:
                file.close()
