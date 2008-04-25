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

def test_no_trailing_whitespace():
    message = "File contains trailing whitespace: %s, line %s."
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
                    if line.endswith(" \n"):
                        assert False, message % (fname, idx+1)
            finally:
                file.close()
