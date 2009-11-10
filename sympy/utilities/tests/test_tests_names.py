from os import walk, sep, pardir
from os.path import split, join, abspath, basename
from glob import glob

def test_no_duplicate_names():
    message = "Test '%s' has a duplicate name, py.test will not run it!"
    base_path = split(__file__)[0]
    base_path += sep + pardir + sep + pardir # go to sympy/
    base_path = abspath(base_path)
    tests = {}
    for root, dirs, files in walk(base_path):
        for fname in glob(join(root, "test_*.py")):
            tname = basename(fname)
            if tname in tests:
                tests[tname].append(fname)
            else:
                tests[tname] = [fname]
    for fname in tests:
        if len(tests[fname]) > 1:
            print tests[fname]
            assert False, message % (fname)
