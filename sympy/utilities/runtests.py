"""
This is our testing framework.

Goals:

* it should be compatible with py.test and operate very similarly (or
identically)
* doesn't require any external dependencies
* preferably all the functionality should be in this file only
* no magic, just import the test file and execute the test functions, that's it
* portable

"""
import os
import sys
import inspect
import traceback
import pdb
from glob import glob
from timeit import default_timer as clock
import doctest as pdoctest # avoid clashing with our doctest() function
from sympy.utilities import any

def isgeneratorfunction(object):
    """
    Return true if the object is a user-defined generator function.

    Generator function objects provides same attributes as functions.

    See isfunction.__doc__ for attributes listing.

    Adapted from Python 2.6.
    """
    CO_GENERATOR = 0x20
    if (inspect.isfunction(object) or inspect.ismethod(object)) and \
        object.func_code.co_flags & CO_GENERATOR:
        return True
    return False

def setup_pprint():
    from sympy import pprint_use_unicode
    # force pprint to be in ascii mode in doctests
    pprint_use_unicode(False)

    # hook our nice, hash-stable strprinter
    from sympy.interactive import init_printing
    from sympy.printing import sstrrepr
    init_printing(sstrrepr)

def convert_to_native_paths(lst):
    """
    Converts a list of '/' separated paths into a list of
    native (os.sep separated) paths.
    """
    return [os.path.join(*x.split("/")) for x in lst]

def test(*args, **kwargs):
    """
    Run all tests containing any of the given strings in their path.

    If sort=False, run them in random order (not default).

    Warning: Tests in *very* deeply nested directories are not found.

    Examples:

    >> import sympy

    Run all tests:
    >> sympy.test()

    Run one file:
    >> sympy.test("sympy/core/tests/test_basic.py")

    Run all tests in sympy/functions/ and some particular file:
    >> sympy.test("sympy/core/tests/test_basic.py", "sympy/functions")

    Run all tests in sympy/core and sympy/utilities:
    >> sympy.test("core", "util")
    """
    from glob import glob
    verbose = kwargs.get("verbose", False)
    tb = kwargs.get("tb", "short")
    kw = kwargs.get("kw", "")
    post_mortem = kwargs.get("pdb", False)
    colors = kwargs.get("colors", True)
    sort = kwargs.get("sort", True)
    r = PyTestReporter(verbose, tb, colors)
    t = SymPyTests(r, kw, post_mortem)
    if len(args) == 0:
        t.add_paths(["sympy"])
    else:
        mypaths = []
        for p in t.get_paths(dir='sympy'):
            mypaths.extend(glob(p))
        mypaths = set(mypaths)
        t.add_paths([p for p in mypaths if any(a in p for a in args)])
    return t.test(sort=sort)

def doctest(*paths, **kwargs):
    """
    Runs the doctests specified by paths, or all tests if paths=[].

    Note: paths are specified relative to the sympy root directory in a unix
    format (on all platforms including windows).

    Examples:

    Run all tests:
    >> import sympy
    >> sympy.doctest()

    Run one file:
    >> import sympy
    >> sympy.doctest("sympy/core/tests/test_basic.py")

    Run all tests in sympy/functions/ and some particular file:
    >> import sympy
    >> sympy.doctest("sympy/core/tests/test_basic.py", "sympy/functions")
    """
    verbose = kwargs.get("verbose", False)
    blacklist = kwargs.get("blacklist", [])
    blacklist.extend([
        "sympy/thirdparty/pyglet", # segfaults
        "sympy/mpmath", # needs to be fixed upstream
        "sympy/plotting", # generates live plots
        "sympy/utilities/compilef.py", # needs tcc
        "sympy/galgebra/GA.py", # needs numpy
        "sympy/galgebra/latex_ex.py", # needs numpy
        "sympy/conftest.py", # needs py.test
        "sympy/utilities/benchmarking.py", # needs py.test
        ])

    r = PyTestReporter(verbose)
    t = SymPyDocTests(r, blacklist=blacklist)
    if len(paths) > 0:
        t.add_paths(paths)
    else:
        t.add_paths(["sympy"])
    dtest = t.test()
    if not dtest:
        return False

    if sys.version_info[:2] <= (2,4):
        return True

    if len(paths) == 0:
        # test documentation under doc/src/ only if we are running the full
        # test suite and this is python2.5 or newer:
        excluded = convert_to_native_paths(['doc/src/modules/plotting.txt'])
        doc_globs = convert_to_native_paths(['doc/src/*.txt',
                'doc/src/modules/*.txt'])
        doc_files = sum([glob(x) for x in doc_globs], [])
        for ex in excluded:
            doc_files.remove(ex)
    else:
        doc_files = paths

    setup_pprint()
    doc_tests_succeeded = True
    for doc_file in doc_files:
        if not os.path.isfile(doc_file):
            continue
        old_displayhook = sys.displayhook
        try:
            out = pdoctest.testfile(doc_file, module_relative=False,
                    optionflags=pdoctest.ELLIPSIS | \
                    pdoctest.NORMALIZE_WHITESPACE)
        finally:
            # make sure we return to the original displayhook in case some
            # doctest has changed that
            sys.displayhook = old_displayhook
        print "Testing ", doc_file
        print "Failed %s, tested %s" % out
        if out[0] != 0:
            doc_tests_succeeded = False
    if not doc_tests_succeeded:
        print("DO *NOT* COMMIT!")
    return doc_tests_succeeded

class SymPyTests(object):

    def __init__(self, reporter, kw="", post_mortem=False):
        self._post_mortem = post_mortem
        self._kw = kw
        self._count = 0
        self._root_dir = self.get_sympy_dir()
        self._reporter = reporter
        self._reporter.root_dir(self._root_dir)
        self._tests = []

    def add_paths(self, paths):
        for path in paths:
            path2 = os.path.join(self._root_dir, *path.split("/"))
            if path2.endswith(".py"):
                self._tests.append(path2)
            else:
                self._tests.extend(self.get_tests(path2))

    def test(self, sort=False):
        """
        Runs the tests.

        If sort=False run tests in random order.

        Returns True if all tests pass, otherwise False.
        """
        if sort:
            self._tests.sort()
        else:
            from random import shuffle
            shuffle(self._tests)
        self._reporter.start()
        for f in self._tests:
            try:
                self.test_file(f)
            except KeyboardInterrupt:
                print " interrupted by user"
                break
        return self._reporter.finish()

    def test_file(self, filename):
        name = "test%d" % self._count
        name = os.path.splitext(os.path.basename(filename))[0]
        self._count += 1
        gl = {'__file__':filename}
        try:
            execfile(filename, gl)
        except (ImportError, SyntaxError):
            self._reporter.import_error(filename, sys.exc_info())
            return
        pytestfile = ""
        if "XFAIL" in gl:
            pytestfile = inspect.getsourcefile(gl["XFAIL"])
        disabled = gl.get("disabled", False)
        if disabled:
            funcs = []
        else:
            # we need to filter only those functions that begin with 'test_'
            # that are defined in the testing file or in the file where
            # is defined the XFAIL decorator
            funcs = [gl[f] for f in gl.keys() if f.startswith("test_") and
                                                 (inspect.isfunction(gl[f])
                                                    or inspect.ismethod(gl[f])) and
                                                 (inspect.getsourcefile(gl[f]) == filename or
                                                   inspect.getsourcefile(gl[f]) == pytestfile)]
            # Sorting of XFAILed functions isn't fixed yet :-(
            funcs.sort(key=lambda x: inspect.getsourcelines(x)[1])
            i = 0
            while i is not len(funcs):
                if isgeneratorfunction(funcs[i]):
                # some tests can be generators, that return the actual
                # test functions. We unpack it below:
                    f = funcs.pop(i)
                    for fg in f():
                        func = fg[0]
                        args = fg[1:]
                        fgw = lambda: func(*args)
                        funcs.insert(i, fgw)
                        i += 1
                else:
                    i += 1
            # drop functions that are not selected with the keyword expression:
            funcs = [x for x in funcs if self.matches(x)]

        self._reporter.entering_filename(filename, len(funcs))
        for f in funcs:
            self._reporter.entering_test(f)
            try:
                f()
            except KeyboardInterrupt:
                raise
            except:
                t, v, tr = sys.exc_info()
                if t is AssertionError:
                    self._reporter.test_fail((t, v, tr))
                    if self._post_mortem:
                        pdb.post_mortem(tr)
                elif t.__name__ == "Skipped":
                    self._reporter.test_skip()
                elif t.__name__ == "XFail":
                    self._reporter.test_xfail()
                elif t.__name__ == "XPass":
                    self._reporter.test_xpass(v)
                else:
                    self._reporter.test_exception((t, v, tr))
                    if self._post_mortem:
                        pdb.post_mortem(tr)
            else:
                self._reporter.test_pass()
        self._reporter.leaving_filename()

    def get_sympy_dir(self):
        """
        Returns the root sympy directory.
        """
        this_file = os.path.abspath(__file__)
        sympy_dir = os.path.join(os.path.dirname(this_file), "..", "..")
        sympy_dir = os.path.normpath(sympy_dir)
        return sympy_dir

    def matches(self, x):
        """
        Does the keyword expression self._kw match "x"? Returns True/False.

        Always returns True if self._kw is "".
        """
        if self._kw == "":
            return True
        return x.__name__.find(self._kw) != -1

    def get_paths(self, dir="", level=15):
        """
        Generates a set of paths for testfiles searching.

        Example:
        >> get_paths(2)
        ['sympy/test_*.py', 'sympy/*/test_*.py', 'sympy/*/*/test_*.py']
        >> get_paths(6)
        ['sympy/test_*.py', 'sympy/*/test_*.py', 'sympy/*/*/test_*.py',
        'sympy/*/*/*/test_*.py', 'sympy/*/*/*/*/test_*.py',
        'sympy/*/*/*/*/*/test_*.py', 'sympy/*/*/*/*/*/*/test_*.py']
        """
        wildcards = [dir]
        for i in range(level):
            wildcards.append(os.path.join(wildcards[-1], "*"))
        p = [os.path.join(x, "test_*.py") for x in wildcards]
        return p

    def get_tests(self, dir):
        """
        Returns the list of tests.
        """
        g = []
        for x in self.get_paths(dir):
            g.extend(glob(x))
        g = list(set(g))
        g.sort()
        return g

class SymPyDocTests(object):

    def __init__(self, reporter, blacklist=[]):
        self._count = 0
        self._root_dir = self.get_sympy_dir()
        self._reporter = reporter
        self._reporter.root_dir(self._root_dir)
        self._tests = []
        self._blacklist = convert_to_native_paths(blacklist)

    def add_paths(self, paths):
        for path in paths:
            path2 = os.path.join(self._root_dir, *path.split("/"))
            if path2.endswith(".py"):
                self._tests.append(path2)
            else:
                self._tests.extend(self.get_tests(path2))

    def test(self):
        """
        Runs the tests.

        Returns True if all tests pass, otherwise False.
        """
        self._reporter.start()
        for f in self._tests:
            try:
                self.test_file(f)
            except KeyboardInterrupt:
                print " interrupted by user"
                break
        return self._reporter.finish()

    def test_file(self, filename):

        import doctest
        import unittest
        from StringIO import StringIO

        rel_name = filename[len(self._root_dir)+1:]
        module = rel_name.replace(os.sep, '.')[:-3]
        setup_pprint()
        try:
            module = doctest._normalize_module(module)
            tests = doctest.DocTestFinder().find(module)
        except:
            self._reporter.import_error(filename, sys.exc_info())
            return

        tests.sort()
        tests = [test for test in tests if len(test.examples) > 0]
        self._reporter.entering_filename(filename, len(tests))
        for test in tests:
            assert len(test.examples) != 0
            runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS | \
                    doctest.NORMALIZE_WHITESPACE)
            old = sys.stdout
            new = StringIO()
            sys.stdout = new
            try:
                f, t = runner.run(test, out=new.write, clear_globs=False)
            finally:
                sys.stdout = old
            if f > 0:
                self._reporter.doctest_fail(test.name, new.getvalue())
            else:
                self._reporter.test_pass()
        self._reporter.leaving_filename()

    def get_sympy_dir(self):
        """
        Returns the root sympy directory.
        """
        this_file = os.path.abspath(__file__)
        sympy_dir = os.path.join(os.path.dirname(this_file), "..", "..")
        sympy_dir = os.path.normpath(sympy_dir)
        return sympy_dir


    def get_paths(self, dir="", level=15):
        """
        Generates a set of paths for testfiles searching.

        Example:
        >> get_paths(2)
        ['sympy/test_*.py', 'sympy/*/test_*.py', 'sympy/*/*/test_*.py']
        >> get_paths(6)
        ['sympy/test_*.py', 'sympy/*/test_*.py', 'sympy/*/*/test_*.py',
        'sympy/*/*/*/test_*.py', 'sympy/*/*/*/*/test_*.py',
        'sympy/*/*/*/*/*/test_*.py', 'sympy/*/*/*/*/*/*/test_*.py']
        """
        wildcards = [dir]
        for i in range(level):
            wildcards.append(os.path.join(wildcards[-1], "*"))
        p = [os.path.join(x, "*.py") for x in wildcards]
        return p

    def is_on_blacklist(self, x):
        """
        Returns True if "x" is on the blacklist. Otherwise False.
        """
        for p in self._blacklist:
            if x.find(p) != -1:
                return True
        return False

    def get_tests(self, dir):
        """
        Returns the list of tests.
        """
        def importable(x):
            """
            Checks if given pathname x is an importable module by checking for
            __init__.py file.

            Returns True/False.

            Currently we only test if the __init__.py file exists in the
            directory with the file "x" (in theory we should also test all the
            parent dirs) and if "x" is not on self._blacklist.
            """
            if self.is_on_blacklist(x):
                return False
            init_py = os.path.dirname(x) + os.path.sep + "__init__.py"
            return os.path.exists(init_py)

        g = []
        for x in self.get_paths(dir):
            g.extend(glob(x))
        g = list(set(g))
        g.sort()
        # skip files that are not importable (i.e. missing __init__.py)
        g = [x for x in g if importable(x)]
        return g

class Reporter(object):
    """
    Parent class for all reporters.
    """
    pass

class PyTestReporter(Reporter):
    """
    Py.test like reporter. Should produce output identical to py.test.
    """

    def __init__(self, verbose=False, tb="short", colors=True):
        self._verbose = verbose
        self._tb_style = tb
        self._colors = colors
        self._xfailed = 0
        self._xpassed = []
        self._failed = []
        self._failed_doctest = []
        self._passed = 0
        self._skipped = 0
        self._exceptions = []

        # this tracks the x-position of the cursor (useful for positioning
        # things on the screen), without the need for any readline library:
        self._write_pos = 0
        self._line_wrap = False

    def root_dir(self, dir):
        self._root_dir = dir

    def write(self, text, color="", align="left", width=80):
        """
        Prints a text on the screen.

        It uses sys.stdout.write(), so no readline library is necessary.

        color ... choose from the colors below, "" means default color
        align ... left/right, left is a normal print, right is aligned on the
                  right hand side of the screen, filled with " " if necessary
        width ... the screen width
        """
        color_templates = (
            ("Black"       , "0;30"),
            ("Red"         , "0;31"),
            ("Green"       , "0;32"),
            ("Brown"       , "0;33"),
            ("Blue"        , "0;34"),
            ("Purple"      , "0;35"),
            ("Cyan"        , "0;36"),
            ("LightGray"   , "0;37"),
            ("DarkGray"    , "1;30"),
            ("LightRed"    , "1;31"),
            ("LightGreen"  , "1;32"),
            ("Yellow"      , "1;33"),
            ("LightBlue"   , "1;34"),
            ("LightPurple" , "1;35"),
            ("LightCyan"   , "1;36"),
            ("White"       , "1;37"),  )

        colors = {}

        for name, value in color_templates:
            colors[name] = value
        c_normal = '\033[0m'
        c_color = '\033[%sm'

        if align == "right":
            if self._write_pos+len(text) > width:
                # we don't fit on the current line, create a new line
                self.write("\n")
            self.write(" "*(width-self._write_pos-len(text)))

        if not sys.stdout.isatty():
            # the stdout is not a terminal, this for example happens if the
            # output is piped to less, e.g. "bin/test | less". In this case,
            # the terminal control sequences would be printed verbatim, so
            # don't use any colors.
            color = ""
        if sys.platform == "win32":
            # Windows consoles don't support ANSI escape sequences
            color = ""

        if self._line_wrap:
            if text[0] != "\n":
                sys.stdout.write("\n")

        if color == "":
            sys.stdout.write(text)
        else:
            sys.stdout.write("%s%s%s" % (c_color % colors[color], text, c_normal))
        sys.stdout.flush()
        l = text.rfind("\n")
        if l == -1:
            self._write_pos += len(text)
        else:
            self._write_pos = len(text)-l-1
        self._line_wrap = self._write_pos >= width
        self._write_pos %= width

    def write_center(self, text, delim="="):
        width = 80
        if text != "":
            text = " %s " % text
        idx = (width-len(text)) // 2
        t = delim*idx + text + delim*(width-idx-len(text))
        self.write(t+"\n")

    def write_exception(self, e, val, tb):
        t = traceback.extract_tb(tb)
        # remove the first item, as that is always runtests.py
        t = t[1:]
        t = traceback.format_list(t)
        self.write("".join(t))
        t = traceback.format_exception_only(e, val)
        self.write("".join(t))

    def start(self):
        self.write_center("test process starts")
        executable = sys.executable
        v = sys.version_info
        python_version = "%s.%s.%s-%s-%s" % v
        self.write("executable:   %s  (%s)\n\n" % (executable, python_version))
        self._t_start = clock()

    def finish(self):
        self._t_end = clock()
        self.write("\n")
        text = "tests finished: %d passed" % self._passed
        if len(self._failed) > 0:
            text += ", %d failed" % len(self._failed)
        if len(self._failed_doctest) > 0:
            text += ", %d failed" % len(self._failed_doctest)
        if self._skipped > 0:
            text += ", %d skipped" % self._skipped
        if self._xfailed > 0:
            text += ", %d xfailed" % self._xfailed
        if len(self._xpassed) > 0:
            text += ", %d xpassed" % len(self._xpassed)
        if len(self._exceptions) > 0:
            text += ", %d exceptions" % len(self._exceptions)
        text += " in %.2f seconds" % (self._t_end - self._t_start)


        if len(self._xpassed) > 0:
            self.write_center("xpassed tests", "_")
            for e in self._xpassed:
                self.write("%s:%s\n" % (e[0], e[1]))
            self.write("\n")

        if self._tb_style != "no" and len(self._exceptions) > 0:
            #self.write_center("These tests raised an exception", "_")
            for e in self._exceptions:
                filename, f, (t, val, tb) = e
                self.write_center("", "_")
                if f is None:
                    s = "%s" % filename
                else:
                    s = "%s:%s" % (filename, f.__name__)
                self.write_center(s, "_")
                self.write_exception(t, val, tb)
            self.write("\n")

        if self._tb_style != "no" and len(self._failed) > 0:
            #self.write_center("Failed", "_")
            for e in self._failed:
                filename, f, (t, val, tb) = e
                self.write_center("", "_")
                self.write_center("%s:%s" % (filename, f.__name__), "_")
                self.write_exception(t, val, tb)
            self.write("\n")

        if self._tb_style != "no" and len(self._failed_doctest) > 0:
            #self.write_center("Failed", "_")
            for e in self._failed_doctest:
                filename, msg = e
                self.write_center("", "_")
                self.write_center("%s" % filename, "_")
                self.write(msg)
            self.write("\n")

        self.write_center(text)
        ok = len(self._failed) == 0 and len(self._exceptions) == 0 and \
                len(self._failed_doctest) == 0
        if not ok:
            self.write("DO *NOT* COMMIT!\n")
        return ok

    def entering_filename(self, filename, n):
        rel_name = filename[len(self._root_dir)+1:]
        self._active_file = rel_name
        self._active_file_error = False
        self.write(rel_name)
        self.write("[%d] " % n)

    def leaving_filename(self):
        if self._colors:
            self.write(" ")
            if self._active_file_error:
                self.write("[FAIL]", "Red", align="right")
            else:
                self.write("[OK]", "Green", align="right")
        self.write("\n")
        if self._verbose:
            self.write("\n")

    def entering_test(self, f):
        self._active_f = f
        if self._verbose:
            self.write("\n"+f.__name__+" ")

    def test_xfail(self):
        self._xfailed += 1
        self.write("f")

    def test_xpass(self, fname):
        self._xpassed.append((self._active_file, fname))
        self.write("X")

    def test_fail(self, exc_info):
        self._failed.append((self._active_file, self._active_f, exc_info))
        self.write("F")
        self._active_file_error = True

    def doctest_fail(self, name, error_msg):
        # the first line contains "******", remove it:
        error_msg = "\n".join(error_msg.split("\n")[1:])
        self._failed_doctest.append((name, error_msg))
        self.write("F")
        self._active_file_error = True

    def test_pass(self):
        self._passed += 1
        if self._verbose:
            self.write("ok")
        else:
            self.write(".")

    def test_skip(self):
        self._skipped += 1
        self.write("s")

    def test_exception(self, exc_info):
        self._exceptions.append((self._active_file, self._active_f, exc_info))
        self.write("E")
        self._active_file_error = True

    def import_error(self, filename, exc_info):
        self._exceptions.append((filename, None, exc_info))
        rel_name = filename[len(self._root_dir)+1:]
        self.write(rel_name)
        self.write("[?]   Failed to import")
        if self._colors:
            self.write(" ")
            self.write("[FAIL]", "Red", align="right")
        self.write("\n")
