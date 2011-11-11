"""
This is our testing framework.

Goals:
======

* it should be compatible with py.test and operate very similarly (or
identically)
* doesn't require any external dependencies
* preferably all the functionality should be in this file only
* no magic, just import the test file and execute the test functions, that's it
* portable


The code structure.
===================

Modules loading:
----------------

* Load standard modules
* Tune standard doctest module (utf8 encoding, indent)

Entering points:
----------------

* test()      - Run all tests in `test_*.py` files, test code quality and so on.
* doctest()   - Run all doctests in `*.py` files, and run tests in `doc/src/*/*.txt` files.
* wikitest()  - Run test in wiki pages.

Tester classes:
---------------

* SymPyTesterBase                 - base tester class
* SymPyTester(SymPyTesterBase)    - find and test 'test_*.py' files
* SymPyDocTester(SymPyTesterBase) - find and test docstrings and *.txt files
* SymPyWikiTester(SymPyDocTester) - find and test wiki pages


Finder, Parser, Runner classes:
-------------------------------
* SymPyTextTestParser(DocTestParser) - Parser to find tests in the *txt files
    and in the wiki page
    It extends "get_doctest" to group examples to blocks (tests)

* SymPyDocTestFinder(DocTestFinder)

* SymPyDocTestRunner(DocTestRunner)
    Modified from the doctest version to not reset the sys.displayhook

Reporter Class:
---------------

* PyTestReporter

"""

#
# Load standard modules
#
import os
import sys
import inspect
import traceback
import pdb
import re
import linecache
from fnmatch import fnmatch
from timeit import default_timer as clock
import doctest as pdoctest # avoid clashing with our doctest() function
from doctest import DocTestFinder, DocTestRunner, DocTestParser
import re as pre
import random
import subprocess
import imp

from StringIO import StringIO
import string
import codecs


#
# Tune standard doctest module (utf8 encoding, indent)
#

# Use sys.stdout encoding for ouput.
# This was only added to Python's doctest in Python 2.6, so we must duplicate
# it here to make utf8 files work in Python 2.5.
pdoctest._encoding = getattr(sys.__stdout__, 'encoding', None) or 'utf-8'

# There is the same function in pdoctest already,
# but we must sure to be independed from python version
def _indent(s, indent=4):
    """
    Add the given number of space characters to the beginning of
    every non-blank line in `s`, and return the result.
    If the string `s` is Unicode, it is encoded using the stdout
    encoding and the `backslashreplace` error handler.
    """
    # After a 2to3 run the below code is bogus, so wrap it with a version check
    if sys.version_info[0] < 3:
        if isinstance(s, unicode):
            s = s.encode(pdoctest._encoding, 'backslashreplace')
    # This regexp matches the start of non-blank lines:
    return re.sub('(?m)^(?!$)', indent*' ', s)

pdoctest._indent = _indent

# Take some classes from the pdoctest

_exception_traceback = pdoctest._exception_traceback
Example = pdoctest.Example

# The Python 2.5 doctest runner uses a tuple, but in 2.6+, it uses a namedtuple
# (which doesn't exist in 2.5-)
if sys.version_info[:2] > (2,5):
    from collections import namedtuple
    SymPyTestResults = namedtuple('TestResults', 'failed attempted')
else:
    SymPyTestResults = lambda a, b: (a, b)

# Take some constants from the pdoctest. And define the new ones.

REPORT_ONLY_FIRST_FAILURE = pdoctest.REPORT_ONLY_FIRST_FAILURE
SKIP = pdoctest.SKIP
IGNORE_EXCEPTION_DETAIL = pdoctest.IGNORE_EXCEPTION_DETAIL
PRETTY = pdoctest.register_optionflag('PRETTY')
pdoctest.PRETTY = PRETTY

################################################################################
####                           Utilities                                    ####
################################################################################


def sys_normcase(f):
    if sys_case_insensitive:
        return f.lower()
    return f

def convert_to_native_paths(lst):
    """
    Converts a list of '/' separated paths into a list of
    native (os.sep separated) paths and converts to lowercase
    if the system is case insensitive.
    """
    newlst = []
    for i, rv in enumerate(lst):
        rv = os.path.join(*rv.split("/"))
        # on windows the slash after the colon is dropped
        if sys.platform == "win32":
            pos = rv.find(':')
            if pos != -1:
                if rv[pos+1] != '\\':
                    rv = rv[:pos+1] + '\\' + rv[pos+1:]
        newlst.append(sys_normcase(rv))
    return newlst

def get_sympy_dir():
    """
    Returns the root sympy directory and set the global value
    indicating whether the system is case sensitive or not.
    """
    global sys_case_insensitive

    this_file = os.path.abspath(__file__)
    sympy_dir = os.path.join(os.path.dirname(this_file), "..", "..")
    sympy_dir = os.path.normpath(sympy_dir)
    sys_case_insensitive = (os.path.isdir(sympy_dir) and
                        os.path.isdir(sympy_dir.lower()) and
                        os.path.isdir(sympy_dir.upper()))
    return sys_normcase(sympy_dir)

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

class PrintManager(object):
    """
    Class for manage test's output printer methods (pprint or strprinter).
    """
    def __init__(self):
        self._pprint_status = None

    def save_displayhook(self):
        self._displayhook = sys.displayhook

    def restore_displayhook(self):
        sys.displayhook = self._displayhook

    def setup_pprint(self, pretty_print = False, hard=False):
        """
        Load and set printer according the options.
        """
        if (self._pprint_status <> pretty_print) or hard:
            self._setup_pprint(pretty_print)
            self._pprint_status = pretty_print

    def _setup_pprint(self, pretty_print = False):

        from sympy import pprint_use_unicode, init_printing

        # force pprint to be in ascii mode in doctests
        pprint_use_unicode(False)

        # hook our nice, hash-stable strprinter
        init_printing(pretty_print=pretty_print)

print_manager = PrintManager()

sympy_dir = get_sympy_dir()

################################################################################
####                           Entering points                              ####
################################################################################

def test(*paths, **kwargs):
    """
    Run all tests in test_*.py files which match any of the given
    strings in `paths` or all tests if paths=[].

    Notes:
       o if sort=False, tests are run in random order (not default).
       o paths can be entered in native system format or in unix,
         forward-slash format.

    Examples:

    >> import sympy

    Run all tests:
    >> sympy.test()

    Run one file:
    >> sympy.test("sympy/core/tests/test_basic.py")
    >> sympy.test("_basic")

    Run all tests in sympy/functions/ and some particular file:
    >> sympy.test("sympy/core/tests/test_basic.py", "sympy/functions")

    Run all tests in sympy/core and sympy/utilities:
    >> sympy.test("/core", "/util")

    Run specific test from a file:
    >> sympy.test("sympy/core/tests/test_basic.py", kw="test_equality")

    Run specific test from any file:
    >> sympy.test(kw="subs")

    Run the tests with verbose mode on:
    >> sympy.test(verbose=True)

    Don't sort the test output:
    >> sympy.test(sort=False)

    Turn on post-mortem pdb:
    >> sympy.test(pdb=True)

    Turn off colors:
    >> sympy.test(colors=False)

    The traceback verboseness can be set to "short" or "no" (default is "short")
    >> sympy.test(tb='no')

    """
    verbose = kwargs.get("verbose", False)
    tb = kwargs.get("tb", "short")
    kw = kwargs.get("kw", "")
    post_mortem = kwargs.get("pdb", False)
    colors = kwargs.get("colors", True)
    sort = kwargs.get("sort", True)
    seed = kwargs.get("seed", None)
    if seed is None:
        seed = random.randrange(100000000)

    r = PyTestReporter(verbose, tb, colors)
    t = SymPyTester(r, kw, post_mortem, seed)
    t.load_sympy_modules()

    test_files = t.get_test_files('sympy')
    if len(paths) == 0:
        t._testfiles.extend(test_files)
    else:
        paths = convert_to_native_paths(paths)
        matched = []
        for f in test_files:
            basename = os.path.basename(f)
            for p in paths:
                if p in f or fnmatch(basename, p):
                    matched.append(f)
                    break
        t._testfiles.extend(matched)

    return t.test(sort=sort)


def doctest(*paths, **kwargs):
    """
    Runs doctests in all *.py files in the sympy directory which match
    any of the given strings in `paths` or all tests if paths=[].

    Also run doctests in all *.txt files in the sympy/doc/src directory.

    Note:
       o paths can be entered in native system format or in unix,
         forward-slash format.
       o files that are on the blacklist can be tested by providing
         their path; they are only excluded if no paths are given.

    Examples:

    >> import sympy

    Run all tests:
    >> sympy.doctest()

    Run one file:
    >> sympy.doctest("sympy/core/basic.py")
    >> sympy.doctest("polynomial.txt")

    Run all tests in sympy/functions/ and some particular file:
    >> sympy.doctest("/functions", "basic.py")

    Run any file having polynomial in its name, doc/src/modules/polynomial.txt,
    sympy\functions\special\polynomials.py, and sympy\polys\polynomial.py:
    >> sympy.doctest("polynomial")
    """
    normal = kwargs.get("normal", False)
    verbose = kwargs.get("verbose", False)
    colors = kwargs.get("colors", True)
    first_only = kwargs.get("first_only", False)
    blacklist = kwargs.get("blacklist", [])
    blacklist.extend([
                    "doc/src/modules/mpmath", # needs to be fixed upstream
                    "sympy/mpmath", # needs to be fixed upstream
                    "doc/src/modules/plotting.txt", # generates live plots
                    "sympy/plotting", # generates live plots
                    "sympy/utilities/compilef.py", # needs tcc
                    "sympy/utilities/autowrap.py", # needs installed compiler
                    "sympy/galgebra/GA.py", # needs numpy
                    "sympy/galgebra/latex_ex.py", # needs numpy
                    "sympy/conftest.py", # needs py.test
                    "sympy/utilities/benchmarking.py", # needs py.test
                    "examples/advanced/autowrap_integrators.py", # needs numpy
                    "examples/advanced/autowrap_ufuncify.py", # needs numpy
                    "examples/intermediate/mplot2d.py", # needs numpy and matplotlib
                    "examples/intermediate/mplot3d.py", # needs numpy and matplotlib
                    "examples/intermediate/sample.py", # needs numpy
                    ])
    blacklist = convert_to_native_paths(blacklist)

    r = PyTestReporter(verbose, colors=colors, first_only=first_only)
    t = SymPyDocTester(r, normal)
    t.load_sympy_modules()
    t.parser = SymPyTextTestParser()
    if first_only:
        t.optionflags |= pdoctest.REPORT_ONLY_FIRST_FAILURE

    r.start()
    # Step 1: Here we test doctests in all *py files
    test_files = t.get_test_files('sympy')
    test_files.extend(t.get_test_files('examples', init_only=False))
    matched = t.filter_files(test_files, blacklist, paths)
    t._testfiles.extend(matched)

    # setup prtinter
    print_manager.setup_pprint(False, hard=True)

    # run the tests and record the result for this *py portion of the tests
    if t._testfiles:
        failed = not t.test_docstrings()
    else:
        failed = False

    # Step 2: Here we test *.txt files at or below doc/src.
    #
    # Code from these must
    # be self supporting in terms of imports since there is no importing
    # of necessary modules by doctest.testfile. If you try to pass *.py
    # files through this they might fail because they will lack the needed
    # imports and smarter parsing that can be done with source code.
    #
    failed_2 = False
    test_files = t.get_test_files('doc/src', '*.txt', init_only=False)
    test_files.sort()
    matched = t.filter_files(test_files, blacklist, paths)
    t._testfiles = matched
    t.pprint_status = None
    print_manager.setup_pprint(False, hard=True)

    if t._testfiles:
        failed_2 = not t.test_text_files()

    r.finish()
    ok = not (failed or failed_2)
    assert r.ok == ok
    return ok

def wikitest(*paths, **kwargs):
    """
    Runs wikitest in files in the wiki_dir directory.

    The files are to be passing throu the test contain the special directives.
    For reStructedText comments is shaped as:
        .. wikitest release
    For markdown:
        <!-- wikitest release -->

    "wikitest release" directive means that tests for current wiki-page have
    to be passed against release version only.
    "wikitest master" - against master.
    "wikitest master,release" - against master and release version

    To run tests  only for a few files, `paths` is used. If paths=[] then all
    files with directives tests.

    Note:
       o paths can be entered in native system format or in unix,
         forward-slash format.
       o files that are on the blacklist can be tested by providing
         their path; they are only excluded if no paths are given.

    Examples:

    >> from sympy.utilities.runtests import wikitest

    Run all tests:

    >> wikitest()

    """
    normal = kwargs.get("normal", False)
    verbose = kwargs.get("verbose", False)
    colors = kwargs.get("colors", True)
    first_only = kwargs.get("first_only", False)

    againstlist = kwargs.get("againstlist", ["master"])
    blacklist = kwargs.get("blacklist", [])
    blacklist.extend([])
    blacklist = convert_to_native_paths(blacklist)
    wiki_dir = kwargs["wiki_dir"]

    r = PyTestReporter(verbose, colors=colors, first_only=first_only)
    t = SymPyWikiTester(r, normal, wiki_dir)
    t.parser = SymPyTextTestParser()
    if first_only:
        t.optionflags |= pdoctest.REPORT_ONLY_FIRST_FAILURE

    failed_total = False

    for against in againstlist:
        t.set_against_dir(against, kwargs["%s_dir" % against])
        t.load_sympy_version()
        t.load_sympy_modules()

        test_files = t.get_test_files(wiki_dir, pat=r"^.*\.(md|rest|mediawiki)$",init_only=False)

        matched = t.filter_files(test_files, blacklist, paths)
        t._testfiles = matched
        t.pprint_status = None
        print_manager.setup_pprint(False, hard=True)
        if t._testfiles:
            r.start()
            failed = not t.test_text_files()
            r.finish()

        if failed:
            failed_total = True
    return not failed_total



################################################################################
####                           Tester classes                               ####
################################################################################



#    Hierarchy:
#
#
#    SymPyTesterBase
#              // base tester class
#        methods:
#         - filter_files()      // Filter files by blacklist
#         - load_sympy_modules() // load sympy modules
#         - clear_cache()
#
#
#    SymPyTester(SymPyTesterBase)
#               // find and test 'test_*.py' files
#        methods:
#        - get_test_files()     // find 'test_*.py'
#        - test()               // test them
#        - test_python_file()   // test a 'test_*.py' file
#
#
#
#    SymPyDocTester(SymPyTesterBase)
#              // find and test docstrings and *.txt files
#        methods:
#        - get_test_files()          // find `*.py' or `*.txt` files
#        - test_docstrings()         // test docstrings in the files founded out
#        - test_docstrings_in_file() // test docstrings in the '*.py' file
#        - test_text_files()         // test examples in the files founded out
#        - test_text_file()          // test examples in the '*.txt' file and in the wiki page
#
#
#
#    SymPyWikiTester(SymPyDocTester)
#               // find and test wiki pages
#       methods:
#        - get_test_files()       // find `*.rst` and `*.md` files
#        - load_sympy_version()   // load right SymPy's version.
#

class SymPyTesterBase(object):
    """
    Base tester class.

    """
    def __init__(self):
        self.pprint_status = None

    def load_sympy_modules(self):
        """
        Manually load SymPy modules.

        E.g. clear_cache() function, without loading of the sympy itself because another
        release can be tested) It is needed for wikitest, wehere test passed against
        another release of SymPy so the modules only of that release must be loaded.
        """
        #fn_core_cache  = os.path.join(os.path.dirname(__file__), "../core/cache.py")
        #module_core_cache = imp.load_source("cache", fn_core_cache)
        #clear_cache = module_core_cache.clear_cache

        from sympy.core.cache import clear_cache
        self._clear_cache = clear_cache

        # Disable warnings for external modules
        import sympy.external
        sympy.external.importtools.WARN_OLD_VERSION = False
        sympy.external.importtools.WARN_NOT_INSTALLED = False

        global ARCH
        from sympy.utilities.misc import ARCH
        global USE_CACHE
        try:
            from sympy.core.cache import USE_CACHE
        except:
            USE_CACHE = "unknown"
        global GROUND_TYPES
        from sympy.polys.domains import GROUND_TYPES

    def clear_cache(self):
        self._clear_cache()

    def filter_files(self, test_files, blacklist, argument_paths):
        """
        Filter files by blacklist and/or by argument of command.
        """
        not_blacklisted = [f for f in test_files
                             if not any(b in f for b in blacklist)]
        if len(argument_paths) == 0:
            matched = not_blacklisted
        else:
            # Take only what was requested as long as it's not on the blacklist.
            # Paths were already made native in *py tests so don't repeat here.
            # There's no chance of having a *py file slip through since we
            # only have *txt files in test_files.
            matched =  []
            for f in not_blacklisted:
                basename = os.path.basename(f)
                for p in argument_paths:
                    if p in f or fnmatch(basename, p):
                        matched.append(f)
                        break
        return matched

class SymPyTester(SymPyTesterBase):
    """
    Class for the testing if `test_*.py` files.
    """

    def __init__(self, reporter, kw="", post_mortem=False,
                seed=random.random()):
        SymPyTesterBase.__init__(self)

        self._post_mortem = post_mortem
        self._kw = kw
        self._count = 0
        self._root_dir = sympy_dir
        self._reporter = reporter
        self._reporter.root_dir(self._root_dir)
        self._testfiles = []
        self._seed = seed

    def get_test_files(self, dir, pat = 'test_*.py'):
        """
        Returns the list of test_*.py (default) files at or below directory
        `dir` relative to the sympy home directory.
        """
        dir = os.path.join(self._root_dir, convert_to_native_paths([dir])[0])

        g = []
        for path, folders, files in os.walk(dir):
            g.extend([os.path.join(path, f) for f in files if fnmatch(f, pat)])

        return [sys_normcase(gi) for gi in g]

    def test(self, sort=False):
        """
        Runs the tests returning True if all tests pass, otherwise False.

        If sort=False run tests in random order.
        """
        if sort:
            self._testfiles.sort()
        else:
            from random import shuffle
            random.seed(self._seed)
            shuffle(self._testfiles)
        self._reporter.start(self._seed)
        for f in self._testfiles:
            try:
                self.test_python_file(f)
            except KeyboardInterrupt:
                print " interrupted by user"
                break
        return self._reporter.finish()

    def test_python_file(self, filename):
        """
        Test the tests in the "test___.py" file.
        """
        # Clear cache for every file tested
        self.clear_cache()

        self._count += 1
        gl = {'__file__':filename}
        random.seed(self._seed)

        try:
            execfile(filename, gl)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
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
            while i < len(funcs):
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

        if not funcs:
            return
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
                    self._reporter.test_skip(v)
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

    def matches(self, x):
        """
        Does the keyword expression self._kw match "x"? Returns True/False.

        Always returns True if self._kw is "".
        """
        if self._kw == "":
            return True
        return x.__name__.find(self._kw) != -1


class SymPyDocTester(SymPyTesterBase):
    """
    Class for the testing of docstring and examples in the documentation.

    It tests examples in the docstrings of the `*.py` files, also it tests
    the examples in the `doc/*.txt` files.
    """
    def __init__(self, reporter, normal):
        SymPyTesterBase.__init__(self)

        self._count = 0
        self._root_dir = sympy_dir
        self._reporter = reporter
        self._reporter.root_dir(self._root_dir)
        self._normal = normal

        self._testfiles = []
        self.optionflags = pdoctest.ELLIPSIS | pdoctest.NORMALIZE_WHITESPACE | \
            pdoctest.DONT_ACCEPT_TRUE_FOR_1 | pdoctest.IGNORE_EXCEPTION_DETAIL

    def get_test_files(self, dir, pat='*.py', init_only=True):
        """
        Returns the list of *py files (default) from which docstrings
        will be tested which are at or below directory `dir`. By default,
        only those that have an __init__.py in their parent directory
        and do not start with `test_` will be included.
        """
        def importable(x):
            """
            Checks if given pathname x is an importable module by checking for
            __init__.py file.

            Returns True/False.

            Currently we only test if the __init__.py file exists in the
            directory with the file "x" (in theory we should also test all the
            parent dirs).
            """
            init_py = os.path.join(os.path.dirname(x), "__init__.py")
            return os.path.exists(init_py)

        dir = os.path.join(self._root_dir, convert_to_native_paths([dir])[0])

        g = []
        for path, folders, files in os.walk(dir):
            g.extend([os.path.join(path, f) for f in files
                      if not f.startswith('test_') and fnmatch(f, pat)])
        if init_only:
            # skip files that are not importable (i.e. missing __init__.py)
            g = [x for x in g if importable(x)]

        return [sys_normcase(gi) for gi in g]

    def test_docstrings(self):
        """
        Runs the tests and returns True if all tests pass, otherwise False.
        """
        for f in self._testfiles:
            try:
                self.test_docstrings_in_file(f)
            except KeyboardInterrupt:
                print " interrupted by user"
                break
        return self._reporter.ok

    def test_docstrings_in_file(self, filename):
        """
        Test docstrings in the *.py file.
        """
        # Clear cache for every file tested
        self.clear_cache()
        rel_name = filename[len(self._root_dir)+1:]
        dirname, file = os.path.split(filename)
        module = rel_name.replace(os.sep, '.')[:-3]

        if rel_name.startswith("examples"):
            # Example files do not have __init__.py files,
            # So we have to temporarily extend sys.path to import them
            sys.path.insert(0, dirname)
            module = file[:-3] # remove ".py"

        print_manager.setup_pprint(False, hard=True)

        try:
            module = pdoctest._normalize_module(module)
            tests = SymPyDocTestFinder().find(module)
        except (SystemExit, KeyboardInterrupt):
            raise
        except:
            self._reporter.import_error(filename, sys.exc_info())
            return

        finally:
            if rel_name.startswith("examples"):
                del sys.path[0]

        tests = [test for test in tests if len(test.examples) > 0]
        # By default tests are sorted by alphabetical order by function name.
        # We sort by line number so one can edit the file sequentially from
        # bottom to top. However, if there are decorated functions, their line
        # numbers will be too large and for now one must just search for these
        # by text and function name.
        tests.sort(key=lambda x: -x.lineno)

        if not tests:
            return
        self._reporter.entering_filename(filename, len(tests))
        for test in tests:
            assert len(test.examples) != 0
            runner = SymPyDocTestRunner(optionflags=pdoctest.ELLIPSIS | \
                    pdoctest.NORMALIZE_WHITESPACE | \
                    pdoctest.IGNORE_EXCEPTION_DETAIL)

            old = sys.stdout
            new = StringIO()
            sys.stdout = new
            # If the testing is normal, the doctests get importing magic to
            # provide the global namespace. If not normal (the default) then
            # then must run on their own; all imports must be explicit within
            # a function's docstring. Once imported that import will be
            # available to the rest of the tests in a given function's
            # docstring (unless clear_globs=True below).
            if not self._normal:
                test.globs = {}
                # if this is uncommented then all the test would get is what
                # comes by default with a "from sympy import *"
                #exec('from sympy import *') in test.globs
            try:
                f, t = runner.run(test, out=new.write, clear_globs=False)
            finally:
                sys.stdout = old
            if f > 0:
                self._reporter.doctest_fail(test.name, new.getvalue())
            else:
                self._reporter.test_pass()
        self._reporter.leaving_filename()

    def test_text_files(self):
        """
        Runs the tests and returns True if all tests pass, otherwise False.
        """
        for f in self._testfiles:
            # make sure we return to the original displayhook in case some
            # doctest has changed that
            old_displayhook = sys.displayhook
            old_optionflags = self.optionflags
            self.pre_options(f)

            try:
                self.test_text_file(f, module_relative=False,
                    optionflags=self.optionflags,
                    parser=self.parser, encoding='utf-8')
                sys.displayhook = old_displayhook
            except KeyboardInterrupt:
                sys.displayhook = old_displayhook
                print " interrupted by user"
                break
            finally:
                sys.displayhook = old_displayhook
                self.optionflags = old_optionflags

        return self._reporter.ok

    def test_text_file(self, filename, module_relative=True, name=None, package=None,
                 globs=None, verbose=None, report=True, optionflags=0,
                 extraglobs=None, raise_on_error=False,
                 parser=DocTestParser, encoding=None):
        """
        Test examples in the given (*.txt) file.

        Return (#failures, #tests).

        Optional keyword arg "module_relative" specifies how filenames
        should be interpreted:

          - If "module_relative" is True (the default), then "filename"
             specifies a module-relative path.  By default, this path is
             relative to the calling module's directory; but if the
             "package" argument is specified, then it is relative to that
             package.  To ensure os-independence, "filename" should use
             "/" characters to separate path segments, and should not
             be an absolute path (i.e., it may not begin with "/").

          - If "module_relative" is False, then "filename" specifies an
            os-specific path.  The path may be absolute or relative (to
            the current working directory).

        Optional keyword arg "name" gives the name of the test; by default
        use the file's basename.

        Optional keyword argument "package" is a Python package or the
        name of a Python package whose directory should be used as the
        base directory for a module relative filename.  If no package is
        specified, then the calling module's directory is used as the base
        directory for module relative filenames.  It is an error to
        specify "package" if "module_relative" is False.

        Optional keyword arg "globs" gives a dict to be used as the globals
        when executing examples; by default, use {}.  A copy of this dict
        is actually used for each docstring, so that each docstring's
        examples start with a clean slate.

        Optional keyword arg "extraglobs" gives a dictionary that should be
        merged into the globals that are used to execute examples.  By
        default, no extra globals are used.

        Optional keyword arg "verbose" prints lots of stuff if true, prints
        only failures if false; by default, it's true iff "-v" is in sys.argv.

        Optional keyword arg "report" prints a summary at the end when true,
        else prints nothing at the end.  In verbose mode, the summary is
        detailed, else very brief (in fact, empty if all tests passed).

        Optional keyword arg "optionflags" or's together module constants,
        and defaults to 0.  Possible values (see the docs for details):

            DONT_ACCEPT_TRUE_FOR_1
            DONT_ACCEPT_BLANKLINE
            NORMALIZE_WHITESPACE
            ELLIPSIS
            SKIP
            IGNORE_EXCEPTION_DETAIL
            REPORT_UDIFF
            REPORT_CDIFF
            REPORT_NDIFF
            REPORT_ONLY_FIRST_FAILURE

        Optional keyword arg "raise_on_error" raises an exception on the
        first unexpected exception or failure. This allows failures to be
        post-mortem debugged.

        Optional keyword arg "parser" specifies a DocTestParser (or
        subclass) that should be used to extract tests from the files.

        Optional keyword arg "encoding" specifies an encoding that should
        be used to convert the file to unicode.

        Advanced tomfoolery:  testmod runs methods of a local instance of
        class doctest.Tester, then merges the results into (or creates)
        global Tester instance doctest.master.  Methods of doctest.master
        can be called directly too, if you want to do something unusual.
        Passing report=0 to testmod is especially useful then, to delay
        displaying a summary.  Invoke doctest.master.summarize(verbose)
        when you're done fiddling.
        """
        if package and not module_relative:
            raise ValueError("Package may only be specified for module-"
                             "relative paths.")

        # Relativize the path
        if sys.version_info[0] < 3:
            # Fix doctest runner for Python 3
            text, filename = pdoctest._load_testfile(filename, package, module_relative)
        else:
            encoding = None
            text, filename = pdoctest._load_testfile(filename, package, module_relative, encoding)

        # If no name was given, then use the file's name.
        if name is None:
            name = os.path.basename(filename)

        # Assemble the globals.
        if globs is None:
            globs = {}
        else:
            globs = globs.copy()
        if extraglobs is not None:
            globs.update(extraglobs)
        if '__name__' not in globs:
            globs['__name__'] = '__main__'

        if raise_on_error:
            runner = pdoctest.DebugRunner(verbose=verbose, optionflags=optionflags)
        else:
            runner = SymPyDocTestRunner(verbose=verbose, optionflags=optionflags)

        if encoding is not None:
            text = text.decode(encoding)

        # Read the file, convert it to a test, and run it.
        tests = parser.get_doctests(text, globs, name, filename, 0)

        if len(tests):
            self._reporter.entering_filename(filename, len(tests))
            for test in tests:
                test.globs = globs
                old = sys.stdout
                new = StringIO()
                sys.stdout = new
                try:
                    f, t = runner.run(test, out=new.write, clear_globs=False)
                finally:
                    sys.stdout = old
                globs.update(test.globs)
                if f > 0:
                    self._reporter.doctest_fail(test.name, new.getvalue())
                else:
                    self._reporter.test_pass()
            self._reporter.leaving_filename()

    def pre_options(self, fn):
        pass


class SymPyWikiTester(SymPyDocTester):
    """
    Class for the testing examples in the wiki pages.
    """

    def __init__(self, reporter, normal, root_dir):
        SymPyDocTester.__init__(self, reporter, normal)

        self._count = 0
        self._root_dir = root_dir
        self._reporter = reporter
        self._reporter.root_dir(self._root_dir)
        self._normal = normal

        self._testfiles = []
        re_directives = []
        re_directives.append(re.compile(r"<!--\s+wikitest\s+(?P<against>[\w,]+)(\s+\((?P<options>[\w,]+)\))?\s+-->", re.M))
        re_directives.append(re.compile(r"\.\.\s+wikitest\s+(?P<against>[\w,]+)(\s+\((?P<options>[\w,]+)\))?\s+$", re.M))
        self.re_directives = re_directives

    def set_against_dir(self, against, against_dir):
        """Set against dir.

        While testing this dir will be injected into the test for catching right
        version of SymPy (release, or master)
        """
        self.against = against
        self.against_dir = against_dir

    def get_test_files(self, dir, pat=r'^.*\.md$', init_only=True):
        """
        Returns the list of wiki-pages (default) from which docstrings
        will be tested which are at or below directory `dir`.
        """
        re_filename = re.compile(pat)

        def _fnmatch(fn):
            if re_filename.match(fn):
                return True
            return False
        #dir = os.path.join(self._root_dir, convert_to_native_paths([dir])[0])
        dir = os.path.join("/", convert_to_native_paths([dir])[0])

        g = []
        for path, folders, files in os.walk(dir):
            g.extend([os.path.join(path, f) for f in files
                      if _fnmatch(f) and self.has_directive(path, f)])
        return [sys_normcase(gi) for gi in g]

    def has_directive(self, pathdir, fn):
        """
        Checks if file contain the directive.
        """
        fn = os.path.join(pathdir, fn)
        f = codecs.open(fn, mode="r", encoding="utf8")
        s = f.read()
        f.close()
        for re_directive in self.re_directives:
            m = re_directive.search(s)
            if m:
                against_permitted = m.group("against").split(",")
                if self.against in against_permitted:
                    return True
        return False

    def get_directive_options(self, filename):
        res = {}
        f = codecs.open(filename, mode="r", encoding="utf8")
        s = f.read()
        f.close()
        for re_directive in self.re_directives:
            m = re_directive.search(s)
            if m:
                against_permitted = m.group("against").split(",")
                if self.against in against_permitted:
                    options = m.group("options")
                    if options:
                        options = options.split(",")
                        for op in options:
                            res[op] = True
        return res

    def pre_options(self, fn):
        directive_options = self.get_directive_options(fn)
        pretty_print = directive_options.get("pretty_print", False)
        if pretty_print:
            self.optionflags |= PRETTY

        print_manager.setup_pprint(pretty_print, hard=True)

    def load_sympy_version(self):
        """
        Load right SymPy's version.
        And unload previous if it's needed
        """
        for _key in sys.modules.keys():
            if len(_key.split("sympy")) > 1:
                del sys.modules[_key]
        sys.path.insert(0, self.against_dir)


################################################################################
####           Finder, Parser, Runner classes                                ####
################################################################################

class SymPyTextTestParser(DocTestParser):
    """
    Parser to find tests in the wiki page.

    It extend standard DocTestParser to catch Markdown code-block syntax with
    pigment like ```py   >>> x = Symbol('x')```
    """
    _EXAMPLE_RE = re.compile(r'''
        # Source consists of a PS1 line followed by zero or more PS2 lines.
        (?P<source>
            (?:^(?P<indent> [ ]*) >>>    .*)    # PS1 line
            (?:\n           [ ]*  \.\.\. .*)*)  # PS2 lines
        \n?
        # Want consists of any non-blank lines that do not start with PS1.
        (?P<want> (?:(?![ ]*$)    # Not a blank line
                     (?![ ]*>>>)  # Not a line starting with PS1
                     (?![ ]*```)  # Not a ```
                     .*$\n?       # But any other line
                  )*)
        ''', re.MULTILINE | re.VERBOSE)

    def get_doctests(self, text, globs, name, filename, lineno):
        """
        Parse text, and extract examples grouped by tests.

        It extends "get_doctest" to group examples to blocks (tests).
        """

        doctests = []
        splitted_examples = self.parse_ext(text, name)
        examples = []
        examples_fullsource = []
        for x in splitted_examples:
            if isinstance(x, Example):
                examples.append(x)
                examples_fullsource.append(x._fullsource)
            elif self._IS_BLANK_OR_COMMENT(x):
                examples_fullsource.append(x)
            else:
                # form DocTest
                if len(examples):
                    doctest = pdoctest.DocTest(examples, globs, name, filename,
                        0,
                        "".join(examples_fullsource))
                    doctests.append(doctest)
                    examples = []
                    examples_fullsource = []
        return doctests

    def parse_ext(self, string, name='<string>'):
        """
        Parse text, and extract examples grouped by tests.

        It extends "get_doctest" to group examples to blocks (tests).
        The code is based on the DocTestParser.parse() method and extend it
        by only adding aditional info for examples ("_fullsource").
        """
        string = string.expandtabs()
        # If all lines begin with the same indentation, then strip it.
        min_indent = self._min_indent(string)
        if min_indent > 0:
            string = '\n'.join([l[min_indent:] for l in string.split('\n')])

        output = []
        charno, lineno = 0, 0
        # Find all doctest examples in the string:
        for m in self._EXAMPLE_RE.finditer(string):
            # Add the pre-example text to `output`.
            output.append(string[charno:m.start()])
            # Update lineno (lines before this example)
            lineno += string.count('\n', charno, m.start())
            # Extract info from the regexp match.
            (source, options, want, exc_msg) = \
                     self._parse_example(m, name, lineno)
            # Create an Example, and add it to the list.
            if not self._IS_BLANK_OR_COMMENT(source):
                example = Example(source, want, exc_msg,
                                    lineno=lineno,
                                    indent=min_indent+len(m.group('indent')),
                                    options=options)
                example._fullsource = string[m.start():m.end()]
                output.append(example)

            # Update lineno (lines inside this example)
            lineno += string.count('\n', m.start(), m.end())
            # Update charno.
            charno = m.end()
        # Add any remaining post-example text to `output`.
        output.append(string[charno:])
        return output


class SymPyDocTestFinder(DocTestFinder):
    """
    A class used to extract the DocTests that are relevant to a given
    object, from its docstring and the docstrings of its contained
    objects.  Doctests can currently be extracted from the following
    object types: modules, functions, classes, methods, staticmethods,
    classmethods, and properties.

    Modified from doctest's version by looking harder for code in the
    case that it looks like the the code comes from a different module.
    In the case of decorated functions (e.g. @vectorize) they appear
    to come from a different module (e.g. multidemensional) even though
    their code is not there.
    """

    def _find(self, tests, obj, name, module, source_lines, globs, seen):
        """
        Find tests for the given object and any contained objects, and
        add them to `tests`.
        """
        if self._verbose:
            print 'Finding tests in %s' % name

        # If we've already processed this object, then ignore it.
        if id(obj) in seen:
            return
        seen[id(obj)] = 1

        # Make sure we don't run doctests for classes outside of sympy, such
        # as in numpy or scipy.
        if inspect.isclass(obj):
            if obj.__module__.split('.')[0] != 'sympy':
                return

        # Find a test for this object, and add it to the list of tests.
        test = self._get_test(obj, name, module, globs, source_lines)
        if test is not None:
            tests.append(test)

        # Look for tests in a module's contained objects.
        if inspect.ismodule(obj) and self._recurse:
            for rawname, val in obj.__dict__.items():
                # Recurse to functions & classes.
                if inspect.isfunction(val) or inspect.isclass(val):
                    in_module = self._from_module(module, val)
                    if not in_module:
                        # double check in case this function is decorated
                        # and just appears to come from a different module.
                        pat = r'\s*(def|class)\s+%s\s*\(' % rawname
                        PAT = pre.compile(pat)
                        in_module = any(PAT.match(line) for line in source_lines)
                    if in_module:
                        try:
                            valname = '%s.%s' % (name, rawname)
                            self._find(tests, val, valname, module, source_lines, globs, seen)
                        except ValueError, msg:
                            raise
                        except:
                            pass

        # Look for tests in a module's __test__ dictionary.
        if inspect.ismodule(obj) and self._recurse:
            for valname, val in getattr(obj, '__test__', {}).items():
                if not isinstance(valname, basestring):
                    raise ValueError("SymPyDocTestFinder.find: __test__ keys "
                                     "must be strings: %r" %
                                     (type(valname),))
                if not (inspect.isfunction(val) or inspect.isclass(val) or
                        inspect.ismethod(val) or inspect.ismodule(val) or
                        isinstance(val, basestring)):
                    raise ValueError("SymPyDocTestFinder.find: __test__ values "
                                     "must be strings, functions, methods, "
                                     "classes, or modules: %r" %
                                     (type(val),))
                valname = '%s.__test__.%s' % (name, valname)
                self._find(tests, val, valname, module, source_lines,
                           globs, seen)

        # Look for tests in a class's contained objects.
        if inspect.isclass(obj) and self._recurse:
            for valname, val in obj.__dict__.items():
                # Special handling for staticmethod/classmethod.
                if isinstance(val, staticmethod):
                    val = getattr(obj, valname)
                if isinstance(val, classmethod):
                    val = getattr(obj, valname).im_func

                # Recurse to methods, properties, and nested classes.
                if (inspect.isfunction(val) or
                     inspect.isclass(val) or
                     isinstance(val, property)):
                    in_module = self._from_module(module, val)
                    if not in_module:
                        # "double check" again
                        pat = r'\s*(def|class)\s+%s\s*\(' % valname
                        PAT = pre.compile(pat)
                        in_module = any(PAT.match(line) for line in source_lines)
                    if in_module:
                        valname = '%s.%s' % (name, valname)
                        self._find(tests, val, valname, module, source_lines,
                                   globs, seen)

    def _get_test(self, obj, name, module, globs, source_lines):
        """
        Return a DocTest for the given object, if it defines a docstring;
        otherwise, return None.
        """
        # Extract the object's docstring.  If it doesn't have one,
        # then return None (no test for this object).
        if isinstance(obj, basestring):
            docstring = obj
        else:
            try:
                if obj.__doc__ is None:
                    docstring = ''
                else:
                    docstring = obj.__doc__
                    if not isinstance(docstring, basestring):
                        docstring = str(docstring)
            except (TypeError, AttributeError):
                docstring = ''

        # Find the docstring's location in the file.
        lineno = self._find_lineno(obj, source_lines)

        if lineno is None:
            # if None, then _find_lineno couldn't find the docstring.
            # But IT IS STILL THERE.  Likely it was decorated or something
            # (i.e., @property docstrings have lineno == None)
            # TODO: Write our own _find_lineno that is smarter in this regard
            # Until then, just give it a dummy lineno.  This is just used for
            # sorting the tests, so the only bad effect is that they will appear
            # last instead of the order that they really are in the file.
            # lineno is also used to report the offending line of a failing
            # doctest, which is another reason to fix this.  See issue 1947.
            lineno = 0

        # Don't bother if the docstring is empty.
        if self._exclude_empty and not docstring:
            return None

        # Return a DocTest for this object.
        if module is None:
            filename = None
        else:
            filename = getattr(module, '__file__', module.__name__)
            if filename[-4:] in (".pyc", ".pyo"):
                filename = filename[:-1]
        return self._parser.get_doctest(docstring, globs, name,
                                        filename, lineno)


class SymPyDocTestRunner(DocTestRunner):
    """
    A class used to run DocTest test cases, and accumulate statistics.
    The `run` method is used to process a single DocTest case.  It
    returns a tuple `(f, t)`, where `t` is the number of test cases
    tried, and `f` is the number of test cases that failed.

    Modified from the doctest version to not reset the sys.displayhook (see
    issue 2041).

    See the docstring of the original DocTestRunner for more information.
    """

    def run(self, test, compileflags=None, out=None, clear_globs=True):
        """
        Run the examples in `test`, and display the results using the
        writer function `out`.

        The examples are run in the namespace `test.globs`.  If
        `clear_globs` is true (the default), then this namespace will
        be cleared after the test runs, to help with garbage
        collection.  If you would like to examine the namespace after
        the test completes, then use `clear_globs=False`.

        `compileflags` gives the set of flags that should be used by
        the Python compiler when running the examples.  If not
        specified, then it will default to the set of future-import
        flags that apply to `globs`.

        The output of each example is checked using
        `SymPyDocTestRunner.check_output`, and the results are formatted by
        the `SymPyDocTestRunner.report_*` methods.
        """
        self.test = test

        if compileflags is None:
            compileflags = pdoctest._extract_future_flags(test.globs)

        save_stdout = sys.stdout
        if out is None:
            out = save_stdout.write
        sys.stdout = self._fakeout

        # Patch pdb.set_trace to restore sys.stdout during interactive
        # debugging (so it's not still redirected to self._fakeout).
        # Note that the interactive output will go to *our*
        # save_stdout, even if that's not the real sys.stdout; this
        # allows us to write test cases for the set_trace behavior.
        save_set_trace = pdb.set_trace
        self.debugger = pdoctest._OutputRedirectingPdb(save_stdout)
        self.debugger.reset()
        pdb.set_trace = self.debugger.set_trace

        # Patch linecache.getlines, so we can see the example's source
        # when we're inside the debugger.
        self.save_linecache_getlines = pdoctest.linecache.getlines
        linecache.getlines = self.__patched_linecache_getlines

        try:
            return self.__run(test, compileflags, out)
        finally:
            sys.stdout = save_stdout
            pdb.set_trace = save_set_trace
            linecache.getlines = self.save_linecache_getlines
            if clear_globs:
                test.globs.clear()

    # this method have gotten from original doctest module
    # except
    #    (a) the checking of `PRETTY` option for every example.
    #    (b) updating compileflags for every example
    #        to handle properly  with "from __future__ import division".
    def __run(self, test, compileflags, out):
        """
        Run the examples in `test`.  Write the outcome of each example
        with one of the `DocTestRunner.report_*` methods, using the
        writer function `out`.  `compileflags` is the set of compiler
        flags that should be used to execute examples.  Return a tuple
        `(f, t)`, where `t` is the number of examples tried, and `f`
        is the number of examples that failed.  The examples are run
        in the namespace `test.globs`.
        """
        # Keep track of the number of failures and tries.
        failures = tries = 0

        # Save the option flags (since option directives can be used
        # to modify them).
        original_optionflags = self.optionflags

        SUCCESS, FAILURE, BOOM = range(3) # `outcome` state

        check = self._checker.check_output

        # Process each example.
        for examplenum, example in enumerate(test.examples):

            # (b) updating compileflags for every example.
            if compileflags ==0:
                compileflags = pdoctest._extract_future_flags(test.globs)

            # If REPORT_ONLY_FIRST_FAILURE is set, then suppress
            # reporting after the first failure.
            quiet = (self.optionflags & REPORT_ONLY_FIRST_FAILURE and
                     failures > 0)

            # Merge in the example's options.
            self.optionflags = original_optionflags
            if example.options:
                for (optionflag, val) in example.options.items():
                    if val:
                        self.optionflags |= optionflag
                    else:
                        self.optionflags &= ~optionflag

            # If 'SKIP' is set, then skip this example.
            if self.optionflags & SKIP:
                continue

            # Record that we started this example.
            tries += 1
            if not quiet:
                self.report_start(out, test, example)

            # Use a special filename for compile(), so we can retrieve
            # the source code during interactive debugging (see
            # __patched_linecache_getlines).
            filename = '<doctest %s[%d]>' % (test.name, examplenum)

            # (a) the checking of `PRETTY` option for every example.
            # Set printer
            if self.optionflags & PRETTY:
                print_manager.setup_pprint(True)
            else:
                print_manager.setup_pprint(False)


            # Run the example in the given context (globs), and record
            # any exception that gets raised.  (But don't intercept
            # keyboard interrupts.)
            try:
                # Don't blink!  This is where the user's code gets run.
                exec compile(example.source, filename, "single",
                             compileflags, 1) in test.globs
                self.debugger.set_continue() # ==== Example Finished ====
                exception = None
            except KeyboardInterrupt:
                raise
            except:
                exception = sys.exc_info()
                self.debugger.set_continue() # ==== Example Finished ====

            got = self._fakeout.getvalue()  # the actual output
            self._fakeout.truncate(0)
            outcome = FAILURE   # guilty until proved innocent or insane

            # If the example executed without raising any exceptions,
            # verify its output.
            if exception is None:
                if check(example.want, got, self.optionflags):
                    outcome = SUCCESS

            # The example raised an exception:  check if it was expected.
            else:
                exc_msg = traceback.format_exception_only(*exception[:2])[-1]
                if not quiet:
                    got += _exception_traceback(exception)

                # If `example.exc_msg` is None, then we weren't expecting
                # an exception.
                if example.exc_msg is None:
                    outcome = BOOM

                # We expected an exception:  see whether it matches.
                elif check(example.exc_msg, exc_msg, self.optionflags):
                    outcome = SUCCESS

                # Another chance if they didn't care about the detail.
                elif self.optionflags & IGNORE_EXCEPTION_DETAIL:
                    m1 = re.match(r'[^:]*:', example.exc_msg)
                    m2 = re.match(r'[^:]*:', exc_msg)
                    if m1 and m2 and check(m1.group(0), m2.group(0),
                                           self.optionflags):
                        outcome = SUCCESS

            # Report the outcome.
            if outcome is SUCCESS:
                if not quiet:
                    self.report_success(out, test, example, got)
            elif outcome is FAILURE:
                if not quiet:
                    self.report_failure(out, test, example, got)
                failures += 1
            elif outcome is BOOM:
                if not quiet:
                    self.report_unexpected_exception(out, test, example,
                                                     exception)
                failures += 1
            else:
                assert False, ("unknown outcome", outcome)

        # Restore the option flags (in case they were modified)
        self.optionflags = original_optionflags

        # Record and return the number of failures and tries.
        self.__record_outcome(test, failures, tries)
        return SymPyTestResults(failures, tries)

# We have to override the name mangled methods.
SymPyDocTestRunner._SymPyDocTestRunner__patched_linecache_getlines = \
    DocTestRunner._DocTestRunner__patched_linecache_getlines
#SymPyDocTestRunner._SymPyDocTestRunner__run = DocTestRunner._DocTestRunner__run
SymPyDocTestRunner._SymPyDocTestRunner__record_outcome = \
    DocTestRunner._DocTestRunner__record_outcome

################################################################################
####                           Reporter Class                               ####
################################################################################


class Reporter(object):
    """
    Parent class for all reporters.
    """
    pass

class PyTestReporter(Reporter):
    """
    Py.test like reporter. Should produce output identical to py.test.
    """

    def __init__(self, verbose=False, tb="short", colors=True, first_only=False):
        self._verbose = verbose
        self._tb_style = tb
        self._colors = colors
        self._first_only = first_only
        self._xfailed = 0
        self._xpassed = []
        self._failed = []
        self._failed_doctest = []
        self._passed = 0
        self._skipped = 0
        self._exceptions = []
        self._terminal_width = None
        self._default_width = 80

        # this tracks the x-position of the cursor (useful for positioning
        # things on the screen), without the need for any readline library:
        self._write_pos = 0
        self._line_wrap = False


    def root_dir(self, dir):
        self._root_dir = dir

    @property
    def terminal_width(self):
        if self._terminal_width is not None:
            return self._terminal_width

        def findout_terminal_width():
            if sys.platform == "win32":
                # Windows support is based on:
                #
                #  http://code.activestate.com/recipes/440694-determine-size-of-console-window-on-windows/

                from ctypes import windll, create_string_buffer

                h = windll.kernel32.GetStdHandle(-12)
                csbi = create_string_buffer(22)
                res = windll.kernel32.GetConsoleScreenBufferInfo(h, csbi)

                if res:
                    import struct
                    (_, _, _, _, _, left, _, right, _, _, _) = struct.unpack("hhhhHhhhhhh", csbi.raw)
                    return right - left
                else:
                    return self._default_width

            if hasattr(sys.stdout, 'isatty') and not sys.stdout.isatty():
                return self._default_width # leave PIPEs alone

            try:
                process = subprocess.Popen(['stty', '-a'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout = process.stdout.read()
                if sys.version_info[0] > 2:
                    stdout = stdout.decode("utf-8")
            except (OSError, IOError):
                pass
            else:
                # We support the following output formats from stty:
                #
                # 1) Linux   -> columns 80
                # 2) OS X    -> 80 columns
                # 3) Solaris -> columns = 80

                re_linux   = r"columns\s+(?P<columns>\d+);"
                re_osx     = r"(?P<columns>\d+)\s*columns;"
                re_solaris = r"columns\s+=\s+(?P<columns>\d+);"

                for regex in (re_linux, re_osx, re_solaris):
                    match = re.search(regex, stdout)

                    if match is not None:
                        columns = match.group('columns')

                        try:
                            return int(columns)
                        except ValueError:
                            pass

            return self._default_width

        width = findout_terminal_width()
        self._terminal_width = width

        return width


    def write(self, text, color="", align="left", width=None):
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

        if width is None:
            width = self.terminal_width

        if align == "right":
            if self._write_pos+len(text) > width:
                # we don't fit on the current line, create a new line
                self.write("\n")
            self.write(" "*(width-self._write_pos-len(text)))

        if hasattr(sys.stdout, 'isatty') and not sys.stdout.isatty():
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
        if not self._colors or (color == ""):
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
        width = self.terminal_width
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

    def start(self, seed=None):
        self.write_center("test process starts")
        executable = sys.executable
        v = tuple(sys.version_info)
        python_version = "%s.%s.%s-%s-%s" % v
        self.write("executable:   %s  (%s)\n" % (executable, python_version))
        self.write("architecture: %s\n" % ARCH)
        self.write("cache:        %s\n" % USE_CACHE)
        self.write("ground types: %s\n" % GROUND_TYPES)
        if seed is not None:
            self.write("random seed:  %d\n" % seed)
        self.write('\n')
        self._t_start = clock()

    def finish(self):
        self._t_end = clock()
        self.write("\n")
        global text, linelen
        text = "tests finished: %d passed, " % self._passed
        linelen = len(text)
        def add_text(mytext):
            global text, linelen
            """Break new text if too long."""
            if linelen + len(mytext) > self.terminal_width:
                text += '\n'
                linelen = 0
            text += mytext
            linelen += len(mytext)

        if len(self._failed) > 0:
            add_text("%d failed, " % len(self._failed))
        if len(self._failed_doctest) > 0:
            add_text("%d failed, " % len(self._failed_doctest))
        if self._skipped > 0:
            add_text("%d skipped, " % self._skipped)
        if self._xfailed > 0:
            add_text("%d expected to fail, " % self._xfailed)
        if len(self._xpassed) > 0:
            add_text("%d expected to fail but passed, " % len(self._xpassed))
        if len(self._exceptions) > 0:
            add_text("%d exceptions, " % len(self._exceptions))
        add_text("in %.2f seconds" % (self._t_end - self._t_start))


        if len(self._xpassed) > 0:
            self.write_center("xpassed tests", "_")
            for e in self._xpassed:
                self.write("%s: %s\n" % (e[0], e[1]))
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
                if self._first_only:
                    break
            self.write("\n")

        self.write_center(text)
        ok = len(self._failed) == 0 and len(self._exceptions) == 0 and \
                len(self._failed_doctest) == 0
        if not self.ok:
            self.write("DO *NOT* COMMIT!\n", "Red")
        return ok

    @property
    def ok(self):
        ok = len(self._failed) == 0 and len(self._exceptions) == 0 and \
                len(self._failed_doctest) == 0
        return ok

    def entering_filename(self, filename, n):
        rel_name = filename[len(self._root_dir)+1:]
        self._active_file = rel_name
        self._active_file_error = False
        self.write(rel_name)
        self.write("[%d] " % n)

    def leaving_filename(self):
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
        self.write("f", "Green")

    def test_xpass(self, v):
        if sys.version_info[:2] < (2, 6):
            message = getattr(v, 'message', '')
        else:
            message = str(v)
        self._xpassed.append((self._active_file, message))
        self.write("X", "Green")

    def test_fail(self, exc_info):
        self._failed.append((self._active_file, self._active_f, exc_info))
        self.write("F", "Red")
        self._active_file_error = True

    def doctest_fail(self, name, error_msg):
        # the first line contains "******", remove it:
        error_msg = "\n".join(error_msg.split("\n")[1:])
        self._failed_doctest.append((name, error_msg))
        self.write("F", "Red")
        self._active_file_error = True

    def test_pass(self):
        self._passed += 1
        if self._verbose:
            self.write("ok", "Green")
        else:
            self.write(".", "Green")

    def test_skip(self, v):
        if sys.version_info[:2] < (2, 6):
            message = getattr(v, 'message', '')
        else:
            message = str(v)
        self._skipped += 1
        self.write("s", "Green")
        if self._verbose:
            self.write(" - ", "Green")
            self.write(message, "Green")

    def test_exception(self, exc_info):
        self._exceptions.append((self._active_file, self._active_f, exc_info))
        self.write("E", "Red")
        self._active_file_error = True

    def import_error(self, filename, exc_info):
        self._exceptions.append((filename, None, exc_info))
        rel_name = filename[len(self._root_dir)+1:]
        self.write(rel_name)
        self.write("[?]   Failed to import", "Red")
        self.write(" ")
        self.write("[FAIL]", "Red", align="right")
        self.write("\n")
