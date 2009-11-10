#!/usr/bin/env python
"""Distutils based setup script for Sympy.

This uses Distutils (http://python.org/sigs/distutils-sig/) the standard
python mechanism for installing packages. For the easiest installation
just type the command (you'll probably need root privileges for that):

    python setup.py install

This will install the library in the default location. For instructions on
how to customize the install procedure read the output of:

    python setup.py --help install

In addition, there are some other commands:

    python setup.py clean -> will clean all trash (*.pyc and stuff)
    python setup.py test  -> will run the complete test suite
    python setup.py test_core -> will run only tests concerning core features
    python setup.py test_doc -> will run tests on the examples of the documentation
    python setup.py bench   -> will run the complete benchmark suite

To get a full list of avaiable commands, read the output of:

    python setup.py --help-commands

Or, if all else fails, feel free to write to the sympy list at
sympy@googlegroups.com and ask for help.
"""

from distutils.core import setup
from distutils.core import Command
import sys

import sympy

# Make sure I have the right Python version.
if sys.version_info[:2] < (2,4):
    print "Sympy requires Python 2.4 or newer. Python %d.%d detected" % \
          sys.version_info[:2]
    sys.exit(-1)

# Check that this list is uptodate against the result of the command (you can
# omit the thirdparty/ dir):
# $ find * -name __init__.py |sort
modules = [
    'sympy.assumptions',
    'sympy.concrete',
    'sympy.core',
    'sympy.functions',
    'sympy.functions.combinatorial',
    'sympy.functions.elementary',
    'sympy.functions.special',
    'sympy.galgebra',
    'sympy.geometry',
    'sympy.integrals',
    'sympy.interactive',
    'sympy.matrices',
    'sympy.ntheory',
    'sympy.parsing',
    'sympy.physics',
    'sympy.plotting',
    'sympy.queries',
    'sympy.thirdparty',
    'sympy.logic',
    'sympy.logic.algorithms',
    'sympy.logic.utilities',
    'sympy.mpmath',
    'sympy.polys',
    'sympy.printing',
    'sympy.printing.pretty',
    'sympy.refine',
    'sympy.series',
    'sympy.simplify',
    'sympy.solvers',
    'sympy.statistics',
    'sympy.utilities',
    'sympy.utilities.mathml',
    ]

class clean(Command):
    """Cleans *.pyc and debian trashs, so you should get the same copy as
    is in the svn.
    """

    description = "Clean everything"
    user_options = [("all","a","the same")]

    def initialize_options(self):
        self.all = None

    def finalize_options(self):
        pass

    def run(self):
        import os
        os.system("py.cleanup")
        os.system("rm -f python-build-stamp-2.4")
        os.system("rm -f MANIFEST")
        os.system("rm -rf build")
        os.system("rm -rf dist")

class gen_doc(Command):
    """Generate the (html) api documentation using epydoc

    output is sent to the directory ../api/
    """

    description = "generate the api doc"
    user_options = []

    target_dir = "../api/"

    def initialize_options(self):
        self.all = None

    def finalize_options(self):
        pass

    def run(self):
        import os
        os.system("epydoc --no-frames -o %s sympy" % self.target_dir)


class test_sympy_core(Command):
    """Run only the tests concerning features of sympy.core.
    It's a lot faster than running the complete test suite.
    """

    description = "Automatically run the core test suite for Sympy."
    user_options = []  # distutils complains if this is not here.

    def initialize_options(self):  # distutils wants this
        pass

    def finalize_options(self):    # this too
        pass


    def run(self):
        try:
            import py
        except ImportError:
            print """In order to run the tests, you need codespeak's py.lib
            web page: http://codespeak.net/py/dist/
            If you are on debian systems, the package is named python-codespeak-lib
            """
            sys.exit(-1)
        py.test.cmdline.main(args=["sympy/core/tests"])


class test_sympy(Command):
    """Runs all tests under the sympy/ folder
    """

    description = "Automatically run the test suite for Sympy."
    user_options = []  # distutils complains if this is not here.

    def __init__(self, *args):
        self.args = args[0] # so we can pass it to other classes
        Command.__init__(self, *args)

    def initialize_options(self):  # distutils wants this
        pass

    def finalize_options(self):    # this too
        pass

    def run(self):
        if sympy.test():
            # all regular tests run successfuly, so let's also run doctests
            # (if some regular test fails, the doctests are not run)
            sympy.doctest()

class test_sympy_doc(Command):

    description = "Run the tests for the examples in the documentation"
    user_options = []  # distutils complains if this is not here.

    def initialize_options(self):  # distutils wants this
        pass

    def finalize_options(self):    # this too
        pass

    def run(self):
        sympy.doctest()


class run_benchmarks(Command):
    """Runs all SymPy benchmarks"""

    description = "Automatically run the test suite for Sympy."
    user_options = []  # distutils complains if this is not here.

    def __init__(self, *args):
        self.args = args[0] # so we can pass it to other classes
        Command.__init__(self, *args)

    def initialize_options(self):  # distutils wants this
        pass

    def finalize_options(self):    # this too
        pass

    # we use py.test like architecture:
    #
    # o collector   -- collects benchmarks
    # o runner      -- executes benchmarks
    # o presenter   -- displays benchmarks results
    #
    # this is done in sympy.utilities.benchmarking on top of py.test
    def run(self):
        from sympy.utilities import benchmarking
        benchmarking.main(['sympy'])


# Check that this list is uptodate against the result of the command:
# $ python bin/generate_test_list.py
tests = [
    'sympy.concrete.tests',
    'sympy.core.tests',
    'sympy.functions.combinatorial.tests',
    'sympy.functions.elementary.tests',
    'sympy.functions.special.tests',
    'sympy.galgebra.tests',
    'sympy.geometry.tests',
    'sympy.integrals.tests',
    'sympy.logic.tests',
    'sympy.matrices.tests',
    'sympy.mpmath.tests',
    'sympy.ntheory.tests',
    'sympy.parsing.tests',
    'sympy.physics.tests',
    'sympy.plotting.tests',
    'sympy.polys.tests',
    'sympy.printing.pretty.tests',
    'sympy.printing.tests',
    'sympy.series.tests',
    'sympy.simplify.tests',
    'sympy.slow_tests',
    'sympy.solvers.tests',
    'sympy.statistics.tests',
    'sympy.test_external',
    'sympy.utilities.tests',
    ]

# update the following list from:
# http://pyglet.googlecode.com/svn/trunk/setup.py
# (whenever we update pyglet in sympy)
pyglet_packages=[
    'pyglet',
    'pyglet.gl',
    'pyglet.font',
    'pyglet.image',
    'pyglet.image.codecs',
    'pyglet.media',
    'pyglet.media.drivers',
    'pyglet.media.drivers.alsa',
    'pyglet.media.drivers.directsound',
    'pyglet.media.drivers.openal',
    'pyglet.window',
    'pyglet.window.carbon',
    'pyglet.window.win32',
    'pyglet.window.xlib',
]
pyglet_packages = ["sympy.thirdparty.pyglet." + s for s in pyglet_packages]

setup(
      name = 'sympy',
      version = sympy.__version__,
      description = 'Computer algebra system (CAS) in Python',
      license = 'BSD',
      url = 'http://code.google.com/p/sympy',
      packages = ['sympy'] + modules + tests + pyglet_packages,
      scripts = ['bin/isympy'],
      ext_modules = [],
      package_data = { 'sympy.utilities.mathml' : ['data/*.xsl'] },
      data_files = [('share/man/man1', ['doc/man/isympy.1'])],
      cmdclass    = {'test': test_sympy,
                     'test_core' : test_sympy_core,
                     'test_doc' : test_sympy_doc,
                     'bench'    : run_benchmarks,
                     'gen_doc' : gen_doc,
                     'clean' : clean,
                     },
      )

