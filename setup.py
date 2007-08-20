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
if sys.version_info[1] < 4:
    print "Sympy requires Python 2.4 or newer.  Python %d.%d detected" % \
          sys.version_info[:2]
    sys.exit(-1)


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
        try:
            import py as pylib
        except ImportError:
            print """In order to run the tests, you need codespeak's py.lib
            web page: http://codespeak.net/py/dist/
            If you are on debian systems, the package is named python-codespeak-lib
            """
            sys.exit(-1)
        pylib.test.cmdline.main(args=["sympy"])
        tdoc = test_sympy_doc(self.args)
        tdoc.run() # run also the doc test suite

class test_sympy_doc(Command):

    description = "Run the tests for the examples in the documentation"
    user_options = []  # distutils complains if this is not here.

    def initialize_options(self):  # distutils wants this
        pass

    def finalize_options(self):    # this too
        pass

    def run(self):

        import unittest
        import doctest

        import glob

        print "Testing docstrings."

        files = glob.glob('sympy/*/*.py')# + glob.glob('sympy/modules/*/*.py')
        files = [ f.replace("\\","/") for f in files ]  # make it work on Windows too

        # files without doctests or that don't work
        files.remove('sympy/core/add.py')
        files.remove('sympy/core/relational.py')
        files.remove('sympy/core/interval.py')

        files.remove('sympy/concrete/__init__.py')
        files.remove('sympy/integrals/__init__.py')
        files.remove('sympy/matrices/__init__.py')
        files.remove('sympy/physics/__init__.py')
        files.remove('sympy/printing/__init__.py')
        files.remove('sympy/series/__init__.py')
        files.remove('sympy/simplify/__init__.py')
        files.remove('sympy/solvers/__init__.py')
        files.remove('sympy/utilities/__init__.py')

        files.remove('sympy/specfun/zeta_functions.py')
        files.remove('sympy/specfun/orthogonal_polynomials.py')
        files.remove('sympy/specfun/factorials.py')

        # fix line 164 above

        files.remove('sympy/plotting/__init__.py')
        files.remove('sympy/plotting/plot_camera.py')
        files.remove('sympy/plotting/plot_curve.py')
        files.remove('sympy/plotting/plot_surface.py')
        files.remove('sympy/plotting/color_scheme.py')
        files.remove('sympy/plotting/managed_window.py')
        files.remove('sympy/plotting/plot_controller.py')
        files.remove('sympy/plotting/plot_rotation.py')
        files.remove('sympy/plotting/plot_window.py')
        files.remove('sympy/plotting/plot_mode.py')
        files.remove('sympy/plotting/plot_modes.py')
        files.remove('sympy/plotting/plot_mode_base.py')
        files.remove('sympy/plotting/plot_axes.py')
        files.remove('sympy/plotting/plot.py')

        try:
            #testing for optional libraries
            import libxslt
        except ImportError:
            pass

        modules = []

        for x in files:
            if len(x) > 12 and x[-11:] == '__init__.py':
                x = x.replace('/__init__', '')
                print x

            #put . as separator and strip the extension (.py)
            modules.append(x.replace('/', '.')[:-3])

        suite = unittest.TestSuite()

        for mod in modules:
            suite.addTest(doctest.DocTestSuite(mod))

        runner = unittest.TextTestRunner()
        runner.run(suite)

import sympy

tests = [
    'sympy.core.tests',
    'sympy.concrete.tests',
    'sympy.geometry.tests',
    'sympy.integrals.tests',
    'sympy.matrices.tests',
    'sympy.numerics.tests',
    'sympy.ntheory.tests',
    'sympy.polynomials.tests',
    'sympy.printing.tests',
    'sympy.series.tests',
    'sympy.simplify.tests',
    'sympy.solvers.tests',
    'sympy.specfun.tests',
    'sympy.utilities.tests',
    'sympy.plotting.tests',
        ]

setup(
      name = 'sympy',
      version = sympy.__version__,
      description = 'Computer algebra system (CAS) in Python',
      license = 'BSD',
      url = 'http://code.google.com/p/sympy',
      packages = ['sympy',
                    'sympy.core', 
                    'sympy.concrete', 
                    'sympy.geometry', 
                    'sympy.integrals', 
                    'sympy.matrices', 
                    'sympy.numerics', 
                    'sympy.ntheory', 
                    'sympy.physics', 
                    'sympy.polynomials', 
                    'sympy.polynomials.fast', 
                    'sympy.printing', 
                    'sympy.series', 
                    'sympy.simplify', 
                    'sympy.solvers', 
                    'sympy.specfun', 
                    'sympy.utilities', 

                    'sympy.plotting', 
                    'sympy.plotting.pyglet',
                    'sympy.plotting.pyglet.ext',
                    'sympy.plotting.pyglet.font',
                    'sympy.plotting.pyglet.gl',
                    'sympy.plotting.pyglet.image',
                    'sympy.plotting.pyglet.image.codecs',
                    'sympy.plotting.pyglet.media',
                    'sympy.plotting.pyglet.window',
                    'sympy.plotting.pyglet.window.carbon',
                    'sympy.plotting.pyglet.window.win32',
                    'sympy.plotting.pyglet.window.xlib',
                  ]+tests,
      scripts = ['bin/isympy'],
      ext_modules = [],
      data_files = [('share/man/man1', ['doc/man/isympy.1'])],
      cmdclass    = {'test': test_sympy,
                     'test_core' : test_sympy_core,
                     'test_doc' : test_sympy_doc,
                     'gen_doc' : gen_doc,
                     'clean' : clean,
                     },
      )
