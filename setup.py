#!/usr/bin/env python
"""Distutils based setup script for SymPy.

This uses Distutils (http://python.org/sigs/distutils-sig/) the standard
python mechanism for installing packages.  Optionally, you can use
Setuptools (http://pythonhosted.org/setuptools/setuptools.html) to automatically
handle dependencies.  For the easiest installation
just type the command (you'll probably need root privileges for that):

    python setup.py install

This will install the library in the default location. For instructions on
how to customize the install procedure read the output of:

    python setup.py --help install

In addition, there are some other commands:

    python setup.py clean -> will clean all trash (*.pyc and stuff)
    python setup.py test  -> will run the complete test suite
    python setup.py bench -> will run the complete benchmark suite
    python setup.py audit -> will run pyflakes checker on source code

To get a full list of avaiable commands, read the output of:

    python setup.py --help-commands

Or, if all else fails, feel free to write to the sympy list at
sympy@googlegroups.com and ask for help.
"""

import sys
import subprocess
import os
import shutil
import glob

mpmath_version = '0.19'

try:
    from setuptools import setup, Command
except ImportError:
    from distutils.core import setup, Command

    # handle mpmath deps in the hard way:
    from distutils.version import LooseVersion
    try:
        import mpmath
        if mpmath.__version__ < LooseVersion(mpmath_version):
            raise ImportError
    except ImportError:
        print("Please install the mpmath package with a version >= %s" % mpmath_version)
        sys.exit(-1)

PY3 = sys.version_info[0] > 2

# Make sure I have the right Python version.
if sys.version_info[:2] < (2, 6):
    print("SymPy requires Python 2.6 or newer. Python %d.%d detected" % sys.version_info[:2])
    sys.exit(-1)

# Check that this list is uptodate against the result of the command:
# for i in `find sympy -name __init__.py | rev | cut -f 2- -d '/' | rev | egrep -v "^sympy$" | egrep -v "tests$" `; do echo "'${i//\//.}',"; done | sort
modules = [
    'sympy.assumptions',
    'sympy.assumptions.handlers',
    'sympy.benchmarks',
    'sympy.calculus',
    'sympy.categories',
    'sympy.combinatorics',
    'sympy.concrete',
    'sympy.core',
    'sympy.core.benchmarks',
    'sympy.crypto',
    'sympy.deprecated',
    'sympy.diffgeom',
    'sympy.external',
    'sympy.functions',
    'sympy.functions.combinatorial',
    'sympy.functions.elementary',
    'sympy.functions.elementary.benchmarks',
    'sympy.functions.special',
    'sympy.functions.special.benchmarks',
    'sympy.galgebra',
    'sympy.geometry',
    'sympy.integrals',
    'sympy.integrals.benchmarks',
    'sympy.interactive',
    'sympy.liealgebras',
    'sympy.logic',
    'sympy.logic.algorithms',
    'sympy.logic.utilities',
    'sympy.matrices',
    'sympy.matrices.benchmarks',
    'sympy.matrices.expressions',
    'sympy.ntheory',
    'sympy.parsing',
    'sympy.physics',
    'sympy.physics.hep',
    'sympy.physics.mechanics',
    'sympy.physics.optics',
    'sympy.physics.quantum',
    'sympy.physics.unitsystems',
    'sympy.physics.unitsystems.systems',
    'sympy.physics.vector',
    'sympy.plotting',
    'sympy.plotting.intervalmath',
    'sympy.plotting.pygletplot',
    'sympy.polys',
    'sympy.polys.agca',
    'sympy.polys.benchmarks',
    'sympy.polys.domains',
    'sympy.printing',
    'sympy.printing.pretty',
    'sympy.series',
    'sympy.series.benchmarks',
    'sympy.sets',
    'sympy.simplify',
    'sympy.solvers',
    'sympy.solvers.benchmarks',
    'sympy.stats',
    'sympy.strategies',
    'sympy.strategies.branch',
    'sympy.tensor',
    'sympy.unify',
    'sympy.utilities',
    'sympy.utilities.mathml',
    'sympy.vector',
]

class audit(Command):
    """Audits SymPy's source code for following issues:
        - Names which are used but not defined or used before they are defined.
        - Names which are redefined without having been used.
    """

    description = "Audit SymPy source with PyFlakes"
    user_options = []

    def initialize_options(self):
        self.all = None

    def finalize_options(self):
        pass

    def run(self):
        import os
        try:
            import pyflakes.scripts.pyflakes as flakes
        except ImportError:
            print("In order to run the audit, you need to have PyFlakes installed.")
            sys.exit(-1)
        dirs = (os.path.join(*d) for d in (m.split('.') for m in modules))
        warns = 0
        for dir in dirs:
            for filename in os.listdir(dir):
                if filename.endswith('.py') and filename != '__init__.py':
                    warns += flakes.checkPath(os.path.join(dir, filename))
        if warns > 0:
            print("Audit finished with total %d warnings" % warns)


class clean(Command):
    """Cleans *.pyc and debian trashs, so you should get the same copy as
    is in the VCS.
    """

    description = "remove build files"
    user_options = [("all", "a", "the same")]

    def initialize_options(self):
        self.all = None

    def finalize_options(self):
        pass

    def run(self):
        dir_setup = os.path.dirname(os.path.realpath(__file__))
        curr_dir = os.getcwd()
        for root, dirs, files in os.walk(dir_setup):
            for file in files:
                if file.endswith('.pyc') and os.path.isfile:
                    os.remove(os.path.join(root, file))

        os.chdir(dir_setup)
        names = ["python-build-stamp-2.4", "MANIFEST", "build", "dist", "doc/_build", "sample.tex"]

        for f in names:
            if os.path.isfile(f):
                os.remove(f)
            elif os.path.isdir(f):
                shutil.rmtree(f)

        for name in glob.glob(os.path.join(dir_setup, "doc", "src", "modules", \
                                           "physics", "vector", "*.pdf")):
            if os.path.isfile(name):
                os.remove(name)

        os.chdir(curr_dir)


class test_sympy(Command):
    """Runs all tests under the sympy/ folder
    """

    description = "run all tests and doctests; also see bin/test and bin/doctest"
    user_options = []  # distutils complains if this is not here.

    def __init__(self, *args):
        self.args = args[0]  # so we can pass it to other classes
        Command.__init__(self, *args)

    def initialize_options(self):  # distutils wants this
        pass

    def finalize_options(self):    # this too
        pass

    def run(self):
        from sympy.utilities import runtests
        runtests.run_all_tests()


class run_benchmarks(Command):
    """Runs all SymPy benchmarks"""

    description = "run all benchmarks"
    user_options = []  # distutils complains if this is not here.

    def __init__(self, *args):
        self.args = args[0]  # so we can pass it to other classes
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
    'sympy.assumptions.tests',
    'sympy.calculus.tests',
    'sympy.categories.tests',
    'sympy.combinatorics.tests',
    'sympy.concrete.tests',
    'sympy.core.tests',
    'sympy.crypto.tests',
    'sympy.deprecated.tests',
    'sympy.diffgeom.tests',
    'sympy.external.tests',
    'sympy.functions.combinatorial.tests',
    'sympy.functions.elementary.tests',
    'sympy.functions.special.tests',
    'sympy.galgebra.tests',
    'sympy.geometry.tests',
    'sympy.integrals.tests',
    'sympy.interactive.tests',
    'sympy.liealgebras.tests',
    'sympy.logic.tests',
    'sympy.matrices.expressions.tests',
    'sympy.matrices.tests',
    'sympy.ntheory.tests',
    'sympy.parsing.tests',
    'sympy.physics.hep.tests',
    'sympy.physics.mechanics.tests',
    'sympy.physics.optics.tests',
    'sympy.physics.quantum.tests',
    'sympy.physics.tests',
    'sympy.physics.unitsystems.tests',
    'sympy.physics.vector.tests',
    'sympy.plotting.intervalmath.tests',
    'sympy.plotting.pygletplot.tests',
    'sympy.plotting.tests',
    'sympy.polys.agca.tests',
    'sympy.polys.domains.tests',
    'sympy.polys.tests',
    'sympy.printing.pretty.tests',
    'sympy.printing.tests',
    'sympy.series.tests',
    'sympy.sets.tests',
    'sympy.simplify.tests',
    'sympy.solvers.tests',
    'sympy.stats.tests',
    'sympy.strategies.branch.tests',
    'sympy.strategies.tests',
    'sympy.tensor.tests',
    'sympy.unify.tests',
    'sympy.utilities.tests',
    'sympy.vector.tests',
    ]

long_description = '''SymPy is a Python library for symbolic mathematics. It aims
to become a full-featured computer algebra system (CAS) while keeping the code
as simple as possible in order to be comprehensible and easily extensible.
SymPy is written entirely in Python and does not require any external libraries.'''

exec(open('sympy/release.py').read())
with open('sympy/__init__.py') as f:
    long_description = f.read().split('"""')[1]

setup(name='sympy',
      version=__version__,
      description='Computer algebra system (CAS) in Python',
      long_description=long_description,
      author='SymPy development team',
      author_email='sympy@googlegroups.com',
      license='BSD',
      keywords="Math CAS",
      url='http://sympy.org',
      packages=['sympy'] + modules + tests,
      scripts=['bin/isympy'],
      ext_modules=[],
      package_data={
          'sympy.utilities.mathml': ['data/*.xsl'],
          'sympy.logic.benchmarks': ['input/*.cnf'],
          },
      data_files=[('share/man/man1', ['doc/man/isympy.1'])],
      cmdclass={'test': test_sympy,
                'bench': run_benchmarks,
                'clean': clean,
                'audit': audit},
      classifiers=[
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        ],
      install_requires=['mpmath>=%s' % mpmath_version]
      )
