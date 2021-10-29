#!/usr/bin/env python
"""Distutils based setup script for SymPy.

This uses Distutils (https://python.org/sigs/distutils-sig/) the standard
python mechanism for installing packages. Optionally, you can use
Setuptools (https://setuptools.readthedocs.io/en/latest/)
to automatically handle dependencies. For the easiest installation
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

To get a full list of available commands, read the output of:

    python setup.py --help-commands

Or, if all else fails, feel free to write to the sympy list at
sympy@googlegroups.com and ask for help.
"""

import sys
import os
import shutil
import glob
import subprocess

from distutils.command.sdist import sdist


min_mpmath_version = '0.19'

# This directory
dir_setup = os.path.dirname(os.path.realpath(__file__))

extra_kwargs = {}

try:
    from setuptools import setup, Command
    extra_kwargs['zip_safe'] = False
    extra_kwargs['entry_points'] = {
        'console_scripts': [
            'isympy = isympy:main',
        ]
    }
except ImportError:
    from distutils.core import setup, Command

    extra_kwargs['scripts'] = ['bin/isympy']

    # handle mpmath deps in the hard way:
    from sympy.external.importtools import version_tuple
    try:
        import mpmath
        if version_tuple(mpmath.__version__) < version_tuple(min_mpmath_version):
            raise ImportError
    except ImportError:
        print("Please install the mpmath package with a version >= %s"
              % min_mpmath_version)
        sys.exit(-1)

if sys.version_info < (3, 6):
    print("SymPy requires Python 3.6 or newer. Python %d.%d detected"
          % sys.version_info[:2])
    sys.exit(-1)

# Check that this list is uptodate against the result of the command:
# python bin/generate_module_list.py
modules = [
    'sympy.algebras',
    'sympy.assumptions',
    'sympy.assumptions.handlers',
    'sympy.assumptions.predicates',
    'sympy.assumptions.relation',
    'sympy.benchmarks',
    'sympy.calculus',
    'sympy.categories',
    'sympy.codegen',
    'sympy.combinatorics',
    'sympy.concrete',
    'sympy.core',
    'sympy.core.benchmarks',
    'sympy.crypto',
    'sympy.diffgeom',
    'sympy.discrete',
    'sympy.external',
    'sympy.functions',
    'sympy.functions.combinatorial',
    'sympy.functions.elementary',
    'sympy.functions.elementary.benchmarks',
    'sympy.functions.special',
    'sympy.functions.special.benchmarks',
    'sympy.geometry',
    'sympy.holonomic',
    'sympy.integrals',
    'sympy.integrals.benchmarks',
    'sympy.integrals.rubi',
    'sympy.integrals.rubi.parsetools',
    'sympy.integrals.rubi.rubi_tests',
    'sympy.integrals.rubi.rules',
    'sympy.interactive',
    'sympy.liealgebras',
    'sympy.logic',
    'sympy.logic.algorithms',
    'sympy.logic.utilities',
    'sympy.matrices',
    'sympy.matrices.benchmarks',
    'sympy.matrices.expressions',
    'sympy.multipledispatch',
    'sympy.ntheory',
    'sympy.parsing',
    'sympy.parsing.autolev',
    'sympy.parsing.autolev._antlr',
    'sympy.parsing.c',
    'sympy.parsing.fortran',
    'sympy.parsing.latex',
    'sympy.parsing.latex._antlr',
    'sympy.physics',
    'sympy.physics.continuum_mechanics',
    'sympy.physics.control',
    'sympy.physics.hep',
    'sympy.physics.mechanics',
    'sympy.physics.optics',
    'sympy.physics.quantum',
    'sympy.physics.units',
    'sympy.physics.units.definitions',
    'sympy.physics.units.systems',
    'sympy.physics.vector',
    'sympy.plotting',
    'sympy.plotting.intervalmath',
    'sympy.plotting.pygletplot',
    'sympy.polys',
    'sympy.polys.agca',
    'sympy.polys.benchmarks',
    'sympy.polys.domains',
    'sympy.polys.matrices',
    'sympy.polys.numberfields',
    'sympy.printing',
    'sympy.printing.pretty',
    'sympy.sandbox',
    'sympy.series',
    'sympy.series.benchmarks',
    'sympy.sets',
    'sympy.sets.handlers',
    'sympy.simplify',
    'sympy.solvers',
    'sympy.solvers.benchmarks',
    'sympy.solvers.diophantine',
    'sympy.solvers.ode',
    'sympy.stats',
    'sympy.stats.sampling',
    'sympy.strategies',
    'sympy.strategies.branch',
    'sympy.tensor',
    'sympy.tensor.array',
    'sympy.tensor.array.expressions',
    'sympy.testing',
    'sympy.unify',
    'sympy.utilities',
    'sympy.utilities._compilation',
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
        curr_dir = os.getcwd()
        for root, dirs, files in os.walk(dir_setup):
            for file in files:
                if file.endswith('.pyc') and os.path.isfile:
                    os.remove(os.path.join(root, file))

        os.chdir(dir_setup)
        names = ["python-build-stamp-2.4", "MANIFEST", "build",
                 "dist", "doc/_build", "sample.tex"]

        for f in names:
            if os.path.isfile(f):
                os.remove(f)
            elif os.path.isdir(f):
                shutil.rmtree(f)

        for name in glob.glob(os.path.join(dir_setup, "doc", "src", "modules",
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


class antlr(Command):
    """Generate code with antlr4"""
    description = "generate parser code from antlr grammars"
    user_options = []  # distutils complains if this is not here.

    def __init__(self, *args):
        self.args = args[0]  # so we can pass it to other classes
        Command.__init__(self, *args)

    def initialize_options(self):  # distutils wants this
        pass

    def finalize_options(self):    # this too
        pass

    def run(self):
        from sympy.parsing.latex._build_latex_antlr import build_parser
        if not build_parser():
            sys.exit(-1)


class sdist_sympy(sdist):
    def run(self):
        # Fetch git commit hash and write down to commit_hash.txt before
        # shipped in tarball.
        commit_hash = None
        commit_hash_filepath = 'doc/commit_hash.txt'
        try:
            commit_hash = \
                subprocess.check_output(['git', 'rev-parse', 'HEAD'])
            commit_hash = commit_hash.decode('ascii')
            commit_hash = commit_hash.rstrip()
            print('Commit hash found : {}.'.format(commit_hash))
            print('Writing it to {}.'.format(commit_hash_filepath))
        except:
            pass

        if commit_hash:
            with open(commit_hash_filepath, 'w') as f:
                f.write(commit_hash)

        super(sdist_sympy, self).run()

        try:
            os.remove(commit_hash_filepath)
            print(
                'Successfully removed temporary file {}.'
                .format(commit_hash_filepath))
        except OSError as e:
            print("Error deleting %s - %s." % (e.filename, e.strerror))


# Check that this list is uptodate against the result of the command:
# python bin/generate_test_list.py
tests = [
    'sympy.algebras.tests',
    'sympy.assumptions.tests',
    'sympy.calculus.tests',
    'sympy.categories.tests',
    'sympy.codegen.tests',
    'sympy.combinatorics.tests',
    'sympy.concrete.tests',
    'sympy.core.tests',
    'sympy.crypto.tests',
    'sympy.diffgeom.tests',
    'sympy.discrete.tests',
    'sympy.external.tests',
    'sympy.functions.combinatorial.tests',
    'sympy.functions.elementary.tests',
    'sympy.functions.special.tests',
    'sympy.geometry.tests',
    'sympy.holonomic.tests',
    'sympy.integrals.rubi.parsetools.tests',
    'sympy.integrals.rubi.rubi_tests.tests',
    'sympy.integrals.rubi.tests',
    'sympy.integrals.tests',
    'sympy.interactive.tests',
    'sympy.liealgebras.tests',
    'sympy.logic.tests',
    'sympy.matrices.expressions.tests',
    'sympy.matrices.tests',
    'sympy.multipledispatch.tests',
    'sympy.ntheory.tests',
    'sympy.parsing.tests',
    'sympy.physics.continuum_mechanics.tests',
    'sympy.physics.control.tests',
    'sympy.physics.hep.tests',
    'sympy.physics.mechanics.tests',
    'sympy.physics.optics.tests',
    'sympy.physics.quantum.tests',
    'sympy.physics.tests',
    'sympy.physics.units.tests',
    'sympy.physics.vector.tests',
    'sympy.plotting.intervalmath.tests',
    'sympy.plotting.pygletplot.tests',
    'sympy.plotting.tests',
    'sympy.polys.agca.tests',
    'sympy.polys.domains.tests',
    'sympy.polys.matrices.tests',
    'sympy.polys.numberfields.tests',
    'sympy.polys.tests',
    'sympy.printing.pretty.tests',
    'sympy.printing.tests',
    'sympy.sandbox.tests',
    'sympy.series.tests',
    'sympy.sets.tests',
    'sympy.simplify.tests',
    'sympy.solvers.diophantine.tests',
    'sympy.solvers.ode.tests',
    'sympy.solvers.tests',
    'sympy.stats.sampling.tests',
    'sympy.stats.tests',
    'sympy.strategies.branch.tests',
    'sympy.strategies.tests',
    'sympy.tensor.array.expressions.tests',
    'sympy.tensor.array.tests',
    'sympy.tensor.tests',
    'sympy.testing.tests',
    'sympy.unify.tests',
    'sympy.utilities._compilation.tests',
    'sympy.utilities.tests',
    'sympy.vector.tests',
]


with open(os.path.join(dir_setup, 'sympy', 'release.py')) as f:
    # Defines __version__
    exec(f.read())


if __name__ == '__main__':
    setup(name='sympy',
          version=__version__,
          description='Computer algebra system (CAS) in Python',
          author='SymPy development team',
          author_email='sympy@googlegroups.com',
          license='BSD',
          keywords="Math CAS",
          url='https://sympy.org',
          py_modules=['isympy'],
          packages=['sympy'] + modules + tests,
          ext_modules=[],
          package_data={
              'sympy.utilities.mathml': ['data/*.xsl'],
              'sympy.logic.benchmarks': ['input/*.cnf'],
              'sympy.parsing.autolev': [
                  '*.g4', 'test-examples/*.al', 'test-examples/*.py',
                  'test-examples/pydy-example-repo/*.al',
                  'test-examples/pydy-example-repo/*.py',
                  'test-examples/README.txt',
                  ],
              'sympy.parsing.latex': ['*.txt', '*.g4'],
              'sympy.integrals.rubi.parsetools': ['header.py.txt'],
              'sympy.plotting.tests': ['test_region_*.png'],
              },
          data_files=[('share/man/man1', ['doc/man/isympy.1'])],
          cmdclass={'test': test_sympy,
                    'bench': run_benchmarks,
                    'clean': clean,
                    'audit': audit,
                    'antlr': antlr,
                    'sdist': sdist_sympy,
                    },
          python_requires='>=3.6',
          classifiers=[
            'License :: OSI Approved :: BSD License',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Topic :: Scientific/Engineering :: Physics',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3 :: Only',
            'Programming Language :: Python :: Implementation :: CPython',
            'Programming Language :: Python :: Implementation :: PyPy',
            ],
          install_requires=[
            'mpmath>=%s' % min_mpmath_version,
            ],
          **extra_kwargs
          )
