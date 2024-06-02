#!/usr/bin/env python
"""Setup script for SymPy.

This uses Setuptools (https://setuptools.pypa.io/en/latest/) the standard
python mechanism for installing packages.
For the easiest installation just type the command (you'll probably need
root privileges for that):

    pip install .

This will install the library in the default location. For instructions on
how to customize the installation procedure read the output of:

    pip install --help

In addition, there are some other commands:

    python setup.py test  -> will run the complete test suite

To get a full list of available commands, read the output of:

    python setup.py --help-commands

Or, if all else fails, feel free to write to the sympy list at
sympy@googlegroups.com and ask for help.
"""

import sys
import os
import subprocess
from pathlib import Path

from setuptools import setup, Command
from setuptools.command.sdist import sdist


# This directory
dir_setup = os.path.dirname(os.path.realpath(__file__))

extra_kwargs = {
    'zip_safe': False,
    'entry_points': {
        'console_scripts': [
            'isympy = isympy:main',
        ]
    }
}

if sys.version_info < (3, 8):
    print("SymPy requires Python 3.8 or newer. Python %d.%d detected"
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
    'sympy.parsing.latex.lark',
    'sympy.physics',
    'sympy.physics.biomechanics',
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
    'sympy.plotting.backends',
    'sympy.plotting.backends.matplotlibbackend',
    'sympy.plotting.backends.textbackend',
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
    'sympy.utilities.mathml.data',
    'sympy.vector',
]


class test_sympy(Command):
    """Runs all tests under the sympy/ folder
    """

    description = "run all tests and doctests; also see bin/test and bin/doctest"
    user_options = []  # setuptools complains if this is not here.

    def __init__(self, *args):
        self.args = args[0]  # so we can pass it to other classes
        Command.__init__(self, *args)

    def initialize_options(self):  # setuptools wants this
        pass

    def finalize_options(self):    # this too
        pass

    def run(self):
        from sympy.testing import runtests
        runtests.run_all_tests()


class antlr(Command):
    """Generate code with antlr4"""
    description = "generate parser code from antlr grammars"
    user_options = []  # setuptools complains if this is not here.

    def __init__(self, *args):
        self.args = args[0]  # so we can pass it to other classes
        Command.__init__(self, *args)

    def initialize_options(self):  # setuptools wants this
        pass

    def finalize_options(self):    # this too
        pass

    def run(self):
        from sympy.parsing.latex._build_latex_antlr import build_parser as build_latex_parser
        if not build_latex_parser():
            sys.exit(-1)

        from sympy.parsing.autolev._build_autolev_antlr import build_parser as build_autolev_parser
        if not build_autolev_parser():
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

        super().run()

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
    'sympy.integrals.tests',
    'sympy.interactive.tests',
    'sympy.liealgebras.tests',
    'sympy.logic.tests',
    'sympy.matrices.expressions.tests',
    'sympy.matrices.tests',
    'sympy.multipledispatch.tests',
    'sympy.ntheory.tests',
    'sympy.parsing.tests',
    'sympy.physics.biomechanics.tests',
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
          long_description=(Path(__file__).parent / 'README.md').read_text("UTF-8"),
          long_description_content_type='text/markdown',
          author='SymPy development team',
          author_email='sympy@googlegroups.com',
          license='BSD',
          keywords="Math CAS",
          url='https://sympy.org',
          project_urls={
              'Source': 'https://github.com/sympy/sympy',
          },
          # Set upper bound when making the release branch.
          install_requires=[
              'mpmath >= 1.1.0',
          ],
          py_modules=['isympy'],
          packages=['sympy'] + modules + tests,
          ext_modules=[],
          package_data={
              'sympy.utilities.mathml.data': ['*.xsl'],
              'sympy.logic.benchmarks': ['input/*.cnf'],
              'sympy.parsing.autolev': [
                  '*.g4', 'test-examples/*.al', 'test-examples/*.py',
                  'test-examples/pydy-example-repo/*.al',
                  'test-examples/pydy-example-repo/*.py',
                  'test-examples/README.txt',
                  ],
              'sympy.parsing.latex': ['*.txt', '*.g4', 'lark/grammar/*.lark'],
              'sympy.plotting.tests': ['test_region_*.png'],
              'sympy': ['py.typed']
              },
          data_files=[('share/man/man1', ['doc/man/isympy.1'])],
          cmdclass={'test': test_sympy,
                    'antlr': antlr,
                    'sdist': sdist_sympy,
                    },
          python_requires='>=3.8',
          classifiers=[
            'License :: OSI Approved :: BSD License',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Topic :: Scientific/Engineering :: Physics',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: 3.10',
            'Programming Language :: Python :: 3.11',
            'Programming Language :: Python :: 3 :: Only',
            'Programming Language :: Python :: Implementation :: CPython',
            'Programming Language :: Python :: Implementation :: PyPy',
            ],
          extras_require={
              "dev": ["pytest>=7.1.0", "hypothesis>=6.70.0"],
            },
          **extra_kwargs
          )
