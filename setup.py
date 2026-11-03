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

from setuptools import setup, Command


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
    'sympy.parsing.smtlib',
    'sympy.parsing.smtlib.lark',
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
    'sympy.polys.series',
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
    'sympy.polys.series.tests',
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


if __name__ == '__main__':
    setup(packages=['sympy'] + modules + tests,
          cmdclass={'test': test_sympy,
                    'antlr': antlr,
                    },
          )
