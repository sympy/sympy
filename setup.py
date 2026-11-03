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


if __name__ == '__main__':
    setup(cmdclass={'test': test_sympy, 'antlr': antlr})
