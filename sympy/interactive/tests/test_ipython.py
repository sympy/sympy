"""Tests of tools for setting up interactive IPython sessions. """

from sympy.interactive.session import init_ipython_session, enable_automatic_symbols

from sympy.core import Symbol
from sympy.external import import_module
from sympy.utilities.pytest import raises

# TODO: The code below could be made more granular with something like:
#
# @requires('IPython', version=">=0.11")
# def test_automatic_symbols(ipython):

ipython = import_module("IPython", min_module_version="0.11")

if not ipython:
    #bin/test will not execute any tests now
    disabled = True

# TODO: Add tests that would verify that enable_automatic_symbols() doesn't
# break anything. For example typing `factorial` or `all` in an interpreter
# shouldn't result in a new symbol.
def test_automatic_symbols():
    app = init_ipython_session()
    app.run_cell("from sympy import *")

    enable_automatic_symbols(app)

    symbol = "verylongsymbolname"
    assert symbol not in app.user_ns
    app.run_cell(symbol, False)
    assert symbol in app.user_ns
    assert isinstance(app.user_ns[symbol], Symbol)
