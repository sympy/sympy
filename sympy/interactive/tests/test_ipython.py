"""Tests of tools for setting up interactive IPython sessions. """

from sympy.interactive.session import (init_ipython_session,
    enable_automatic_symbols, enable_automatic_int_sympification)

from sympy.core import Symbol, Rational, Integer
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

    # Check that built-in names aren't overridden
    app.run_cell("a = all == __builtin__.all", False)
    assert "all" not in app.user_ns
    assert app.user_ns['a'] == True

    # Check that sympy names aren't overridden
    app.run_cell("import sympy")
    app.run_cell("a = factorial == sympy.factorial")
    assert app.user_ns['a'] == True

def test_int_to_Integer():
    # XXX: Warning, don't test with == here.  0.5 == Rational(1, 2) is True!
    app = init_ipython_session()
    app.run_cell("a = 1")
    assert isinstance(app.user_ns['a'], int)

    enable_automatic_int_sympification(app)
    app.run_cell("a = 1/2")
    assert isinstance(app.user_ns['a'], Rational)
    app.run_cell("a = 1")
    assert isinstance(app.user_ns['a'], Integer)
    app.run_cell("a = int(1)")
    assert isinstance(app.user_ns['a'], int)
    app.run_cell("a = (1/\n2)")
    assert app.user_ns['a'] == Rational(1, 2)
    # TODO: How can we test that the output of a SyntaxError is the original
    # input, not the transformed input?
