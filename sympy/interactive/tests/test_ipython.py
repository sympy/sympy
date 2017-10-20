"""Tests of tools for setting up interactive IPython sessions. """

from sympy.interactive.session import (init_ipython_session,
    enable_automatic_symbols, enable_automatic_int_sympification,
    enable_automatic_rationalize)

from sympy.core import Symbol, Rational, Integer, Float
from sympy.external import import_module

# TODO: The code below could be made more granular with something like:
#
# @requires('IPython', version=">=0.11")
# def test_automatic_symbols(ipython):

ipython = import_module("IPython", min_module_version="0.11")

if not ipython:
    #bin/test will not execute any tests now
    disabled = True


def test_automatic_symbols():
    # NOTE: Because of the way the hook works, you have to use run_cell(code,
    # True).  This means that the code must have no Out, or it will be printed
    # during the tests.
    app = init_ipython_session()
    app.run_cell("from sympy import *")

    enable_automatic_symbols(app)

    symbol = "verylongsymbolname"
    assert symbol not in app.user_ns
    app.run_cell("a = %s" % symbol)
    assert symbol in app.user_ns
    app.run_cell("a = type(%s)" % symbol)
    assert app.user_ns['a'] == Symbol
    app.run_cell("%s = Symbol('%s')" % (symbol, symbol))
    assert symbol in app.user_ns

    # Check that built-in names aren't overridden
    app.run_cell("a = all == __builtin__.all")
    assert "all" not in app.user_ns
    assert app.user_ns['a'] is True

    # Check that sympy names aren't overridden
    app.run_cell("import sympy")
    app.run_cell("a = factorial == sympy.factorial")
    assert app.user_ns['a'] is True


def test_int_to_Integer():
    # XXX: Warning, don't test with == here.  0.5 == Rational(1, 2) is True!
    app = init_ipython_session()
    app.run_cell("from __future__ import division")
    app.run_cell("from sympy import Integer")
    app.run_cell("a = 1")
    assert isinstance(app.user_ns['a'], int)

    enable_automatic_int_sympification(app)
    app.run_cell("a = 1/2")
    assert isinstance(app.user_ns['a'], Rational)
    app.run_cell("a = 1")
    assert isinstance(app.user_ns['a'], Integer)
    app.run_cell("a = int(1)")
    assert isinstance(app.user_ns['a'], int)
    app.run_cell("a = int(1)/int(2)")
    assert isinstance(app.user_ns['a'], float)
    app.run_cell("a = (1/\n2)")
    assert isinstance(app.user_ns['a'], Rational)


def test_rationalize():
    app = init_ipython_session()
    app.run_cell("from sympy import Integer")
    app.run_cell("a = 0.5")
    assert isinstance(app.user_ns['a'], float)

    enable_automatic_rationalize(app)
    app.run_cell("a = 0.5")
    assert isinstance(app.user_ns['a'], Rational)
    app.run_cell("a = float(0.5)")
    assert isinstance(app.user_ns['a'], float)
    app.run_cell("a = Float(1.234567890123456)")
    assert isinstance(app.user_ns['a'], Float)
    assert app.user_ns['a'] == Float(1.234567890123456)
