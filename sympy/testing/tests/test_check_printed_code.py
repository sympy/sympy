from __future__ import annotations

import os
import tempfile
from contextlib import contextmanager

import pytest

from sympy.codegen.ast import float80, real
from sympy.core.numbers import I
from sympy.core.relational import Gt
from sympy.core.symbol import Symbol, symbols
from sympy.functions.elementary.trigonometric import sin
from sympy.logic.boolalg import And, ITE
from sympy.matrices.expressions.matexpr import MatrixSymbol
from sympy.printing.c import C89CodePrinter, C99CodePrinter, ccode
from sympy.printing.codeprinter import CodePrinter, rust_code
from sympy.printing.octave import octave_code
from sympy.printing.rust import RustCodePrinter
from sympy.testing.check_printed_code import (PrintedCodeError, check_code, check_compiles,
                                              check_printed_code, compiler_problem,
                                              compiling_enabled,
                                              declarations, declare, wrap_in_function)
from sympy.testing.pytest import raises, skip

x = Symbol('x', real=True)
n = Symbol('n', integer=True)
b, z = symbols('b z')


def test_declare():
    assert declare(C99CodePrinter(), x, False) == ('x', 'real')
    assert declare(C99CodePrinter(), n, False) == ('n', 'int')
    assert declare(C99CodePrinter(), b, True) == ('b', 'bool')
    e = raises(PrintedCodeError, lambda: declare(C99CodePrinter(), z, False))
    assert str(e.value).startswith('The type of z is unknown')


def test_declarations():
    assert declarations(C99CodePrinter(), sin(x) + n, None) == (
        {n: ('n', 'int'), x: ('x', 'real')}, None)
    assert declarations(C99CodePrinter(), x + n, 'x') == (
        {n: ('n', 'int'), x: ('x', 'real')}, ('x', 'real'))
    assert declarations(C99CodePrinter(), And(b, Gt(x, 0)), Symbol('c')) == (
        {b: ('b', 'bool'), x: ('x', 'real')}, ('c', 'bool'))
    assert declarations(C99CodePrinter(), x + n, 'z') == (
        {n: ('n', 'int'), x: ('x', 'real')}, ('z', 'real'))
    assert declarations(C99CodePrinter(), 2*n, 'z')[1] == ('z', 'int')
    raises(PrintedCodeError, lambda: declarations(C99CodePrinter(), I*x, 'z'))
    m = Symbol('m', integer=True)
    assert declarations(RustCodePrinter(), ITE(Gt(n, 1), m, n), Symbol('r', integer=True)) == (
        {m: ('m', 'int'), n: ('n', 'int')}, ('r', 'int'))
    nested = ITE(Gt(n, 0), ITE(Gt(n, 1), m, n), ITE(Gt(n, 2), n, m))
    assert declarations(RustCodePrinter(), nested, Symbol('r', integer=True)) == (
        {m: ('m', 'int'), n: ('n', 'int')}, ('r', 'int'))
    raises(PrintedCodeError, lambda: declarations(RustCodePrinter(), nested, 'r'))
    assert declarations(RustCodePrinter(), ITE(b, b, Gt(x, 0)), 'c') == (
        {b: ('b', 'bool'), x: ('x', 'real')}, ('c', 'bool'))


def _problem(printer, expr, assign_to=None):
    """The message of the PrintedCodeError that check_code() raises, or None."""
    try:
        check_code(printer, expr, assign_to, printer.doprint(expr, assign_to))
    except PrintedCodeError as problem:
        return str(problem)
    return None


def test_check_code():
    assert _problem(C99CodePrinter(), sin(x) + n) is None
    assert _problem(RustCodePrinter(), And(b, Gt(x, 0)), Symbol('y')) is None
    assert _problem(C99CodePrinter(), x, 'x') is None
    assert _problem(C99CodePrinter(), 5) is None
    assert _problem(C99CodePrinter(), z + x).startswith('The type of z is unknown')
    assert _problem(C99CodePrinter(), x, 'z') is None
    assert _problem(C99CodePrinter(), I*x, 'z').startswith('The type of z is unknown')
    assert 'since it uses arrays' in _problem(C99CodePrinter(), MatrixSymbol('A', 2, 2)[0, 1])
    assert 'since it uses arrays' in _problem(C99CodePrinter(), x, MatrixSymbol('A', 1, 1)[0, 0])
    assert "since it's a list" in _problem(C99CodePrinter(), [x, x])
    assert 'setting user_functions' in _problem(
        C99CodePrinter({'user_functions': {'sin': 'mysin'}}), sin(x)).replace('\n', ' ')
    assert 'setting type_aliases' in _problem(
        C99CodePrinter({'type_aliases': {real: float80}}), x).replace('\n', ' ')
    assert 'setting human' in _problem(C99CodePrinter({'human': False}), z)


def test_wrap_in_function():
    assert wrap_in_function(C89CodePrinter, 'y = M_PI*x;', [('x', 'real')], ('y', 'real')) == (
        '#define _XOPEN_SOURCE 700\n#include <math.h>\n#include <stdlib.h>\n'
        '\nvoid f(double x)\n{\n    double y;\n    y = M_PI*x;\n}\n')
    assert wrap_in_function(RustCodePrinter, 'b && x > 0.0', [('b', 'bool'), ('x', 'real')], None
                  ).endswith('\npub fn f(b: bool, x: f64)\n{\n    let _ = {\n'
                             '        b && x > 0.0\n    };\n}\n')


@pytest.mark.thread_unsafe(reason='runs tests that replace CodePrinter.doprint')
def test_check_printed_code():
    doprint = CodePrinter.doprint

    class SubclassPrinter(C99CodePrinter):
        pass

    @check_printed_code
    def test_declared():
        assert rust_code(sin(x) + n) == 'n as f64 + x.sin()'
        assert octave_code(z) == 'z'  # not checked
        assert SubclassPrinter().doprint(z) == 'z'  # not checked
    test_declared()
    assert test_declared.__name__ == 'test_declared'
    assert [mark.name for mark in test_declared.pytestmark] == ['thread_unsafe']

    @check_printed_code
    def test_undeclared():
        assert ccode(z + 1) == 'z + 1'
    e = raises(PrintedCodeError, test_undeclared)
    assert str(e.value).startswith('C99CodePrinter printed\n\n    z + 1\n\nThe type of z')

    @check_printed_code
    def test_fails():
        assert ccode(x + 1) == 'wrong'
    e = raises(AssertionError, test_fails)
    assert type(e.value) is AssertionError

    @check_printed_code
    def test_prints_nothing():
        pass
    raises(PrintedCodeError, test_prints_nothing)

    @check_printed_code
    def test_unchecked():
        assert ccode(sin(x)) == 'sin(x)'
        assert ccode(sin(x), user_functions={'sin': 'sin'}) == 'sin(x)'
        assert ccode(z, human=False) == (set(), set(), 'z')
        ccode([x, x])
    e = raises(PrintedCodeError, test_unchecked)
    assert [problem in str(e.value).replace('\n', ' ') for problem in (
        'setting user_functions', 'setting human', "since it's a list")] == [True]*3
    assert CodePrinter.doprint is doprint


@pytest.mark.thread_unsafe(reason='changes PATH and SYMPY_CHECK_PRINTED_CODE')
def test_compiling_enabled():
    with _environment('SYMPY_CHECK_PRINTED_CODE', '0'):
        assert not compiling_enabled(C99CodePrinter)
    compiler_problem.cache_clear()
    try:
        with _environment('PATH', ''):
            with _environment('SYMPY_CHECK_PRINTED_CODE', ''):
                assert not compiling_enabled(C99CodePrinter)
            with _environment('SYMPY_CHECK_PRINTED_CODE', 'yes'):
                assert not compiling_enabled(C99CodePrinter)
            with _environment('SYMPY_CHECK_PRINTED_CODE', '1'):
                e = raises(RuntimeError, lambda: compiling_enabled(C99CodePrinter))
                assert str(e.value) == 'gcc is needed to check printed code, but it is not found'
    finally:
        compiler_problem.cache_clear()


@pytest.mark.thread_unsafe(reason='changes PATH')
def test_compiler_problem_cannot_run():
    if os.name == 'nt':
        skip('needs a shebang line')
    with tempfile.TemporaryDirectory() as directory:
        gcc = os.path.join(directory, 'gcc')
        with open(gcc, 'w') as file:
            file.write('#!/nonexistent/interpreter\n')
        os.chmod(gcc, 0o755)
        compiler_problem.cache_clear()
        try:
            with _environment('PATH', directory):
                assert "gcc is needed to check printed code, but it can't be run" in (
                    compiler_problem(C99CodePrinter))
        finally:
            compiler_problem.cache_clear()


@contextmanager
def _environment(name, value):
    old = os.environ.get(name)
    os.environ[name] = value
    try:
        yield
    finally:
        if old is None:
            del os.environ[name]
        else:
            os.environ[name] = old


def _wrap_and_compile(printer_class, code):
    check_compiles(printer_class, wrap_in_function(printer_class, code, [('x', 'real')], None))


def test_check_compiles_c():
    if not compiling_enabled(C99CodePrinter):
        skip("gcc isn't found or doesn't work")
    e = raises(PrintedCodeError, lambda: _wrap_and_compile(C99CodePrinter, 'x % x'))
    assert str(e.value).startswith("This doesn't compile.")
    _wrap_and_compile(C99CodePrinter, 'fmod(x, x)')
    _wrap_and_compile(C99CodePrinter, 'x + 0*1.0  /* \u03bb */')


def test_check_compiles_rust():
    if not compiling_enabled(RustCodePrinter):
        skip("rustc isn't found or doesn't work")
    e = raises(PrintedCodeError, lambda: _wrap_and_compile(RustCodePrinter, 'x + 42'))
    assert str(e.value).startswith("This doesn't compile.")
    _wrap_and_compile(RustCodePrinter, 'x + 42.0')

    # https://github.com/sympy/sympy/pull/30867 first asserted this code, which
    # rustc parses as the start of generic arguments for f64
    parameters = [('i', 'int'), ('j', 'int'), ('w', 'real')]
    raises(PrintedCodeError, lambda: check_compiles(RustCodePrinter, wrap_in_function(
        RustCodePrinter, '(i + j) as f64 < w', parameters, None)))
    check_compiles(RustCodePrinter, wrap_in_function(
        RustCodePrinter, '((i + j) as f64) < w', parameters, None))
