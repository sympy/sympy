from __future__ import annotations
import os
import logging
from pathlib import Path

from sympy.external import import_module

_lark = import_module("lark")

__doctest_requires__ = {('parse_smtlib',): ['lark']}


class SMTLibSyntaxError(Exception):
    """Raised when the input is not syntactically valid SMT-LIB."""
    def __init__(self, message, line=None, col=None):
        super().__init__(message)
        self.line = line
        self.col = col


class UnknownSMTLibCommandError(Exception):
    """Raised for SMT-LIB commands that the parser does not recognize."""


class UnknownSMTLibOperatorError(Exception):
    """Raised for operators or functions in a term that the parser does not
    recognize."""


class LarkSMTLibParser:
    r"""Class for converting input SMT-LIB strings into SymPy Expressions.
    It holds all the necessary internal data for doing so, and exposes hooks for
    customizing its behavior.
    """
    def __init__(self, print_debug_output=False, transform=True, grammar_file=None):
        grammar_dir_path = os.path.join(os.path.dirname(__file__), "grammar/")

        if grammar_file is None:
            grammar = Path(os.path.join(grammar_dir_path, "smtlib.lark")).read_text(encoding="utf-8")
        else:
            grammar = Path(grammar_file).read_text(encoding="utf-8")

        self.parser = _lark.Lark(
            grammar,
            source_path=grammar_dir_path,
            parser="lalr",
            start="start",
            lexer="auto",
            propagate_positions=False,
            maybe_placeholders=False,
            keep_all_tokens=False)

        self.print_debug_output = print_debug_output
        self.transform_expr = transform

    def doparse(self, s: str):
        prev_log_level = _lark.logger.level
        if self.print_debug_output:
            _lark.logger.setLevel(logging.DEBUG)

        try:
            try:
                parse_tree = self.parser.parse(s)
            except _lark.exceptions.LarkError as e:
                # str(e) already includes the line/column information
                raise SMTLibSyntaxError(f"Syntax error: {str(e)}",
                                        getattr(e, 'line', None),
                                        getattr(e, 'column', None)) from e

            if not self.transform_expr:
                return parse_tree

            try:
                from sympy.parsing.smtlib.lark.transformer import SMTLibTransformer
                transformer = SMTLibTransformer()
                transformer.transform(parse_tree)
                return transformer.get_result()
            except _lark.exceptions.VisitError as e:
                if isinstance(e.orig_exc, (UnknownSMTLibCommandError,
                                           UnknownSMTLibOperatorError,
                                           NotImplementedError)):
                    raise e.orig_exc from e
                raise
        finally:
            _lark.logger.setLevel(prev_log_level)


_lark_smtlib_parser = None


def parse_smtlib(s):
    """
    Parses a string in the SMT-LIB format and returns the declared symbols
    and the list of assertions.

    Parameters
    ==========

    s : str
        The SMT-LIB source text to parse. This must be the text itself; to
        parse a file, read its contents first, e.g. with
        ``Path(filename).read_text()``.

    Returns
    =======

    tuple
        A tuple ``(symbols, assertions)`` where ``symbols`` is a dictionary
        mapping symbol names to SymPy symbols, and ``assertions`` is a list
        of SymPy expressions representing the assertions.

    Examples
    ========

    >>> from sympy.parsing.smtlib import parse_smtlib
    >>> source = '''
    ... (declare-const x Int)
    ... (declare-const y Int)
    ... (assert (= (+ x y) 10))
    ... '''
    >>> symbols, assertions = parse_smtlib(source)
    >>> symbols['x']
    x
    >>> assertions
    [Q.integer(x), Q.integer(y), Eq(x + y, 10)]

    Notes
    =====

    Declaring a constant of sort ``Int`` or ``Real`` appends a corresponding
    assumption predicate (``Q.integer`` or ``Q.real``) to the returned
    ``assertions`` list, as shown in the example above. This is intentional:
    the sort information is part of what the input asserts about each symbol,
    and including it in the assertions makes it available to functions that
    consume them, such as :func:`~.satisfiable`. Constants of other sorts
    (including uninterpreted sorts) produce no assumption predicate.
    """
    global _lark_smtlib_parser
    if _lark is None:
        raise ImportError("Lark is not installed")
    if _lark_smtlib_parser is None:
        _lark_smtlib_parser = LarkSMTLibParser()

    return _lark_smtlib_parser.doparse(s)
