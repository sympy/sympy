from __future__ import annotations
import os
import logging
from pathlib import Path

from sympy.external import import_module

_lark = import_module("lark")

class SMTLibSyntaxError(Exception):
    def __init__(self, message, line, col):
        super().__init__(f"{message} at line {line}, column {col}")
        self.line = line
        self.col = col

class UnknownSMTLibCommandError(Exception):
    pass

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
        if self.print_debug_output:
            _lark.logger.setLevel(logging.DEBUG)

        try:
            parse_tree = self.parser.parse(s)
        except _lark.exceptions.LarkError as e:
            line = getattr(e, 'line', 0)
            col = getattr(e, 'column', 0)
            raise SMTLibSyntaxError("Syntax error", line, col)

        if not self.transform_expr:
            return parse_tree

        try:
            from sympy.parsing.smtlib.lark.transformer import SMTLibTransformer
            transformer = SMTLibTransformer()
            transformer.transform(parse_tree)
            return transformer.get_result()
        except _lark.exceptions.VisitError as e:
            if isinstance(e.orig_exc, UnknownSMTLibCommandError):
                raise e.orig_exc
            if isinstance(e.orig_exc, NotImplementedError):
                raise e.orig_exc
            raise

if _lark is not None:
    _lark_smtlib_parser = LarkSMTLibParser()

def parse_smtlib(s: str):
    """
    Experimental SMT-LIB parser using Lark.
    """
    if _lark is None:
        raise ImportError("Lark is probably not installed")
    return _lark_smtlib_parser.doparse(s)
