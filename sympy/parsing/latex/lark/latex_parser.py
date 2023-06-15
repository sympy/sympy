import operator
import os
import functools

import sympy
from sympy.external import import_module


class LaTeXParsingError(Exception):
    pass


def parse_latex_lark(s: str):
    # default behavior is to regenerate grammar from `LaTeX.lark` if
    # `lark` is not importable it falls back to using the standalone
    # precompiled lark latex grammar if it exists

    # TODO: move the import_module call outside of the function so that all that
    # code isn't run every time the function is called.

    # TODO: if lark isn't found, just define a dummy transformer class, and don't do
    # any type of parsing.
    _lark = import_module('lark')
    if _lark is None:
        # TODO: Emit appropriate error message if Lark module not found.
        raise ImportError("Could not load 'lark'")

    # TODO: should we use pkg_resource to get grammar file?  I
    # think this would make sympy depend on setuptools which we
    # would not like
    with open(os.path.join(os.path.dirname(__file__), 'latex.lark')) as f:
        latex_grammar = f.read()

    parser = _lark.Lark(latex_grammar, parser='earley', start='string',
                        lexer='auto',
                        ambiguity='explicit',
                        debug=True,
                        propagate_positions=False,
                        maybe_placeholders=False,
                        keep_all_tokens=True)

    string = parser.parse(s)
    print(string)
    print(string.pretty())

if __name__ == "__main__":
    parse_latex_lark(r"a \;\thickspace * b")

# with open('latex.lark', 'r') as f:
#     latex_parser = lark.Lark(f, start='string')
#
#     string = latex_parser.parse("a * b")
#     print(string)
#     print(string.pretty())
