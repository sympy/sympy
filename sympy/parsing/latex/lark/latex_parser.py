import os

from sympy.external import import_module


class LaTeXParsingError(Exception):
    pass


# TODO: if lark isn't found, just define a dummy transformer class, and don't do
# any type of parsing.
_lark = import_module('lark')

if _lark is None:
    raise ImportError("Could not load 'lark'")


def parse_latex_lark(s: str):
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

    print("expression =", s)
    string = parser.parse(s)
    print(string)
    print(string.pretty())


if __name__ == "__main__":
    parse_latex_lark(r"a \;\thickspace * b")
