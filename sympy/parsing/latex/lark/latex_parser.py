import os
import logging
import re

from sympy.external import import_module


class LaTeXParsingError(Exception):
    pass


# TODO: if lark isn't found, just define a dummy transformer class, and don't do
# any type of parsing.
_lark = import_module('lark')

if _lark is None:
    raise ImportError("Could not load 'lark'")


def parse_latex_lark(s: str, *, logger=False):
    # TODO: should we use pkg_resource to get grammar file?  I
    # think this would make sympy depend on setuptools which we
    # would not like
    with open(os.path.join(os.path.dirname(__file__), 'latex.lark')) as f:
        latex_grammar = f.read()

    parser = _lark.Lark(latex_grammar, parser='earley', start='latex_string',
                        lexer='auto',
                        ambiguity='explicit',
                        debug=True,
                        propagate_positions=False,
                        maybe_placeholders=False,
                        keep_all_tokens=True)

    if logger:
        _lark.logger.setLevel(logging.DEBUG)

    print("expression =", s)
    string = parser.parse(s)
    # print(string)
    # print(string.pretty())
    return string


def pretty_print_lark_trees(tree):
    if isinstance(tree, _lark.Token):
        return tree.value
    data = tree.data
    data = str(data)
    if data.startswith("expression"):
        data = re.sub(r"^expression", "E", data)
    output = str(data) + "(" + ", ".join([pretty_print_lark_trees(i) for i in tree.children]) + ")"
    return output


if __name__ == "__main__":
    # temporary, for sanity testing and catching errors in the lark grammar.
    parse_latex_lark(r"\frac{1 + 1}{7\cdot 6} + 7")
