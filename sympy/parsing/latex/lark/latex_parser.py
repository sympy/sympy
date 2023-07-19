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


def pretty_print_lark_trees(tree, indent=0, show_expr=True):
    if isinstance(tree, _lark.Token):
        return tree.value
    data = tree.data
    data = str(data)
    is_expr = data.startswith("expression")
    if is_expr:
        data = re.sub(r"^expression", "E", data)
    ambig = data == "_ambig"
    if ambig:
        new_indent = indent + 2
    else:
        new_indent = indent
    output = ""
    show_node = not is_expr or show_expr
    if show_node:
        output += str(data) + "("
    if ambig:
        output += "\n" + "\n".join([" "*new_indent + pretty_print_lark_trees(i, new_indent, show_expr) for i in tree.children])
    else:
        output += ",".join([pretty_print_lark_trees(i, new_indent, show_expr) for i in tree.children])
    if show_node:
        output += ")"
    return output


if __name__ == "__main__":
    # temporary, for sanity testing and catching errors in the lark grammar.
    parse_latex_lark(r"\frac{1 + 1}{7\cdot 6} + 7")
