import os
import logging
import re

from sympy.external import import_module
from sympy.parsing.latex.lark.transformer import TransformToSymPyExpr

_lark = import_module("lark")


class LarkLatexParser:
    def __init__(self, logger=False, print_debug_output=False, transform=True):
        self.parser = _lark.Lark.open("latex.lark",
                                      rel_to=os.path.join(os.path.dirname(__file__), "grammar/"),
                                      parser="earley",
                                      start="latex_string",
                                      lexer="auto",
                                      ambiguity="explicit",
                                      debug=True,
                                      propagate_positions=False,
                                      maybe_placeholders=False,
                                      keep_all_tokens=True)

        self.logger = logger
        self.print_debug_output = print_debug_output
        self.transform_expr = transform

        self.transformer = TransformToSymPyExpr()

    def doparse(self, s: str):
        if self.logger:
            _lark.logger.setLevel(logging.DEBUG)

        parse_tree = self.parser.parse(s)

        if not self.transform_expr:
            # exit early and return the parse tree
            _lark.logger.debug("expression =", s)
            _lark.logger.debug(parse_tree)
            _lark.logger.debug(parse_tree.pretty())
            return parse_tree

        if self.print_debug_output:
            # print this stuff before attempting to run the transformer
            _lark.logger.debug("expression =", s)
            # print the `parse_tree` variable
            _lark.logger.debug(parse_tree.pretty())

        sympy_expression = self.transformer.transform(parse_tree)

        if self.print_debug_output:
            _lark.logger.debug("SymPy expression =", sympy_expression)

        return sympy_expression


if _lark is not None:
    _lark_latex_parser = LarkLatexParser()


def parse_latex_lark(s: str):
    """
    Experimental LaTeX parser using Lark.

    This function is still under development and its API may change with the
    next releases of SymPy.
    """
    if _lark is None:
        raise ImportError("Lark is probably not installed")
    return _lark_latex_parser.doparse(s)


def _pretty_print_lark_trees(tree, indent=0, show_expr=True):
    if isinstance(tree, _lark.Token):
        return tree.value

    data = str(tree.data)

    is_expr = data.startswith("expression")

    if is_expr:
        data = re.sub(r"^expression", "E", data)

    is_ambig = (data == "_ambig")

    if is_ambig:
        new_indent = indent + 2
    else:
        new_indent = indent

    output = ""
    show_node = not is_expr or show_expr

    if show_node:
        output += str(data) + "("

    if is_ambig:
        output += "\n" + "\n".join([" " * new_indent + _pretty_print_lark_trees(i, new_indent, show_expr) for i in tree.children])
    else:
        output += ",".join([_pretty_print_lark_trees(i, new_indent, show_expr) for i in tree.children])

    if show_node:
        output += ")"

    return output
