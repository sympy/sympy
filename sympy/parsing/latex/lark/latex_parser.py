import os
import logging
import re
import sympy

from sympy.external import import_module


class LaTeXParsingError(Exception):
    pass


class DummyTransformer:
    # This class is needed to properly handle the case when Lark could not be found,
    # because we need our custom TransformToSymPyExpr class to inherit from lark's
    # Transformer class. TODO: Write more descriptive comment.
    pass


_lark = import_module('lark')

if _lark is None:
    Transformer = DummyTransformer
    raise ImportError("Could not load 'lark'")
else:
    Transformer = _lark.Transformer


# noinspection PyPep8Naming,PyMethodMayBeStatic
class TransformToSymPyExpr(Transformer):
    INT = int
    FLOAT = float
    # TODO: Decide whether to use Python floats or SymPy floats (sympy.core.numbers.Float)
    SYMBOL = sympy.Symbol

    def SUBSCRIPTED_SYMBOL(self, token):
        symbol, sub = token.value.split('_')
        if sub.startswith('{'):
            return sympy.Symbol('%s_{%s}' % (symbol, sub[1:-1]))
        else:
            return sympy.Symbol('%s_{%s}' % (symbol, sub))

    def number(self, token):
        if "." not in token[0]:
            return int(token[0])
        else:
            pass # TODO: Handle this case later

    def latex_string(self, token):
        return token[0]

    def infinity(self, token):
        return sympy.oo

    def add(self, token):
        return sympy.Add(token[0], token[2], evaluate=False)

    def mul(self, token):
        return sympy.Mul(token[0], token[2], evaluate=False)

def parse_latex_lark(s: str, *, logger=False, print_debug_output=False, transform=True):
    # last option is temporary, for quick prototyping
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

    parse_tree = parser.parse(s)

    if not transform:
        # exit early and return the parse tree
        print("expression =", s)
        # print(parse_tree)
        print(parse_tree.pretty())
        return parse_tree

    sympy_expression = ""
    try:
        sympy_expression = TransformToSymPyExpr().transform(parse_tree)
    except Exception as e:
        raise LaTeXParsingError(str(e))

    if print_debug_output:
        print("expression =", s)
        print("SymPy expression =", sympy_expression)
        # print(parse_tree)
        print(parse_tree.pretty())

    return sympy_expression


def pretty_print_lark_trees(tree, indent=0, show_expr=True):
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
        output += "\n" + "\n".join([" "*new_indent + pretty_print_lark_trees(i, new_indent, show_expr) for i in tree.children])
    else:
        output += ",".join([pretty_print_lark_trees(i, new_indent, show_expr) for i in tree.children])

    if show_node:
        output += ")"

    return output


if __name__ == "__main__":
    # temporary, for sanity testing and catching errors in the lark grammar.
    # parse_latex_lark(r"\frac{1}{7\cdot 6} + 7", print_debug_output=True)
    # parse_latex_lark(r"1 + 1", print_debug_output=True)
    parse_latex_lark(r"1 * 1", print_debug_output=True)
