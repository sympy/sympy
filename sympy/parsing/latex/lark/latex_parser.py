import os
import logging
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
        symbol, sub = token.value.replace('\\', '').split('_')
        if sub.startswith('{'):
            return sympy.Symbol('%s_{%s}' % (symbol, sub[1:-1]))
        else:
            return sympy.Symbol('%s_{%s}' % (symbol, sub))


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
    print(string)
    print(string.pretty())
    return string


def mypretty(tree):
    data = tree.data
    if isinstance(data, _lark.Token):
        data = data.value
    output = str(data) + "(" + ", ".join([i for i in tree.children]) + ")"
    return output


if __name__ == "__main__":
    # temporary, for sanity testing and catching errors in the lark grammar.
    parse_latex_lark(r"\frac{1}{7\cdot 6} + 7")


# try:
#     s = transform_string(s)
#     print(s)
#     tree = parser.parse(s)
#     # this could be done within the parse step of lark however we
#     # would like to keep it separate to allow for the possibility of
#     # multiple backends for the generated expression
#     sympy_expression = TreeToSympy().transform(tree)
# except Exception as e:
#     raise LaTeXParsingError(str(e))
