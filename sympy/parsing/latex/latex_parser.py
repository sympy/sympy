import operator

import sympy
from sympy.external import import_module


class LaTeXParsingError(Exception):
    pass


def _assert_nargs_func(args, n):
    if len(args) != 3 + n:
        raise LaTeXParsingError('latex function=%s expected %d args got %d' % (args[0].value, n, len(args)-3))


def parse_latex(s):
    r"""Converts the string ``s`` to a SymPy ``Expr``

    Parameters
    ==========

    s : str
        The LaTeX string to parse. In Python source containing LaTeX,
        *raw strings* (denoted with ``r"``, like this one) are preferred,
        as LaTeX makes liberal use of the ``\`` character, which would
        trigger escaping in normal Python strings.

    Examples
    ========

    >>> from sympy.parsing.latex import parse_latex
    >>> expr = parse_latex(r"\frac {1 + \sqrt {\a}} {\b}")
    >>> expr
    (sqrt(a) + 1)/b
    >>> expr.evalf(4, subs=dict(a=5, b=2))
    1.618
    """

    # default behavior is to regenerate grammar from `LaTeX.lark` if
    # `lark` is not importable it falls back to using the standalone
    # precompiled lark latex grammar if it exists
    _lark = import_module('lark')
    if _lark is None:
        # generated via lark standalone
        # python -m lark.tools.standalone latex.lark > latex_grammar.py
        standalone_lark = import_module('sympy.parsing.latex.latex_grammar')
        if standalone_lark is None:
            # TODO: throw proper error
            raise ValueError('precompiled lark grammar "latex_grammar" not avialable')
        parser = standalone_lark.Lark_StandAlone()
        Transformer = standalone_lark.Transformer
    else:
        # TODO: proper pkg_resource get text
        with open('latex.lark') as f:
            latex_grammar = f.read()

        parser = _lark.Lark(latex_grammar, parser='lalr',
                            lexer='standard',
                            propagate_positions=False,
                            maybe_placeholders=False)
        Transformer = _lark.Transformer

    class TreeToSympy(Transformer):
        INTEGER = int
        FLOAT = float
        SYMBOL = sympy.Symbol

        def group(self, args):
            return args[1]

        def abs_group(self, args):
             return sympy.Abs(args[1])

        def func(self, args):
            unary_funcs = {
                'FUNC_LOG': sympy.log,
                'FUNC_LN': sympy.log,
                'FUNC_SIN': sympy.sin,
                'FUNC_COS': sympy.cos,
                'FUNC_TAN': sympy.tan,
                'FUNC_CSC': sympy.csc,
                'FUNC_SEC': sympy.sec,
                'FUNC_COT': sympy.cot,
                'FUNC_ARCSIN': sympy.asin,
                'FUNC_ARCCOS': sympy.acos,
                'FUNC_ARCTAN': sympy.atan,
                'FUNC_ARCCSC': sympy.acsc,
                'FUNC_ARCSEC': sympy.asec,
                'FUNC_ARCCOT': sympy.acot,
                'FUNC_SINH': sympy.sinh,
                'FUNC_COSH': sympy.cosh,
                'FUNC_TANH': sympy.tanh,
                'FUNC_ARSINH': sympy.asinh,
                'FUNC_ARCOSH': sympy.acosh,
                'FUNC_ARTANH': sympy.atanh
            }
            if args[0].type in unary_funcs:
                _assert_nargs_func(args, 1)
                return unary_funcs[args[0].type](args[2])

        def expr(self, args):
            op_map = {
                'ADD': operator.add,
                'SUB': operator.sub,
                'MUL': operator.mul,
                'DIV': operator.truediv,
            }
            return op_map[args[1].type](args[0], args[2])

        def relation(self, args):
            op_map = {
                'EQUAL': sympy.Eq,
                'LT': sympy.StictLessThan,
                'LTE': sympy.LessThan,
                'GT': sympy.StictGreaterThan,
                'GTE': sympy.GreaterThan,
            }
            return op_map[args[1].type](args[0], args[2])

    tree = parser.parse(s)
    # this could be done within the parse step of lark however we
    # would like to keep it seperate to allow for the possiblity of
    # multiple backends for the generated expression
    sympy_expression = TreeToSympy().transform(tree)
    return sympy_expression
