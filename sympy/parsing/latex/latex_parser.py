import operator
import os
import functools

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
        # TODO: should we use pkg_resource to get grammar file?  I
        # think this would make sympy depend on setuptools which we
        # would not like
        with open(os.path.join(os.path.dirname(__file__), 'latex.lark')) as f:
            latex_grammar = f.read()

        parser = _lark.Lark(latex_grammar, parser='lalr',
                            lexer='standard',
                            propagate_positions=False,
                            maybe_placeholders=False)
        Transformer = _lark.Transformer

    class TreeToSympy(Transformer):
        INTEGER = int
        FLOAT = float
        LETTER = sympy.Symbol

        def SYMBOL(self, token):
            symbol = token.value.replace('\\', '')
            if symbol == 'infty':
                return sympy.oo
            elif symbol == 'pi':
                return sympy.pi
            return sympy.Symbol(symbol)

        def CMD_MATHIT(self, token):
            return sympy.Symbol(token.value[8:-1])

        def group(self, args):
            print(args)
            if len(args) == 5: # e.g. \\left ( expr \right )
                return args[2]
            else: # e.g. ( expr )
                return args[1]

        def abs_group(self, args):
            print(args)
            if len(args) == 5: # e.g. \\left | expr \\right |
                return sympy.Abs(args[2])
            else: # e.g. | expr |
                return sympy.Abs(args[1])

        def normal_function(self, args):
            unary_func_map = {
                'FUNC_LOG': functools.partial(sympy.log, base=10),
                'FUNC_LN': functools.partial(sympy.log, base=sympy.E),
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
            if args[0].type in unary_func_map:
                _assert_nargs_func(args, 1)
                return unary_func_map[args[0].type](args[-2])

        def named_function(self, args):
            if hasattr(args[2], 'children'):
                return sympy.Function(str(args[0]))(*args[2].children)
            return sympy.Function(str(args[0]))(args[2])

        def implicit_mul(self, args):
            print('implicit_mul', args)
            result = args[0] * args[1]
            for i in range(len(args) - 2):
                result = result * args[i + 2]
            return result

        def unary_op(self, args):
            op_map = {
                'SUB': lambda expr: -expr,
                'ADD': lambda expr: expr
            }
            op, left = args[0].type, args[1]
            print('unary', op, left)
            return op_map[op](left)

        def factorial(self, args):
            # ensure we handle x!!!!
            result = sympy.factorial(args[0])
            for _ in range(len(args) - 2):
                result = sympy.factorial(result)
            return result

        def limit(self, args):
            sub = args[1].children
            if len(sub) == 9: # e.g. \\lim_{x \\to 3^{+}} a
                direction = '+' if sub[7].type == 'ADD' else '-'
                return sympy.Limit(args[2], sub[2], sub[4], dir=direction)
            else: # e.g. \\lim_{x \\to 3} a
                return sympy.Limit(args[2], sub[2], sub[4])

        def binary_op_1(self, args):
            print('op_1', args)
            return self._binary_op(args)

        def binary_op_2(self, args):
            print('op_2', args)
            return self._binary_op(args)

        def binary_op_3(self, args):
            left, op, right = args[0], args[1][0], args[1][1]
            print('op_3', op, left, right)
            return self._binary_op([left, op, right])

        def sup_expression(self, args):
            if len(args) == 2: # ^ atom form
                return args
            else: # ^ { atom } form
                return [args[0], args[2]]

        def _binary_op(self, args):
            op_map = {
                'ADD': operator.add,
                'SUB': operator.sub,
                'MUL': operator.mul,
                'CMD_TIMES': operator.mul,
                'CMD_CDOT': operator.mul,
                'DIV': operator.truediv,
                'CMD_DIV': operator.truediv,
                'CARET': pow
            }
            left, op, right = args[0], args[1].type, args[2]
            result = op_map[op](left, right)
            for i in range((len(args) - 3) // 2):
                left, op, right = result, args[i*2 + 3].type, args[i*2 + 4]
                result = op_map[op](left, right)
            return result

        def fraction(self, args):
            return args[2] * sympy.Pow(args[5], -1)

        def summation(self, args):
            if hasattr(args[1], 'children'):
                sub, sup = args[1].children, args[2]
            else:
                sup, sub = args[1], args[2].children
            return sympy.Sum(args[3], (sub[2], sub[4], sup[1]))

        def product(self, args):
            if hasattr(args[1], 'children'):
                sub, sup = args[1], args[2]
            else:
                sup, sub = args[1], args[2]
            return sympy.Product(args[3], (sub.children[2], sub.children[4], sup[1]))

        def sqrt(self, args):
            if len(args) == 7: # e.g \\sqrt [ 2 ] { a }
                return sympy.root(args[5], args[2])
            else: # e.g. \\sqrt { a }
                return sympy.sqrt(args[2])

        def binomial(self, args):
            return sympy.binomial(args[2], args[5])

        def relation(self, args):
            if len(args) == 1: # relation: expr
                return args[0]

            op_map = {
                'EQUAL': sympy.Eq,
                'LT': sympy.StrictLessThan,
                'LTE': sympy.LessThan,
                'GT': sympy.StrictGreaterThan,
                'GTE': sympy.GreaterThan,
            }
            return op_map[args[1].type](args[0], args[2])

    try:
        tree = parser.parse(s)
        # this could be done within the parse step of lark however we
        # would like to keep it seperate to allow for the possiblity of
        # multiple backends for the generated expression
        sympy_expression = TreeToSympy().transform(tree)
    except Exception as e:
        raise LaTeXParsingError(str(e))
    return sympy_expression
