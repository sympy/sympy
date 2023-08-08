import os
import logging
import re
import sympy

from sympy.external import import_module


class LaTeXParsingError(Exception):
    pass


_lark = import_module("lark")

if _lark is not None:
    from lark import Transformer  # type: ignore
else:
    class Transformer:  # type: ignore
        def transform(self, *args):
            pass


# noinspection PyPep8Naming,PyMethodMayBeStatic
class TransformToSymPyExpr(Transformer):
    SYMBOL = sympy.Symbol
    DIGIT = int

    def GREEK_SYMBOL(self, tokens):
        # we omit the first character because it is a backslash. Also, if the variable name has "var" in it,
        # like "varphi" or "varepsilon", we remove that
        variable_name = re.sub("var", "", tokens[1:])
        if variable_name == "lambda":
            # we do the same name change as sympy.abc because lambda is a Python keyword
            return sympy.Symbol("lamda")
        else:
            return sympy.Symbol(variable_name)

    def BASIC_SUBSCRIPTED_SYMBOL(self, tokens):
        symbol, sub = tokens.value.split("_")
        if sub.startswith("{"):
            return sympy.Symbol("%s_{%s}" % (symbol, sub[1:-1]))
        else:
            return sympy.Symbol("%s_{%s}" % (symbol, sub))

    def GREEK_SUBSCRIPTED_SYMBOL(self, tokens):
        greek_letter, sub = tokens.value.split("_")

        greek_letter = re.sub("var", "", greek_letter[1:])
        if greek_letter == "lambda":
            # we do the same name change as sympy.abc because lambda is a Python keyword
            greek_letter = "lamda"

        if sub.startswith("{"):
            return sympy.Symbol("%s_{%s}" % (greek_letter, sub[1:-1]))
        else:
            return sympy.Symbol("%s_{%s}" % (greek_letter, sub))

    def SYMBOL_WITH_GREEK_SUBSCRIPT(self, tokens):
        symbol, sub = tokens.value.split("_")
        if sub.startswith("{"):
            greek_letter = sub[2:-1]
            greek_letter = re.sub("var", "", greek_letter)
            if greek_letter == "lambda":
                # we do the same name change as sympy.abc because lambda is a Python keyword
                greek_letter = "lamda"
            return sympy.Symbol("%s_{%s}" % (symbol, greek_letter))
        else:
            greek_letter = sub[1:]
            greek_letter = re.sub("var", "", greek_letter)
            if greek_letter == "lambda":
                # we do the same name change as sympy.abc because lambda is a Python keyword
                greek_letter = "lamda"
            return sympy.Symbol("%s_{%s}" % (symbol, greek_letter))

    def multiletter_symbol(self, tokens):
        return sympy.Symbol(tokens[2])

    def number(self, tokens):
        if "." in tokens[0]:
            return sympy.core.numbers.Float(tokens[0])
        else:
            return int(tokens[0])

    def latex_string(self, tokens):
        return tokens[0]

    def infinity(self, tokens):
        return sympy.oo

    def group_round_parentheses(self, tokens):
        return tokens[1]

    def group_square_brackets(self, tokens):
        return tokens[1]

    def group_curly_parentheses(self, tokens):
        return tokens[1]

    def relation(self, tokens):
        relation_type = tokens[1].type
        if relation_type == "EQUAL":
            return sympy.Eq(tokens[0], tokens[2])
        if relation_type == "NOT_EQUAL":
            return sympy.Ne(tokens[0], tokens[2])
        if relation_type == "LT":
            return sympy.Lt(tokens[0], tokens[2])
        if relation_type == "LTE":
            return sympy.Le(tokens[0], tokens[2])
        if relation_type == "GT":
            return sympy.Gt(tokens[0], tokens[2])
        if relation_type == "GTE":
            return sympy.Ge(tokens[0], tokens[2])

    def add(self, tokens):
        return sympy.Add(tokens[0], tokens[2], evaluate=False)

    def sub(self, tokens):
        if len(tokens) == 2:
            return -tokens[1]
        elif len(tokens) == 3:
            return sympy.Add(tokens[0], -tokens[2], evaluate=False)

    def mul(self, tokens):
        if len(tokens) == 2:
            return sympy.Mul(tokens[0], tokens[1], evaluate=False)
        elif len(tokens) == 3:
            return sympy.Mul(tokens[0], tokens[2], evaluate=False)
        else:
            raise LaTeXParsingError() # TODO: fill out descriptive error message

    def div(self, tokens):
        return sympy.Mul(tokens[0], sympy.Pow(tokens[2], -1, evaluate=False), evaluate=False)

    def superscript(self, tokens):
        return sympy.Pow(tokens[0], tokens[2], evaluate=False)

    def fraction(self, tokens):
        return sympy.Mul(tokens[1], sympy.Pow(tokens[2], -1, evaluate=False), evaluate=False)

    def binomial(self, tokens):
        return sympy.binomial(tokens[1], tokens[2], evaluate=False)

    def integral(self, tokens):
        underscore_index = None
        caret_index = None

        if "_" in tokens:
            # we need to know the index because the next item in the list is the
            # arguments for the lower bound of the integral
            underscore_index = tokens.index("_")

        if "^" in tokens:
            # we need to know the index because the next item in the list is the
            # arguments for the upper bound of the integral
            caret_index = tokens.index("^")

        lower_bound = tokens[underscore_index + 1] if underscore_index else None
        upper_bound = tokens[caret_index + 1] if caret_index else None

        if "d" in tokens:
            differential_symbol_index = tokens.index("d")
        elif r"\text{d}" in tokens:
            differential_symbol_index = tokens.index(r"\text{d}")
        elif r"\mathrm{d}" in tokens:
            differential_symbol_index = tokens.index(r"\mathrm{d}")
        else:
            # differential symbol was not found
            raise LaTeXParsingError() # TODO: fill out descriptive error message

        differential_symbol = tokens[differential_symbol_index + 1]

        if (lower_bound is not None and upper_bound is None) or (upper_bound is not None and lower_bound is None):
            # then one was given and the other wasn't
            
            # we can"t simply do something like `if (lower_bound and not upper_bound) ...` because this would evaluate
            # to True if the lower_bound is 0
            raise LaTeXParsingError() # TODO: fill out descriptive error message

        # check if any expression was given or not. If it wasn't, then set the integrand to 1.
        if underscore_index is not None and underscore_index == differential_symbol_index - 2:
            integrand = 1
        elif caret_index is not None and caret_index == differential_symbol_index - 2:
            integrand = 1
        elif differential_symbol_index == 1:
            # this means we have something like \int dx, because the \int symbol will always be
            # at index 0 in `tokens`
            integrand = 1
        else:
            integrand = tokens[differential_symbol_index - 1]

        if lower_bound is not None:
            # we can assume that either both the lower and upper bounds are given, or
            # neither of them are
            return sympy.Integral(integrand, (differential_symbol, lower_bound, upper_bound))

        return sympy.Integral(integrand, differential_symbol)

    def group_curly_parentheses_special(self, tokens):
        underscore_index = tokens.index("_")
        caret_index = tokens.index("^")

        # given the type of expressions we are parsing, we can assume that the lower limit
        # will always use braces around its arguments. This is because we don't support
        # converting unconstrained sums into SymPy expressions.

        # first we isolate the bottom limit
        left_brace_index = tokens.index("{", underscore_index)
        right_brace_index = tokens.index("}", underscore_index)

        bottom_limit = tokens[left_brace_index + 1: right_brace_index]

        # next, we isolate the upper limit
        top_limit = tokens[caret_index + 1:]

        # the code below will be useful for supporting things like `\sum_{n = 0}^{n = 5} n^2`
        # if "{" in top_limit:
        #     left_brace_index = tokens.index("{", caret_index)
        #     if left_brace_index != -1:
        #         # then there's a left brace in the string, and we need to find the closing right brace
        #         right_brace_index = tokens.index("}", caret_index)
        #         top_limit = tokens[left_brace_index + 1: right_brace_index]

        # print(f"top  limit = {top_limit}")

        index_variable = bottom_limit[0]
        lower_limit = bottom_limit[-1]
        upper_limit = top_limit[0]  # for now, the index will always be 0

        # print(f"return value = ({index_variable}, {lower_limit}, {upper_limit})")

        return index_variable, lower_limit, upper_limit

    def summation(self, tokens):
        return sympy.Sum(tokens[2], tokens[1])

    def product(self, tokens):
        return sympy.Product(tokens[2], tokens[1])

    def limit_dir_expr(self, tokens):
        caret_index = tokens.index("^")

        if "{" in tokens:
            left_curly_brace_index = tokens.index("{", caret_index)
            direction = tokens[left_curly_brace_index + 1]
        else:
            direction = tokens[caret_index + 1]

        if direction == "+":
            return tokens[0], "+"
        elif direction == "-":
            return tokens[0], "-"
        else:
            return tokens[0], "+-"

    def group_curly_parentheses_lim(self, tokens):
        limit_variable = tokens[1]
        if isinstance(tokens[3], tuple):
            destination, direction = tokens[3]
        else:
            destination = tokens[3]
            direction = "+-"

        return limit_variable, destination, direction

    def limit(self, tokens):
        limit_variable, destination, direction = tokens[2]

        return sympy.Limit(tokens[-1], limit_variable, destination, direction)

    def list_of_expressions(self, tokens):
        if len(tokens) == 1:
            # we return it like verbatim because the function_applied node expects
            # a list
            return tokens
        else:
            def remove_tokens(args):
                if isinstance(args, _lark.lexer.Token):
                    if args.type != "COMMA":
                        # an unexpected token was encountered
                        raise LaTeXParsingError()  # TODO: write descriptive error message
                    return False
                return True

            return filter(remove_tokens, tokens)

    # Function-related stuff
    def function_applied(self, tokens):
        return sympy.Function(tokens[0])(*tokens[2])

    def sin(self, tokens):
        return sympy.sin(tokens[1], evaluate=False)

    def cos(self, tokens):
        return sympy.cos(tokens[1], evaluate=False)

    def tan(self, tokens):
        return sympy.tan(tokens[1], evaluate=False)

    def csc(self, tokens):
        return sympy.csc(tokens[1], evaluate=False)

    def sec(self, tokens):
        return sympy.sec(tokens[1], evaluate=False)

    def cot(self, tokens):
        return sympy.cot(tokens[1], evaluate=False)

    def sin_power(self, tokens):
        exponent = tokens[2]
        if exponent == -1:
            return sympy.asin(tokens[-1], evaluate=False)
        else:
            return sympy.Pow(sympy.sin(tokens[-1], evaluate=False), exponent, evaluate=False)

    def cos_power(self, tokens):
        exponent = tokens[2]
        if exponent == -1:
            return sympy.acos(tokens[-1], evaluate=False)
        else:
            return sympy.Pow(sympy.cos(tokens[-1], evaluate=False), exponent, evaluate=False)

    def tan_power(self, tokens):
        exponent = tokens[2]
        if exponent == -1:
            return sympy.atan(tokens[-1], evaluate=False)
        else:
            return sympy.Pow(sympy.tan(tokens[-1], evaluate=False), exponent, evaluate=False)

    def csc_power(self, tokens):
        exponent = tokens[2]
        if exponent == -1:
            return sympy.acsc(tokens[-1], evaluate=False)
        else:
            return sympy.Pow(sympy.csc(tokens[-1], evaluate=False), exponent, evaluate=False)

    def sec_power(self, tokens):
        exponent = tokens[2]
        if exponent == -1:
            return sympy.asec(tokens[-1], evaluate=False)
        else:
            return sympy.Pow(sympy.sec(tokens[-1], evaluate=False), exponent, evaluate=False)

    def cot_power(self, tokens):
        exponent = tokens[2]
        if exponent == -1:
            return sympy.acot(tokens[-1], evaluate=False)
        else:
            return sympy.Pow(sympy.cot(tokens[-1], evaluate=False), exponent, evaluate=False)

    def arcsin(self, tokens):
        return sympy.asin(tokens[1], evaluate=False)

    def arccos(self, tokens):
        return sympy.acos(tokens[1], evaluate=False)

    def arctan(self, tokens):
        return sympy.atan(tokens[1], evaluate=False)

    def arccsc(self, tokens):
        return sympy.acsc(tokens[1], evaluate=False)

    def arcsec(self, tokens):
        return sympy.asec(tokens[1], evaluate=False)

    def arccot(self, tokens):
        return sympy.acot(tokens[1], evaluate=False)

    def sinh(self, tokens):
        return sympy.sinh(tokens[1], evaluate=False)

    def cosh(self, tokens):
        return sympy.cosh(tokens[1], evaluate=False)

    def tanh(self, tokens):
        return sympy.tanh(tokens[1], evaluate=False)

    def asinh(self, tokens):
        return sympy.asinh(tokens[1], evaluate=False)

    def acosh(self, tokens):
        return sympy.acosh(tokens[1], evaluate=False)

    def atanh(self, tokens):
        return sympy.atanh(tokens[1], evaluate=False)

    def abs(self, tokens):
        return sympy.Abs(tokens[1], evaluate=False)

    def floor(self, tokens):
        return sympy.floor(tokens[1], evaluate=False)

    def ceil(self, tokens):
        return sympy.ceiling(tokens[1], evaluate=False)

    def factorial(self, tokens):
        return sympy.factorial(tokens[0], evaluate=False)

    def conjugate(self, tokens):
        return sympy.conjugate(tokens[1], evaluate=False)

    def square_root(self, tokens):
        if len(tokens) == 2:
            # then there was no square bracket argument
            return sympy.sqrt(tokens[1], evaluate=False)
        elif len(tokens) == 3:
            # then there _was_ a square bracket argument
            return sympy.root(tokens[2], tokens[1], evaluate=False)

    def exponential(self, tokens):
        return sympy.exp(tokens[1], evaluate=False)

    def log(self, tokens):
        if tokens[0].type == "FUNC_LG":
            # we don't need to check if there's an underscore or not because having one
            # in this case would be meaningless
            # TODO: ANTLR refers to ISO 80000-2:2019. should we keep base 10 or base 2?
            return sympy.log(tokens[1], 10, evaluate=False)
        elif tokens[0].type == "FUNC_LN":
            return sympy.log(tokens[1], evaluate=False)
        elif tokens[0].type == "FUNC_LOG":
            # we check if a base was passed in or not
            if "_" in tokens:
                # then a base was specified
                return sympy.log(tokens[3], tokens[2], evaluate=False) # fix the arguments
            else:
                # a base was not specified
                return sympy.log(tokens[1], evaluate=False)


class LarkLatexParser:
    def __init__(self, logger=False, print_debug_output=False, transform=True):
        import lark

        with open(os.path.join(os.path.dirname(__file__), "latex.lark"), encoding="utf-8") as f:
            latex_grammar = f.read()

        self.parser = lark.Lark(
            latex_grammar, parser="earley", start="latex_string",
            lexer="auto",
            ambiguity="explicit",
            debug=True,
            propagate_positions=False,
            maybe_placeholders=False,
            keep_all_tokens=True)

        self.logger = logger
        self.print_debug_output = print_debug_output
        self.transform = transform

        self.transform_to_sympy_expr = TransformToSymPyExpr()

    def doparse(self, s: str):
        import lark

        if self.logger:
            lark.logger.setLevel(logging.DEBUG)

        parse_tree = self.parser.parse(s)

        if not self.transform:
            # exit early and return the parse tree
            lark.logger.debug("expression =", s)
            lark.logger.debug(parse_tree)
            lark.logger.debug(parse_tree.pretty())
            return parse_tree

        if self.print_debug_output:
            # print this stuff before attempting to run the transformer
            lark.logger.debug("expression =", s)
            # print(parse_tree)
            lark.logger.debug(parse_tree.pretty())

        sympy_expression = self.transform_to_sympy_expr.transform(parse_tree)

        if self.print_debug_output:
            lark.logger.debug("SymPy expression =", sympy_expression)

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


def pr_ltx(s: str):
    LarkLatexParser(print_debug_output=True, transform=False).doparse(s)


def trfm_ltx(s: str):
    LarkLatexParser(print_debug_output=True).doparse(s)


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
