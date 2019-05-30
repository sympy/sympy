from sympy.printing.codeprinter import CodePrinter

import lfortran.asr.asr as asr
import lfortran.asr.builder as builder


class ASRConverter(CodePrinter):

    def __init__(self):

        # The translation unit that will be modified during the conversion
        # process. Since this object is stateful, the intention is that each
        # instance of this class is used only once.
        self.translation_unit = builder.make_translation_unit()
        self.type_integer = builder.make_type_integer()

        # Variables that will be used in the body of a function declaration.
        self.variables = []

    def _print(self, expr, **kwargs):
        convertmethod = '_print_' + type(expr).__name__
        if hasattr(self, convertmethod):
            return getattr(self, convertmethod)(expr, **kwargs)
        else:
            raise NotImplementedError("Conversion method %s not implemented." %
                                      convertmethod)

    def _print_arith(self, expr, op):
        "A generalized converter for binary arithmetic expressions."
        terms = list(expr.args)

        first_left = self._print(terms[0])
        first_right = self._print(terms[1])

        lf_node = builder.make_binop(first_left, op, first_right)

        for term in terms[2:]:
            converted_term = self._print(term)
            lf_node = builder.make_binop(lf_node, op, converted_term)

        return lf_node

    # Only works with Add and Mul so far
    def _print_Add(self, expr):
        return self._print_arith(expr, asr.Add())

    def _print_Mul(self, expr):
        return self._print_arith(expr, asr.Mul())

    def _print_Integer(self, expr):
        return asr.Num(expr.p, type=self.type_integer)

    _print_One = _print_Integer
    _print_Zero = _print_Integer

    def _print_Symbol(self, expr):
        sym_name = expr.name

        # TODO: Symbols can also be assignments, in which case the intent
        # should not be "in"
        var = asr.Variable(name=sym_name, intent="in", type=self.type_integer)
        self.variables.append(var)
        return var

    def wrap_function_definition(self, expr):
        """Convert a sympy expression, wrapping it in a function definition."""

        converted_expr = self._print(expr)
        f = builder.scope_add_function(self.translation_unit.global_scope, "f",
                                       return_var="ret")
        f.args = self.variables
        ret = f.return_var

        f.body.extend([asr.Assignment(ret, converted_expr)])
        return self.translation_unit


def sympy_to_asr(expr):
    """Convert a Sympy expression to an Lfortran expression, either standalone or
       wrapped in a function definition

    """

    converter = ASRConverter()
    return converter.wrap_function_definition(expr)
