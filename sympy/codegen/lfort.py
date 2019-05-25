import lfortran.asr.asr as asr
import lfortran.asr.builder as builder


class ASRConverter():

    def __init__(self):

        # The translation unit that will be modified during the conversion
        # process. Since this object is stateful, the intention is that each
        # instance of this class is used only once.
        self.translation_unit = builder.make_translation_unit()
        self.type_integer = builder.make_type_integer()

        # Variables that will be used in the body of a function declaration.
        self.variables = []

    def _convert(self, expr, **kwargs):
        convertmethod = '_convert_' + type(expr).__name__
        if hasattr(self, convertmethod):
            return getattr(self, convertmethod)(expr, **kwargs)
        else:
            raise NotImplementedError("Conversion method %s not implemented." %
                                      convertmethod)

    def _convert_arith(self, expr, op):
        "A generalized converter for binary arithmetic expressions."
        terms = list(expr.args)

        first_left = self._convert(terms[0])
        first_right = self._convert(terms[1])

        lf_node = builder.make_binop(first_left, op, first_right)

        for term in terms[2:]:
            converted_term = self._convert(term)
            lf_node = builder.make_binop(lf_node, op, converted_term)

        return lf_node

    # Only works with Add and Mul so far
    def _convert_Add(self, expr):
        return self._convert_arith(expr, asr.Add())

    def _convert_Mul(self, expr):
        return self._convert_arith(expr, asr.Mul())

    def _convert_Integer(self, expr):
        return asr.Num(expr.p, type=self.type_integer)

    _convert_One = _convert_Integer
    _convert_Zero = _convert_Integer

    def _convert_Symbol(self, expr):
        sym_name = expr.name

        # TODO: Symbols can also be assignments, in which case the intent
        # should not be "in"
        var = asr.Variable(name=sym_name, intent="in", type=self.type_integer)
        self.variables.append(var)
        return var

    def wrap_function_definition(self, expr):
        """Convert a sympy expression, wrapping it in a function definition."""

        converted_expr = self._convert(expr)
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
