from sympy.printing.codeprinter import CodePrinter

from sympy.external import import_module

lfortran = import_module("lfortran")
asr = lfortran.asr.asr
builder = lfortran.asr.builder


class ASRConverter(CodePrinter):

    def __init__(self):

        # The translation unit that will be modified during the conversion
        # process. Since this object is stateful, the intention is that each
        # instance of this class is used only once.
        self.translation_unit = builder.make_translation_unit()

        # Build the types that should be used in conversion of reals
        self.type_integer = builder.make_type_integer()
        self.type_real = builder.make_type_real()

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

    def _print_Float(self, expr):
        # Floating point numbers in LFortran are represented as strings.
        # LFortran supports the _dp suffix, so we default to double precision
        # for everything
        printed = CodePrinter._print_Float(self, expr)
        return asr.Num("%s_dp" % printed, type=self.type_real)

    def _print_Rational(self, expr):
        p = asr.Num("%d_dp" % expr.p, type=self.type_real)
        q = asr.Num("%d_dp" % expr.q, type=self.type_real)
        return builder.make_binop(p, asr.Div(), q)

    def _print_Pow(self, expr):
        # FIXME: Results in a type mismatch when the base contains a variable,
        # since variables default to integers
        if expr.exp == -1:
            one = asr.Num("1_dp", type=self.type_real)
            base = self._print(expr.base)
            return builder.make_binop(one, asr.Div(), base)
        else:
            exp = self._print(expr.exp)
            base = self._print(expr.base)
            return builder.make_binop(base, asr.Pow(), exp)

    def _print_Function(self, expr):
        args = [self._print(x) for x in expr.args]
        return asr.FuncCall(args=args)

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



def sympy_to_lfortran(expr):
    """Convert a SymPy expression to an LFortran ASR expression."""
    converter = ASRConverter()
    return converter._print(expr)


def sympy_to_lfortran_wrapped(expr):
    """Convert a Sympy expression to an Lfortran expression, wrapping it in a
       function definition."""
    converter = ASRConverter()
    return converter.wrap_function_definition(expr)
