from str import StrPrinter

class LambdaPrinter(StrPrinter):
    """
    This printer converts expressions into strings that can be used by
    lambdify.
    """

    def _print_Matrix(self, expr):
        return "Matrix([%s])"%expr._format_str(self._print, ",")


def lambdarepr(expr):
    """
    Returns a string usable for lambdifying.
    """
    return LambdaPrinter().doprint(expr)
