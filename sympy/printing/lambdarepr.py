from printer import Printer

class LambdaPrinter(Printer):
    """
    This printer converts expressions into strings that can be used by
    lambdify.
    """

    tostr = lambda self, args, sep: sep.join([self._print(arg) for arg in args])

    def _print_Add(self, expr):
        return "(%s)"%self.tostr(expr.args, "+")

    def _print_Mul(self, expr):
        return "(%s)"%self.tostr(expr.args, "*")

    def _print_Pow(self, expr):
        return "(%s)"%self.tostr(expr.args, ")**(")

    def _print_Function(self, expr, exp=None):
        return expr.func.__name__ + "(%s)"%self.tostr(expr.args, ", ")

    def _print_list(self, expr):
        return "[%s]"%self.tostr(expr, ", ")

    def _print_tuple(self, expr):
        return "(%s,)"%self.tostr(expr, ", ")

    def _print_dict(self, expr):
        items = ["%s:%s"%(self._print(arg), self._print(arg)) for
                 key, value in expr.iteritems]
        return "{%s}"%", ".join(items)

    def _print_Matrix(self, expr):
        return "Matrix([%s])"%expr._format_str(self._print, ",")

def lambdarepr(expr):
    """
    Returns a string usable for lambdifying.
    """
    return LambdaPrinter().doprint(expr)