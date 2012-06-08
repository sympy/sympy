from sympy import Expr, Add, Mul, Matrix, Pow, sympify, Matrix
from expr import Expr

def _is_scalar(e):
    """Convert from a sympy scalar to a Python scalar."""
    if isinstance(e, Expr):
        if (e.is_Integer or  e.is_Float or e.is_Rational or e.is_Number or
            e.is_NumberSymbol or e == I):
            return true
    return false

class Tr(Expr):
    """ Generic Trace operation than can trace over:

    a) sympy matrix
    b) operators
    c) outer products

    Parameters
    ==========
    o : operator, matrix, expr

    Examples
    ========

    """

    def __new__(cls, *args):
        """ Construct a Trace object. Return the following expr.

        a) Trace(A+B) = Tr(A) + Tr(B)
        b) Trace(scalar*Operator) = scalar*Trace(Operator)


        """
        expr = args[0]
        indices = args[1] if len(args) == 2 else -1 #-1 indicates full trace
        if isinstance(expr, Matrix):
            return expr.trace()
        elif isinstance(expr, Add):
            return Add(*[Tr(arg,indices) for arg in expr.args])
        elif isinstance(expr, Mul):
            c_part, nc_part = expr.args_cnc()
            if len(nc_part) == 0:
                return Mul(*c_part)
            else:
                return Mul(*c_part)*Expr.__new__(cls, Mul(*nc_part), indices)
        else:
            inst = Expr.__new__(cls, expr, indices)
            return inst

    def doit(self,**kwargs):
        """ Perform the trace operation.

        #TODO: Current version ignores the indices set for partial trace.

        """
        if hasattr(self.args[0], '_eval_trace'):
            return self.args[0]._eval_trace()

        return self
