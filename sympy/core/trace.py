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
    i : indices (optional)

    Examples
    ========

    #TODO: Need to handle printing

    a) Trace(A+B) = Tr(A) + Tr(B)
    b) Trace(scalar*Operator) = scalar*Trace(Operator)

    >>> from sympy.core.trace import Tr
    >>> from sympy import symbols, Matrix
    >>> a, b = symbols('a b', commutative=True)
    >>> A, B = symbols('A B', commutative=False)
    >>> Tr(a*A,2)
    a*Tr(A, 2)
    >>> m = Matrix([[1,2],[1,1]])
    >>> Tr(m)
    2

    """

    def __new__(cls, *args):
        """ Construct a Trace object. Return the following expr.

        """
        expr = args[0]
        indices = args[1] if len(args) == 2 else -1 #-1 indicates full trace
        if isinstance(expr, Matrix):
            return expr.trace()
        elif hasattr(expr, 'trace') and callable(t.x):
            #for any objects that have trace() defined e.g numpy
            return expr.trace()
        elif isinstance(expr, Add):
            return Add(*[Tr(arg, indices) for arg in expr.args])
        elif isinstance(expr, Mul):
            c_part, nc_part = expr.args_cnc()
            if len(nc_part) == 0:
                return Mul(*c_part)
            else:
                # cyclic permute nc_part for canonical ordering
                idx = nc_part.index(min(nc_part))
                nc_part_ordered = nc_part[idx:]
                nc_part_ordered.extend(nc_part[:idx])

                return Mul(*c_part) * Expr.__new__(cls, Mul(*nc_part_ordered),
                                                   indices)
        else:
            inst = Expr.__new__(cls, expr, indices)
            return inst

    def doit(self,**kwargs):
        """ Perform the trace operation.

        #TODO: Current version ignores the indices set for partial trace.

        >>> from sympy.core.trace import Tr
        >>> from sympy.physics.quantum.operator import OuterProduct
        >>> from sympy.physics.quantum.spin import JzKet, JzBra
        >>> t = Tr(OuterProduct(JzKet(1,1), JzBra(1,1)))
        >>> t.doit()
        1

        """
        if hasattr(self.args[0], '_eval_trace'):
            return self.args[0]._eval_trace()

        return self

    @property
    def is_number(self):
        #TODO : This function to be reviewed
        # and implementation improved.

        return True
