from sympy.core.basic import Basic
from sympy.core.function import diff, Function, Subs, Derivative
from sympy.core.symbol import var, Symbol
from sympy.functions.elementary.exponential import exp, log
from sympy.functions.elementary.trigonometric import sin
from sympy.tensor.tensor import tensor_indices
from sympy.tensor.vtensor import VTensExpr, VTensorHead, vtensorhead, \
    VTensorIndexType
from sympy.core import sympify
from sympy.tensor.multiarray import MultiArray
from sympy.core.expr import Expr
from sympy.core.add import Add


class _DiffOpAdd(Expr):
    def _eval_diff(self):
        return Add(*[_._eval_diff() for _ in self.args])


class _DiffOp(Expr):

    # _op_priority = 3.0

    def __new__(cls, dvar, coeff=1):
        cls = Expr.__new__(cls, dvar, sympify(coeff))
        if len(dvar.args) == 2:
            coeff *= dvar.args[0]
            dvar = dvar.args[1]
        cls._args = (dvar, coeff,)
        return cls

    def _eval_diff(self):
        return diff(self.args[1], self.args[0])

    def __mul__(self, other):
        # if isinstance(other, MultiArray):
        #    return self.args[1] * other.applyfunc(lambda xv: diff(xv, self._args[0]))
        return _DiffOp(self.args[0], self.args[1] * other)

    def __rmul__(self, other):
        return _DiffOp(self.args[0], other * self.args[1])

    def __div__(self, other):
        return _DiffOp(self.args[0], self.args[1] / other)

    def __rdiv__(self, other):
        return _DiffOp(self.args[0], other / self.args[1])

    def __add__(self, other):
        return _DiffOpAdd(self, other)

    def __radd__(self, other):
        return _DiffOpAdd(other, self)

    def __sub__(self, other):
        return _DiffOpAdd(self, -other)

    def __rsub__(self, other):
        return _DiffOpAdd(other, -self)

    def __str__(self):
        return "diff(*, %s)" % (self.args[0])

    def __neg__(self):
        return _DiffOp(self.args[0], -self.args[1])


def tdiff(expr, dvar):
    """
    Tensor differentiation.

    Notes
    =====

    Keep in mind that the differentiating tensor index will be the first one of the resulting derivative,
    unless it is contracted.

    Tensor differentiation will care about contraction metric.

    Examples
    ========

    >>> from sympy import symbols, exp, log
    >>> from sympy.tensor.vtensor import VTensorIndexType, vtensorhead, tensor_indices
    >>> from sympy.tensor.tdiff import tdiff
    >>> L = VTensorIndexType('L', [1, -1])
    >>> i1, i2 = tensor_indices('i1, i2', L)
    >>> x, y = symbols('x, y')
    >>> A = vtensorhead('A', [L], [[1]], [x/y, log(x)*exp(y)])
    >>> dv = vtensorhead('d', [L], [[1]], [x, y])
    >>> tdiff(A(i1), dv(-i1))
    -exp(y)*log(x) + 1/y
    """
    if isinstance(dvar, VTensExpr):
        # only rank 0 or rank 1 tensor as deriving variable
        assert dvar.rank == 1

        dvar_index = dvar.abstract.free[0][0]
        dvar_dim = dvar._multiarray.dimensions[0]

        if isinstance(expr, Function):
            # TODO: identify scalar functions of tensors and use the chain rule.
            raise NotImplementedError("cannot derive a function.")
            for i, arg in enumerate(expr.args):
                dummy_var = Symbol("_xi_%i" % (i))
                return tdiff(arg, dvar) * Subs(Derivative(expr.func, dummy_var), (dummy_var,), (dvar,))
        else:
            _diffth = vtensorhead('partial', [dvar_index.tensortype], [[1]], values=[_DiffOp(dvar._multiarray[_i]) for _i in xrange(dvar_dim)])
            retu_expr = (_diffth(dvar_index) * expr)
            if isinstance(retu_expr, VTensExpr):
                return retu_expr.applyfunc(lambda x: x._eval_diff())
            else:
                return Add(*[_._eval_diff() for _ in retu_expr.atoms(_DiffOp)])

    else:
        return diff(expr, dvar)
