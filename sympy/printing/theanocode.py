from sympy.utilities import default_sort_key
from sympy.external import import_module

from sympy.printing.printer import Printer
import sympy
from functools import partial

theano = import_module('theano')
if theano:
    ts = theano.scalar
    tt = theano.tensor
    from theano import sandbox
    from theano.sandbox import linalg as tlinalg

    mapping = {
            sympy.Add: tt.add,
            sympy.Mul: tt.mul,
            sympy.Abs: tt.abs_,
            sympy.sign: tt.sgn,
            sympy.ceiling: tt.ceil,
            sympy.floor: tt.floor,
            sympy.log: tt.log,
            sympy.exp: tt.exp,
            sympy.sqrt: tt.sqrt,
            sympy.cos: tt.cos,
            sympy.acos: tt.arccos,
            sympy.sin: tt.sin,
            sympy.asin: tt.arcsin,
            sympy.tan: tt.tan,
            sympy.atan: tt.arctan,
            sympy.atan2: tt.arctan2,
            sympy.cosh: tt.cosh,
            sympy.acosh: tt.arccosh,
            sympy.sinh: tt.sinh,
            sympy.asinh: tt.arcsinh,
            sympy.tanh: tt.tanh,
            sympy.atanh: tt.arctanh,
            sympy.re: tt.real,
            sympy.im: tt.imag,
            sympy.arg: tt.angle,
            sympy.erf: tt.erf,
            sympy.gamma: tt.gamma,
            sympy.loggamma: tt.gammaln,
            sympy.Pow: tt.pow,
            sympy.Eq: tt.eq,
            sympy.Gt: tt.gt,
            sympy.Lt: tt.lt,
            sympy.Le: tt.le,
            sympy.Ge: tt.ge,
            sympy.Max: tt.maximum,  # Sympy accept >2 inputs, Theano only 2
            sympy.Min: tt.minimum,  # Sympy accept >2 inputs, Theano only 2

            # Matrices
            sympy.MatAdd: tt.Elemwise(ts.add),
            sympy.HadamardProduct: tt.Elemwise(ts.mul),
            sympy.Trace: tlinalg.trace,
            sympy.Determinant : tlinalg.det,
            sympy.Inverse: tlinalg.matrix_inverse,
            sympy.Transpose: tt.DimShuffle((False, False), [1, 0]),
    }

class TheanoPrinter(Printer):
    """ Code printer for Theano computations """
    printmethod = "_theano"

    cache = dict()

    def _print_Symbol(self, s, dtypes={}, broadcastables={}):
        dtype = dtypes.get(s, 'floatX')
        broadcastable = broadcastables.get(s, ())
        key = (s.name, dtype, broadcastable, type(s))
        if key in self.cache:
            return self.cache[key]
        else:
            value = tt.tensor(name=s.name, dtype=dtype, broadcastable=broadcastable)
            self.cache[key] = value
            return value

    def _print_Basic(self, expr, **kwargs):
        op = mapping[type(expr)]
        children = [self._print(arg, **kwargs) for arg in expr.args]
        return op(*children)

    def _print_Number(self, n, **kwargs):
        return eval(str(n))

    def _print_MatrixSymbol(self, X, dtypes={}, **kwargs):
        dtype = dtypes.get(X, 'floatX')
        # shape = [self._print(d, dtypes) for d in X.shape]
        key = (X.name, dtype, type(X))
        if key in self.cache:
            return self.cache[key]
        else:
            value = tt.Tensor(dtype, (False, False))(X.name)
            self.cache[key] = value
            return value

    def _print_MatMul(self, expr, **kwargs):
        children = [self._print(arg, **kwargs) for arg in expr.args]
        result = children[0]
        for child in children[1:]:
            result = tt.dot(result, child)
        return result

    def _print_MatrixSlice(self, expr, **kwargs):
        parent = self._print(expr.parent, **kwargs)
        rowslice = self._print(slice(*expr.rowslice), **kwargs)
        colslice = self._print(slice(*expr.colslice), **kwargs)
        return parent[rowslice, colslice]

    def _print_BlockMatrix(self, expr, **kwargs):
        nrows, ncols = expr.blocks.shape
        blocks = [[self._print(expr.blocks[r, c], **kwargs)
                        for c in range(ncols)]
                        for r in range(nrows)]
        return tt.join(0, *[tt.join(1, *row) for row in blocks])


    def _print_slice(self, expr, **kwargs):
        return slice(*[self._print(i, **kwargs)
                        if isinstance(i, sympy.Basic) else i
                        for i in (expr.start, expr.stop, expr.step)])

    def _print_Pi(self, expr, **kwargs):
        return 3.141592653589793

    def _print_Rational(self, expr, **kwargs):
        return tt.true_div(self._print(expr.p, **kwargs),
                           self._print(expr.q, **kwargs))

    def _print_Integer(self, expr, **kwargs):
        return expr.p

    def _print_factorial(self, expr, **kwargs):
        return self._print(sympy.gamma(expr.args[0] + 1), **kwargs)

    def _print_Derivative(self, deriv, **kwargs):
        rv = self._print(deriv.expr, **kwargs)
        for var in deriv.variables:
            var = self._print(var, **kwargs)
            rv = tt.Rop(rv, var, tt.ones_like(var))
        return rv

    def emptyPrinter(self, expr):
        return expr

    def doprint(self, expr, **kwargs):
        """Returns printer's representation for expr (as a string)"""
        return self._print(expr, **kwargs)


def theano_code(expr, **kwargs):
    return TheanoPrinter({}).doprint(expr, **kwargs)


def dim_handling(inputs, dim=None, dims={}, broadcastables={}, keys=()):
    """ Handle various input types for dimensions in tensor_wrap

    See Also:
        tensor_wrap
        theano_funciton
    """
    if dim:
        dims = dict(zip(inputs, [dim]*len(inputs)))
    if dims:
        maxdim = max(dims.values())
        broadcastables = dict((i, (False,)*dims[i] + (True,)*(maxdim-dims[i]))
                         for i in inputs)
    return broadcastables


def theano_function(inputs, outputs, dtypes={}, **kwargs):
    """ Create Theano function from SymPy expressions """
    broadcastables = dim_handling(inputs, **kwargs)
    code = partial(theano_code, dtypes=dtypes, broadcastables=broadcastables)
    tinputs  = map(code, inputs)
    toutputs = map(code, outputs)
    toutputs = toutputs[0] if len(toutputs) == 1 else toutputs
    return theano.function(tinputs, toutputs)
