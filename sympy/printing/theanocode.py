from sympy.utilities import default_sort_key

from sympy.printing.printer import Printer
import theano
import sympy

Scalar = theano.scalar.Scalar
ts = theano.scalar
tt = theano.tensor
from theano import sandbox
from theano.sandbox import linalg as tlinalg

mapping = {sympy.Add: ts.add,
           sympy.Mul: ts.mul,
#           sympy.Sub: ts.sub,
#           sympy.Mul(x, Pow(y, -1)): ts.true_div,  # int_div,
           #sympy.Mod: ts.mod, # master branch
           #clip,  # second,
           #identity,
           #cast,
           sympy.Abs: ts.abs_,
           sympy.sign: ts.sgn,
           sympy.ceiling: ts.ceil,
           sympy.floor: ts.floor,
#           round_half_to_even, round_half_away_from_zero,
           #lambda x: sympy.Mul(x, -1): ts.neg,
           #lambda x: sympy.Pow(x, -1): ts.inv,
           sympy.log: ts.log,
           #lambda x: log(1 + p): ts.log1p,
           #log2, log10,
           sympy.exp: ts.exp,
           #exp2,
           #lambda x: sympy.Mul(x, x): ts.sqr,
           sympy.sqrt: ts.sqrt,
           sympy.cos: ts.cos,
           sympy.acos: ts.arccos,
           sympy.sin: ts.sin,
           sympy.asin: ts.arcsin,
           sympy.tan: ts.tan,
           sympy.atan: ts.arctan,
           sympy.atan2: ts.arctan2,
           sympy.cosh: ts.cosh,
           sympy.acosh: ts.arccosh,
           sympy.sinh: ts.sinh,
           sympy.asinh: ts.arcsinh,
           sympy.tanh: ts.tanh,
           sympy.atanh: ts.arctanh,
           sympy.re: ts.real,
           sympy.im: ts.imag,
           sympy.arg: ts.angle,
           sympy.erf: ts.erf,
           sympy.gamma: ts.gamma,
           sympy.loggamma: ts.gammaln,
           sympy.Pow: ts.pow,
           sympy.Eq: ts.eq,
           sympy.Gt: ts.gt,
           sympy.Lt: ts.lt,
           sympy.Le: ts.le,
           sympy.Ge: ts.ge,
#           switch,
           sympy.Max: ts.maximum,  # Sympy accept >2 inputs, Theano only 2
           sympy.Min: ts.minimum,  # Sympy accept >2 inputs, Theano only 2

           sympy.MatAdd: tt.Elemwise(ts.add),
           sympy.HadamardProduct: tt.Elemwise(ts.mul),
           sympy.Trace: tlinalg.trace,
           sympy.Inverse: tlinalg.matrix_inverse,
}

class TheanoPrinter(Printer):
    """ Code printer for Theano computations """
    printmethod = "_theano"

    cache = dict()

    def _print_Symbol(self, s, dtypes={}):
        dtype = dtypes.get(s, 'floatX')
        key = (s.name, dtype, type(s))
        if key in self.cache:
            return self.cache[key]
        else:
            value = Scalar(dtype)(s.name)
            self.cache[key] = value
            return value

    def _print_Basic(self, expr, dtypes={}):
        op = mapping[type(expr)]
        children = [self._print(arg, dtypes) for arg in expr.args]
        return op(*children)

    def _print_Number(self, n, dtypes={}):
        return eval(str(n))

    def _print_MatrixSymbol(self, X, dtypes={}):
        dtype = dtypes.get(X, 'floatX')
        # shape = [self._print(d, dtypes) for d in X.shape]
        key = (X.name, dtype, type(X))
        if key in self.cache:
            return self.cache[key]
        else:
            value = tt.Tensor(dtype, (False, False))(X.name)
            self.cache[key] = value
            return value

    def _print_MatMul(self, expr, dtypes):
        children = [self._print(arg, dtypes) for arg in expr.args]
        result = children[0]
        for child in children[1:]:
            result = tt.dot(result, child)
        return result

    def _print_Pi(self, expr, dtypes):
        return 3.141592653589793

    def _print_Rational(self, expr, dtypes):
        return ts.true_div(self._print(expr.p, dtypes),
                           self._print(expr.q, dtypes))

    def _print_Integer(self, expr, dtypes):
        return expr.p

    def _print_factorial(self, expr, dtypes):
        return self._print(sympy.gamma(expr.args[0] + 1), dtypes)

    def emptyPrinter(self, expr):
        return expr

    def doprint(self, expr, dtypes={}):
        """Returns printer's representation for expr (as a string)"""
        return self._print(expr, dtypes)

def theano_code(expr, dtypes={}, **settings):
    return TheanoPrinter(settings).doprint(expr, dtypes)

def unpack(coll):
    if len(coll) == 1:
        return coll[0]
    else:
        return coll

def dim_handling(inputs, dim=None, dims={}, broadcastable={}, keys=()):
    """ Handle various input types for dimensions in tensor_wrap

    See Also:
        tensor_wrap
        theano_funciton
    """
    if keys:
        dims = {i: dims[oi] for i, oi in zip(inputs, keys) if oi in dims}
        broadcastable = {i: broadcastable[oi]
                for i, oi in zip(inputs, keys) if oi in broadcastable}
    if dim:
        dims = dict(zip(inputs, [dim]*len(inputs)))
    if dims:
        maxdim = max(dims.values())
        broadcastable = {i: (False,)*dims[i] + (True,)*(maxdim-dims[i])
                         for i in inputs}
    return broadcastable


def tensor_wrap(inputs, outputs, **kwargs):
    """ Convert scalar io-graph to tensor io-graph """
    broadcastable = dim_handling(inputs, **kwargs)

    Tinputs = [tt.Tensor(i.type.dtype, broadcastable[i])(i.name)
                            for i in inputs]
    Toutputs = tt.Elemwise(ts.Composite(inputs, outputs))(*Tinputs)
    return Tinputs, Toutputs

def theano_function(inputs, outputs, dtypes={}, **kwargs):
    """ Create Theano function from SymPy expressions """
    tinputs  = [theano_code(i, dtypes) for i in inputs]
    toutputs = [theano_code(o, dtypes) for o in outputs]
    if not kwargs:
        toutputs = unpack(toutputs)
        return theano.function(tinputs, toutputs)
    else:
        Tinputs, Toutputs = tensor_wrap(tinputs, toutputs, keys=inputs, **kwargs)
        return theano.function(Tinputs, Toutputs)
