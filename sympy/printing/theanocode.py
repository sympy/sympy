from sympy.utilities import default_sort_key

from sympy.printing.printer import Printer
import theano
import sympy

Scalar = theano.scalar.Scalar
ts = theano.scalar
tt = theano.tensor


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

}# Implement factorial

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

    def doprint(self, expr, dtypes={}):
        """Returns printer's representation for expr (as a string)"""
        return self._print(expr, dtypes)

def theano_code(expr, dtypes={}, **settings):
    return TheanoPrinter(settings).doprint(expr, dtypes)
