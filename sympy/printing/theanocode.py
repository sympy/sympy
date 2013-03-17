from sympy.utilities import default_sort_key

from sympy.printing.printer import Printer
import theano
import sympy

Scalar = theano.scalar.Scalar
s = theano.scalar
tt = theano.tensor


mapping = {sympy.Add: s.add,
           sympy.Mul: s.mul,
#           sympy.Sub: s.sub,
#           sympy.Mul(x, Pow(y, -1)): s.true_div,  # int_div,
           #sympy.Mod: s.mod, # master branch
           #clip,  # second,
           #identity,
           #cast,
           sympy.Abs: s.abs_,
           sympy.sign: s.sgn,
           sympy.ceiling: s.ceil,
           sympy.floor: s.floor,
#           round_half_to_even, round_half_away_from_zero,
           #lambda x: sympy.Mul(x, -1): s.neg,
           #lambda x: sympy.Pow(x, -1): s.inv,
           sympy.log: s.log,
           #lambda x: log(1 + p): s.log1p,
           #log2, log10,
           sympy.exp: s.exp,
           #exp2,
           #lambda x: sympy.Mul(x, x): s.sqr,
           sympy.sqrt: s.sqrt,
           sympy.cos: s.cos,
           sympy.acos: s.arccos,
           sympy.sin: s.sin,
           sympy.asin: s.arcsin,
           sympy.tan: s.tan,
           sympy.atan: s.arctan,
           sympy.atan2: s.arctan2,
           sympy.cosh: s.cosh,
           sympy.acosh: s.arccosh,
           sympy.sinh: s.sinh,
           sympy.asinh: s.arcsinh,
           sympy.tanh: s.tanh,
           sympy.atanh: s.arctanh,
           sympy.re: s.real,
           sympy.im: s.imag,
           sympy.arg: s.angle,
           sympy.erf: s.erf,
           sympy.gamma: s.gamma,
           sympy.loggamma: s.gammaln,
           sympy.Pow: s.pow,
           sympy.Eq: s.eq,
           sympy.Gt: s.gt,
           sympy.Lt: s.lt,
           sympy.Le: s.le,
           sympy.Ge: s.ge,
#           switch,
           sympy.Max: s.maximum,  # Sympy accept >2 inputs, Theano only 2
           sympy.Min: s.minimum,  # Sympy accept >2 inputs, Theano only 2

           sympy.MatAdd: tt.Elemwise(s.add),
           sympy.HadamardProduct: tt.Elemwise(s.mul),

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

    def _print_Expr(self, expr, dtypes={}):
        op = mapping[type(expr)]
        children = [self._print(arg, dtypes) for arg in expr.args]
        return op(*children)

    def _print_Number(self, n, dtypes={}):
        return eval(str(n))

    def _print_MatrixSymbol(self, X, dtypes={}):
        dtype = dtypes.get(s, 'floatX')
        # shape = [self._print(d, dtypes) for d in X.shape]
        key = (X.name, dtype, type(s))
        if key in self.cache:
            return self.cache[key]
        else:
            value = tt.Tensor(dtype, (False, False))(X.name)
            self.cache[key] = value
            return value

    def _print_MatrixExpr(self, expr, dtypes={}):
        op = mapping[type(expr)]
        children = [self._print(arg, dtypes) for arg in expr.args]
        return op(*children)

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
