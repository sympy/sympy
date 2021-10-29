from sympy.core import S
from .pycode import PythonCodePrinter, _known_functions_math, _print_known_const, _print_known_func, _unpack_integral_limits
from .codeprinter import CodePrinter


_not_in_numpy = 'erf erfc factorial gamma loggamma'.split()
_in_numpy = [(k, v) for k, v in _known_functions_math.items() if k not in _not_in_numpy]
_known_functions_numpy = dict(_in_numpy, **{
    'acos': 'arccos',
    'acosh': 'arccosh',
    'asin': 'arcsin',
    'asinh': 'arcsinh',
    'atan': 'arctan',
    'atan2': 'arctan2',
    'atanh': 'arctanh',
    'exp2': 'exp2',
    'sign': 'sign',
    'logaddexp': 'logaddexp',
    'logaddexp2': 'logaddexp2',
})
_known_constants_numpy = {
    'Exp1': 'e',
    'Pi': 'pi',
    'EulerGamma': 'euler_gamma',
    'NaN': 'nan',
    'Infinity': 'PINF',
    'NegativeInfinity': 'NINF'
}

_numpy_known_functions = {k: 'numpy.' + v for k, v in _known_functions_numpy.items()}
_numpy_known_constants = {k: 'numpy.' + v for k, v in _known_constants_numpy.items()}

class NumPyPrinter(PythonCodePrinter):
    """
    Numpy printer which handles vectorized piecewise functions,
    logical operators, etc.
    """

    _module = 'numpy'
    _kf = _numpy_known_functions
    _kc = _numpy_known_constants

    def __init__(self, settings=None):
        """
        `settings` is passed to CodePrinter.__init__()
        `module` specifies the array module to use, currently 'NumPy' or 'CuPy'
        """
        self.language = "Python with {}".format(self._module)
        self.printmethod = "_{}code".format(self._module)

        self._kf = {**PythonCodePrinter._kf, **self._kf}

        super().__init__(settings=settings)


    def _print_seq(self, seq):
        "General sequence printer: converts to tuple"
        # Print tuples here instead of lists because numba supports
        #     tuples in nopython mode.
        delimiter=', '
        return '({},)'.format(delimiter.join(self._print(item) for item in seq))

    def _print_MatMul(self, expr):
        "Matrix multiplication printer"
        if expr.as_coeff_matrices()[0] is not S.One:
            expr_list = expr.as_coeff_matrices()[1]+[(expr.as_coeff_matrices()[0])]
            return '({})'.format(').dot('.join(self._print(i) for i in expr_list))
        return '({})'.format(').dot('.join(self._print(i) for i in expr.args))

    def _print_MatPow(self, expr):
        "Matrix power printer"
        return '{}({}, {})'.format(self._module_format(self._module + '.linalg.matrix_power'),
            self._print(expr.args[0]), self._print(expr.args[1]))

    def _print_Inverse(self, expr):
        "Matrix inverse printer"
        return '{}({})'.format(self._module_format(self._module + '.linalg.inv'),
            self._print(expr.args[0]))

    def _print_DotProduct(self, expr):
        # DotProduct allows any shape order, but numpy.dot does matrix
        # multiplication, so we have to make sure it gets 1 x n by n x 1.
        arg1, arg2 = expr.args
        if arg1.shape[0] != 1:
            arg1 = arg1.T
        if arg2.shape[1] != 1:
            arg2 = arg2.T

        return "%s(%s, %s)" % (self._module_format(self._module + '.dot'),
                               self._print(arg1),
                               self._print(arg2))

    def _print_MatrixSolve(self, expr):
        return "%s(%s, %s)" % (self._module_format(self._module + '.linalg.solve'),
                               self._print(expr.matrix),
                               self._print(expr.vector))

    def _print_ZeroMatrix(self, expr):
        return '{}({})'.format(self._module_format(self._module + '.zeros'),
            self._print(expr.shape))

    def _print_OneMatrix(self, expr):
        return '{}({})'.format(self._module_format(self._module + '.ones'),
            self._print(expr.shape))

    def _print_FunctionMatrix(self, expr):
        from sympy.core.function import Lambda
        from sympy.abc import i, j
        lamda = expr.lamda
        if not isinstance(lamda, Lambda):
            lamda = Lambda((i, j), lamda(i, j))
        return '{}(lambda {}: {}, {})'.format(self._module_format(self._module + '.fromfunction'),
            ', '.join(self._print(arg) for arg in lamda.args[0]),
            self._print(lamda.args[1]), self._print(expr.shape))

    def _print_HadamardProduct(self, expr):
        func = self._module_format(self._module + '.multiply')
        return ''.join('{}({}, '.format(func, self._print(arg)) \
            for arg in expr.args[:-1]) + "{}{}".format(self._print(expr.args[-1]),
            ')' * (len(expr.args) - 1))

    def _print_KroneckerProduct(self, expr):
        func = self._module_format(self._module + '.kron')
        return ''.join('{}({}, '.format(func, self._print(arg)) \
            for arg in expr.args[:-1]) + "{}{}".format(self._print(expr.args[-1]),
            ')' * (len(expr.args) - 1))

    def _print_Adjoint(self, expr):
        return '{}({}({}))'.format(
            self._module_format(self._module + '.conjugate'),
            self._module_format(self._module + '.transpose'),
            self._print(expr.args[0]))

    def _print_DiagonalOf(self, expr):
        vect = '{}({})'.format(
            self._module_format(self._module + '.diag'),
            self._print(expr.arg))
        return '{}({}, (-1, 1))'.format(
            self._module_format(self._module + '.reshape'), vect)

    def _print_DiagMatrix(self, expr):
        return '{}({})'.format(self._module_format(self._module + '.diagflat'),
            self._print(expr.args[0]))

    def _print_DiagonalMatrix(self, expr):
        return '{}({}, {}({}, {}))'.format(self._module_format(self._module + '.multiply'),
            self._print(expr.arg), self._module_format(self._module + '.eye'),
            self._print(expr.shape[0]), self._print(expr.shape[1]))

    def _print_Piecewise(self, expr):
        "Piecewise function printer"
        exprs = '[{}]'.format(','.join(self._print(arg.expr) for arg in expr.args))
        conds = '[{}]'.format(','.join(self._print(arg.cond) for arg in expr.args))
        # If [default_value, True] is a (expr, cond) sequence in a Piecewise object
        #     it will behave the same as passing the 'default' kwarg to select()
        #     *as long as* it is the last element in expr.args.
        # If this is not the case, it may be triggered prematurely.
        return '{}({}, {}, default={})'.format(
            self._module_format(self._module + '.select'), conds, exprs,
            self._print(S.NaN))

    def _print_Relational(self, expr):
        "Relational printer for Equality and Unequality"
        op = {
            '==' :'equal',
            '!=' :'not_equal',
            '<'  :'less',
            '<=' :'less_equal',
            '>'  :'greater',
            '>=' :'greater_equal',
        }
        if expr.rel_op in op:
            lhs = self._print(expr.lhs)
            rhs = self._print(expr.rhs)
            return '{op}({lhs}, {rhs})'.format(op=self._module_format(self._module + '.'+op[expr.rel_op]),
                                               lhs=lhs, rhs=rhs)
        return super()._print_Relational(expr)

    def _print_And(self, expr):
        "Logical And printer"
        # We have to override LambdaPrinter because it uses Python 'and' keyword.
        # If LambdaPrinter didn't define it, we could use StrPrinter's
        # version of the function and add 'logical_and' to NUMPY_TRANSLATIONS.
        return '{}.reduce(({}))'.format(self._module_format(self._module + '.logical_and'), ','.join(self._print(i) for i in expr.args))

    def _print_Or(self, expr):
        "Logical Or printer"
        # We have to override LambdaPrinter because it uses Python 'or' keyword.
        # If LambdaPrinter didn't define it, we could use StrPrinter's
        # version of the function and add 'logical_or' to NUMPY_TRANSLATIONS.
        return '{}.reduce(({}))'.format(self._module_format(self._module + '.logical_or'), ','.join(self._print(i) for i in expr.args))

    def _print_Not(self, expr):
        "Logical Not printer"
        # We have to override LambdaPrinter because it uses Python 'not' keyword.
        # If LambdaPrinter didn't define it, we would still have to define our
        #     own because StrPrinter doesn't define it.
        return '{}({})'.format(self._module_format(self._module + '.logical_not'), ','.join(self._print(i) for i in expr.args))

    def _print_Pow(self, expr, rational=False):
        # XXX Workaround for negative integer power error
        from sympy.core.power import Pow
        if expr.exp.is_integer and expr.exp.is_negative:
            expr = Pow(expr.base, expr.exp.evalf(), evaluate=False)
        return self._hprint_Pow(expr, rational=rational, sqrt=self._module + '.sqrt')

    def _print_Min(self, expr):
        return '{}(({}), axis=0)'.format(self._module_format(self._module + '.amin'), ','.join(self._print(i) for i in expr.args))

    def _print_Max(self, expr):
        return '{}(({}), axis=0)'.format(self._module_format(self._module + '.amax'), ','.join(self._print(i) for i in expr.args))

    def _print_arg(self, expr):
        return "%s(%s)" % (self._module_format(self._module + '.angle'), self._print(expr.args[0]))

    def _print_im(self, expr):
        return "%s(%s)" % (self._module_format(self._module + '.imag'), self._print(expr.args[0]))

    def _print_Mod(self, expr):
        return "%s(%s)" % (self._module_format(self._module + '.mod'), ', '.join(
            map(lambda arg: self._print(arg), expr.args)))

    def _print_re(self, expr):
        return "%s(%s)" % (self._module_format(self._module + '.real'), self._print(expr.args[0]))

    def _print_sinc(self, expr):
        return "%s(%s)" % (self._module_format(self._module + '.sinc'), self._print(expr.args[0]/S.Pi))

    def _print_MatrixBase(self, expr):
        func = self.known_functions.get(expr.__class__.__name__, None)
        if func is None:
            func = self._module_format(self._module + '.array')
        return "%s(%s)" % (func, self._print(expr.tolist()))

    def _print_Identity(self, expr):
        shape = expr.shape
        if all(dim.is_Integer for dim in shape):
            return "%s(%s)" % (self._module_format(self._module + '.eye'), self._print(expr.shape[0]))
        else:
            raise NotImplementedError("Symbolic matrix dimensions are not yet supported for identity matrices")

    def _print_BlockMatrix(self, expr):
        return '{}({})'.format(self._module_format(self._module + '.block'),
                                 self._print(expr.args[0].tolist()))

    def _print_ArrayTensorProduct(self, expr):
        array_list = [j for i, arg in enumerate(expr.args) for j in
                (self._print(arg), "[%i, %i]" % (2*i, 2*i+1))]
        return "%s(%s)" % (self._module_format(self._module + '.einsum'), ", ".join(array_list))

    def _print_ArrayContraction(self, expr):
        from ..tensor.array.expressions.array_expressions import ArrayTensorProduct
        base = expr.expr
        contraction_indices = expr.contraction_indices
        if not contraction_indices:
            return self._print(base)
        if isinstance(base, ArrayTensorProduct):
            counter = 0
            d = {j: min(i) for i in contraction_indices for j in i}
            indices = []
            for rank_arg in base.subranks:
                lindices = []
                for i in range(rank_arg):
                    if counter in d:
                        lindices.append(d[counter])
                    else:
                        lindices.append(counter)
                    counter += 1
                indices.append(lindices)
            elems = ["%s, %s" % (self._print(arg), ind) for arg, ind in zip(base.args, indices)]
            return "%s(%s)" % (
                self._module_format(self._module + '.einsum'),
                ", ".join(elems)
            )
        raise NotImplementedError()

    def _print_ArrayDiagonal(self, expr):
        diagonal_indices = list(expr.diagonal_indices)
        if len(diagonal_indices) > 1:
            # TODO: this should be handled in sympy.codegen.array_utils,
            # possibly by creating the possibility of unfolding the
            # ArrayDiagonal object into nested ones. Same reasoning for
            # the array contraction.
            raise NotImplementedError
        if len(diagonal_indices[0]) != 2:
            raise NotImplementedError
        return "%s(%s, 0, axis1=%s, axis2=%s)" % (
            self._module_format("numpy.diagonal"),
            self._print(expr.expr),
            diagonal_indices[0][0],
            diagonal_indices[0][1],
        )

    def _print_PermuteDims(self, expr):
        return "%s(%s, %s)" % (
            self._module_format("numpy.transpose"),
            self._print(expr.expr),
            self._print(expr.permutation.array_form),
        )

    def _print_ArrayAdd(self, expr):
        return self._expand_fold_binary_op(self._module + '.add', expr.args)

    _print_lowergamma = CodePrinter._print_not_supported
    _print_uppergamma = CodePrinter._print_not_supported
    _print_fresnelc = CodePrinter._print_not_supported
    _print_fresnels = CodePrinter._print_not_supported

for func in _numpy_known_functions:
    setattr(NumPyPrinter, f'_print_{func}', _print_known_func)

for const in _numpy_known_constants:
    setattr(NumPyPrinter, f'_print_{const}', _print_known_const)


_known_functions_scipy_special = {
    'erf': 'erf',
    'erfc': 'erfc',
    'besselj': 'jv',
    'bessely': 'yv',
    'besseli': 'iv',
    'besselk': 'kv',
    'cosm1': 'cosm1',
    'factorial': 'factorial',
    'gamma': 'gamma',
    'loggamma': 'gammaln',
    'digamma': 'psi',
    'RisingFactorial': 'poch',
    'jacobi': 'eval_jacobi',
    'gegenbauer': 'eval_gegenbauer',
    'chebyshevt': 'eval_chebyt',
    'chebyshevu': 'eval_chebyu',
    'legendre': 'eval_legendre',
    'hermite': 'eval_hermite',
    'laguerre': 'eval_laguerre',
    'assoc_laguerre': 'eval_genlaguerre',
    'beta': 'beta',
    'LambertW' : 'lambertw',
}

_known_constants_scipy_constants = {
    'GoldenRatio': 'golden_ratio',
    'Pi': 'pi',
}
_scipy_known_functions = {k : "scipy.special." + v for k, v in _known_functions_scipy_special.items()}
_scipy_known_constants = {k : "scipy.constants." + v for k, v in _known_constants_scipy_constants.items()}

class SciPyPrinter(NumPyPrinter):

    _kf = {**NumPyPrinter._kf, **_scipy_known_functions}
    _kc = {**NumPyPrinter._kc, **_scipy_known_constants}

    def __init__(self, settings=None):
        super().__init__(settings=settings)
        self.language = "Python with SciPy and NumPy"

    def _print_SparseRepMatrix(self, expr):
        i, j, data = [], [], []
        for (r, c), v in expr.todok().items():
            i.append(r)
            j.append(c)
            data.append(v)

        return "{name}(({data}, ({i}, {j})), shape={shape})".format(
            name=self._module_format('scipy.sparse.coo_matrix'),
            data=data, i=i, j=j, shape=expr.shape
        )

    _print_ImmutableSparseMatrix = _print_SparseRepMatrix

    # SciPy's lpmv has a different order of arguments from assoc_legendre
    def _print_assoc_legendre(self, expr):
        return "{0}({2}, {1}, {3})".format(
            self._module_format('scipy.special.lpmv'),
            self._print(expr.args[0]),
            self._print(expr.args[1]),
            self._print(expr.args[2]))

    def _print_lowergamma(self, expr):
        return "{0}({2})*{1}({2}, {3})".format(
            self._module_format('scipy.special.gamma'),
            self._module_format('scipy.special.gammainc'),
            self._print(expr.args[0]),
            self._print(expr.args[1]))

    def _print_uppergamma(self, expr):
        return "{0}({2})*{1}({2}, {3})".format(
            self._module_format('scipy.special.gamma'),
            self._module_format('scipy.special.gammaincc'),
            self._print(expr.args[0]),
            self._print(expr.args[1]))

    def _print_betainc(self, expr):
        betainc = self._module_format('scipy.special.betainc')
        beta = self._module_format('scipy.special.beta')
        args = [self._print(arg) for arg in expr.args]
        return f"({betainc}({args[0]}, {args[1]}, {args[3]}) - {betainc}({args[0]}, {args[1]}, {args[2]})) \
            * {beta}({args[0]}, {args[1]})"

    def _print_betainc_regularized(self, expr):
        return "{0}({1}, {2}, {4}) - {0}({1}, {2}, {3})".format(
            self._module_format('scipy.special.betainc'),
            self._print(expr.args[0]),
            self._print(expr.args[1]),
            self._print(expr.args[2]),
            self._print(expr.args[3]))

    def _print_fresnels(self, expr):
        return "{}({})[0]".format(
                self._module_format("scipy.special.fresnel"),
                self._print(expr.args[0]))

    def _print_fresnelc(self, expr):
        return "{}({})[1]".format(
                self._module_format("scipy.special.fresnel"),
                self._print(expr.args[0]))

    def _print_airyai(self, expr):
        return "{}({})[0]".format(
                self._module_format("scipy.special.airy"),
                self._print(expr.args[0]))

    def _print_airyaiprime(self, expr):
        return "{}({})[1]".format(
                self._module_format("scipy.special.airy"),
                self._print(expr.args[0]))

    def _print_airybi(self, expr):
        return "{}({})[2]".format(
                self._module_format("scipy.special.airy"),
                self._print(expr.args[0]))

    def _print_airybiprime(self, expr):
        return "{}({})[3]".format(
                self._module_format("scipy.special.airy"),
                self._print(expr.args[0]))

    def _print_Integral(self, e):
        integration_vars, limits = _unpack_integral_limits(e)

        if len(limits) == 1:
            # nicer (but not necessary) to prefer quad over nquad for 1D case
            module_str = self._module_format("scipy.integrate.quad")
            limit_str = "%s, %s" % tuple(map(self._print, limits[0]))
        else:
            module_str = self._module_format("scipy.integrate.nquad")
            limit_str = "({})".format(", ".join(
                "(%s, %s)" % tuple(map(self._print, l)) for l in limits))

        return "{}(lambda {}: {}, {})[0]".format(
                module_str,
                ", ".join(map(self._print, integration_vars)),
                self._print(e.args[0]),
                limit_str)


for func in _scipy_known_functions:
    setattr(SciPyPrinter, f'_print_{func}', _print_known_func)

for const in _scipy_known_constants:
    setattr(SciPyPrinter, f'_print_{const}', _print_known_const)


_cupy_known_functions = {k : "cupy." + v for k, v in _known_functions_numpy.items()}
_cupy_known_constants = {k : "cupy." + v for k, v in _known_constants_numpy.items()}

class CuPyPrinter(NumPyPrinter):
    """
    CuPy printer which handles vectorized piecewise functions,
    logical operators, etc.
    """

    _module = 'cupy'
    _kf = _cupy_known_functions
    _kc = _cupy_known_constants

    def __init__(self, settings=None):
        super().__init__(settings=settings)

for func in _cupy_known_functions:
    setattr(CuPyPrinter, f'_print_{func}', _print_known_func)

for const in _cupy_known_constants:
    setattr(CuPyPrinter, f'_print_{const}', _print_known_const)
