from sympy.external.importtools import version_tuple
from collections.abc import Iterable

from sympy.core.mul import Mul
from sympy.core.singleton import S
from sympy.codegen.cfunctions import Sqrt
from sympy.external import import_module
from sympy.printing.precedence import PRECEDENCE
from sympy.printing.pycode import AbstractPythonCodePrinter
import sympy

tensorflow = import_module('tensorflow')

class TensorflowPrinter(AbstractPythonCodePrinter):
    """
    Tensorflow printer which handles vectorized piecewise functions,
    logical operators, max/min, and relational operators.
    """
    printmethod = "_tensorflowcode"

    mapping = {
        sympy.Abs: "tensorflow.math.abs",
        sympy.sign: "tensorflow.math.sign",

        # XXX May raise error for ints.
        sympy.ceiling: "tensorflow.math.ceil",
        sympy.floor: "tensorflow.math.floor",
        sympy.log: "tensorflow.math.log",
        sympy.exp: "tensorflow.math.exp",
        Sqrt: "tensorflow.math.sqrt",
        sympy.cos: "tensorflow.math.cos",
        sympy.acos: "tensorflow.math.acos",
        sympy.sin: "tensorflow.math.sin",
        sympy.asin: "tensorflow.math.asin",
        sympy.tan: "tensorflow.math.tan",
        sympy.atan: "tensorflow.math.atan",
        sympy.atan2: "tensorflow.math.atan2",
        # XXX Also may give NaN for complex results.
        sympy.cosh: "tensorflow.math.cosh",
        sympy.acosh: "tensorflow.math.acosh",
        sympy.sinh: "tensorflow.math.sinh",
        sympy.asinh: "tensorflow.math.asinh",
        sympy.tanh: "tensorflow.math.tanh",
        sympy.atanh: "tensorflow.math.atanh",

        sympy.re: "tensorflow.math.real",
        sympy.im: "tensorflow.math.imag",
        sympy.arg: "tensorflow.math.angle",

        # XXX May raise error for ints and complexes
        sympy.erf: "tensorflow.math.erf",
        sympy.loggamma: "tensorflow.math.lgamma",

        sympy.Eq: "tensorflow.math.equal",
        sympy.Ne: "tensorflow.math.not_equal",
        sympy.StrictGreaterThan: "tensorflow.math.greater",
        sympy.StrictLessThan: "tensorflow.math.less",
        sympy.LessThan: "tensorflow.math.less_equal",
        sympy.GreaterThan: "tensorflow.math.greater_equal",

        sympy.And: "tensorflow.math.logical_and",
        sympy.Or: "tensorflow.math.logical_or",
        sympy.Not: "tensorflow.math.logical_not",
        sympy.Max: "tensorflow.math.maximum",
        sympy.Min: "tensorflow.math.minimum",

        # Matrices
        sympy.MatAdd: "tensorflow.math.add",
        sympy.HadamardProduct: "tensorflow.math.multiply",
        sympy.Trace: "tensorflow.linalg.trace",

        # XXX May raise error for integer matrices.
        sympy.Determinant : "tensorflow.linalg.det",
    }

    _default_settings = dict(
        AbstractPythonCodePrinter._default_settings,
        tensorflow_version=None
    )

    def __init__(self, settings=None):
        super().__init__(settings)

        version = self._settings['tensorflow_version']
        if version is None and tensorflow:
            version = tensorflow.__version__
        self.tensorflow_version = version

    def _print_Function(self, expr):
        op = self.mapping.get(type(expr), None)
        if op is None:
            return super()._print_Basic(expr)
        children = [self._print(arg) for arg in expr.args]
        if len(children) == 1:
            return "%s(%s)" % (
                self._module_format(op),
                children[0]
            )
        else:
            return self._expand_fold_binary_op(op, children)

    _print_Expr = _print_Function
    _print_Application = _print_Function
    _print_MatrixExpr = _print_Function
    # TODO: a better class structure would avoid this mess:
    _print_Relational = _print_Function
    _print_Not = _print_Function
    _print_And = _print_Function
    _print_Or = _print_Function
    _print_HadamardProduct = _print_Function
    _print_Trace = _print_Function
    _print_Determinant = _print_Function

    def _print_Inverse(self, expr):
        op = self._module_format('tensorflow.linalg.inv')
        return "{}({})".format(op, self._print(expr.arg))

    def _print_Transpose(self, expr):
        version = self.tensorflow_version
        if version and version_tuple(version) < version_tuple('1.14'):
            op = self._module_format('tensorflow.matrix_transpose')
        else:
            op = self._module_format('tensorflow.linalg.matrix_transpose')
        return "{}({})".format(op, self._print(expr.arg))

    def _print_Derivative(self, expr):
        variables = expr.variables
        if any(isinstance(i, Iterable) for i in variables):
            raise NotImplementedError("derivation by multiple variables is not supported")
        def unfold(expr, args):
            if not args:
                return self._print(expr)
            return "%s(%s, %s)[0]" % (
                    self._module_format("tensorflow.gradients"),
                    unfold(expr, args[:-1]),
                    self._print(args[-1]),
                )
        return unfold(expr.expr, variables)

    def _print_Piecewise(self, expr):
        version = self.tensorflow_version
        if version and version_tuple(version) < version_tuple('1.0'):
            tensorflow_piecewise = "tensorflow.select"
        else:
            tensorflow_piecewise = "tensorflow.where"

        from sympy.functions.elementary.piecewise import Piecewise
        e, cond = expr.args[0].args
        if len(expr.args) == 1:
            return '{}({}, {}, {})'.format(
                self._module_format(tensorflow_piecewise),
                self._print(cond),
                self._print(e),
                0)

        return '{}({}, {}, {})'.format(
            self._module_format(tensorflow_piecewise),
            self._print(cond),
            self._print(e),
            self._print(Piecewise(*expr.args[1:])))

    def _print_Pow(self, expr):
        # XXX May raise error for
        # int**float or int**complex or float**complex
        base, exp = expr.args
        if expr.exp == S.Half:
            return "{}({})".format(
                self._module_format("tensorflow.math.sqrt"), self._print(base))
        return "{}({}, {})".format(
            self._module_format("tensorflow.math.pow"),
            self._print(base), self._print(exp))

    def _print_MatrixBase(self, expr):
        tensorflow_f = "tensorflow.Variable" if expr.free_symbols else "tensorflow.constant"
        data = "["+", ".join(["["+", ".join([self._print(j) for j in i])+"]" for i in expr.tolist()])+"]"
        return "%s(%s)" % (
            self._module_format(tensorflow_f),
            data,
        )

    def _print_MatMul(self, expr):
        from sympy.matrices.expressions import MatrixExpr
        mat_args = [arg for arg in expr.args if isinstance(arg, MatrixExpr)]
        args = [arg for arg in expr.args if arg not in mat_args]
        if args:
            return "%s*%s" % (
                self.parenthesize(Mul.fromiter(args), PRECEDENCE["Mul"]),
                self._expand_fold_binary_op(
                    "tensorflow.linalg.matmul", mat_args)
            )
        else:
            return self._expand_fold_binary_op(
                "tensorflow.linalg.matmul", mat_args)

    def _print_MatPow(self, expr):
        return self._expand_fold_binary_op(
            "tensorflow.linalg.matmul", [expr.base]*expr.exp)

    def _print_Assignment(self, expr):
        # TODO: is this necessary?
        return "%s = %s" % (
            self._print(expr.lhs),
            self._print(expr.rhs),
        )

    def _print_CodeBlock(self, expr):
        # TODO: is this necessary?
        ret = []
        for subexpr in expr.args:
            ret.append(self._print(subexpr))
        return "\n".join(ret)

    def _get_letter_generator_for_einsum(self):
        for i in range(97, 123):
            yield chr(i)
        for i in range(65, 91):
            yield chr(i)
        raise ValueError("out of letters")

    def _print_ArrayTensorProduct(self, expr):
        letters = self._get_letter_generator_for_einsum()
        contraction_string = ",".join(["".join([next(letters) for j in range(i)]) for i in expr.subranks])
        return '%s("%s", %s)' % (
                self._module_format('tensorflow.linalg.einsum'),
                contraction_string,
                ", ".join([self._print(arg) for arg in expr.args])
        )

    def _print_ArrayContraction(self, expr):
        from sympy.tensor.array.expressions.array_expressions import ArrayTensorProduct
        base = expr.expr
        contraction_indices = expr.contraction_indices
        contraction_string, letters_free, letters_dum = self._get_einsum_string(base.subranks, contraction_indices)

        if not contraction_indices:
            return self._print(base)
        if isinstance(base, ArrayTensorProduct):
            elems = ["%s" % (self._print(arg)) for arg in base.args]
            return "%s(\"%s\", %s)" % (
                self._module_format("tensorflow.linalg.einsum"),
                contraction_string,
                ", ".join(elems)
            )
        raise NotImplementedError()

    def _print_ArrayDiagonal(self, expr):
        from sympy.tensor.array.expressions.array_expressions import ArrayTensorProduct
        diagonal_indices = list(expr.diagonal_indices)
        if len(diagonal_indices) > 1:
            # TODO: this should be handled in sympy.codegen.array_utils,
            # possibly by creating the possibility of unfolding the
            # ArrayDiagonal object into nested ones. Same reasoning for
            # the array contraction.
            raise NotImplementedError
        if len(diagonal_indices[0]) != 2:
            raise NotImplementedError
        if isinstance(expr.expr, ArrayTensorProduct):
            subranks = expr.expr.subranks
            elems = expr.expr.args
        else:
            subranks = expr.subranks
            elems = [expr.expr]
        diagonal_string, letters_free, letters_dum = self._get_einsum_string(subranks, diagonal_indices)
        elems = [self._print(i) for i in elems]
        return '%s("%s", %s)' % (
            self._module_format("tensorflow.linalg.einsum"),
            "{}->{}{}".format(diagonal_string, "".join(letters_free), "".join(letters_dum)),
            ", ".join(elems)
        )

    def _print_PermuteDims(self, expr):
        return "%s(%s, %s)" % (
            self._module_format("tensorflow.transpose"),
            self._print(expr.expr),
            self._print(expr.permutation.array_form),
        )

    def _print_ArrayAdd(self, expr):
        return self._expand_fold_binary_op('tensorflow.math.add', expr.args)


def tensorflow_code(expr, **settings):
    printer = TensorflowPrinter(settings)
    return printer.doprint(expr)
