import functools as ft

from sympy.printing.pycode import AbstractPythonCodePrinter, ArrayPrinter
from sympy.tensor.array.expressions.array_expressions import ArrayTensorProduct
from sympy.matrices.expressions.fourier import DFT
from sympy.combinatorics.permutations import Permutation
from sympy.matrices.expressions import MatrixExpr
from sympy.core.mul import Mul
from sympy.printing.precedence import PRECEDENCE
from sympy.external import import_module
from sympy.codegen.cfunctions import Sqrt
from sympy import S, MutableDenseMatrix
from sympy import Integer

import sympy


torch = import_module('torch')


def _reduce(fn):
    def fn_(*args):
        return ft.reduce(fn, args)
    return fn_


def piecewise(*expr_conds):
    output = None
    already_used = None
    for expr, cond in expr_conds:
        if output is None:
            already_used = cond
            output = torch.where(cond, expr, torch.zeros_like(expr))
        else:
            output += torch.where(cond.bool() & ~already_used, expr, torch.zeros_like(expr))
            already_used = already_used | cond.bool()
    return output


def expr_cond_pair(expr, cond):
    return expr, cond


def if_then_else(cond, if_true, if_false):
    return torch.where(cond, if_true, if_false)


def _get_einsum_spec(subranks, contraction_indices=None, diagonal_indices=None):
    letters = [chr(i) for i in range(97, 123)] + [chr(i) for i in range(65, 91)]
    if len(subranks) * max(subranks, default=0) > len(letters):
        raise ValueError("Too many indices for available letters")

    index_map = {}
    counter = 0
    input_specs = []
    for rank in subranks:
        arg_indices = []
        for _ in range(rank):
            if counter not in index_map:
                index_map[counter] = letters[len(index_map)]
            arg_indices.append(index_map[counter])
            counter += 1
        input_specs.append("".join(arg_indices))

    if diagonal_indices:
        if all(isinstance(i, (int, sympy.Integer)) for i in diagonal_indices):
            diagonal_indices = (tuple(diagonal_indices),)
        d = {j: min(group) for group in diagonal_indices for j in group}
        indices = []
        counter = 0
        for rank in subranks:
            lindices = []
            for _ in range(rank):
                if counter in d:
                    lindices.append(d[counter])
                else:
                    lindices.append(counter)
                counter += 1
            indices.append(lindices)
        mapping = {}
        letters_free = []
        letters_dum = []
        input_specs = []
        for i in indices:
            spec = ""
            for j in i:
                if j not in mapping:
                    l = letters[len(mapping)]
                    mapping[j] = l
                else:
                    l = mapping[j]
                spec += l
                if j in d:
                    if l not in letters_dum:
                        letters_dum.append(l)
                else:
                    letters_free.append(l)
            input_specs.append(spec)

    output_spec = "".join(sorted(letters_free + letters_dum))
    input_spec = ",".join(input_specs)
    return f"{input_spec}->{output_spec}"


def array_tensor_product(*args):
    if not args:
        raise ValueError("At least one argument is required for ArrayTensorProduct")
    total_dims = sum(len(arg.shape) for arg in args)
    letters = [chr(97 + i) for i in range(total_dims)]  # 'a', 'b', 'c', ...
    input_specs = []
    start = 0
    for arg in args:
        rank = len(arg.shape)
        input_specs.append("".join(letters[start:start + rank]))
        start += rank
    input_spec = ",".join(input_specs)
    output_spec = "".join(letters)
    return torch.einsum(f"{input_spec}->{output_spec}", *args)


def array_contraction(expr, contraction_indices, *args):
    subranks = expr.expr.subranks if isinstance(expr.expr, ArrayTensorProduct) else [len(expr.expr.shape)]
    operands = expr.expr.args if isinstance(expr.expr, ArrayTensorProduct) else [expr.expr]
    spec = _get_einsum_spec(subranks, contraction_indices=contraction_indices)
    return torch.einsum(spec, *operands)


def array_diagonal(expr, diagonal_indices, *args):
    if isinstance(expr, torch.Tensor):
        # Single tensor case: adjust to match ArrayPrinter's multi-tensor behavior
        spec = _get_einsum_spec([expr.dim()], diagonal_indices=diagonal_indices)
        return torch.einsum(spec, expr)
    else:
        subranks = expr.expr.subranks if isinstance(expr.expr, ArrayTensorProduct) else [len(expr.expr.shape)]
        operands = expr.expr.args if isinstance(expr.expr, ArrayTensorProduct) else [expr.expr]
        spec = _get_einsum_spec(subranks, diagonal_indices=diagonal_indices)
        return torch.einsum(spec, *operands)


def permute_dims(tensor, perm_indices, *args):
    return tensor.permute(*perm_indices)


def array_add(*args):
    """Add multiple arrays using PyTorch's broadcasting."""
    if not args:
        raise ValueError("At least one argument is required for ArrayAdd")
    if len(args) == 1:
        return args[0]
    result = args[0]
    for arg in args[1:]:
        result = torch.add(result, arg)
    return result


_TORCH_FUNCTION_MAP = {
    sympy.Mul: ("torch.mul", _reduce(torch.mul)),
    sympy.Add: ("torch.add", _reduce(torch.add)),
    sympy.div: ("torch.div", torch.div),
    sympy.Abs: ("torch.abs", torch.abs),
    sympy.sign: ("torch.sign", torch.sign),
    sympy.ceiling: ("torch.ceil", torch.ceil),
    sympy.floor: ("torch.floor", torch.floor),
    sympy.log: ("torch.log", torch.log),
    sympy.exp: ("torch.exp", torch.exp),
    Sqrt: ("torch.sqrt", torch.sqrt),
    sympy.cos: ("torch.cos", torch.cos),
    sympy.acos: ("torch.acos", torch.acos),
    sympy.sin: ("torch.sin", torch.sin),
    sympy.asin: ("torch.asin", torch.asin),
    sympy.tan: ("torch.tan", torch.tan),
    sympy.atan: ("torch.atan", torch.atan),
    sympy.atan2: ("torch.atan2", torch.atan2),
    sympy.cosh: ("torch.cosh", torch.cosh),
    sympy.acosh: ("torch.acosh", torch.acosh),
    sympy.sinh: ("torch.sinh", torch.sinh),
    sympy.asinh: ("torch.asinh", torch.asinh),
    sympy.tanh: ("torch.tanh", torch.tanh),
    sympy.atanh: ("torch.atanh", torch.atanh),
    sympy.Pow: ("torch.pow", torch.pow),
    sympy.re: ("torch.real", torch.real),
    sympy.im: ("torch.imag", torch.imag),
    sympy.arg: ("torch.angle", torch.angle),
    sympy.conjugate: ("torch.conj", torch.conj),
    sympy.core.numbers.ImaginaryUnit: ("1j", lambda *args: torch.tensor(1j)),
    sympy.erf: ("torch.erf", torch.erf),
    sympy.loggamma: ("torch.lgamma", torch.lgamma),
    sympy.Eq: ("torch.eq", torch.eq),
    sympy.Ne: ("torch.ne", torch.ne),
    sympy.StrictGreaterThan: ("torch.gt", torch.gt),
    sympy.StrictLessThan: ("torch.lt", torch.lt),
    sympy.LessThan: ("torch.le", torch.le),
    sympy.GreaterThan: ("torch.ge", torch.ge),
    sympy.And: ("torch.logical_and", torch.logical_and),
    sympy.Or: ("torch.logical_or", torch.logical_or),
    sympy.Not: ("torch.logical_not", torch.logical_not),
    sympy.Max: ("torch.max", torch.max),
    sympy.Min: ("torch.min", torch.min),
    sympy.MatAdd: ("torch.add", _reduce(torch.add)),
    sympy.HadamardProduct: ("torch.mul", torch.mul),
    sympy.Trace: ("torch.trace", torch.trace),
    sympy.Determinant: ("torch.det", torch.det),
    sympy.MatMul: ("torch.matmul", _reduce(torch.matmul)),
    sympy.MatPow: ("torch.linalg.matrix_power", torch.linalg.matrix_power),
    sympy.Inverse: ("torch.linalg.inv", lambda x, *args: torch.linalg.inv(x)),
    sympy.Piecewise: ("torch.where", piecewise),
    sympy.functions.elementary.piecewise.ExprCondPair: ("expr_cond_pair", expr_cond_pair),
    sympy.Mod: ("torch.remainder", torch.remainder),
    sympy.Heaviside: ("torch.heaviside", lambda x, h=None: torch.heaviside(x, torch.tensor(
        float(h) if h is not None else 0.5, dtype=x.dtype))),
    sympy.core.numbers.Half: ("0.5", lambda *args: torch.tensor(0.5)),
    sympy.core.numbers.One: ("1.0", lambda *args: torch.tensor(1.0)),
    sympy.logic.boolalg.ITE: ("torch.where", if_then_else),
}


_TORCHMODULE_ONLY_FUNCTION_MAP = {
    # Core matrix and symbol operations
    sympy.matrices.expressions.matexpr.MatrixSymbol: lambda value, *args: value,
    sympy.matrices.expressions.special.Identity: lambda n: torch.eye(n, dtype=torch.float64),
    sympy.matrices.expressions.special.ZeroMatrix: lambda m, n: torch.zeros((m, n), dtype=torch.float64),
    sympy.matrices.expressions.special.OneMatrix: lambda m, n: torch.ones((m, n), dtype=torch.float64),
    # Array operations
    sympy.tensor.array.expressions.array_expressions.ArrayAdd: array_add,
    sympy.tensor.array.expressions.array_expressions.PermuteDims: permute_dims,
    sympy.tensor.array.expressions.array_expressions.ArrayTensorProduct: array_tensor_product,
    sympy.tensor.array.expressions.array_expressions.ArrayContraction: array_contraction,
    sympy.tensor.array.expressions.array_expressions.ArrayDiagonal: array_diagonal,
    # Other operations
    sympy.Transpose: lambda x: x.transpose(-2, -1),
    sympy.matrices.expressions.diagonal.DiagonalMatrix: torch.diag,
    sympy.matrices.expressions.hadamard.HadamardProduct: torch.mul,
    sympy.matrices.expressions.hadamard.HadamardPower: torch.pow,
    sympy.matrices.expressions.determinant.Determinant: torch.linalg.det,
    DFT: torch.fft.fft,
    sympy.core.symbol.Str: lambda value, *args: value
}


number_symbols = [cls for cls in sympy.NumberSymbol.__subclasses__()]


_TORCH_FUNCTION_MAP.update({
    s: (s.__name__, ft.partial(lambda val, *args: torch.tensor(float(val)), s()))
    for s in number_symbols
})


_TORCHMODULE_ONLY_FUNCTION_MAP .update(
    {s: ft.partial(lambda val, *args: torch.tensor(float(val)), s())
     for s in number_symbols}
)


class TorchPrinter(ArrayPrinter, AbstractPythonCodePrinter):
    printmethod = "_torchcode"

    _default_settings = dict(
        AbstractPythonCodePrinter._default_settings,
        torch_version=None,
        requires_grad=False,
        dtype="torch.float64",
    )

    def __init__(self, settings=None):
        super().__init__(settings)

        version = self._settings['torch_version']
        self.requires_grad = self._settings['requires_grad']
        self.dtype = self._settings['dtype']
        if version is None and torch:
            version = torch.__version__
        self.torch_version = version

    def _get_op(self, expr_type):
        return _TORCH_FUNCTION_MAP.get(expr_type, (None, None))

    def _print_Function(self, expr):
        op_str, _ = self._get_op(expr.func)
        if op_str is None:
            return super()._print_Basic(expr)
        children = [self._print(arg) for arg in expr.args]
        if len(children) == 1:
            return f"{self._module_format(op_str)}({children[0]})"
        return self._expand_fold_binary_op(op_str, children)

    # mirrors the tensorflow version
    _print_Expr = _print_Function
    _print_Application = _print_Function
    _print_MatrixExpr = _print_Function
    _print_Relational = _print_Function
    _print_Not = _print_Function
    _print_And = _print_Function
    _print_Or = _print_Function
    _print_HadamardProduct = _print_Function
    _print_Trace = _print_Function
    _print_Determinant = _print_Function

    def _print_Inverse(self, expr):
        return '{}({})'.format(self._module_format("torch.linalg.inv"),
                               self._print(expr.args[0]))

    def _print_Transpose(self, expr):
        if expr.arg.is_Matrix and expr.arg.shape[0] == expr.arg.shape[1]:
            # For square matrices, we can use the .t() method
            return "{}({}).t()".format("torch.transpose", self._print(expr.arg))
        else:
            # For non-square matrices or more general cases
            # transpose first and second dimensions (typical matrix transpose)
            return "{}.permute({})".format(
                self._print(expr.arg),
                ", ".join([str(i) for i in range(len(expr.arg.shape))])[::-1]
            )

    def _print_PermuteDims(self, expr):
        return "%s.permute(%s)" % (
            self._print(expr.expr),
            ", ".join(str(i) for i in expr.permutation.array_form)
        )

    def _print_Derivative(self, expr):
        # this version handles multi-variable and mixed partial derivatives. The tensorflow version does not.
        variables = expr.variables
        expr_arg = expr.expr

        # Handle multi-variable or repeated derivatives
        if len(variables) > 1 or (
            len(variables) == 1 and not isinstance(variables[0], tuple) and variables.count(variables[0]) > 1):
            result = self._print(expr_arg)
            var_groups = {}

            # Group variables by base symbol
            for var in variables:
                if isinstance(var, tuple):
                    base_var, order = var
                    var_groups[base_var] = var_groups.get(base_var, 0) + order
                else:
                    var_groups[var] = var_groups.get(var, 0) + 1

            # Apply gradients in sequence
            for var, order in var_groups.items():
                for _ in range(order):
                    result = "torch.autograd.grad({}, {}, create_graph=True)[0]".format(result, self._print(var))
            return result

        # Handle single variable case
        if len(variables) == 1:
            variable = variables[0]
            if isinstance(variable, tuple) and len(variable) == 2:
                base_var, order = variable
                if not isinstance(order, Integer): raise NotImplementedError("Only integer orders are supported")
                result = self._print(expr_arg)
                for _ in range(order):
                    result = "torch.autograd.grad({}, {}, create_graph=True)[0]".format(result, self._print(base_var))
                return result
            return "torch.autograd.grad({}, {})[0]".format(self._print(expr_arg), self._print(variable))

        return self._print(expr_arg)  # Empty variables case

    def _print_Piecewise(self, expr):
        from sympy import Piecewise
        e, cond = expr.args[0].args
        if len(expr.args) == 1:
            return '{}({}, {}, {})'.format(
                self._module_format("torch.where"),
                self._print(cond),
                self._print(e),
                0)

        return '{}({}, {}, {})'.format(
            self._module_format("torch.where"),
            self._print(cond),
            self._print(e),
            self._print(Piecewise(*expr.args[1:])))

    def _print_Pow(self, expr):
        # XXX May raise error for
        # int**float or int**complex or float**complex
        base, exp = expr.args
        if expr.exp == S.Half:
            return "{}({})".format(
                self._module_format("torch.sqrt"), self._print(base))
        return "{}({}, {})".format(
            self._module_format("torch.pow"),
            self._print(base), self._print(exp))

    def _print_MatMul(self, expr):
        # Separate matrix and scalar arguments
        mat_args = [arg for arg in expr.args if isinstance(arg, MatrixExpr)]
        args = [arg for arg in expr.args if arg not in mat_args]
        # Handle scalar multipliers if present
        if args:
            return "%s*%s" % (
                self.parenthesize(Mul.fromiter(args), PRECEDENCE["Mul"]),
                self._expand_fold_binary_op("torch.matmul", mat_args)
            )
        else:
            return self._expand_fold_binary_op("torch.matmul", mat_args)

    def _print_MatPow(self, expr):
        return self._expand_fold_binary_op("torch.mm", [expr.base]*expr.exp)

    def _print_MatrixBase(self, expr):
        data = "[" + ", ".join(["[" + ", ".join([self._print(j) for j in i]) + "]" for i in expr.tolist()]) + "]"
        params = [str(data)]
        params.append(f"dtype={self.dtype}")
        if self.requires_grad:
            params.append("requires_grad=True")

        return "{}({})".format(
            self._module_format("torch.tensor"),
            ", ".join(params)
        )

    def _print_isnan(self, expr):
        return f'torch.isnan({self._print(expr.args[0])})'

    def _print_isinf(self, expr):
        return f'torch.isinf({self._print(expr.args[0])})'

    def _print_Identity(self, expr):
        if all(dim.is_Integer for dim in expr.shape):
            return "{}({})".format(
                self._module_format("torch.eye"),
                self._print(expr.shape[0])
            )
        else:
            # For symbolic dimensions, fall back to a more general approach
            return "{}({}, {})".format(
                self._module_format("torch.eye"),
                self._print(expr.shape[0]),
                self._print(expr.shape[1])
            )

    def _print_ZeroMatrix(self, expr):
        return "{}({})".format(
            self._module_format("torch.zeros"),
            self._print(expr.shape)
        )

    def _print_OneMatrix(self, expr):
        return "{}({})".format(
            self._module_format("torch.ones"),
            self._print(expr.shape)
        )

    def _print_conjugate(self, expr):
        return f"{self._module_format('torch.conj')}({self._print(expr.args[0])})"

    def _print_ImaginaryUnit(self, expr):
        return "1j"  # uses the Python built-in 1j notation for the imaginary unit

    def _print_Heaviside(self, expr):
        args = [self._print(expr.args[0]), "0.5"]
        if len(expr.args) > 1:
            args[1] = self._print(expr.args[1])
        return f"{self._module_format('torch.heaviside')}({args[0]}, {args[1]})"

    def _print_gamma(self, expr):
        return f"{self._module_format('torch.special.gamma')}({self._print(expr.args[0])})"

    def _print_polygamma(self, expr):
        if expr.args[0] == S.Zero:
            return f"{self._module_format('torch.special.digamma')}({self._print(expr.args[1])})"
        else:
            raise NotImplementedError("PyTorch only supports digamma (0th order polygamma)")

    _module = "torch"
    _einsum = "einsum"
    _add = "add"
    _transpose = "t"
    _ones = "ones"
    _zeros = "zeros"


class _Node(torch.nn.Module):
    def __init__(self, expr, _memodict, _func_lookup, **kwargs):
        super().__init__(**kwargs)
        self._sympy_func = expr.func

        if issubclass(expr.func, sympy.Float):
            self._value = torch.nn.Parameter(torch.tensor(float(expr)))
            self._torch_func = lambda: self._value
            self._args = ()
        elif issubclass(expr.func, sympy.Integer):
            self.register_buffer("_value", torch.tensor(int(expr)))
            self._torch_func = lambda: self._value
            self._args = ()
        elif issubclass(expr.func, sympy.core.containers.Tuple):
            self._value = tuple(expr.args)
            self._torch_func = lambda: self._value
            self._args = ()
        elif issubclass(expr.func, sympy.Rational):
            self.register_buffer("_numerator", torch.tensor(expr.p))
            self.register_buffer("_denominator", torch.tensor(expr.q))
            self._torch_func = lambda: self._numerator / self._denominator
            self._args = ()
        elif issubclass(expr.func, sympy.combinatorics.permutations.Permutation):
            self._value = list(expr.array_form)
            self._torch_func = lambda: self._value
            self._args = ()
        elif issubclass(expr.func, sympy.UnevaluatedExpr):
            if len(expr.args) != 1 or not issubclass(expr.args[0].func, sympy.Float):
                raise ValueError("UnevaluatedExpr should only be used to wrap floats.")
            self.register_buffer("_value", torch.tensor(float(expr.args[0])))
            self._torch_func = lambda: self._value
            self._args = ()
        elif issubclass(expr.func, sympy.Symbol):
            self._name = expr.name
            self._torch_func = lambda value, *args: value
            self._args = ((lambda memodict: memodict[expr.name]),)
        elif issubclass(expr.func, sympy.matrices.expressions.matexpr.MatrixSymbol):
            self._name = expr.name
            self._torch_func = lambda value, *args: value
            self._args = ((lambda memodict: memodict[expr.name]),)
        elif issubclass(expr.func, sympy.core.numbers.ImaginaryUnit):
            self._torch_func = lambda: torch.tensor(1j)
            self._args = ()
        else:
            try:
                self._torch_func = _func_lookup[self._sympy_func]
            except KeyError:
                raise ValueError(f"Unsupported SymPy function: {self._sympy_func.__name__}")
            args = []
            for arg in expr.args:
                try:
                    arg_ = _memodict[arg]
                except KeyError:
                    arg_ = _Node(expr=arg, _memodict=_memodict, _func_lookup=_func_lookup, **kwargs)
                    _memodict[arg] = arg_
                args.append(arg_)
            self._args = torch.nn.ModuleList(args)

    def forward(self, memodict) -> torch.Tensor:
        args = []
        for arg in self._args:
            try:
                arg_ = memodict[arg]
            except KeyError:
                arg_ = arg(memodict)
                memodict[arg] = arg_
            args.append(arg_)
        return self._torch_func(*args)

    def to_sympy(self, _memodict):
        if issubclass(self._sympy_func, sympy.Float):
            return self._sympy_func(self._value.item())
        elif issubclass(self._sympy_func, sympy.UnevaluatedExpr):
            return self._sympy_func(self._value.item())
        elif issubclass(self._sympy_func, (type(sympy.S.NegativeOne), type(sympy.S.One), type(sympy.S.Zero))):
            return self._sympy_func()
        elif issubclass(self._sympy_func, sympy.Integer):
            return self._sympy_func(self._value.item())
        elif issubclass(self._sympy_func, sympy.Rational):
            if issubclass(self._sympy_func, type(sympy.S.Half)):
                return sympy.S.Half
            else:
                return self._sympy_func(self._numerator.item(), self._denominator.item())
        elif issubclass(self._sympy_func, sympy.Symbol):
            return self._sympy_func(self._name)
        elif issubclass(self._sympy_func, sympy.core.numbers.ImaginaryUnit):
            return sympy.I
        elif issubclass(self._sympy_func, sympy.core.numbers.NumberSymbol):
            return self._sympy_func()
        else:
            if issubclass(self._sympy_func, (sympy.Min, sympy.Max)):
                evaluate = False
            else:
                evaluate = True
            args = []
            for arg in self._args:
                try:
                    arg_ = _memodict[arg]
                except KeyError:
                    arg_ = arg.to_sympy(_memodict)
                    _memodict[arg] = arg_
                args.append(arg_)
            return self._sympy_func(*args, evaluate=evaluate)


class SympyModule(torch.nn.Module):
    _default_settings = {
        'torch_version': None,
    }

    def __init__(self, args, expr, extra_funcs=None, **settings):
        super().__init__()
        self._settings = self._default_settings.copy()
        if settings:
            self._settings.update(settings)

        version = self._settings['torch_version']
        if version is None and torch:
            version = torch.__version__
        self.torch_version = version

        self.args = tuple(args)
        self.arg_names = [arg.name for arg in self.args]
        self.expressions = tuple(expr) if isinstance(expr, (tuple, list, MutableDenseMatrix)) else (expr,)

        func_lookup = {k: v[1] if isinstance(v, tuple) else v for k, v in _TORCH_FUNCTION_MAP.items()}
        func_lookup = {**func_lookup, **_TORCHMODULE_ONLY_FUNCTION_MAP}
        if extra_funcs:
            func_lookup.update(extra_funcs)
        _memodict = {}
        self._nodes = torch.nn.ModuleList([
            _Node(expr, _memodict, func_lookup) for expr in self.expressions
        ])
        self._expressions_string = str(self.expressions)

    def __repr__(self):
        return f"{type(self).__name__}(expressions={self._expressions_string})"

    def forward(self, *inputs, **kwargs):
        memodict = dict(kwargs)

        for i, input_tensor in enumerate(inputs):
            if i < len(self.arg_names):
                memodict[self.arg_names[i]] = input_tensor

        missing = [name for name in self.arg_names if name not in memodict]
        if missing:
            raise ValueError(f"Missing inputs for symbols: {missing}")

        out = [node(memodict) for node in self._nodes]

        if len(out) == 1:
            return out[0]

        out = torch.broadcast_tensors(*out)
        return torch.stack(out, dim=-1)

    def to_sympy(self):
        _memodict = {}
        return [node.to_sympy(_memodict) for node in self._nodes]


def torch_code(expr, requires_grad=False, dtype="torch.float64", **settings):
    printer = TorchPrinter(settings={'requires_grad': requires_grad, 'dtype': dtype})
    return printer.doprint(expr, **settings)


def torch_module(args, expr, extra_funcs=None, **settings):
    if isinstance(args, sympy.Symbol):
        args = (args,)
    return SympyModule(args, expr, extra_funcs, **settings)
