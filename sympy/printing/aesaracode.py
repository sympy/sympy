# sympy/printing/aesaracode.py
from __future__ import annotations

import math
from functools import partial
from typing import Any

import sympy
from sympy.external import import_module
from sympy.printing.printer import Printer
from sympy.utilities.exceptions import sympy_deprecation_warning
from sympy.utilities.iterables import is_sequence

# ---------------------------------------------------------------------------
# Optional dependency handling
# ---------------------------------------------------------------------------
aesara = import_module("aesara")  # -> module or None

# Keep for doctest runners that honor it (nose, xdoctest, etc.)
__doctest_requires__ = {
    "aesara_function": ["aesara"],
    "aesara_code": ["aesara"],
    "sympy.printing.aesaracode.aesara_function": ["aesara"],
    "sympy.printing.aesaracode.aesara_code": ["aesara"],
}

# ---------------------------------------------------------------------------
# Aesara-dependent bindings
# ---------------------------------------------------------------------------
if aesara:
    aes = aesara.scalar
    aet = aesara.tensor
    from aesara.tensor import nlinalg
    from aesara.tensor.elemwise import Elemwise, DimShuffle

    # Aesara 2.8.11 renamed true_div -> true_divide
    true_divide = getattr(aet, "true_divide", None) or aet.true_div

    _MAPPING = {
        sympy.Add:        aet.add,
        sympy.Mul:        aet.mul,
        sympy.Abs:        aet.abs,
        sympy.sign:       aet.sgn,
        sympy.ceiling:    aet.ceil,
        sympy.floor:      aet.floor,
        sympy.log:        aet.log,
        sympy.exp:        aet.exp,
        sympy.sqrt:       aet.sqrt,
        sympy.cos:        aet.cos,
        sympy.acos:       aet.arccos,
        sympy.sin:        aet.sin,
        sympy.asin:       aet.arcsin,
        sympy.tan:        aet.tan,
        sympy.atan:       aet.arctan,
        sympy.atan2:      aet.arctan2,
        sympy.cosh:       aet.cosh,
        sympy.acosh:      aet.arccosh,
        sympy.sinh:       aet.sinh,
        sympy.asinh:      aet.arcsinh,
        sympy.tanh:       aet.tanh,
        sympy.atanh:      aet.arctanh,
        sympy.re:         aet.real,
        sympy.im:         aet.imag,
        sympy.arg:        aet.angle,
        sympy.erf:        aet.erf,
        sympy.gamma:      aet.gamma,
        sympy.loggamma:   aet.gammaln,
        sympy.Pow:        aet.pow,
        sympy.Eq:         aet.eq,
        sympy.StrictGreaterThan: aet.gt,
        sympy.StrictLessThan:    aet.lt,
        sympy.LessThan:          aet.le,
        sympy.GreaterThan:       aet.ge,
        sympy.And:        aet.bitwise_and,
        sympy.Or:         aet.bitwise_or,
        sympy.Not:        aet.invert,
        sympy.Xor:        aet.bitwise_xor,
        sympy.Max:        aet.maximum,
        sympy.Min:        aet.minimum,
        sympy.conjugate:  aet.conj,
        sympy.core.numbers.ImaginaryUnit: lambda: aet.complex(0, 1),
        # Matrices
        sympy.MatAdd:         Elemwise(aes.add),
        sympy.HadamardProduct: Elemwise(aes.mul),
        sympy.Trace:          nlinalg.trace,
        sympy.Determinant:    nlinalg.det,
        sympy.Inverse:        nlinalg.matrix_inverse,
        sympy.Transpose:      DimShuffle((False, False), [1, 0]),
    }

else:
    # Stubs to keep module importable; APIs will raise when called.
    aes = aet = nlinalg = None
    true_divide = None
    _MAPPING: dict[type, Any] = {}

# ---------------------------------------------------------------------------
# Printer
# ---------------------------------------------------------------------------
class _AesaraMissingError(ImportError):
    def __init__(self, api: str):
        super().__init__(f"Aesara is required for {api}")

class AesaraPrinter(Printer):
    """
    .. deprecated:: 1.14
       The Aesara code printer is deprecated. See :ref:`deprecated-aesaraprinter`.

    Code printer producing Aesara symbolic graphs.
    """
    printmethod = "_aesara"

    def __init__(self, *args, **kwargs):
        if not aesara:
            raise _AesaraMissingError("AesaraPrinter")
        self.cache = kwargs.pop("cache", {})
        super().__init__(*args, **kwargs)

    # ---- cache helpers ----
    def _get_key(self, s, name=None, dtype=None, broadcastable=None):
        if name is None:
            name = s.name
        return (name, type(s), s.args, dtype, broadcastable)

    def _get_or_create(self, s, name=None, dtype=None, broadcastable=None):
        if not aesara:
            raise _AesaraMissingError("AesaraPrinter")
        if name is None:
            name = getattr(s, "name", str(s))
        if dtype is None:
            dtype = "floatX"
        if broadcastable is None:
            broadcastable = ()
        key = self._get_key(s, name, dtype=dtype, broadcastable=broadcastable)
        if key in self.cache:
            return self.cache[key]
        value = aet.tensor(name=name, dtype=dtype, shape=broadcastable)
        self.cache[key] = value
        return value

    # ---- node printers ----
    def _print_Symbol(self, s, **kw):
        return self._get_or_create(
            s,
            dtype=kw.get("dtypes", {}).get(s),
            broadcastable=kw.get("broadcastables", {}).get(s),
        )

    def _print_AppliedUndef(self, s, **kw):
        return self._get_or_create(
            s,
            name=f"{type(s)}_{s.args[0]}",
            dtype=kw.get("dtypes", {}).get(s),
            broadcastable=kw.get("broadcastables", {}).get(s),
        )

    def _print_Basic(self, expr, **kw):
        op = _MAPPING[type(expr)]
        return op(*[self._print(arg, **kw) for arg in expr.args])

    def _print_Number(self, n, **_):
        return float(n.evalf())

    def _print_MatrixSymbol(self, X, **kw):
        return self._get_or_create(
            X,
            dtype=kw.get("dtypes", {}).get(X),
            broadcastable=(None, None),
        )

    def _print_DenseMatrix(self, X, **kw):
        if not hasattr(aet, "stacklists"):
            raise NotImplementedError(
                "Matrix translation not yet supported in this Aesara version"
            )
        return aet.stacklists([[self._print(a, **kw) for a in row] for row in X.tolist()])

    _print_ImmutableMatrix = _print_ImmutableDenseMatrix = _print_DenseMatrix

    def _print_MatMul(self, expr, **kw):
        children = [self._print(a, **kw) for a in expr.args]
        out = children[0]
        for c in children[1:]:
            out = aet.dot(out, c)
        return out

    def _print_MatPow(self, expr, **kw):
        base, power = [self._print(a, **kw) for a in expr.args]
        if isinstance(power, int) and power >= 0:
            out = 1
            for _ in range(power):
                out = aet.dot(out, base)
            return out
        raise NotImplementedError("Only non-negative integer matrix powers are supported")

    def _print_MatrixSlice(self, expr, **kw):
        parent = self._print(expr.parent, **kw)
        rs = self._print(slice(*expr.rowslice), **kw)
        cs = self._print(slice(*expr.colslice), **kw)
        return parent[rs, cs]

    def _print_BlockMatrix(self, expr, **kw):
        nrows, ncols = expr.blocks.shape
        blocks = [[self._print(expr.blocks[r, c], **kw) for c in range(ncols)] for r in range(nrows)]
        return aet.join(0, *[aet.join(1, *row) for row in blocks])

    def _print_slice(self, expr, **kw):
        return slice(
            *[
                self._print(i, **kw) if isinstance(i, sympy.Basic) else i
                for i in (expr.start, expr.stop, expr.step)
            ]
        )

    def _print_Pi(self, *_):
        return math.pi

    def _print_Piecewise(self, expr, **kw):
        import numpy as np
        e, cond = expr.args[0].args
        pe, pc = self._print(e, **kw), self._print(cond, **kw)
        if len(expr.args) == 1:
            return aet.switch(pc, pe, np.nan)
        rest = self._print(sympy.Piecewise(*expr.args[1:]), **kw)
        return aet.switch(pc, pe, rest)

    def _print_Rational(self, expr, **kw):
        return true_divide(self._print(expr.p, **kw), self._print(expr.q, **kw))

    def _print_Integer(self, expr, **_):
        return expr.p

    def _print_factorial(self, expr, **kw):
        return self._print(sympy.gamma(expr.args[0] + 1), **kw)

    def _print_Derivative(self, deriv, **kw):
        from aesara.gradient import Rop
        rv = self._print(deriv.expr, **kw)
        for var in deriv.variables:
            v = self._print(var, **kw)
            rv = Rop(rv, v, aet.ones_like(v))
        return rv

    def emptyPrinter(self, expr):
        return expr

    def doprint(self, expr, dtypes=None, broadcastables=None):
        if dtypes is None:
            dtypes = {}
        if broadcastables is None:
            broadcastables = {}
        return self._print(expr, dtypes=dtypes, broadcastables=broadcastables)

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
global_cache: dict[Any, Any] = {}

def aesara_code(expr, cache=None, **kwargs):
    """
    Convert a SymPy expression into an Aesara graph variable.

    Examples (guarded so plain doctest doesn't fail if Aesara is absent)
    ------------------------------------------------------------------------------
    >>> from sympy.abc import x
    >>> from sympy.external import import_module
    >>> _aes = import_module('aesara')
    >>> if _aes:
    ...     from sympy.printing.aesaracode import aesara_code
    ...     v = aesara_code(x + 1, cache={})
    ...     type(v).__name__ in {'TensorVariable', 'TensorConstant'}
    True
    """
    sympy_deprecation_warning(
        "The aesara_code function is deprecated.",
        deprecated_since_version="1.14",
        active_deprecations_target="deprecated-aesaraprinter",
    )
    if not aesara:
        raise _AesaraMissingError("aesara_code")
    if cache is None:
        cache = global_cache
    return AesaraPrinter(cache=cache, settings={}).doprint(expr, **kwargs)

def dim_handling(inputs, dim=None, dims=None, broadcastables=None):
    if dim is not None:
        return dict.fromkeys(inputs, (False,) * dim)
    if dims is not None:
        maxdim = max(dims.values())
        return {s: (False,) * d + (True,) * (maxdim - d) for s, d in dims.items()}
    if broadcastables is not None:
        return broadcastables
    return {}

def aesara_function(
    inputs,
    outputs,
    scalar=False,
    *,
    dim=None,
    dims=None,
    broadcastables=None,
    **kwargs,
):
    """
    Create an Aesara function from SymPy expressions.

    Examples (guarded)
    ------------------
    >>> from sympy.abc import x, y, z
    >>> from sympy.external import import_module
    >>> _aes = import_module('aesara')
    >>> if _aes:
    ...     from sympy.printing.aesaracode import aesara_function
    ...     f1 = aesara_function([x], [x**2 - 1], scalar=True)
    ...     f1(3)
    8.0
    >>> if _aes:
    ...     f2 = aesara_function([x, y, z], [(x**z + y**z)**(1/z)], scalar=True)
    ...     f2(3, 4, 2)
    5.0
    >>> if _aes:
    ...     f3 = aesara_function([x, y], [x**2 + y**2, x**2 - y**2], scalar=True)
    ...     f3(2, 3)
    [13.0, -5.0]
    """
    sympy_deprecation_warning(
        "The aesara_function function is deprecated.",
        deprecated_since_version="1.14",
        active_deprecations_target="deprecated-aesaraprinter",
    )
    if not aesara:
        raise _AesaraMissingError("aesara_function")

    cache = kwargs.pop("cache", {})
    dtypes = kwargs.pop("dtypes", {})

    broadcastables = dim_handling(
        inputs, dim=dim, dims=dims, broadcastables=broadcastables
    )

    code = partial(aesara_code, cache=cache, dtypes=dtypes, broadcastables=broadcastables)
    tinputs = list(map(code, inputs))
    toutputs = list(map(code, outputs))

    # Ensure constants are tensors
    from aesara.graph.basic import Variable as AesaraVariable
    toutputs = [
        out if isinstance(out, AesaraVariable) else aet.as_tensor_variable(out)
        for out in toutputs
    ]
    if len(toutputs) == 1:
        toutputs = toutputs[0]

    func = aesara.function(tinputs, toutputs, **kwargs)
    is_0d = [len(o.variable.broadcastable) == 0 for o in func.outputs]

    if not scalar or not any(is_0d):
        func.aesara_function = func
        return func

    def wrapper(*args):
        out = func(*args)
        if is_sequence(out):
            return [o[()] if is_0d[i] else o for i, o in enumerate(out)]
        return out[()]

    wrapper.__wrapped__ = func
    wrapper.__doc__ = func.__doc__
    wrapper.aesara_function = func
    return wrapper

__all__ = [
    "AesaraPrinter",
    "aesara_code",
    "aesara_function",
    "dim_handling",
]
