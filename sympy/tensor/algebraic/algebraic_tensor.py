from __future__ import annotations

from sympy.core.mul import Mul
from sympy.core.numbers import Number
from sympy.core.sympify import sympify


class ZeroTensor:
    """Zero tensor carrying a specific tensor shape.

    A ZeroTensor of shape ``(m, n)`` acts as the additive identity for all
    AlgebraicTensors and PureTensors of the same shape.  ZeroTensors of
    different shapes belong to different tensor spaces and are not summable.
    """

    __slots__ = ("_shape",)
    is_ZeroTensor = True

    def __init__(self, shape):
        self._shape = tuple(shape)

    @property
    def shape(self):
        return self._shape

    @property
    def tensor_shape(self):
        return self._shape

    def __neg__(self):
        return self

    def __add__(self, other):
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        return AlgebraicTensor(self, other)

    def __radd__(self, other):
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        return AlgebraicTensor(other, self)

    def __sub__(self, other):
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        return AlgebraicTensor(self, -other)

    def __rsub__(self, other):
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        return AlgebraicTensor(other, -self)

    def __eq__(self, other):
        if isinstance(other, ZeroTensor):
            return self._shape == other._shape
        return NotImplemented

    def __hash__(self):
        return hash(("ZeroTensor", self._shape))

    def __repr__(self):
        return f"ZeroTensor{self._shape}"

    def __str__(self):
        return f"0_{self._shape}"

    def __bool__(self):
        return False


class PureTensor(Mul):
    """Pure tensor as an unevaluated (non-commutative) tensor product of factors.

    Extends SymPy's non-commutative Mul so that all existing simplification,
    differentiation, and pattern-matching machinery is available out of the
    box.  The tensor-product operator is non-commutative: order of factors
    is always preserved.

    Each factor must carry ``.shape`` (e.g. any ``MatrixExpr`` or a 1-x1
    wrapper around a non-commutative Symbol).
    """

    __slots__ = ()

    is_PureTensor = True

    _eval_is_commutative = lambda self: False

    @property
    def factors(self):
        """Individual tensor-product factors in left-to-right order."""
        return self.args

    @property
    def num_factors(self):
        return len(self.args)

    @property
    def tensor_shape(self):
        """Outer-product shape: (dim_0-of-first, dim_1-of-last)."""
        shapes = [f.shape for f in self.factors]
        return (shapes[0][0], shapes[-1][1])

    def __str__(self):
        return " \u2297 ".join(str(f) for f in self.factors)

    def __repr__(self):
        return f"PureTensor({', '.join(repr(f) for f in self.factors)})"

    def __new__(cls, *args, evaluate=False):
        if not args:
            raise ValueError("PureTensor requires at least one factor")

        processed = []
        for a in args:
            a = sympify(a)
            if isinstance(a, Number):
                raise TypeError(
                    f"Scalar numbers are not valid tensor factors "
                    f"(use a 1x1 MatrixExpr wrapper): {a}"
                )
            if not hasattr(a, "shape"):
                raise TypeError(
                    f"Tensor factor must have a .shape attribute, "
                    f"got {type(a).__name__}"
                )
            processed.append(a)

        if len(processed) == 1:
            return processed[0]

        obj = Mul.__new__(cls, *processed, evaluate=False)
        return obj

    def __neg__(self):
        """Return Mul(-1, self) as a regular Mul, not a PureTensor."""
        from sympy.core.singleton import S
        return Mul(S.NegativeOne, self, evaluate=False)

    def __mul__(self, other):
        """PureTensor * scalar -> Mul(scalar, PureTensor)."""
        from sympy.core.singleton import S
        other = sympify(other)
        if isinstance(other, Number):
            if other is S.One:
                return self
            if other is S.Zero:
                return ZeroTensor(self.tensor_shape)
            return Mul(other, self, evaluate=False)
        return Mul(self, other)

    def __rmul__(self, other):
        """scalar * PureTensor -> Mul(scalar, PureTensor)."""
        from sympy.core.singleton import S
        if other is S.Zero:
            return ZeroTensor(self.tensor_shape)
        if other is S.One:
            return self
        other = sympify(other)
        if isinstance(other, Number):
            if other is S.Zero:
                return ZeroTensor(self.tensor_shape)
            return Mul(other, self, evaluate=False)
        return Mul(other, self)

    def __add__(self, other):
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        return AlgebraicTensor(self, other)

    def __radd__(self, other):
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        return AlgebraicTensor(other, self)

    def __sub__(self, other):
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        return AlgebraicTensor(self, -other)

    def __rsub__(self, other):
        from sympy.tensor.algebraic.algebraic_tensor import AlgebraicTensor
        return AlgebraicTensor(other, -self)


def tensor_product(*args):
    """Convenience constructor for PureTensor.

    >>> from sympy.tensor.algebraic.algebraic_tensor import tensor_product
    >>> from sympy.matrices.expressions import MatrixSymbol
    >>> A = MatrixSymbol("A", 2, 3)
    >>> v = MatrixSymbol("v", 3, 1)
    >>> tensor_product(A, v)
    A \u2297 v
    """
    return PureTensor(*args)


def zero_tensor(shape):
    """Convenience constructor for ZeroTensor.

    Parameters
    ----------
    shape : tuple
        A pair ``(rows, cols)`` describing the tensor shape.

    Returns
    -------
    ZeroTensor
    """
    return ZeroTensor(shape)


# ---- Heavy imports (after all early-defined classes) ----

from sympy.core.add import Add, add
from sympy.core.basic import Basic
from sympy.core.singleton import S


class ShapeMismatchError(TypeError):
    """Raised when attempting to add tensors of incompatible shapes."""
    pass


class AlgebraicTensor(Basic):
    """Sum of PureTensors (and/or ZeroTensors) sharing the same tensor shape.

    Built on top of ``Basic`` (not ``Add``) to avoid the is_commutative
    descriptor conflict in ``AssocOp._from_args``, while still exposing
    ``is_Add = True`` so that the wider SymPy ecosystem recognises this
    as an additive expression.

    Shape enforcement and ZeroTensor handling are done in ``__new__``.
    The flattening logic is kept intentionally simple so that a future
    ``.simplify`` extension can layer noncommutative-polynomial factoring
    on top without fighting against ``Add.flatten`` internals.
    """

    __slots__ = ()

    is_AlgebraicTensor = True
    is_Add = True

    identity = None  # no single identity; use ZeroTensor(shape) instead

    def __new__(cls, *args, evaluate=True, _sympify=True):
        if not args:
            raise ValueError("AlgebraicTensor requires at least one argument")

        if _sympify:
            args = tuple(
                a if isinstance(a, ZeroTensor) else sympify(a)
                for a in args
            )

        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, (PureTensor, ZeroTensor)):
                return arg
            if isinstance(arg, AlgebraicTensor):
                return arg
            if isinstance(arg, (list, tuple)) and len(arg) == 1:
                return arg[0]

        # Flatten nested AlgebraicTensors and collect all leaf terms
        flat, shape, zero_term = cls._flatten_args(args)

        if not flat:
            if zero_term is not None:
                return zero_term
            raise ValueError("AlgebraicTensor resulted in zero terms")

        # Remove identity Numbers (S.Zero) from flat list
        flat = [a for a in flat if a is not S.Zero]

        if not flat:
            if zero_term is not None:
                return zero_term
            raise ValueError("AlgebraicTensor resulted in zero terms")

        # If a ZeroTensor was provided, always keep it as an anchor.
        # Only unwrap to a single term when there is exactly one non-zero
        # term AND no ZeroTensor was explicitly given.
        if len(flat) == 1 and zero_term is None and not isinstance(flat[0], ZeroTensor):
            return flat[0]

        # Append ZeroTensor if present (for shape anchoring)
        if zero_term is not None:
            flat.append(zero_term)

        obj = Basic.__new__(cls, *flat)
        return obj

    @classmethod
    def _flatten_args(cls, args):
        """Flatten nested AlgebraicTensors, validate shapes, collect ZeroTensors.

        Returns
        -------
        (terms : list, shape : tuple | None, zero_term : ZeroTensor | None)
        """
        shape = None
        zero_term = None
        terms = []

        work = list(args)
        while work:
            o = work.pop()

            if isinstance(o, ZeroTensor):
                if shape is None:
                    shape = o.shape
                elif o.shape != shape:
                    raise ShapeMismatchError(
                        f"Cannot add tensors of different shapes: "
                        f"{shape} vs {o.shape}"
                    )
                zero_term = o
                continue

            if isinstance(o, AlgebraicTensor):
                work.extend(o.args)
                continue

            if isinstance(o, PureTensor):
                if shape is None:
                    shape = o.tensor_shape
                elif o.tensor_shape != shape:
                    raise ShapeMismatchError(
                        f"Cannot add tensors of different shapes: "
                        f"{shape} vs {o.tensor_shape}"
                    )
                terms.append(o)
                continue

            if isinstance(o, Mul) and not isinstance(o, PureTensor):
                # coeff * PureTensor pattern
                for arg in o.args:
                    if isinstance(arg, PureTensor):
                        if shape is None:
                            shape = arg.tensor_shape
                        elif arg.tensor_shape != shape:
                            raise ShapeMismatchError(
                                f"Cannot add tensors of different shapes: "
                                f"{shape} vs {arg.tensor_shape}"
                            )
                        break
                terms.append(o)
                continue

            if o.is_Number:
                terms.append(o)
                continue

            # Catch-all: treat as a regular term
            terms.append(o)

        return terms, shape, zero_term

    @property
    def tensor_shape(self):
        """Shape shared by every term in this sum."""
        for arg in self.args:
            if isinstance(arg, ZeroTensor):
                return arg.shape
            if isinstance(arg, PureTensor):
                return arg.tensor_shape
            if isinstance(arg, Mul) and not isinstance(arg, PureTensor):
                for f in arg.args:
                    if isinstance(f, PureTensor):
                        return f.tensor_shape
        for arg in self.args:
            if hasattr(arg, "tensor_shape"):
                return arg.tensor_shape
            if hasattr(arg, "shape"):
                return arg.shape
        raise AttributeError("Cannot determine tensor_shape")

    @property
    def terms(self):
        """Non-zero, non-coefficient terms in this sum."""
        return tuple(
            a for a in self.args
            if not isinstance(a, Number) and not isinstance(a, ZeroTensor)
        )

    def has_zero_term(self):
        """Return True if a ZeroTensor anchors this sum."""
        return any(isinstance(a, ZeroTensor) for a in self.args)

    # ---- Arithmetic delegation ----

    def __neg__(self):
        """Negate each term in the sum."""
        return AlgebraicTensor(*(-a for a in self.args))

    def __add__(self, other):
        return AlgebraicTensor(self, other)

    def __radd__(self, other):
        return AlgebraicTensor(other, self)

    def __sub__(self, other):
        return AlgebraicTensor(self, -other)

    def __rsub__(self, other):
        return AlgebraicTensor(other, -self)

    # ---- Factorization helpers (infrastructure for future .simplify) ----

    def as_common_left(self):
        """Extract common left factors from all PureTensor terms.

        Returns
        -------
        (left_factors, rest, right_factors)
            *left_factors* is a tuple of shared leading factors.
            *rest* is an AlgebraicTensor of the middle parts (or a single
            PureTensor / sum of PureTensors when there is only one term).
            *right_factors* is always () for left extraction.

        If no common left factors exist, returns ((), self, ()).
        """
        pure_terms = [
            a for a in self.args
            if isinstance(a, PureTensor) or (isinstance(a, Mul) and
                any(isinstance(f, PureTensor) for f in a.args))
        ]
        if not pure_terms:
            return ((), self, ())

        def _pure_part(expr):
            if isinstance(expr, PureTensor):
                return expr.factors
            if isinstance(expr, Mul):
                for f in expr.args:
                    if isinstance(f, PureTensor):
                        return f.factors
            return ()

        factor_lists = [_pure_part(t) for t in pure_terms]
        if not factor_lists[0]:
            return ((), self, ())

        common_len = 0
        for pos in range(min(len(fl) for fl in factor_lists)):
            refs = factor_lists[0][pos]
            if all(fl[pos] == refs for fl in factor_lists):
                common_len += 1
            else:
                break

        if common_len == 0:
            return ((), self, ())

        left_factors = factor_lists[0][:common_len]
        middles = []
        for t, fl in zip(pure_terms, factor_lists):
            mid = fl[common_len:]
            if isinstance(t, Mul) and not isinstance(t, PureTensor):
                coeff = S.One
                for f in t.args:
                    if isinstance(f, PureTensor):
                        coeff = t / f
                        break
                if mid:
                    middles.append(coeff * PureTensor(*mid))
                else:
                    middles.append(coeff)
            else:
                if mid:
                    middles.append(PureTensor(*mid))
                else:
                    middles.append(S.One)

        if len(middles) == 1:
            rest = middles[0]
        else:
            rest = AlgebraicTensor(*middles)

        return (left_factors, rest, ())

    def as_common_right(self):
        """Extract common right factors from all PureTensor terms.

        Symmetric to :meth:`as_common_left` but for trailing factors.

        Returns
        -------
        (left_factors, rest, right_factors)
        """
        pure_terms = [
            a for a in self.args
            if isinstance(a, PureTensor) or (isinstance(a, Mul) and
                any(isinstance(f, PureTensor) for f in a.args))
        ]
        if not pure_terms:
            return ((), self, ())

        def _pure_part(expr):
            if isinstance(expr, PureTensor):
                return expr.factors
            if isinstance(expr, Mul):
                for f in expr.args:
                    if isinstance(f, PureTensor):
                        return f.factors
            return ()

        factor_lists = [_pure_part(t) for t in pure_terms]
        if not factor_lists[0]:
            return ((), self, ())

        common_len = 0
        for pos in range(1, min(len(fl) for fl in factor_lists) + 1):
            refs = factor_lists[0][-pos]
            if all(fl[-pos] == refs for fl in factor_lists):
                common_len += 1
            else:
                break

        if common_len == 0:
            return ((), self, ())

        right_factors = factor_lists[0][-common_len:]
        middles = []
        for t, fl in zip(pure_terms, factor_lists):
            mid = fl[:-common_len] if common_len < len(fl) else ()
            if isinstance(t, Mul) and not isinstance(t, PureTensor):
                coeff = S.One
                for f in t.args:
                    if isinstance(f, PureTensor):
                        coeff = t / f
                        break
                if mid:
                    middles.append(coeff * PureTensor(*mid))
                else:
                    middles.append(coeff)
            else:
                if mid:
                    middles.append(PureTensor(*mid))
                else:
                    middles.append(S.One)

        if len(middles) == 1:
            rest = middles[0]
        else:
            rest = AlgebraicTensor(*middles)

        return ((), rest, right_factors)

    def as_common_factors(self):
        """Extract both common left and right factors simultaneously.

        Returns
        -------
        (left_factors, middle_expr, right_factors)
        """
        left_factors, rest, _ = self.as_common_left()
        if isinstance(rest, AlgebraicTensor):
            _, rest2, right_factors = rest.as_common_right()
        else:
            rest2, right_factors = rest, ()
        return (left_factors, rest2, right_factors)

    def __str__(self):
        return " + ".join(str(a) for a in self.args)

    def __repr__(self):
        return f"AlgebraicTensor({', '.join(repr(a) for a in self.args)})"


# Register AlgebraicTensor with the SymPy add dispatcher so that
# expressions like (A + B) where A/B involve PureTensor or AlgebraicTensor
# are routed through AlgebraicTensor.__new__ instead of Add.__new__.
add.register_handlerclass(
    (PureTensor, AlgebraicTensor), AlgebraicTensor
)
