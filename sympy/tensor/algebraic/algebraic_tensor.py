from __future__ import annotations

from sympy.core.add import Add, add
from sympy.core.basic import Basic
from sympy.core.mul import Mul
from sympy.core.numbers import Number
from sympy.core.singleton import S
from sympy.core.sympify import sympify

from sympy.tensor.algebraic.pure_tensor import PureTensor, _factor_shapes
from sympy.tensor.algebraic.zero_tensor import ZeroTensor, zero_tensor


# ---------------------------------------------------------------------------
# Shape helpers
# ---------------------------------------------------------------------------

def _normalize_shape(s):
    """Normalize any shape-like input to a tuple of (rows, cols) tuples.

    ``((3,4), (4,5))``  -> ``((3, 4), (4, 5))``  (unchanged)
    ``(3, 4)``          -> ``((3, 4),)``          (bare pair wrapped)
    ``[(3,4)]``         -> ``((3, 4),)``
    """
    if not s:
        return ()
    if len(s) == 2 and not isinstance(s[0], (tuple, list)):
        # Bare (m, n) -> ((m, n),)
        return (tuple(s),)
    return tuple(tuple(x) for x in s)


def _tensor_shape_of(expr):
    """Return the full tensor shape of *expr* as a tuple of factor shapes.

    Handles PureTensor, ZeroTensor, Mul(coeff, PureTensor), and bare
    matrix-like objects (whose single-factor shape is wrapped).
    """
    if isinstance(expr, PureTensor):
        return expr.tensor_shape
    if isinstance(expr, ZeroTensor):
        return expr.shape
    if isinstance(expr, Mul) and not isinstance(expr, PureTensor):
        for f in expr.args:
            if isinstance(f, PureTensor):
                return f.tensor_shape
            if isinstance(f, ZeroTensor):
                return f.shape
    if hasattr(expr, "shape"):
        return (expr.shape,)
    return None


# ---------------------------------------------------------------------------

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

    Tensor shapes are full sequences of per-factor shapes, e.g.
    ``((3, 4), (4, 5))`` for the product ``A_3x4 ⊗ C_4x5``.
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

        *shape* is a tuple of per-factor (rows, cols) pairs, e.g.
        ``((3, 4), (4, 5))``.
        """
        shape = None
        zero_term = None
        terms = []

        work = list(args)
        while work:
            o = work.pop()

            if isinstance(o, ZeroTensor):
                candidate = o.shape
                if shape is None:
                    shape = candidate
                elif candidate != shape:
                    raise ShapeMismatchError(
                        f"Cannot add tensors of different shapes: "
                        f"{shape} vs {candidate}"
                    )
                zero_term = o
                continue

            if isinstance(o, AlgebraicTensor):
                work.extend(o.args)
                continue

            if isinstance(o, PureTensor):
                candidate = o.tensor_shape
                if shape is None:
                    shape = candidate
                elif candidate != shape:
                    raise ShapeMismatchError(
                        f"Cannot add tensors of different shapes: "
                        f"{shape} vs {candidate}"
                    )
                terms.append(o)
                continue

            if isinstance(o, Mul) and not isinstance(o, PureTensor):
                candidate = _tensor_shape_of(o)
                if candidate is not None:
                    if shape is None:
                        shape = candidate
                    elif candidate != shape:
                        raise ShapeMismatchError(
                            f"Cannot add tensors of different shapes: "
                            f"{shape} vs {candidate}"
                        )
                terms.append(o)
                continue

            if o.is_Number:
                terms.append(o)
                continue

            # Catch-all: bare matrix-like objects or anything else
            candidate = _tensor_shape_of(o)
            if candidate is not None:
                if shape is None:
                    shape = candidate
                elif candidate != shape:
                    raise ShapeMismatchError(
                        f"Cannot add tensors of different shapes: "
                        f"{shape} vs {candidate}"
                    )
            terms.append(o)

        return terms, shape, zero_term

    @property
    def tensor_shape(self):
        """Shape shared by every term in this sum.

        Returns a tuple of per-factor (rows, cols) pairs, e.g.
        ``((3, 4), (4, 5))``.
        """
        for arg in self.args:
            ts = _tensor_shape_of(arg)
            if ts is not None:
                return ts
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

    def __mul__(self, other):
        """Multiply each term by a scalar/symbol, or tensor-product with PureTensor."""
        other = sympify(other)
        if isinstance(other, Number) or (hasattr(other, 'is_commutative') and
                other.is_commutative and not isinstance(other, (PureTensor, AlgebraicTensor))):
            return AlgebraicTensor(*(a * other for a in self.args))
        if isinstance(other, PureTensor):
            return AlgebraicTensor(*(PureTensor(other, a) for a in self.args))
        return Mul(self, other)

    def __rmul__(self, other):
        """Multiply each term by a scalar/symbol, or tensor-product with PureTensor."""
        other = sympify(other)
        if isinstance(other, Number) or (hasattr(other, 'is_commutative') and
                other.is_commutative and not isinstance(other, (PureTensor, AlgebraicTensor))):
            return AlgebraicTensor(*(a * other for a in self.args))
        if isinstance(other, PureTensor):
            return AlgebraicTensor(*(PureTensor(other, a) for a in self.args))
        return Mul(other, self)

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
        if not self.args:
            return ""
        parts = []
        first = True
        for a in self.args:
            s = str(a)
            if first:
                parts.append(s)
                first = False
            elif s.startswith("-"):
                parts.append(f"- {s[1:]}")
            else:
                parts.append(f"+ {s}")
        return " ".join(parts)

    def __repr__(self):
        return f"AlgebraicTensor({', '.join(repr(a) for a in self.args)})"


# Register AlgebraicTensor with the SymPy add dispatcher so that
# expressions like (A + B) where A/B involve PureTensor or AlgebraicTensor
# are routed through AlgebraicTensor.__new__ instead of Add.__new__.
add.register_handlerclass(
    (PureTensor, AlgebraicTensor), AlgebraicTensor
)
