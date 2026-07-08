from __future__ import annotations

from sympy.core.add import Add, add
from sympy.core.basic import Basic
from sympy.core.mul import Mul
from sympy.core.numbers import Number
from sympy.core.singleton import S
from sympy.core.sympify import sympify

from sympy.tensor.algebraic.algebraic_pure_tensor import (
    AlgebraicPureTensor, _factor_shapes, _factor_has_noncommutative
)
from sympy.tensor.algebraic.algebraic_zero_tensor import AlgebraicZeroTensor
from sympy.tensor.algebraic.scalar_mul import ScalarMul


# ---------------------------------------------------------------------------
# Shape helpers
# ---------------------------------------------------------------------------

def _tensor_shape_of(expr):
    """Return the full tensor shape of *expr* as a tuple of factor shapes.

    Handles AlgebraicPureTensor, AlgebraicZeroTensor, ScalarMul,
    Mul(coeff, AlgebraicPureTensor), and bare matrix-like objects
    (whose single-factor shape is wrapped).
    """
    if isinstance(expr, AlgebraicPureTensor):
        return expr.tensor_shape
    if isinstance(expr, AlgebraicZeroTensor):
        return expr.shape
    if isinstance(expr, ScalarMul):
        return expr.tensor_shape
    if isinstance(expr, Mul) and not isinstance(expr, AlgebraicPureTensor):
        for f in expr.args:
            if isinstance(f, AlgebraicPureTensor):
                return f.tensor_shape
            if isinstance(f, AlgebraicZeroTensor):
                return f.shape
    if hasattr(expr, "shape"):
        return (expr.shape,)
    return None


def _commutativity_shape_of(expr):
    """Return the commutativity_shape of *expr*, or None if not determinable.

    Handles AlgebraicPureTensor, AlgebraicZeroTensor, AlgebraicTensor,
    ScalarMul, Mul(coeff, AlgebraicPureTensor), and bare matrix-like objects.
    """
    if isinstance(expr, AlgebraicPureTensor):
        return expr.commutativity_shape
    if isinstance(expr, AlgebraicZeroTensor):
        return tuple(1 for _ in expr.shape)
    if isinstance(expr, AlgebraicTensor):
        return expr.commutativity_shape
    if isinstance(expr, ScalarMul):
        return expr.commutativity_shape
    if isinstance(expr, Mul) and not isinstance(expr, AlgebraicPureTensor):
        for f in expr.args:
            if isinstance(f, AlgebraicPureTensor):
                return f.commutativity_shape
            if isinstance(f, AlgebraicZeroTensor):
                return tuple(1 for _ in f.shape)
    if hasattr(expr, "shape"):
        if _factor_has_noncommutative(expr):
            return (0,)
        return (1,)
    return None


# ---------------------------------------------------------------------------

class ShapeMismatchError(TypeError):
    """Raised when attempting to add tensors of incompatible shapes."""
    pass


def _compose_reassemble(results, tensor_shape):
    """Reassemble a list of composition results into a single expression.

    Strips AlgebraicZeroTensor anchors, handles cancellation, and wraps
    back into an AlgebraicTensor if needed.
    """
    from sympy.core.singleton import S
    from sympy.tensor.algebraic.algebraic_zero_tensor import (
        AlgebraicZeroTensor as _ZT,
    )

    real = []
    zero_shape = None
    for r in results:
        if isinstance(r, _ZT):
            zero_shape = r.shape
        elif r is not S.Zero:
            real.append(r)

    if not real:
        if zero_shape:
            return _ZT(zero_shape)
        return _ZT(tensor_shape)

    if len(real) == 1:
        return real[0]

    return AlgebraicTensor(*real, _sympify=False)


def compose_algebraic_tensors(left, right):
    """Compose two algebraic-tensor expressions by linearity.

    Composes *left* and *right* by factor-wise matrix multiplication,
    extending ``compose_algebraic_pure_tensors`` by linearity to sums.

    Parameters
    ----------
    left : AlgebraicTensor, AlgebraicPureTensor, AlgebraicZeroTensor, or matrix
        The left operand.
    right : AlgebraicTensor, AlgebraicPureTensor, AlgebraicZeroTensor, or matrix
        The right operand.

    Returns
    -------
    AlgebraicTensor, AlgebraicPureTensor, AlgebraicZeroTensor, or matrix
        The composition result.

    Raises
    ------
    TypeError
        If either argument is not a recognised tensor or matrix type.
    ValueError
        If the factor structures are incompatible for composition.
    """
    from sympy.tensor.algebraic.algebraic_pure_tensor import (
        AlgebraicPureTensor as _PT,
        compose_algebraic_pure_tensors,
    )
    from sympy.tensor.algebraic.algebraic_zero_tensor import (
        AlgebraicZeroTensor as _ZT,
    )
    from sympy.tensor.algebraic.scalar_mul import ScalarMul as _SM

    # --- AlgebraicZeroTensor shortcuts ---
    if isinstance(left, _ZT):
        return _ZT(left.shape)
    if isinstance(right, _ZT):
        return _ZT(right.shape)

    # --- ScalarMul shortcuts: compose inner tensor, preserve scalar ---
    if isinstance(left, _SM):
        comp = compose_algebraic_tensors(left.tensor, right)
        if isinstance(comp, _ZT):
            return comp
        return _SM(left.scalar, comp)
    if isinstance(right, _SM):
        comp = compose_algebraic_tensors(left, right.tensor)
        if isinstance(comp, _ZT):
            return comp
        return _SM(right.scalar, comp)

    # --- AlgebraicTensor × AlgebraicTensor (check first before single-side) ---
    if isinstance(left, AlgebraicTensor) and isinstance(right, AlgebraicTensor):
        results = []
        for la in left.args:
            if isinstance(la, _ZT):
                results.append(la)
                continue
            for ra in right.args:
                if isinstance(ra, _ZT):
                    results.append(ra)
                else:
                    comp = compose_algebraic_tensors(la, ra)
                    results.append(comp)
        return _compose_reassemble(results, left.tensor_shape)

    # --- PureTensor / bare matrix × AlgebraicTensor ---
    if isinstance(right, AlgebraicTensor):
        if isinstance(left, (_PT, _SM)) or hasattr(left, "shape"):
            results = []
            for a in right.args:
                if isinstance(a, _ZT):
                    results.append(a)
                else:
                    results.append(
                        compose_algebraic_tensors(left, a)
                    )
            return _compose_reassemble(results, right.tensor_shape)
        raise TypeError(
            f"Expected AlgebraicTensor, AlgebraicPureTensor, or matrix on "
            f"the left, got {type(left).__name__}"
        )

    # --- AlgebraicTensor × PureTensor / bare matrix ---
    if isinstance(left, AlgebraicTensor):
        if isinstance(right, (_PT, _SM)) or hasattr(right, "shape"):
            return left._compose_with_term(right)
        raise TypeError(
            f"Expected AlgebraicTensor, AlgebraicPureTensor, or matrix on "
            f"the right, got {type(right).__name__}"
        )

    # --- PureTensor / ScalarMul / bare matrix × PureTensor / ScalarMul / bare matrix ---
    if isinstance(left, _PT) and isinstance(right, _PT):
        return compose_algebraic_pure_tensors(left, right)
    if isinstance(left, _PT) and isinstance(right, _SM):
        return _SM(right.scalar, compose_algebraic_pure_tensors(left, right.tensor))
    if isinstance(left, _SM) and isinstance(right, _PT):
        return _SM(left.scalar, compose_algebraic_pure_tensors(left.tensor, right))
    if isinstance(left, _SM) and isinstance(right, _SM):
        c = compose_algebraic_pure_tensors(left.tensor, right.tensor)
        return _SM(left.scalar * right.scalar, c)
    if isinstance(left, _PT) and hasattr(right, "shape"):
        return compose_algebraic_pure_tensors(left, right)
    if hasattr(left, "shape") and isinstance(right, _PT):
        return compose_algebraic_pure_tensors(left, right)
    if hasattr(left, "shape") and hasattr(right, "shape"):
        return compose_algebraic_pure_tensors(left, right)

    raise TypeError(
        f"Cannot compose {type(left).__name__} with {type(right).__name__}"
    )


class AlgebraicTensor(Basic):
    """Sum of AlgebraicPureTensors (and/or AlgebraicZeroTensors) sharing the same tensor shape.

    Built on top of ``Basic`` (not ``Add``) to avoid the is_commutative
    descriptor conflict in ``AssocOp._from_args``, while still exposing
    ``is_Add = True`` so that the wider SymPy ecosystem recognises this
    as an additive expression.

    Shape enforcement and AlgebraicZeroTensor handling are done in ``__new__``.

    Tensor shapes are full sequences of per-factor shapes, e.g.
    ``((3, 4), (4, 5))`` for the product ``A_3x4 ⊗ C_4x5``.
    """

    __slots__ = ()

    is_AlgebraicTensor = True
    is_Add = True

    identity = None  # no single identity; use AlgebraicZeroTensor(shape) instead

    def __new__(cls, *args, evaluate=True, _sympify=True):
        if not args:
            raise ValueError("AlgebraicTensor requires at least one argument")

        if _sympify:
            args = tuple(
                a if isinstance(a, AlgebraicZeroTensor) else sympify(a)
                for a in args
            )

        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, (AlgebraicPureTensor, AlgebraicZeroTensor)):
                return arg
            if isinstance(arg, AlgebraicTensor):
                return arg
            if isinstance(arg, (list, tuple)) and len(arg) == 1:
                return arg[0]

        # Flatten nested AlgebraicTensors and collect all leaf terms
        flat, shape, zero_term, comm_cs = cls._flatten_args(args)

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

        # If an AlgebraicZeroTensor was provided, always keep it as an anchor.
        # Only unwrap to a single term when there is exactly one non-zero
        # term AND no AlgebraicZeroTensor was explicitly given.
        if len(flat) == 1 and zero_term is None and not isinstance(flat[0], AlgebraicZeroTensor):
            return flat[0]

        # Append AlgebraicZeroTensor if present (for shape anchoring)
        if zero_term is not None:
            flat.append(zero_term)

        obj = Basic.__new__(cls, *flat)
        return obj

    @classmethod
    def _flatten_args(cls, args):
        """Flatten nested AlgebraicTensors, validate shapes, collect AlgebraicZeroTensors.

        Returns
        -------
        (terms : list, shape : tuple | None, zero_term : AlgebraicZeroTensor | None,
         commutativity_shape : tuple | None)

        *shape* is a tuple of per-factor (rows, cols) pairs, e.g.
        ``((3, 4), (4, 5))``.
        *commutativity_shape* is a tuple of binary entries, same length as shape,
        representing the component-wise AND of all term commutativity_shapes.
        """
        shape = None
        zero_term = None
        comm_cs = None
        terms = []

        work = list(args)
        while work:
            o = work.pop()

            if isinstance(o, AlgebraicZeroTensor):
                candidate = o.shape
                if shape is None:
                    shape = candidate
                    comm_cs = tuple(1 for _ in shape)
                elif candidate != shape:
                    raise ShapeMismatchError(
                        f"Cannot add tensors of different shapes: "
                        f"{shape} vs {candidate}"
                    )
                zero_term = o
                continue

            if isinstance(o, AlgebraicTensor):
                if comm_cs is not None:
                    comm_cs = tuple(r & c for r, c in zip(comm_cs, o.commutativity_shape))
                work.extend(o.args)
                continue

            if isinstance(o, AlgebraicPureTensor):
                candidate = o.tensor_shape
                if shape is None:
                    shape = candidate
                    comm_cs = tuple(1 for _ in shape)
                elif candidate != shape:
                    raise ShapeMismatchError(
                        f"Cannot add tensors of different shapes: "
                        f"{shape} vs {candidate}"
                    )
                if comm_cs is not None:
                    comm_cs = tuple(r & c for r, c in zip(comm_cs, o.commutativity_shape))
                terms.append(o)
                continue

            if isinstance(o, ScalarMul):
                candidate = o.tensor_shape
                if shape is None:
                    shape = candidate
                    comm_cs = tuple(1 for _ in shape)
                elif candidate != shape:
                    raise ShapeMismatchError(
                        f"Cannot add tensors of different shapes: "
                        f"{shape} vs {candidate}"
                    )
                if comm_cs is not None:
                    comm_cs = tuple(r & c for r, c in zip(comm_cs, o.commutativity_shape))
                terms.append(o)
                continue

            if isinstance(o, Mul) and not isinstance(o, AlgebraicPureTensor):
                candidate = _tensor_shape_of(o)
                if candidate is not None:
                    if shape is None:
                        shape = candidate
                        comm_cs = tuple(1 for _ in shape)
                    elif candidate != shape:
                        raise ShapeMismatchError(
                            f"Cannot add tensors of different shapes: "
                            f"{shape} vs {candidate}"
                        )
                cs = _commutativity_shape_of(o)
                if cs is not None and comm_cs is not None:
                    comm_cs = tuple(r & c for r, c in zip(comm_cs, cs))
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
                    comm_cs = tuple(1 for _ in shape)
                elif candidate != shape:
                    raise ShapeMismatchError(
                        f"Cannot add tensors of different shapes: "
                        f"{shape} vs {candidate}"
                    )
            cs = _commutativity_shape_of(o)
            if cs is not None and comm_cs is not None:
                comm_cs = tuple(r & c for r, c in zip(comm_cs, cs))
            terms.append(o)

        return terms, shape, zero_term, comm_cs

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
    def commutativity_shape(self):
        """Component-wise AND of commutativity_shape over all terms in this sum.

        For AlgebraicTensor operands, uses their stored commutativity_shape
        rather than recomputing from individual PureTensor factors.
        """
        ts = self.tensor_shape
        if ts is None:
            raise AttributeError("Cannot determine commutativity_shape without tensor_shape")
        result = tuple(1 for _ in ts)
        for arg in self.args:
            cs = _commutativity_shape_of(arg)
            if cs is not None:
                result = tuple(r & c for r, c in zip(result, cs))
        return result

    @property
    def terms(self):
        """Non-zero, non-coefficient terms in this sum."""
        return tuple(
            a for a in self.args
            if not isinstance(a, Number) and not isinstance(a, AlgebraicZeroTensor)
        )

    def has_zero_term(self):
        """Return True if an AlgebraicZeroTensor anchors this sum."""
        return any(isinstance(a, AlgebraicZeroTensor) for a in self.args)

    # ---- Arithmetic delegation ----

    def __neg__(self):
        """Negate each term in the sum."""
        return AlgebraicTensor(*(-a for a in self.args))

    def _compose_with_term(self, other):
        """Compose this AlgebraicTensor with a single term (AlgebraicPureTensor, ScalarMul or bare matrix).

        Composes each term of self with *other* by factor-wise matrix
        multiplication.

        Parameters
        ----------
        other : AlgebraicPureTensor, ScalarMul or bare matrix-like object
            The right operand.

        Returns
        -------
        AlgebraicTensor, AlgebraicPureTensor, or AlgebraicZeroTensor
            The composition result.
        """
        from sympy.core.singleton import S
        from sympy.tensor.algebraic.algebraic_pure_tensor import (
            AlgebraicPureTensor as _PT,
            compose_algebraic_pure_tensors,
        )
        from sympy.tensor.algebraic.algebraic_zero_tensor import (
            AlgebraicZeroTensor as _ZT,
        )
        from sympy.tensor.algebraic.scalar_mul import ScalarMul as _SM

        # If other is ScalarMul, compose with inner tensor, then apply scalar
        if isinstance(other, _SM):
            comp = self._compose_with_term(other.tensor)
            return other.scalar * comp

        results = []
        for a in self.args:
            if isinstance(a, _ZT):
                results.append(a)
                continue

            if isinstance(a, _PT):
                comp = compose_algebraic_pure_tensors(a, other)
                results.append(comp)
            elif isinstance(a, Mul) and not isinstance(a, _PT):
                # Mul(coeff, AlgebraicPureTensor) — compose the PureTensor part
                for f in a.args:
                    if isinstance(f, _PT):
                        coeff = Mul(*[x for x in a.args if x is not f])
                        comp = compose_algebraic_pure_tensors(f, other)
                        results.append(ScalarMul(coeff, comp))
                        break
                else:
                    results.append(a)
            elif isinstance(a, ScalarMul):
                # ScalarMul(scalar, AlgebraicPureTensor) — compose the tensor part
                comp = compose_algebraic_pure_tensors(a.tensor, other)
                results.append(ScalarMul(a.scalar, comp))
            elif hasattr(a, "shape"):
                comp = compose_algebraic_pure_tensors(a, other)
                results.append(comp)
            else:
                results.append(a)

        return _compose_reassemble(results, self.tensor_shape)

    def __add__(self, other):
        return AlgebraicTensor(self, other)

    def __radd__(self, other):
        return AlgebraicTensor(other, self)

    def __sub__(self, other):
        return AlgebraicTensor(self, -other)

    def __rsub__(self, other):
        return AlgebraicTensor(other, -self)

    def __mul__(self, other):
        """Compose or scale this AlgebraicTensor.

        For commutative scalars/symbols the scalar is absorbed as a
        coefficient into each term.  For AlgebraicPureTensor,
        AlgebraicTensor, or bare matrices the result is the tensor
        composition (factor-wise matrix multiplication).

        Parameters
        ----------
        other : scalar, AlgebraicPureTensor, AlgebraicTensor, or matrix
            The right operand.

        Returns
        -------
        AlgebraicTensor, AlgebraicPureTensor, AlgebraicZeroTensor, or Mul
            The scaled/composed result.
        """
        other = sympify(other)
        if isinstance(other, Number) or (hasattr(other, 'is_commutative') and
                other.is_commutative and not isinstance(other, (AlgebraicPureTensor, AlgebraicTensor))):
            return AlgebraicTensor(*(a * other for a in self.args))
        return compose_algebraic_tensors(self, other)

    def __rmul__(self, other):
        """Compose or scale this AlgebraicTensor from the left.

        For commutative scalars/symbols the scalar is absorbed as a
        coefficient into each term.  For AlgebraicPureTensor,
        AlgebraicTensor, or bare matrices the result is the tensor
        composition (factor-wise matrix multiplication).

        Parameters
        ----------
        other : scalar, AlgebraicPureTensor, AlgebraicTensor, or matrix
            The left operand.

        Returns
        -------
        AlgebraicTensor, AlgebraicPureTensor, AlgebraicZeroTensor, or Mul
            The scaled/composed result.
        """
        other = sympify(other)
        if isinstance(other, Number) or (hasattr(other, 'is_commutative') and
                other.is_commutative and not isinstance(other, (AlgebraicPureTensor, AlgebraicTensor))):
            return AlgebraicTensor(*(other * a for a in self.args))
        return compose_algebraic_tensors(other, self)

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

    def simplify(self):
        """Simplify this AlgebraicTensor using the tensor simplification pipeline.

        Applies proportionality factoring and per-term simplification.
        """
        from sympy.tensor.algebraic.simplify import _simplify_algebraic_tensor
        return _simplify_algebraic_tensor(self)

    def display(self, mode="latex"):
        """Display this tensor using IPython display or fallback to print.

        Parameters
        ----------
        mode : str, default 'latex'
            'latex' for LaTeX rendering, 'text' for plain text.
        """
        try:
            from IPython.display import display, Latex
            if mode == "latex":
                display(Latex(self._repr_latex_()))
            else:
                display(self, plain=True)
        except ImportError:
            if mode == "latex":
                print(self._repr_latex_())
            else:
                print(self)

    def expand(self, deep=True, **hints):
        """Expand this AlgebraicTensor by expanding each term.

        Each term is expanded individually, and the results are reassembled
        into a single AlgebraicTensor (flattening any nested sums).
        """
        expanded_args = []
        for a in self.args:
            if isinstance(a, AlgebraicZeroTensor):
                expanded_args.append(a)
            elif hasattr(a, 'expand'):
                expanded_args.append(a.expand(deep=deep, **hints))
            else:
                expanded_args.append(a)

        if expanded_args == list(self.args):
            return self

        return AlgebraicTensor(*expanded_args)


# Register AlgebraicTensor with the SymPy add dispatcher so that
# expressions like (A + B) where A/B involve PureTensor or AlgebraicTensor
# are routed through AlgebraicTensor.__new__ instead of Add.__new__.
add.register_handlerclass(
    (AlgebraicPureTensor, AlgebraicTensor), AlgebraicTensor
)
add.register_handlerclass(
    (AlgebraicPureTensor, ScalarMul), AlgebraicTensor
)
add.register_handlerclass(
    (AlgebraicTensor, ScalarMul), AlgebraicTensor
)
