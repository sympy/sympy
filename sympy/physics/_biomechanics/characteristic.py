"""Implementations of characteristic curves for musculotendon models."""

from __future__ import annotations

from typing import TYPE_CHECKING

from sympy.core.expr import UnevaluatedExpr
from sympy.core.function import ArgumentIndexError, Function
from sympy.core.numbers import Float, Integer
from sympy.functions.elementary.exponential import exp

if TYPE_CHECKING:
    from typing import Any

    from sympy.core.expr import Expr
    from sympy.printing.printer import Printer


__all__ = ['TendonForceLengthDeGroote2016']


class CharacteristicCurveFunction(Function):
    """Base class for all musculotendon characteristic curve functions."""

    @classmethod
    def eval(cls):
        msg = (
            f'Cannot directly instantiate {cls.__name__!r}, instances of '
            f'characteristic curves must be of a concrete subclass.'

        )
        raise TypeError(msg)

    def _print_code(self, printer: Printer) -> str:
        """Print code for the function defining the curve using a printer.

        Explanation
        ===========

        The order of operations may need to be controlled as constant folding
        the numeric terms within the equations of a musculotendon
        characteristic curve can sometimes results in a numerically-unstable
        expression.

        Parameters
        ==========

        printer : Printer
            The printer to be used to print a string representation of the
            characteristic curve as valid code in the target language.

        """
        return printer.doprint(self.doit(deep=False, evaluate=False))

    _ccode = _print_code
    _cupycode = _print_code
    _cxxcode = _print_code
    _fcode = _print_code
    _jaxcode = _print_code
    _lambdacode = _print_code
    _mpmathcode = _print_code
    _octave = _print_code
    _pythoncode = _print_code
    _numpycode = _print_code
    _scipycode = _print_code


class TendonForceLengthDeGroote2016(CharacteristicCurveFunction):
    r"""Tendon force-length curve based on De Groote et al., 2016 [1].

    Explanation
    ===========

    Gives the normalized tendon force produced as a function of normalized
    tendon length.

    The function is defined by the equation:

    $fl^T = c_0 \exp{c_3 \left( \tilde{l}^T - c_1 \right)} - c_2$

    with constant values of $c_0 = 0.2$, $c_1 = 0.995$, $c_2 = 0.25$, and
    $c_3 = 33.93669377311689$.

    While it is possible to change the constant values, these were carefully
    selected in the original publication to give the characteristic curve
    specific and required properties. For example, the function produces no
    force when the tendon is in an unstrained state. It also produces a force
    of 1 normalized unit when the tendon is under a 5% strain.

    References
    ==========

    .. [1] De Groote, F., Kinney, A. L., Rao, A. V., & Fregly, B. J., Evaluation
           of direct collocation optimal control problem formulations for
           solving the muscle redundancy problem, Annals of biomedical
           engineering, 44(10), (2016) pp. 2922-2936

    """

    @classmethod
    def with_default_constants(cls, l_T_tilde: Any) -> TendonForceLengthDeGroote2016:
        r"""Alternative constructor that will use the recommended constants.

        Explanation
        ===========

        Returns a new instance of the tendon force-length function using the
        four constant values specified in the original publication.

        These have the values:

        $c_0 = 0.2$
        $c_1 = 0.995$
        $c_2 = 0.25$
        $c_3 = 33.93669377311689$

        Parameters
        ==========

        l_T_tilde : Any (sympifiable)
            Normalized tendon length.

        """
        c0 = Float('0.2')
        c1 = Float('0.995')
        c2 = Float('0.25')
        c3 = Float('33.93669377311689')
        return cls(l_T_tilde, c0, c1, c2, c3)

    @classmethod
    def eval(cls, l_T_tilde: Any, c0: Any, c1: Any, c2: Any, c3: Any) -> Any:  # type: ignore
        """Evaluation of basic inputs.

        Parameters
        ==========

        l_T_tilde : Any (sympifiable)
            Normalized tendon length.
        c0 : Any (sympifiable)
            The first constant in the characteristic equation. The published
            value is ``0.2``.
        c1 : Any (sympifiable)
            The second constant in the characteristic equation. The published
            value is ``0.995``.
        c2 : Any (sympifiable)
            The third constant in the characteristic equation. The published
            value is ``0.25``.
        c3 : Any (sympifiable)
            The fourth constant in the characteristic equation. The published
            value is ``33.93669377311689``.

        """
        pass

    def _eval_evalf(self, prec):
        """Evaluate the expression numerically using ``evalf``."""
        return self.doit(deep=False, evaluate=False)._eval_evalf(prec)

    def doit(
        self,
        deep: bool = True,
        evaluate : bool = True,
        **hints: Any,
    ) -> Expr:
        """Evaluate the expression defining the function.

        Parameters
        ==========

        deep : bool
            Whether ``doit`` should be recursively called. Default is ``True``.
        evaluate : bool.
            Whether the SymPy expression should be evaluated as it is
            constructed. If ``False``, then no constant folding will be
            conducted which will leave the expression in a more numerically-
            stable for values of ``l_T_tilde`` that correspond to a sensible
            operating range for a musculotendon. Default is ``True``.
        **kwargs : dict[str, Any]
            Additional keyword argument pairs to be recursively passed to
            ``doit``.

        """
        l_T_tilde, *constants = self.args
        if deep:
            hints['evaluate'] = evaluate
            l_T_tilde = l_T_tilde.doit(deep=deep, **hints)
            c0, c1, c2, c3 = [c.doit(deep=deep, **hints) for c in constants]
        else:
            c0, c1, c2, c3 = constants

        if evaluate:
            return c0 * exp(c3 * (l_T_tilde - c1)) - c2

        return c0 * exp(c3 * UnevaluatedExpr(l_T_tilde - c1)) - c2

    def fdiff(self, argindex: int = 1) -> Expr:
        """Derivative of the function with respect to a single argument.

        Parameters
        ==========

        argindex : int
            The index of the function's arguments with respect to which the
            derivative should be taken. Argument indexes start at ``1``.
            Default is ``1``.

        """
        l_T_tilde, c0, c1, c2, c3 = self.args
        if argindex == 1:
            return c0 * c3 * exp(c3 * UnevaluatedExpr(l_T_tilde - c1))  # type: ignore
        elif argindex == 2:
            return exp(c3 * UnevaluatedExpr(l_T_tilde - c1))  # type: ignore
        elif argindex == 3:
            return -c0 * c3 * exp(c3 * UnevaluatedExpr(l_T_tilde - c1))  # type: ignore
        elif argindex == 4:
            return Integer(-1)
        elif argindex == 5:
            return c0 * (l_T_tilde - c1) * exp(c3 * UnevaluatedExpr(l_T_tilde - c1))  # type: ignore

        raise ArgumentIndexError(self, argindex)

    def _latex(self, printer: Printer) -> str:
        """Print a LaTeX representation of the function defining the curve.

        Parameters
        ==========

        printer : Printer
            The printer to be used to print the LaTeX string representation.

        """
        l_T_tilde = self.args[0]
        _l_T_tilde = printer._print(l_T_tilde)
        return r'\operatorname{fl}^T \left( %s \right)' % _l_T_tilde
