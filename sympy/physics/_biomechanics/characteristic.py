"""Implementations of characteristic curves for musculotendon models."""

from __future__ import annotations

from typing import TYPE_CHECKING

from sympy.core.expr import UnevaluatedExpr
from sympy.core.function import ArgumentIndexError, Function
from sympy.core.numbers import Float, Integer
from sympy.functions.elementary.exponential import exp, log

if TYPE_CHECKING:
    from typing import Any

    from sympy.core.expr import Expr
    from sympy.printing.printer import Printer


__all__ = [
    'TendonForceLengthDeGroote2016',
    'TendonForceLengthInverseDeGroote2016',
]


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

    Examples
    ========

    The preferred way to instantiate ``TendonForceLengthDeGroote2016`` is using
    the ``with_default_constants`` constructor because this will automatically
    populate the constants within the characteristic curve equation with the
    floating point values from the original publication. This constructor takes
    a single argument corresponding to normalized tendon length. We'll create a
    ``Symbol`` called ``l_T_tilde`` to represent this.

    >>> from sympy import Symbol
    >>> from sympy.physics._biomechanics import TendonForceLengthDeGroote2016
    >>> l_T_tilde = Symbol('l_T_tilde')
    >>> fl_T = TendonForceLengthDeGroote2016.with_default_constants(l_T_tilde)
    >>> fl_T
    TendonForceLengthDeGroote2016(l_T_tilde, 0.2, 0.995, 0.25,
    33.93669377311689)

    It's also possible to populate the four constants with your own values too.

    >>> from sympy import symbols
    >>> c0, c1, c2, c3 = symbols('c0 c1 c2 c3')
    >>> fl_T = TendonForceLengthDeGroote2016(l_T_tilde, c0, c1, c2, c3)
    >>> fl_T
    TendonForceLengthDeGroote2016(l_T_tilde, c0, c1, c2, c3)

    You don't just have to use symbols as the arguments, it's also possible to
    use expressions. Let's create a new pair of symbols, ``l_T`` and
    ``l_T_slack``, representing tendon length and tendon slack length
    respectively. We can then represent ``l_T_tilde`` as an expression, the
    ratio of these.

    >>> l_T, l_T_slack = symbols('l_T l_T_slack')
    >>> l_T_tilde = l_T/l_T_slack
    >>> fl_T = TendonForceLengthDeGroote2016.with_default_constants(l_T_tilde)
    >>> fl_T
    TendonForceLengthDeGroote2016(l_T/l_T_slack, 0.2, 0.995, 0.25,
    33.93669377311689)

    To inspect the actual symbolic expression that this function represents,
    we can call the ``doit`` method on an instance. We'll use the keyword
    argument ``evaluate=False`` as this will keep the expression in its
    canonical form and won't simplify any constants.

    >>> fl_T.doit(evaluate=False)
    -0.25 + 0.2*exp(33.93669377311689*(l_T/l_T_slack - 0.995))

    The function can also be differentiated. We'll differentiate with respect
    to l_T using the ``diff`` method on an instance with the single positional
    argument ``l_T``.

    >>> fl_T.diff(l_T)
    6.787338754623378*exp(33.93669377311689*(l_T/l_T_slack - 0.995))/l_T_slack

    References
    ==========

    .. [1] De Groote, F., Kinney, A. L., Rao, A. V., & Fregly, B. J., Evaluation
           of direct collocation optimal control problem formulations for
           solving the muscle redundancy problem, Annals of biomedical
           engineering, 44(10), (2016) pp. 2922-2936

    """

    @classmethod
    def with_default_constants(cls, l_T_tilde: Any) -> TendonForceLengthDeGroote2016:
        r"""Recommended constructor that will use the published constants.

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
            return c0*exp(c3*(l_T_tilde - c1)) - c2

        return c0*exp(c3*UnevaluatedExpr(l_T_tilde - c1)) - c2

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
            return c0*c3*exp(c3*UnevaluatedExpr(l_T_tilde - c1))  # type: ignore
        elif argindex == 2:
            return exp(c3*UnevaluatedExpr(l_T_tilde - c1))  # type: ignore
        elif argindex == 3:
            return -c0*c3*exp(c3*UnevaluatedExpr(l_T_tilde - c1))  # type: ignore
        elif argindex == 4:
            return Integer(-1)
        elif argindex == 5:
            return c0*(l_T_tilde - c1)*exp(c3*UnevaluatedExpr(l_T_tilde - c1))  # type: ignore

        raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex: int = 1) -> Function:
        """Inverse function.

        Parameters
        ==========

        argindex : int
            Value to start indexing the arguments at. Default is ``1``.

        """
        return TendonForceLengthInverseDeGroote2016

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


class TendonForceLengthInverseDeGroote2016(CharacteristicCurveFunction):
        r"""Inverse tendon force-length curve based on De Groote et al., 2016 [1].

        Explanation
        ===========

        Gives the normalized tendon length that produces a specific normalized
        tendon force.

        The function is defined by the equation:

        ${fl^T}^{-1} = frac{\log{\frac{fl^T + c_2}{c_0}}}{c_3} + c_1$

        with constant values of $c_0 = 0.2$, $c_1 = 0.995$, $c_2 = 0.25$, and
        $c_3 = 33.93669377311689$. This function is the exact analytical inverse
        of the related tendon force-length curve ``TendonForceLengthDeGroote2016``.

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
        def with_default_constants(cls, fl_T: Any) -> TendonForceLengthInverseDeGroote2016:
            r"""Recommended constructor that will use the published constants.

            Explanation
            ===========

            Returns a new instance of the inverse tendon force-length function
            using the four constant values specified in the original publication.

            These have the values:

            $c_0 = 0.2$
            $c_1 = 0.995$
            $c_2 = 0.25$
            $c_3 = 33.93669377311689$

            Parameters
            ==========

            fl_T : Any (sympifiable)
                Normalized tendon force as a function of tendon length.

            """
            c0 = Float('0.2')
            c1 = Float('0.995')
            c2 = Float('0.25')
            c3 = Float('33.93669377311689')
            return cls(fl_T, c0, c1, c2, c3)

        @classmethod
        def eval(cls, fl_T: Any, c0: Any, c1: Any, c2: Any, c3: Any) -> Any:  # type: ignore
            """Evaluation of basic inputs.

            Parameters
            ==========

            fl_T : Any (sympifiable)
                Normalized tendon force as a function of tendon length.
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
            fl_T, *constants = self.args
            if deep:
                hints['evaluate'] = evaluate
                fl_T = fl_T.doit(deep=deep, **hints)
                c0, c1, c2, c3 = [c.doit(deep=deep, **hints) for c in constants]
            else:
                c0, c1, c2, c3 = constants

            if evaluate:
                return log((fl_T + c2)/c0)/c3 + c1

            return log(UnevaluatedExpr((fl_T + c2)/c0))/c3 + c1

        def fdiff(self, argindex: int = 1) -> Expr:
            """Derivative of the function with respect to a single argument.

            Parameters
            ==========

            argindex : int
                The index of the function's arguments with respect to which the
                derivative should be taken. Argument indexes start at ``1``.
                Default is ``1``.

            """
            fl_T, c0, c1, c2, c3 = self.args
            if argindex == 1:
                return 1/(c3*(fl_T + c2))  # type: ignore
            elif argindex == 2:
                return -1/(c0*c3)  # type: ignore
            elif argindex == 3:
                return Integer(1)
            elif argindex == 4:
                return 1/(c3*(fl_T + c2))  # type: ignore
            elif argindex == 5:
                return -log(UnevaluatedExpr((fl_T + c2)/c0))/c3**2  # type: ignore

            raise ArgumentIndexError(self, argindex)

        def inverse(self, argindex: int = 1) -> Function:
            """Inverse function.

            Parameters
            ==========

            argindex : int
                Value to start indexing the arguments at. Default is ``1``.

            """
            return TendonForceLengthDeGroote2016

        def _latex(self, printer: Printer) -> str:
            """Print a LaTeX representation of the function defining the curve.

            Parameters
            ==========

            printer : Printer
                The printer to be used to print the LaTeX string representation.

            """
            fl_T = self.args[0]
            _fl_T = printer._print(fl_T)
            return r'\left( \operatorname{fl}^T \right)^{-1} \left( %s \right)' % _fl_T


class FiberForceLengthPassiveDeGroote2016(CharacteristicCurveFunction):
    pass
