"""Implementations of characteristic curves for musculotendon models."""

from sympy.core.function import Function


__all__ = ['fl_T_de_groote_2016']


class CharacteristicCurveFunction(Function):
    """Base class for all musculotendon characteristic curve functions."""


class fl_T_de_groote_2016(CharacteristicCurveFunction):
    r"""Tendon force-length curve based on De Groote et al., 2016 [1].

    Explanation
    ===========

    Gives the normalized tendon force produced as a function of normalized
    tendon length.

    The function is defined by the equation:

    $fl^T = c_0 \exp{c_3 \tilde{l}^T - c_1} - c_2$

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
