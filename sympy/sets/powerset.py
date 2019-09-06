from __future__ import print_function, division

from sympy.core.decorators import _sympifyit
from sympy.core.evaluate import global_evaluate
from sympy.core.logic import fuzzy_bool
from sympy.core.singleton import S
from sympy.core.sympify import _sympify

from .sets import Set, tfn


class PowerSet(Set):
    r"""A symbolic object representing a power set.

    Parameters
    ==========

    arg : Set
        The set to take power of.

    evaluate : bool
        The flag to control evaluation.

        If the evaluation is disabled for finite sets, it can take
        advantage of using subset test as a membership test.

    Notes
    =====

    Power set `\mathcal{P}(S)` is defined as a set containing all the
    subsets of `S`.

    If the set `S` is a finite set, its power set would have
    `2^{\left| S \right|}` elements, where `\left| S \right|` denotes
    the cardinality of `S`.

    Examples
    ========

    >>> from sympy.sets.powerset import PowerSet
    >>> from sympy import S, FiniteSet

    A power set of a finite set:

    >>> PowerSet(FiniteSet(1, 2, 3))
    {EmptySet(), {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}, {1, 2, 3}}

    A power set of an empty set:

    >>> PowerSet(S.EmptySet)
    {EmptySet()}
    >>> PowerSet(PowerSet(S.EmptySet))
    {EmptySet(), {EmptySet()}}

    A power set of an infinite set:

    >>> PowerSet(S.Reals)
    PowerSet(Reals)

    An unevaluated power set:

    >>> PowerSet(FiniteSet(1, 2, 3), evaluate=False)
    PowerSet({1, 2, 3})

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Power_set

    .. [2] https://en.wikipedia.org/wiki/Axiom_of_power_set
    """
    def __new__(cls, arg, evaluate=global_evaluate[0]):
        arg = _sympify(arg)

        if not isinstance(arg, Set):
            raise ValueError('{} must be a set.'.format(arg))

        if evaluate:
            ret = arg._eval_powerset()

            if ret is not None and not isinstance(arg, PowerSet):
                return ret

        return super(PowerSet, cls).__new__(cls, arg)

    @property
    def arg(self):
        return self.args[0]

    @_sympifyit('other', NotImplemented)
    def _contains(self, other):
        if not isinstance(other, Set):
            return None

        arg = self.arg
        ret = fuzzy_bool(arg.is_superset(other))
        if ret is not None:
            return ret
        return None

    def __len__(self):
        return 2 ** len(self.arg)
