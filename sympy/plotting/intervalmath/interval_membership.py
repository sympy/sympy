from sympy.core.logic import fuzzy_and, fuzzy_or, fuzzy_not, fuzzy_xor
from collections.abc import Iterator
from typing import Any, Literal


class intervalMembership:
    """Represents a boolean expression returned by the comparison of
    the interval object.

    Parameters
    ==========

    (a, b) : (bool, bool)
        The first value determines the comparison as follows:
        - True: If the comparison is True throughout the intervals.
        - False: If the comparison is False throughout the intervals.
        - None: If the comparison is True for some part of the intervals.

        The second value is determined as follows:
        - True: If both the intervals in comparison are valid.
        - False: If at least one of the intervals is False, else
        - None
    """
    def __init__(self, a, b) -> None:
        self._wrapped = (a, b)

    def __getitem__(self, i):
        try:
            return self._wrapped[i]
        except IndexError:
            raise IndexError(
                "{} must be a valid indexing for the 2-tuple."
                .format(i))

    def __len__(self) -> Literal[2]:
        return 2

    def __iter__(self) -> Iterator[Any]:
        return iter(self._wrapped)

    def __str__(self):
        return "intervalMembership({}, {})".format(*self)
    __repr__ = __str__

    def __and__(self, other) -> "intervalMembership":
        if not isinstance(other, intervalMembership):
            raise ValueError(
                "The comparison is not supported for {}.".format(other))

        a1, b1 = self
        a2, b2 = other
        return intervalMembership(fuzzy_and([a1, a2]), fuzzy_and([b1, b2]))

    def __or__(self, other) -> "intervalMembership":
        if not isinstance(other, intervalMembership):
            raise ValueError(
                "The comparison is not supported for {}.".format(other))

        a1, b1 = self
        a2, b2 = other
        return intervalMembership(fuzzy_or([a1, a2]), fuzzy_and([b1, b2]))

    def __invert__(self) -> "intervalMembership":
        a, b = self
        return intervalMembership(fuzzy_not(a), b)

    def __xor__(self, other) -> "intervalMembership":
        if not isinstance(other, intervalMembership):
            raise ValueError(
                "The comparison is not supported for {}.".format(other))

        a1, b1 = self
        a2, b2 = other
        return intervalMembership(fuzzy_xor([a1, a2]), fuzzy_and([b1, b2]))

    def __eq__(self, other) -> bool:
        return self._wrapped == other

    def __ne__(self, other) -> bool:
        return self._wrapped != other
