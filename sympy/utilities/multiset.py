# -*- coding: utf-8 -*-
"""An implementation of a multiset."""

from collections import defaultdict
try:
    from collections.abc import Iterable, Mapping, MutableMapping, Set, Sized, Container
except ImportError:
    from collections import Iterable, Mapping, MutableMapping, Set, Sized, Container
from itertools import chain, repeat, starmap


class BaseMultiset:
    """A multiset implementation.

    A multiset is similar to the builtin :class:`set`, but elements can occur multiple times in the multiset.
    It is also similar to a :class:`list` without ordering of the values and hence no index-based operations.

    The multiset internally uses a :class:`dict` for storage where the key is the element and the value its
    multiplicity. It supports all operations that the :class:`set` supports.

    In contrast to the builtin :class:`collections.Counter`, no negative counts are allowed, elements with
    zero counts are removed from the :class:`dict`, and set operations are supported.

    The multiset comes in two variants, `Multiset` and `FrozenMultiset` which correspond to the `set` and
    `frozenset` classes, respectively.

    .. warning::

        You cannot instantiate this class directly. Use one of its variants instead.

    :see: https://en.wikipedia.org/wiki/Multiset
    """

    __slots__ = ('_elements', '_total')

    def __init__(self, iterable=None):
        r"""Create a new, empty Multiset object.

        And if given, initialize with elements from input iterable.
        Or, initialize from a mapping of elements to their multiplicity.

        Example:

        >>> ms = Multiset()                 # a new, empty multiset
        >>> ms = Multiset('abc')            # a new multiset from an iterable
        >>> ms = Multiset({'a': 4, 'b': 2}) # a new multiset from a mapping

        Args:
            iterable:
                An optional iterable of elements or mapping of elements to multiplicity to
                initialize the multiset from.
        """
        if isinstance(iterable, BaseMultiset):
            self._elements = iterable._elements.copy()
            self._total = iterable._total
        else:
            self._elements = _elements = defaultdict(int)
            _total = 0
            if iterable is not None:
                if isinstance(iterable, Mapping):
                    for element, multiplicity in iterable.items():
                        if multiplicity > 0:
                            _elements[element] = multiplicity
                            _total += multiplicity
                elif isinstance(iterable, Sized):
                    for element in iterable:
                        _elements[element] += 1
                    _total = len(iterable)
                else:
                    for element in iterable:
                        _elements[element] += 1
                        _total += 1
            self._total = _total

    def __new__(cls, iterable=None):
        if cls is BaseMultiset:
            raise TypeError("Cannot instantiate BaseMultiset directly, use either Multiset or FrozenMultiset.")
        return super(BaseMultiset, cls).__new__(cls)

    def __contains__(self, element):
        return element in self._elements

    def __getitem__(self, element):
        """The multiplicity of an element or zero if it is not in the multiset."""
        return self._elements.get(element, 0)

    def __str__(self):
        return '{%s}' % ', '.join(map(str, self))

    def __repr__(self):
        items = ', '.join('%r: %r' % item for item in self._elements.items())
        return '%s({%s})' % (type(self).__name__, items)

    def __len__(self):
        """Returns the total number of elements in the multiset.

        Note that this is equivalent to the sum of the multiplicities:

        >>> ms = Multiset('aab')
        >>> len(ms)
        3
        >>> sum(ms.multiplicities())
        3

        If you need the total number of distinct elements, use either the :meth:`distinct_elements` method:
        >>> len(ms.distinct_elements())
        2

        or convert to a :class:`set`:
        >>> len(set(ms))
        2
        """
        return self._total

    def __bool__(self):
        return self._total > 0

    def __iter__(self):
        return chain.from_iterable(starmap(repeat, self._elements.items()))

    def isdisjoint(self, other):
        r"""Return True if the set has no elements in common with other.

        Sets are disjoint iff their intersection is the empty set.

        >>> ms = Multiset('aab')
        >>> ms.isdisjoint('bc')
        False
        >>> ms.isdisjoint(Multiset('ccd'))
        True

        Args:
            other: The other set to check disjointedness. Can also be an :class:`~typing.Iterable`\[~T]
                or :class:`~typing.Mapping`\[~T, :class:`int`] which are then converted to :class:`Multiset`\[~T].
        """
        if not isinstance(other, Container):
            other = self._as_multiset(other)
        return all(element not in other for element in self._elements.keys())

    def difference(self, *others):
        r"""Return a new multiset with all elements from the others removed.

        >>> ms = Multiset('aab')
        >>> sorted(ms.difference('bc'))
        ['a', 'a']

        You can also use the ``-`` operator for the same effect. However, the operator version
        will only accept a set as other operator, not any iterable, to avoid errors.

        >>> ms = Multiset('aabbbc')
        >>> sorted(ms - Multiset('abd'))
        ['a', 'b', 'b', 'c']

        For a variant of the operation which modifies the multiset in place see
        :meth:`difference_update`.

        Args:
            others: The other sets to remove from the multiset. Can also be any :class:`~typing.Iterable`\[~T]
                or :class:`~typing.Mapping`\[~T, :class:`int`] which are then converted to :class:`Multiset`\[~T].

        Returns:
            The resulting difference multiset.
        """
        result = self.__copy__()
        _elements = result._elements
        _total = result._total
        for other in map(self._as_multiset, others):
            for element, multiplicity in other.items():
                if element in _elements:
                    old_multiplicity = _elements[element]
                    new_multiplicity = old_multiplicity - multiplicity
                    if new_multiplicity > 0:
                        _elements[element] = new_multiplicity
                        _total -= multiplicity
                    else:
                        del _elements[element]
                        _total -= old_multiplicity
        result._total = _total
        return result

    def __sub__(self, other):
        if not isinstance(other, (Set, BaseMultiset)):
            return NotImplemented
        return self.difference(other)

    def union(self, *others):
        r"""Return a new multiset with all elements from the multiset and the others with maximal multiplicities.

        >>> ms = Multiset('aab')
        >>> sorted(ms.union('bc'))
        ['a', 'a', 'b', 'c']

        You can also use the ``|`` operator for the same effect. However, the operator version
        will only accept a set as other operator, not any iterable, to avoid errors.

        >>> ms = Multiset('aab')
        >>> sorted(ms | Multiset('aaa'))
        ['a', 'a', 'a', 'b']

        For a variant of the operation which modifies the multiset in place see
        :meth:`union`.

        Args:
            *others: The other sets to union the multiset with. Can also be any :class:`~typing.Iterable`\[~T]
                or :class:`~typing.Mapping`\[~T, :class:`int`] which are then converted to :class:`Multiset`\[~T].

        Returns:
            The multiset resulting from the union.
        """
        result = self.__copy__()
        _elements = result._elements
        _total = result._total
        for other in map(self._as_mapping, others):
            for element, multiplicity in other.items():
                old_multiplicity = _elements.get(element, 0)
                if multiplicity > old_multiplicity:
                    _elements[element] = multiplicity
                    _total += multiplicity - old_multiplicity
        result._total = _total
        return result

    def __or__(self, other):
        if not isinstance(other, (Set, BaseMultiset)):
            return NotImplemented
        return self.union(other)

    __ror__ = __or__

    def combine(self, *others):
        r"""Return a new multiset with all elements from the multiset and the others with their multiplicities summed up.

        >>> ms = Multiset('aab')
        >>> sorted(ms.combine('bc'))
        ['a', 'a', 'b', 'b', 'c']

        You can also use the ``+`` operator for the same effect. However, the operator version
        will only accept a set as other operator, not any iterable, to avoid errors.

        >>> ms = Multiset('aab')
        >>> sorted(ms + Multiset('a'))
        ['a', 'a', 'a', 'b']

        For a variant of the operation which modifies the multiset in place see
        :meth:`update`.

        Args:
            others: The other sets to add to the multiset. Can also be any :class:`~typing.Iterable`\[~T]
                or :class:`~typing.Mapping`\[~T, :class:`int`] which are then converted to :class:`Multiset`\[~T].

        Returns:
            The multiset resulting from the addition of the sets.
        """
        result = self.__copy__()
        _elements = result._elements
        _total = result._total
        for other in map(self._as_mapping, others):
            for element, multiplicity in other.items():
                old_multiplicity = _elements.get(element, 0)
                new_multiplicity = old_multiplicity + multiplicity
                if old_multiplicity > 0 and new_multiplicity <= 0:
                    del _elements[element]
                    _total -= old_multiplicity
                elif new_multiplicity > 0:
                    _elements[element] = new_multiplicity
                    _total += multiplicity
        result._total = _total
        return result

    def __add__(self, other):
        if not isinstance(other, (Set, BaseMultiset)):
            return NotImplemented
        return self.combine(other)

    __radd__ = __add__

    def intersection(self, *others):
        r"""Return a new multiset with elements common to the multiset and all others.

        >>> ms = Multiset('aab')
        >>> sorted(ms.intersection('abc'))
        ['a', 'b']

        You can also use the ``&`` operator for the same effect. However, the operator version
        will only accept a set as other operator, not any iterable, to avoid errors.

        >>> ms = Multiset('aab')
        >>> sorted(ms & Multiset('aaac'))
        ['a', 'a']

        For a variant of the operation which modifies the multiset in place see
        :meth:`intersection_update`.

        Args:
            others: The other sets intersect with the multiset. Can also be any :class:`~typing.Iterable`\[~T]
                or :class:`~typing.Mapping`\[~T, :class:`int`] which are then converted to :class:`Multiset`\[~T].

        Returns:
            The multiset resulting from the intersection of the sets.
        """
        result = self.__copy__()
        _elements = result._elements
        _total = result._total
        for other in map(self._as_mapping, others):
            for element, multiplicity in list(_elements.items()):
                new_multiplicity = other.get(element, 0)
                if new_multiplicity < multiplicity:
                    if new_multiplicity > 0:
                        _elements[element] = new_multiplicity
                        _total -= multiplicity - new_multiplicity
                    else:
                        del _elements[element]
                        _total -= multiplicity
        result._total = _total
        return result

    def __and__(self, other):
        if not isinstance(other, (Set, BaseMultiset)):
            return NotImplemented
        return self.intersection(other)

    __rand__ = __and__

    def symmetric_difference(self, other):
        r"""Return a new set with elements in either the set or other but not both.

        >>> ms = Multiset('aab')
        >>> sorted(ms.symmetric_difference('abc'))
        ['a', 'c']

        You can also use the ``^`` operator for the same effect. However, the operator version
        will only accept a set as other operator, not any iterable, to avoid errors.

        >>> ms = Multiset('aab')
        >>> sorted(ms ^ Multiset('aaac'))
        ['a', 'b', 'c']

        For a variant of the operation which modifies the multiset in place see
        :meth:`symmetric_difference_update`.

        Args:
            other: The other set to take the symmetric difference with. Can also be any :class:`~typing.Iterable`\[~T]
                or :class:`~typing.Mapping`\[~T, :class:`int`] which are then converted to :class:`Multiset`\[~T].

        Returns:
            The resulting symmetric difference multiset.
        """
        other = self._as_multiset(other)
        result = self.__class__()
        _total = 0
        _elements = result._elements
        self_elements = self._elements
        other_elements = other._elements
        dist_elements = set(self_elements.keys()) | set(other_elements.keys())
        for element in dist_elements:
            multiplicity = self_elements.get(element, 0)
            other_multiplicity = other_elements.get(element, 0)
            new_multiplicity = (multiplicity - other_multiplicity
                                if multiplicity > other_multiplicity else other_multiplicity - multiplicity)
            _total += new_multiplicity
            if new_multiplicity > 0:
                _elements[element] = new_multiplicity
        result._total = _total
        return result

    def __xor__(self, other):
        if not isinstance(other, (Set, BaseMultiset)):
            return NotImplemented
        return self.symmetric_difference(other)

    __rxor__ = __xor__

    def times(self, factor):
        """Return a new set with each element's multiplicity multiplied with the given scalar factor.

        >>> ms = Multiset('aab')
        >>> sorted(ms.times(2))
        ['a', 'a', 'a', 'a', 'b', 'b']

        You can also use the ``*`` operator for the same effect:

        >>> sorted(ms * 3)
        ['a', 'a', 'a', 'a', 'a', 'a', 'b', 'b', 'b']

        For a variant of the operation which modifies the multiset in place see
        :meth:`times_update`.

        Args:
            factor: The factor to multiply each multiplicity with.
        """
        result = self.__copy__()
        _elements = result._elements
        for element in _elements:
            _elements[element] *= factor
        result._total *= factor
        return result

    def __mul__(self, factor):
        if not isinstance(factor, int):
            return NotImplemented
        return self.times(factor)

    __rmul__ = __mul__

    def _issubset(self, other, strict):
        other = self._as_multiset(other)
        self_len = self._total
        other_len = len(other)
        if self_len > other_len:
            return False
        if self_len == other_len and strict:
            return False
        return all(multiplicity <= other[element] for element, multiplicity in self.items())

    def issubset(self, other):
        """Return True iff this set is a subset of the other.

        >>> Multiset('ab').issubset('aabc')
        True
        >>> Multiset('aabb').issubset(Multiset('aabc'))
        False

        You can also use the ``<=`` operator for this comparison:

        >>> Multiset('ab') <= Multiset('ab')
        True

        When using the ``<`` operator for comparison, the sets are checked
        to be unequal in addition:

        >>> Multiset('ab') < Multiset('ab')
        False

        Args:
            other: The potential superset of the multiset to be checked.

        Returns:
            True iff this set is a subset of the other.
        """
        return self._issubset(other, False)

    def __le__(self, other):
        if not isinstance(other, (Set, BaseMultiset)):
            return NotImplemented
        return self._issubset(other, False)

    def __lt__(self, other):
        if not isinstance(other, (Set, BaseMultiset)):
            return NotImplemented
        return self._issubset(other, True)

    def _issuperset(self, other, strict):
        other = self._as_multiset(other)
        other_len = len(other)
        if len(self) < other_len:
            return False
        if len(self) == other_len and strict:
            return False
        for element, multiplicity in other.items():
            if self[element] < multiplicity:
                return False
        return True

    def issuperset(self, other):
        """Return True iff this multiset is a superset of the other.

        >>> Multiset('aabc').issuperset('ab')
        True
        >>> Multiset('aabc').issuperset(Multiset('abcc'))
        False

        You can also use the ``>=`` operator for this comparison:

        >>> Multiset('ab') >= Multiset('ab')
        True

        When using the ``>`` operator for comparison, the sets are checked
        to be unequal in addition:

        >>> Multiset('ab') > Multiset('ab')
        False

        Args:
            other: The potential subset of the multiset to be checked.

        Returns:
            True iff this set is a subset of the other.
        """
        return self._issuperset(other, False)

    def __ge__(self, other):
        if not isinstance(other, (Set, BaseMultiset)):
            return NotImplemented
        return self._issuperset(other, False)

    def __gt__(self, other):
        if not isinstance(other, (Set, BaseMultiset)):
            return NotImplemented
        return self._issuperset(other, True)

    def __eq__(self, other):
        if not isinstance(other, (Set, BaseMultiset)):
            return NotImplemented
        if isinstance(other, BaseMultiset):
            return self._total == other._total and self._elements == other._elements
        if self._total != len(other):
            return False
        return self._issubset(other, False)

    def __ne__(self, other):
        if not isinstance(other, (Set, BaseMultiset)):
            return NotImplemented
        if isinstance(other, BaseMultiset):
            return self._total != other._total or self._elements != other._elements
        if self._total != len(other):
            return True
        return not self._issubset(other, False)

    def get(self, element, default):
        """Return the multiplicity for *element* if it is in the multiset, else *default*.

        Makes the *default* argument of the original :meth:`dict.get` non-optional.

        Args:
            element: The element of which to get the multiplicity.
            default: The default value to return if the element if not in the multiset.

        Returns:
            The multiplicity for *element* if it is in the multiset, else *default*.
        """
        return self._elements.get(element, default)

    @classmethod
    def from_elements(cls, elements, multiplicity):
        """Create a new multiset with the given *elements* and each multiplicity set to *multiplicity*.

        Uses :meth:`dict.fromkeys` internally.

        Args:
            elements: The element for the new multiset.
            multiplicity: The multiplicity for all elements.

        Returns:
            The new multiset.
        """
        return cls(dict.fromkeys(elements, multiplicity))

    def copy(self):
        """Return a shallow copy of the multiset."""
        return self.__class__(self)

    __copy__ = copy

    def items(self):
        return self._elements.items()

    def distinct_elements(self):
        return self._elements.keys()

    def multiplicities(self):
        return self._elements.values()

    @classmethod
    def _as_multiset(cls, other):
        if not isinstance(other, BaseMultiset):
            if not isinstance(other, Iterable):
                raise TypeError("'%s' object is not iterable" % type(other))
            return cls(other)
        return other

    @staticmethod
    def _as_mapping(iterable):
        if isinstance(iterable, BaseMultiset):
            return iterable._elements
        if isinstance(iterable, Mapping):
            return iterable
        if not isinstance(iterable, Iterable):
            raise TypeError("'%s' object is not iterable" % type(iterable))
        mapping = dict()
        for element in iterable:
            if element in mapping:
                mapping[element] += 1
            else:
                mapping[element] = 1
        return mapping


class Multiset(BaseMultiset):
    """The mutable multiset variant."""
    __slots__ = ()

    def __setitem__(self, element, multiplicity):
        """Set the element's multiplicity.

        This will remove the element if the multiplicity is less than or equal to zero.
        '"""
        if not isinstance(multiplicity, int):
            raise TypeError('multiplicity must be an integer')
        _elements = self._elements
        if element in _elements:
            old_multiplicity = _elements[element]
            if multiplicity > 0:
                _elements[element] = multiplicity
                self._total += multiplicity - old_multiplicity
            else:
                del _elements[element]
                self._total -= old_multiplicity
        elif multiplicity > 0:
            _elements[element] = multiplicity
            self._total += multiplicity

    def __delitem__(self, element):
        del self._elements[element]

    def update(self, *others):
        r"""Like :meth:`dict.update` but add multiplicities instead of replacing them.

        >>> ms = Multiset('aab')
        >>> ms.update('abc')
        >>> sorted(ms)
        ['a', 'a', 'a', 'b', 'b', 'c']

        Note that the operator ``+=`` is equivalent to :meth:`update`, except that the operator will only
        accept sets to avoid accidental errors.

        >>> ms += Multiset('bc')
        >>> sorted(ms)
        ['a', 'a', 'a', 'b', 'b', 'b', 'c', 'c']

        For a variant of the operation which does not modify the multiset, but returns a new
        multiset instead see :meth:`combine`.

        Args:
            others: The other sets to add to this multiset. Can also be any :class:`~typing.Iterable`\[~T]
                or :class:`~typing.Mapping`\[~T, :class:`int`] which are then converted to :class:`Multiset`\[~T].
        """
        _elements = self._elements
        for other in map(self._as_mapping, others):
            for element, multiplicity in other.items():
                self[element] += multiplicity

    def union_update(self, *others):
        r"""Update the multiset, adding elements from all others using the maximum multiplicity.

        >>> ms = Multiset('aab')
        >>> ms.union_update('bc')
        >>> sorted(ms)
        ['a', 'a', 'b', 'c']

        You can also use the ``|=`` operator for the same effect. However, the operator version
        will only accept a set as other operator, not any iterable, to avoid errors.

        >>> ms = Multiset('aab')
        >>> ms |= Multiset('bccd')
        >>> sorted(ms)
        ['a', 'a', 'b', 'c', 'c', 'd']

        For a variant of the operation which does not modify the multiset, but returns a new
        multiset instead see :meth:`union`.

        Args:
            others: The other sets to union this multiset with. Can also be any :class:`~typing.Iterable`\[~T]
                or :class:`~typing.Mapping`\[~T, :class:`int`] which are then converted to :class:`Multiset`\[~T].
        """
        _elements = self._elements
        _total = self._total
        for other in map(self._as_mapping, others):
            for element, multiplicity in other.items():
                old_multiplicity = _elements.get(element, 0)
                if multiplicity > old_multiplicity:
                    _elements[element] = multiplicity
                    _total += multiplicity - old_multiplicity
        self._total = _total

    def __ior__(self, other):
        if not isinstance(other, (Set, BaseMultiset)):
            return NotImplemented
        self.union_update(other)
        return self

    def intersection_update(self, *others):
        r"""Update the multiset, keeping only elements found in it and all others.

        >>> ms = Multiset('aab')
        >>> ms.intersection_update('bc')
        >>> sorted(ms)
        ['b']

        You can also use the ``&=`` operator for the same effect. However, the operator version
        will only accept a set as other operator, not any iterable, to avoid errors.

        >>> ms = Multiset('aabc')
        >>> ms &= Multiset('abbd')
        >>> sorted(ms)
        ['a', 'b']

        For a variant of the operation which does not modify the multiset, but returns a new
        multiset instead see :meth:`intersection`.

        Args:
            others: The other sets to intersect this multiset with. Can also be any :class:`~typing.Iterable`\[~T]
                or :class:`~typing.Mapping`\[~T, :class:`int`] which are then converted to :class:`Multiset`\[~T].
        """
        for other in map(self._as_mapping, others):
            for element, current_count in list(self.items()):
                multiplicity = other.get(element, 0)
                if multiplicity < current_count:
                    self[element] = multiplicity

    def __iand__(self, other):
        if not isinstance(other, (Set, BaseMultiset)):
            return NotImplemented
        self.intersection_update(other)
        return self

    def difference_update(self, *others):
        r"""Remove all elements contained the others from this multiset.

        >>> ms = Multiset('aab')
        >>> ms.difference_update('abc')
        >>> sorted(ms)
        ['a']

        You can also use the ``-=`` operator for the same effect. However, the operator version
        will only accept a set as other operator, not any iterable, to avoid errors.

        >>> ms = Multiset('aabbbc')
        >>> ms -= Multiset('abd')
        >>> sorted(ms)
        ['a', 'b', 'b', 'c']

        For a variant of the operation which does not modify the multiset, but returns a new
        multiset instead see :meth:`difference`.

        Args:
            others: The other sets to remove from this multiset. Can also be any :class:`~typing.Iterable`\[~T]
                or :class:`~typing.Mapping`\[~T, :class:`int`] which are then converted to :class:`Multiset`\[~T].
        """
        for other in map(self._as_multiset, others):
            for element, multiplicity in other.items():
                self.discard(element, multiplicity)

    def __isub__(self, other):
        if not isinstance(other, (Set, BaseMultiset)):
            return NotImplemented
        self.difference_update(other)
        return self

    def symmetric_difference_update(self, other):
        r"""Update the multiset to contain only elements in either this multiset or the other but not both.

        >>> ms = Multiset('aab')
        >>> ms.symmetric_difference_update('abc')
        >>> sorted(ms)
        ['a', 'c']

        You can also use the ``^=`` operator for the same effect. However, the operator version
        will only accept a set as other operator, not any iterable, to avoid errors.

        >>> ms = Multiset('aabbbc')
        >>> ms ^= Multiset('abd')
        >>> sorted(ms)
        ['a', 'b', 'b', 'c', 'd']

        For a variant of the operation which does not modify the multiset, but returns a new
        multiset instead see :meth:`symmetric_difference`.

        Args:
            other: The other set to take the symmetric difference with. Can also be any :class:`~typing.Iterable`\[~T]
                or :class:`~typing.Mapping`\[~T, :class:`int`] which are then converted to :class:`Multiset`\[~T].
        """
        other = self._as_multiset(other)
        elements = set(self.distinct_elements()) | set(other.distinct_elements())
        for element in elements:
            multiplicity = self[element]
            other_count = other[element]
            self[element] = (multiplicity - other_count if multiplicity > other_count else other_count - multiplicity)

    def __ixor__(self, other):
        if not isinstance(other, (Set, BaseMultiset)):
            return NotImplemented
        self.symmetric_difference_update(other)
        return self

    def times_update(self, factor):
        """Update each this multiset by multiplying each element's multiplicity with the given scalar factor.

        >>> ms = Multiset('aab')
        >>> ms.times_update(2)
        >>> sorted(ms)
        ['a', 'a', 'a', 'a', 'b', 'b']

        You can also use the ``*=`` operator for the same effect:

        >>> ms = Multiset('ac')
        >>> ms *= 3
        >>> sorted(ms)
        ['a', 'a', 'a', 'c', 'c', 'c']

        For a variant of the operation which does not modify the multiset, but returns a new
        multiset instead see :meth:`times`.

        Args:
            factor: The factor to multiply each multiplicity with.
        """
        if factor <= 0:
            self.clear()
        else:
            for element in self.distinct_elements():
                self[element] *= factor

    def __imul__(self, factor):
        self.times_update(factor)
        return self

    def add(self, element, multiplicity=1):
        """Adds an element to the multiset.

        >>> ms = Multiset()
        >>> ms.add('a')
        >>> sorted(ms)
        ['a']

        An optional multiplicity can be specified to define how many of the element are added:

        >>> ms.add('b', 2)
        >>> sorted(ms)
        ['a', 'b', 'b']

        This extends the :meth:`MutableSet.add` signature to allow specifying the multiplicity.

        Args:
            element:
                The element to add to the multiset.
            multiplicity:
                The multiplicity i.e. count of elements to add.
        """
        if multiplicity < 1:
            raise ValueError("Multiplicity must be positive")
        self._elements[element] += multiplicity
        self._total += multiplicity

    def remove(self, element, multiplicity=None):
        """Removes an element from the multiset.

        If no multiplicity is specified, the element is completely removed from the multiset:

        >>> ms = Multiset('aabbbc')
        >>> ms.remove('a')
        2
        >>> sorted(ms)
        ['b', 'b', 'b', 'c']

        If the multiplicity is given, it is subtracted from the element's multiplicity in the multiset:

        >>> ms.remove('b', 2)
        3
        >>> sorted(ms)
        ['b', 'c']

        It is not an error to remove more elements than are in the set:

        >>> ms.remove('b', 2)
        1
        >>> sorted(ms)
        ['c']

        This extends the :meth:`MutableSet.remove` signature to allow specifying the multiplicity.

        Args:
            element:
                The element to remove from the multiset.
            multiplicity:
                An optional multiplicity i.e. count of elements to remove.

        Returns:
            The multiplicity of the element in the multiset before
            the removal.

        Raises:
            KeyError: if the element is not contained in the set. Use :meth:`discard` if
                you do not want an exception to be raised.
        """
        _elements = self._elements
        if element not in _elements:
            raise KeyError
        old_multiplicity = _elements[element]
        if multiplicity is None or multiplicity >= old_multiplicity:
            del _elements[element]
            self._total -= old_multiplicity
        elif multiplicity < 0:
            raise ValueError("Multiplicity must be not be negative")
        elif multiplicity > 0:
            _elements[element] -= multiplicity
            self._total -= multiplicity
        return old_multiplicity

    def discard(self, element, multiplicity=None):
        """Removes the `element` from the multiset.

        If multiplicity is ``None``, all occurrences of the element are removed:

        >>> ms = Multiset('aab')
        >>> ms.discard('a')
        2
        >>> sorted(ms)
        ['b']

        Otherwise, the multiplicity is subtracted from the one in the multiset and the
        old multiplicity is removed:

        >>> ms = Multiset('aab')
        >>> ms.discard('a', 1)
        2
        >>> sorted(ms)
        ['a', 'b']

        In contrast to :meth:`remove`, this does not raise an error if the
        element is not in the multiset:

        >>> ms = Multiset('a')
        >>> ms.discard('b')
        0
        >>> sorted(ms)
        ['a']

        It is also not an error to remove more elements than are in the set:

        >>> ms.remove('a', 2)
        1
        >>> sorted(ms)
        []

        Args:
            element:
                The element to remove from the multiset.
            multiplicity:
                An optional multiplicity i.e. count of elements to remove.

        Returns:
            The multiplicity of the element in the multiset before
            the removal.
        """
        _elements = self._elements
        if element in _elements:
            old_multiplicity = _elements[element]
            if multiplicity is None or multiplicity >= old_multiplicity:
                del _elements[element]
                self._total -= old_multiplicity
            elif multiplicity < 0:
                raise ValueError("Multiplicity must not be negative")
            elif multiplicity > 0:
                _elements[element] -= multiplicity
                self._total -= multiplicity
            return old_multiplicity
        else:
            return 0

    def pop(self, element, default):
        """If *element* is in the multiset, remove it and return its multiplicity, else return *default*.

        Makes the *default* argument of the original :meth:`dict.pop` non-optional.

        Args:
            element: The element which is removed.
            default: The default value to return if the element if not in the multiset.

        Returns:
            The multiplicity for *element* if it is in the multiset, else *default*.
        """
        return self._elements.pop(element, default)

    def setdefault(self, element, default):
        """If *element* is in the multiset, return its multiplicity.
        Else add it with a multiplicity of *default* and return *default*.

        Makes the *default* argument of the original :meth:`dict.setdefault` non-optional.

        Args:
            element: The element which is added if not already present.
            default: The default multiplicity to add the element with if not in the multiset.

        Returns:
            The multiplicity for *element* if it is in the multiset, else *default*.
        """
        return self._elements.setdefault(element, default)

    def clear(self):
        """Empty the multiset."""
        self._elements.clear()
        self._total = 0


class FrozenMultiset(BaseMultiset):
    """The frozen multiset variant that is immutable and hashable."""
    __slots__ = ()

    def __hash__(self):
        return hash(tuple(sorted(self._elements.items())))

Mapping.register(BaseMultiset)  # type: ignore
MutableMapping.register(Multiset)  # type: ignore

if __name__ == '__main__':
    import doctest
    doctest.testmod()
