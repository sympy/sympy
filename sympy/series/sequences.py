from sympy.core.expr import Expr
from sympy.core.singleton import (S, Singleton)
from sympy.core.compatibility import (range, integer_types, with_metaclass)
from sympy.utilities.misc import filldedent


class SeqBase(Expr):
    """Base class for sequences"""

    @property
    def gen(self):
        """Returns the generator for the sequence"""
        return self._gen

    @property
    def _gen(self):
        raise NotImplementedError()

    @property
    def interval(self):
        """The interval in which the sequence is defined"""
        return self._interval

    @property
    def _interval(self):
        raise NotImplementedError()

    @property
    def start(self):
        """The starting point of the sequence. This point is included"""
        return self._start

    @property
    def _start(self):
        raise NotImplementedError()

    @property
    def end(self):
        """The ending point of the sequence. This point is included"""
        return self._end

    @property
    def _end(self):
        raise NotImplementedError()

    @property
    def length(self):
        """Length of the sequence"""
        return self._length

    @property
    def _length(self):
        raise NotImplementedError()

    def coeff(self, i):
        """Returns the coefficient at point i"""
        if i < self.start or i > self.end:
            raise IndexError("Index %s out of bounds %s" %(i, self.interval))
        return self._eval_coeff(i)

    def _eval_coeff(self, i):
        raise NotImplementedError(filldedent(""" The _eval_coeff method should\
                                             be added to %s to give\
                                             coefficients so it is available\
                                             when coeff calls it."""\
                                             % self.func))

    def __getitem__(self, index):
        if isinstance(index, integer_types):
            return self.coeff(index + self.start)
        elif isinstance(index, slice):
            start, stop = index.start, index.stop
            if start == None:
                start = 0
            if stop == None:
                stop = self.length
            return [self.coeff(i + self.start) for i in range(start, stop,\
                                                 index.step or 1)]


class EmptySequence(with_metaclass(Singleton, SeqBase)):
    """
    Represents an empty sequence. The empty sequence is available as a
    singleton as S.EmptySequence.

    Examples
    ========

    >>> from sympy import S
    >>> S.EmptySequence
    EmptySequence()
    """

    @property
    def _interval(self):
        return S.EmptySet

    @property
    def _length(self):
        return S.Zero
