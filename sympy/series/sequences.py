from sympy.core.expr import Expr
from sympy.core.singleton import (S, Singleton)
from sympy.core.compatibility import (range, integer_types, with_metaclass)
from sympy.core.sympify import sympify
from sympy.utilities.misc import filldedent
from sympy.sets.sets import Interval


class SeqBase(Expr):
    """Base class for sequences"""

    @property
    def gen(self):
        """Returns the generator for the sequence"""
        return self._gen

    @property
    def _gen(self):
        raise NotImplementedError("(%s)._gen" % self)

    @property
    def interval(self):
        """The interval in which the sequence is defined"""
        return self._interval

    @property
    def _interval(self):
        raise NotImplementedError("(%s)._interval" % self)


    @property
    def start(self):
        """The starting point of the sequence. This point is included"""
        return self._start

    @property
    def _start(self):
        raise NotImplementedError("(%s)._start" % self)

    @property
    def end(self):
        """The ending point of the sequence. This point is included"""
        return self._end

    @property
    def _end(self):
        raise NotImplementedError("(%s)._end" % self)

    @property
    def length(self):
        """Length of the sequence"""
        return self._length

    @property
    def _length(self):
        raise NotImplementedError("(%s)._length" % self)

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


def _parse_interval(interval):
    if len(interval) != 2 or None in interval:
        raise ValueError(filldedent("""Sequence requires values for lower\
                                    and upper bounds."""))
    return Interval(*interval)


class SeqExpr(SeqBase):
    """Sequence expression class
    Various sequences (SeqPer, SeqFormula, SeqFunc...) should inherit from
    this class

    Examples
    ========

    >>> from sympy.series.sequences import SeqExpr
    >>> s = SeqExpr((1, 2, 3), (0, 10))
    >>> s.gen
    (1, 2, 3)
    >>> s.interval
    [0, 10]
    >>> s.length
    11

    """
    def __new__(cls, gen, interval=(0, S.Infinity)):
        interval = _parse_interval(interval)
        if interval is S.EmptySet:
            return S.EmptySequence
        gen = sympify(gen)
        return Expr.__new__(cls, gen, interval)

    @property
    def _gen(self):
        return self.args[0]

    @property
    def _interval(self):
        return self.args[1]

    @property
    def _start(self):
        return self.interval.start

    @property
    def _end(self):
        return self.interval.end

    @property
    def _length(self):
        return self.end - self.start + 1
