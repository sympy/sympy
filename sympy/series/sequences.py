from __future__ import print_function, division

from sympy.core.expr import Expr
from sympy.core.singleton import (S, Singleton)
from sympy.core.compatibility import (range, integer_types, with_metaclass)
from sympy.core.sympify import sympify
from sympy.functions.elementary.integers import ceiling
from sympy.utilities.misc import filldedent
from sympy.sets.sets import Interval


class SeqBase(Expr):
    """Base class for sequences"""

    is_iterable = False

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
    def stop(self):
        """The ending point of the sequence. This point is included"""
        return self._stop

    @property
    def _stop(self):
        raise NotImplementedError("(%s)._stop" % self)

    @property
    def step(self):
        """Increase points by step"""
        return self._step

    @property
    def _step(self):
        raise NotImplementedError("(%s)._step" % self)

    @property
    def length(self):
        """Length of the sequence"""
        return self._length

    @property
    def _length(self):
        raise NotImplementedError("(%s)._length" % self)

    def coeff(self, i):
        """Returns the coefficient at point i"""
        if i < self.start or i > self.stop:
            raise IndexError("Index %s out of bounds %s" %(i, self.interval))
        return self._eval_coeff(i)

    def _eval_coeff(self, i):
        raise NotImplementedError(filldedent(""" The _eval_coeff method should\
                                             be added to %s to return\
                                             coefficient so it is available\
                                             when coeff calls it."""\
                                             % self.func))

    def _ith_point(self, i, step=None):
        """
        Returns the i'th point of a sequence
        If start point is negative infinity, point is returned from the end.
        Assumes the first point to be indexed zero.

        Examples
        =========

        >>> from sympy import oo
        >>> from sympy.series.sequences import SeqExpr

        bounded

        >>> SeqExpr((1, 2, 3), (-10, 10))._ith_point(0)
        -10
        >>> SeqExpr((1, 2, 3), (-10, 10))._ith_point(5)
        -5
        >>> SeqExpr((1, 2, 3), (-10, 10, 2))._ith_point(5)
        0

        End is at infinity

        >>> SeqExpr((1, 2, 3), (0, oo))._ith_point(5)
        5
        >>> SeqExpr((1, 2, 3), (0, oo, 2))._ith_point(5)
        10

        Starts at negative infinity

        >>> SeqExpr((1, 2, 3), (-oo, 0))._ith_point(5)
        -5
        >>> SeqExpr((1, 2, 3), (-oo, 0, 2))._ith_point(5)
        -10

        """
        if self.start is S.NegativeInfinity:
            initial = self.stop
        else:
            initial = self.start

        if step == None:
            if self.start is S.NegativeInfinity:
                step = -self.step
            else:
                step = self.step

        return initial + i*step

    def __getitem__(self, index):
        if isinstance(index, integer_types):
            index = self._ith_point(index)
            return self.coeff(index)
        elif isinstance(index, slice):
            start, stop = index.start, index.stop
            if start == None:
                start = 0
            if stop == None:
                stop = self.length
            return [self.coeff(self._ith_point(i, index.step)) for i in\
                               range(start, stop)]


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

    is_iterable = True

    @property
    def _interval(self):
        return S.EmptySet

    @property
    def _length(self):
        return S.Zero

    def __iter__(self):
        return iter([])


def _parse_interval(interval):
    """
    interval should be of the form (start, step) or (start, stop, step)
    Both start and stop cannot be unbounded
    step cannot be unbounded
    """
    if len(interval) == 2:
        start, stop = interval
        step = 1
    elif len(interval) == 3:
        start, stop, step = interval
        if step == None:
            step = 1
    else:
        raise ValueError(filldedent("""interval should be of the form\
                                    (start, stop) or (start, stop, step)"""))

    if start is S.NegativeInfinity and stop is S.Infinity:
            raise ValueError(filldedent("""Both the start and end value\
                                        cannot be unbounded"""))
    if step in [S.NegativeInfinity, S.Infinity]:
        raise ValueError("step cannot be unbounded")

    return Interval(start, stop), step


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

    change the step size

    >>> SeqExpr((1, 2, 3), (0, 10, 2)).length
    6

    """

    is_iterable = True

    def __new__(cls, gen, interval=(0, S.Infinity, 1)):
        interval, step = _parse_interval(interval)
        if interval is S.EmptySet:
            return S.EmptySequence
        gen = sympify(gen)
        return Expr.__new__(cls, gen, interval, step)

    def __iter__(self):
        i = 0
        while(i < self.length):
            pt = self._ith_point(i)
            yield self.coeff(pt)
            i += 1

    @property
    def _gen(self):
        return self.args[0]

    @property
    def _interval(self):
        return self.args[1]

    @property
    def _start(self):
        return self.interval.inf

    @property
    def _stop(self):
        return self.interval.sup

    @property
    def _step(self):
        return self.args[2]

    @property
    def _length(self):
        return ceiling((self.stop - self.start + 1) / self.step)
