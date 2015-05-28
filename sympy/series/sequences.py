from __future__ import print_function, division

from sympy.core.expr import Expr
from sympy.core.singleton import (S, Singleton)
from sympy.core.compatibility import (range, integer_types, with_metaclass\
                                      , is_sequence)
from sympy.core.sympify import sympify
from sympy.core.containers import Tuple
from sympy.functions.elementary.integers import ceiling
from sympy.utilities.misc import filldedent
from sympy.sets.sets import Interval, Set


class SeqBase(Expr):
    """Base class for sequences"""

    is_iterable = False

    is_EmptySequence = False
    is_Periodic = False
    is_Functional = False
    is_Formula = False

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

    @property
    def dummy(self):
        """Returns a tuple of variables that are bounded"""
        return self._dummy

    @property
    def _dummy(self):
        return ()

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
    is_EmptySequence = True

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

    Allowed:
    * Any instance of set
    * (Set, step)
    * (start, stop)
    * (start, stop, step)

    returns an Interval object and a step value

    Examples
    ========

    >>> from sympy.series.sequences import _parse_interval as pari
    >>> from sympy import Interval
    >>> pari(Interval(0, 5))
    ([0, 5], 1)
    >>> pari((Interval(0, 5), 2))
    ([0, 5], 2)
    >>> pari((0, 5))
    ([0, 5], 1)
    >>> pari((0, 5, 2))
    ([0, 5], 2)
    """
    start, stop, step = None, None, None
    if isinstance(interval, Set):
        start, stop = interval.inf, interval.sup
    elif is_sequence(interval, Tuple):
        if len(interval) == 2:
            if isinstance(interval[0], Set):
                start, stop = interval[0].inf, interval[0].sup
                step = interval[1]
            else:
                start, stop = interval
        elif len(interval) == 3:
            start, stop, step = interval

    if step == None:
        step = 1 # default

    if start == None or stop == None:
        raise ValueError('Invalid limits given: %s' % str(interval))

    if start is S.NegativeInfinity and stop is S.Infinity:
            raise ValueError(filldedent("""Both the start and end value\
                                        cannot be unbounded"""))

    if step in [S.NegativeInfinity, S.Infinity]:
        raise ValueError("step cannot be unbounded")

    return (Interval(start, stop), sympify(step))


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

    changing the step size

    >>> SeqExpr((1, 2, 3), (0, 10, 2)).length
    6

    """

    is_iterable = True

    def __new__(cls, gen, interval=(0, S.Infinity, 1)):
        bounds, step = _parse_interval(interval)
        if bounds is S.EmptySet:
            return S.EmptySequence
        gen = sympify(gen)
        interval = Tuple(bounds, step)
        return Expr.__new__(cls, gen, interval)

    def __iter__(self):
        i = 0
        while(i < self.length):
            pt = self._ith_point(i)
            yield self.coeff(pt)
            i += 1

    @property
    def free_symbols(self):
        """
        This method returns the symbols in the object, excluding those
        that take on a specific value (i.e. the dummy symbols).

        Examples
        ========

        >>> from sympy import SeqFormula
        >>> from sympy.abc import n, m
        >>> SeqFormula((m*n**2, n), (0, 5)).free_symbols
        set([m])
        """
        fsyms = set().union(*[a.free_symbols for a in self.args])
        for d in self.dummy:
            if d in fsyms:
                fsyms.remove(d)
        return fsyms

    @property
    def _gen(self):
        return self.args[0]

    @property
    def _interval(self):
        return self.args[1][0]

    @property
    def _start(self):
        return self.interval.inf

    @property
    def _stop(self):
        return self.interval.sup

    @property
    def _step(self):
        return self.args[1][1]

    @property
    def _length(self):
        return ceiling((self.stop - self.start + 1) / self.step)


class SeqPer(SeqExpr):
    """Represents a periodical sequence

    The elements are repeated after a given period.

    Examples
    ========

    >>> from sympy import SeqPer, oo
    >>> s = SeqPer((1, 2, 3), (0, 5))
    >>> s.periodical
    (1, 2, 3)
    >>> s.period
    3

    For value at a particular point

    >>> s.coeff(3)
    1

    supports slicing

    >>> s[:]
    [1, 2, 3, 1, 2, 3]

    iterable

    >>> list(s)
    [1, 2, 3, 1, 2, 3]

    changing step size

    >>> SeqPer((1, 2, 3), (0, 5, 2))[:]
    [1, 3, 2]

    sequence starts from negative infinity

    >>> SeqPer((1, 2, 3), (-oo, 0))[0:6]
    [1, 2, 3, 1, 2, 3]

    See Also
    ========

    sympy.series.sequences.SeqFormula
    sympy.series.sequences.SeqFunc
    """

    is_Periodic = True

    @property
    def period(self):
        return len(self.gen)

    @property
    def periodical(self):
        return self.gen

    def _eval_coeff(self, i):
        if self.start is S.NegativeInfinity:
            idx = (self.stop - i) % self.period
        else:
            idx = (i - self.start) % self.period
        return self.periodical[idx]


class SeqFormula(SeqExpr):
    """Represents sequence based on a formula

    Elements are generated using a formula

    Examples
    ========

    >>> from sympy import SeqFormula, oo, Symbol
    >>> n = Symbol('n')
    >>> s = SeqFormula((n**2, n), (0, 5))
    >>> s.formula
    n**2

    For value at a particular point

    >>> s.coeff(3)
    9

    supports slicing

    >>> s[:]
    [0, 1, 4, 9, 16, 25]

    iterable

    >>> list(s)
    [0, 1, 4, 9, 16, 25]

    changing step size

    >>> SeqFormula((n**2, n), (0, 5, 2))[:]
    [0, 4, 16]

    sequence starts from negative infinity

    >>> SeqFormula((n**2, n), (-oo, 0))[0:6]
    [0, 1, 4, 9, 16, 25]

    See Also
    ========

    sympy.series.sequences.SeqPer
    sympy.series.sequences.SeqFunc
    """

    is_Formula = True

    def __new__(cls, formula, interval=(0, S.Infinity, 1)):
        # try to find the dummy symbol
        formula = sympify(formula)
        if not is_sequence(formula, Tuple):
            free = formula.free_symbols
            if len(free) != 1:
                raise ValueError(filldedent(
                    " specify dummy variables for %s. If the formula contains"
                    " more than one free symbol, a dummy variable should be"
                    " supplied explicitly e.g., SeqFormula((m*n**2, n), (0, 5))"
                    % formula))
            formula = (formula, free.pop())
        return SeqExpr.__new__(cls, formula, interval)

    @property
    def formula(self):
        return self.gen[0]

    @property
    def _dummy(self):
        return (self.gen[1],)

    def _eval_coeff(self, i):
        d = self.dummy[0]
        return self.formula.subs(d, i)


class SeqFunc(SeqExpr):
    """Represents sequence based on a function

    Elements are generated by calling a function.
    Only single argument functions are allowed.

    Examples
    ========

    >>> from sympy import SeqFunc, oo, Lambda, Symbol
    >>> n = Symbol('n')
    >>> s = SeqFunc(Lambda(n, n**2), (0, 5))
    >>> s.function
    Lambda(n, n**2)

    For value at a particular point

    >>> s.coeff(3)
    9

    supports slicing

    >>> s[:]
    [0, 1, 4, 9, 16, 25]

    iterable

    >>> list(s)
    [0, 1, 4, 9, 16, 25]

    changing step size

    >>> SeqFunc(Lambda(n, n**2), (0, 5, 2))[:]
    [0, 4, 16]

    sequence starts from negative infinity

    >>> SeqFunc(Lambda(n, n**2), (-oo, 0))[0:6]
    [0, 1, 4, 9, 16, 25]

    See Also
    ========

    sympy.series.sequences.SeqPer
    sympy.series.sequences.SeqFormula
    """

    is_Functional = True

    def __new__(cls, function, interval=(0, S.Infinity, 1)):
        function = sympify(function)
        if len(function.variables) != 1:
            raise ValueError(filldedent(
                "Only single argument functions are allowed"))
        return SeqExpr.__new__(cls, function, interval)

    @property
    def function(self):
        return self.gen

    def _eval_coeff(self, i):
        return self.function(i)


def sequence(**kwargs):
    """Returns appropriate sequence object.

    Examples
    ========

    >>> from sympy import sequence, SeqPer, SeqFunc, SeqFormula, Lambda
    >>> from sympy.abc import n

    >>> sequence(formula=(n**2, n), interval=(0, 5))
    SeqFormula((n**2, n), ([0, 5], 1))

    >>> sequence(periodical=(1, 2, 3), interval=(0, 5))
    SeqPer((1, 2, 3), ([0, 5], 1))

    >>> sequence(func=Lambda(n, n**2), interval=(0, 5))
    SeqFunc(Lambda(n, n**2), ([0, 5], 1))

    See Also
    ========

    sympy.series.sequences.SeqPer
    sympy.series.sequences.SeqFormula
    sympy.series.sequences.SeqFunc
    """
    interval = kwargs.pop('interval', None)
    if interval == None:
        interval = (0, S.Infinity, 1)

    key = kwargs.pop('periodical', None)
    if not key == None:
        return SeqPer(key, interval)

    key = kwargs.pop('func', None)
    if not key == None:
        return SeqFunc(key, interval)

    key = kwargs.pop('formula', None)
    if not key == None:
        return SeqFormula(key, interval)

    raise ValueError('Invalid Arguments')
