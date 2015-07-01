from __future__ import print_function, division

from sympy.core.basic import Basic
from sympy.core.mul import Mul
from sympy.core.singleton import S, Singleton
from sympy.core.symbol import Dummy, Symbol
from sympy.core.compatibility import (range, integer_types, with_metaclass,
                                      is_sequence, iterable, ordered)
from sympy.core.decorators import call_highest_priority
from sympy.core.cache import cacheit
from sympy.core.sympify import sympify
from sympy.core.containers import Tuple
from sympy.core.evaluate import global_evaluate
from sympy.polys import lcm
from sympy.sets.sets import Interval, Intersection
from sympy.utilities.iterables import flatten
from sympy.tensor.indexed import Idx


###############################################################################
#                            SEQUENCES                                        #
###############################################################################


class SeqBase(Basic):
    """Base class for sequences"""

    is_commutative = True
    _op_priority = 15

    @staticmethod
    def _start_key(expr):
        """
        Return start (if possible) else S.Infinity.
        adapted from Set._infimum_key
        """
        try:
            start = expr.start
        except (NotImplementedError,
                AttributeError, ValueError):
            start = S.Infinity
        return start

    def _intersect_interval(self, other):
        """
        returns the start, stop
        takes intersection over the two intervals
        """
        interval = Intersection(self.interval, other.interval)
        return interval.inf, interval.sup

    @property
    def gen(self):
        """Returns the generator for the sequence"""
        raise NotImplementedError("(%s).gen" % self)

    @property
    def interval(self):
        """The interval on which the sequence is defined"""
        raise NotImplementedError("(%s).interval" % self)

    @property
    def start(self):
        """The starting point of the sequence. This point is included"""
        raise NotImplementedError("(%s).start" % self)

    @property
    def stop(self):
        """The ending point of the sequence. This point is included"""
        raise NotImplementedError("(%s).stop" % self)

    @property
    def length(self):
        """Length of the sequence"""
        raise NotImplementedError("(%s).length" % self)

    @property
    def variables(self):
        """Returns a tuple of variables that are bounded"""
        return ()

    @property
    def free_symbols(self):
        """
        This method returns the symbols in the object, excluding those
        that take on a specific value (i.e. the dummy symbols).

        Examples
        ========

        >>> from sympy import SeqFormula
        >>> from sympy.abc import n, m
        >>> SeqFormula(m*n**2, (n, 0, 5)).free_symbols
        set([m])
        """
        return (set(j for i in self.args for j in i.free_symbols
                   .difference(self.variables)))

    @cacheit
    def coeff(self, pt):
        """Returns the coefficient at point pt"""
        if pt < self.start or pt > self.stop:
            raise IndexError("Index %s out of bounds %s" % (pt, self.interval))
        return self._eval_coeff(pt)

    def _eval_coeff(self, pt):
        raise NotImplementedError("The _eval_coeff method should be added to"
                                  "%s to return coefficient so it is available"
                                  "when coeff calls it."
                                  % self.func)

    def _ith_point(self, i):
        """
        Returns the i'th point of a sequence
        If start point is negative infinity, point is returned from the end.
        Assumes the first point to be indexed zero.

        Examples
        =========

        >>> from sympy import oo
        >>> from sympy.series.sequences import SeqPer

        bounded

        >>> SeqPer((1, 2, 3), (-10, 10))._ith_point(0)
        -10
        >>> SeqPer((1, 2, 3), (-10, 10))._ith_point(5)
        -5

        End is at infinity

        >>> SeqPer((1, 2, 3), (0, oo))._ith_point(5)
        5

        Starts at negative infinity

        >>> SeqPer((1, 2, 3), (-oo, 0))._ith_point(5)
        -5
        """
        if self.start is S.NegativeInfinity:
            initial = self.stop
        else:
            initial = self.start

        if self.start is S.NegativeInfinity:
            step = -1
        else:
            step = 1

        return initial + i*step

    def _add(self, other):
        """
        Should only be used internally

        self._add(other) returns a new, term-wise added sequence if self
        knows how to add with other, otherwise it returns ``None``.

        ``other`` should only be a sequence object.

        Used within :class:`SeqAdd` class
        """
        return None

    def _mul(self, other):
        """
        Should only be used internally

        self._mul(other) returns a new, term-wise multiplied sequence if self
        knows how to multiply with other, otherwise it returns ``None``.

        ``other`` should only be a sequence object.

        Used within :class:`SeqMul` class
        """
        return None

    def coeff_mul(self, other):
        """
        Should be used when ``other`` is not a sequence. Should be
        defined to define custom behaviour.

        Examples
        ========

        >>> from sympy import S, oo, SeqFormula
        >>> from sympy.abc import n

        >>> SeqFormula(n**2).coeff_mul(2)
        SeqFormula(2*n**2, (n, 0, oo))

        Notes
        =====

        '*' defines multiplication of sequences with sequences only.
        For multiplying sequences use ``mul`` method.
        """
        return Mul(self, other)

    def __add__(self, other):
        """
        Returns the term-wise addition of 'self' and 'other'.
        ``other`` should be a sequence.

        Examples
        ========

        >>> from sympy import S, oo, SeqFormula
        >>> from sympy.abc import n
        >>> SeqFormula(n**2) + SeqFormula(n**3)
        SeqFormula(n**3 + n**2, (n, 0, oo))
        """
        if not isinstance(other, SeqBase):
            raise TypeError('cannot add sequence and %s' % type(other))
        return SeqAdd(self, other)

    @call_highest_priority('__add__')
    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        """
        Returns the term-wise subtraction of 'self' and 'other'.
        ``other`` should be a sequence.

        Examples
        ========

        >>> from sympy import S, oo, SeqFormula
        >>> from sympy.abc import n
        >>> SeqFormula(n**2) - (SeqFormula(n))
        SeqFormula(n**2 - n, (n, 0, oo))
        """
        if not isinstance(other, SeqBase):
            raise TypeError('cannot subtract sequence and %s' % type(other))
        return SeqAdd(self, -other)

    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return (-self) + other

    def __neg__(self):
        """
        Negates the sequence

        Examples
        ========

        >>> from sympy import S, oo, SeqFormula
        >>> from sympy.abc import n
        >>> -SeqFormula(n**2)
        SeqFormula(-n**2, (n, 0, oo))
        """
        return self.coeff_mul(-1)

    def __mul__(self, other):
        """
        Returns the term-wise multiplication of 'self' and 'other'.
        ``other`` should be a sequence. For ``other`` not being a
        sequence see ``coeff_mul`` method.

        Examples
        ========

        >>> from sympy import S, oo, SeqFormula
        >>> from sympy.abc import n

        >>> SeqFormula(n**2) * (SeqFormula(n))
        SeqFormula(n**3, (n, 0, oo))
        """
        if not isinstance(other, SeqBase):
            raise TypeError('cannot multiply sequence and %s' % type(other))
        return SeqMul(self, other)

    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return self * other

    def __iter__(self):
        for i in range(self.length):
            pt = self._ith_point(i)
            yield self.coeff(pt)

    def __getitem__(self, index):
        if isinstance(index, integer_types):
            index = self._ith_point(index)
            return self.coeff(index)
        elif isinstance(index, slice):
            start, stop = index.start, index.stop
            if start is None:
                start = 0
            if stop is None:
                stop = self.length
            return [self.coeff(self._ith_point(i)) for i in
                    range(start, stop, index.step or 1)]


class EmptySequence(with_metaclass(Singleton, SeqBase)):
    """
    Represents an empty sequence. The empty sequence is available as a
    singleton as S.EmptySequence.

    Examples
    ========

    >>> from sympy import S, SeqPer, oo
    >>> from sympy.abc import x
    >>> S.EmptySequence
    EmptySequence()

    >>> SeqPer((1, 2), (x, 0, 10)) + S.EmptySequence
    SeqPer((1, 2), (x, 0, 10))
    >>> SeqPer((1, 2)) * S.EmptySequence
    EmptySequence()
    >>> S.EmptySequence.coeff_mul(-1)
    EmptySequence()
    """

    @property
    def interval(self):
        return S.EmptySet

    @property
    def length(self):
        return S.Zero

    def coeff_mul(self, coeff):
        """See docstring of SeqBase.coeff_mul"""
        return self

    def __iter__(self):
        return iter([])


class SeqExpr(SeqBase):
    """Sequence expression class
    Various sequences should inherit from this class

    Examples
    ========

    >>> from sympy.series.sequences import SeqExpr
    >>> from sympy.abc import x
    >>> s = SeqExpr((1, 2, 3), (x, 0, 10))
    >>> s.gen
    (1, 2, 3)
    >>> s.interval
    [0, 10]
    >>> s.length
    11

    See Also
    ========

    sympy.series.sequences.SeqPer
    sympy.series.sequences.SeqFormula
    """

    @property
    def gen(self):
        return self.args[0]

    @property
    def interval(self):
        return Interval(self.args[1][1], self.args[1][2])

    @property
    def start(self):
        return self.interval.inf

    @property
    def stop(self):
        return self.interval.sup

    @property
    def length(self):
        return self.stop - self.start + 1

    @property
    def variables(self):
        return (self.args[1][0],)


class SeqPer(SeqExpr):
    """Represents a periodic sequence

    The elements are repeated after a given period.

    Examples
    ========

    >>> from sympy import SeqPer, oo
    >>> from sympy.abc import k

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

    sequence starts from negative infinity

    >>> SeqPer((1, 2, 3), (-oo, 0))[0:6]
    [1, 2, 3, 1, 2, 3]

    Periodic formulas

    >>> SeqPer((k, k**2, k**3), (k, 0, oo))[0:6]
    [0, 1, 8, 3, 16, 125]

    See Also
    ========

    sympy.series.sequences.SeqFormula
    """

    def __new__(cls, periodical, limits=None):
        periodical = sympify(periodical)

        def _find_x(periodical):
            free = periodical.free_symbols
            if len(periodical.free_symbols) == 1:
                return free.pop()
            else:
                return Dummy('k')

        x, start, stop = None, None, None
        if limits is None:
            x, start, stop = _find_x(periodical), 0, S.Infinity
        if is_sequence(limits, Tuple):
            if len(limits) == 3:
                x, start, stop = limits
            elif len(limits) == 2:
                x = _find_x(periodical)
                start, stop = limits

        if not isinstance(x, (Symbol, Idx)) or start is None or stop is None:
            raise ValueError('Invalid limits given: %s' % str(limits))

        if start is S.NegativeInfinity and stop is S.Infinity:
                raise ValueError("Both the start and end value"
                                 "cannot be unbounded")

        limits = sympify((x, start, stop))

        if is_sequence(periodical, Tuple):
            periodical = sympify(tuple(flatten(periodical)))
        else:
            raise ValueError("invalid period %s should be something "
                             "like e.g (1, 2) " % periodical)

        if Interval(limits[1], limits[2]) is S.EmptySet:
            return S.EmptySequence

        return Basic.__new__(cls, periodical, limits)

    @property
    def period(self):
        return len(self.gen)

    @property
    def periodical(self):
        return self.gen

    def _eval_coeff(self, pt):
        if self.start is S.NegativeInfinity:
            idx = (self.stop - pt) % self.period
        else:
            idx = (pt - self.start) % self.period
        return self.periodical[idx].subs(self.variables[0], pt)

    def _add(self, other):
        """See docstring of SeqBase._add"""
        if isinstance(other, SeqPer):
            per1, lper1 = self.periodical, self.period
            per2, lper2 = other.periodical, other.period

            per_length = lcm(lper1, lper2)

            new_per = []
            for x in range(per_length):
                ele1 = per1[x % lper1]
                ele2 = per2[x % lper2]
                new_per.append(ele1 + ele2)

            start, stop = self._intersect_interval(other)
            return SeqPer(new_per, (self.variables[0], start, stop))

    def _mul(self, other):
        """See docstring of SeqBase._mul"""
        if isinstance(other, SeqPer):
            per1, lper1 = self.periodical, self.period
            per2, lper2 = other.periodical, other.period

            per_length = lcm(lper1, lper2)

            new_per = []
            for x in range(per_length):
                ele1 = per1[x % lper1]
                ele2 = per2[x % lper2]
                new_per.append(ele1 * ele2)

            start, stop = self._intersect_interval(other)
            return SeqPer(new_per, (self.variables[0], start, stop))

    def coeff_mul(self, coeff):
        """See docstring of SeqBase.coeff_mul"""
        coeff = sympify(coeff)
        per = [x * coeff for x in self.periodical]
        return SeqPer(per, self.args[1])


class SeqFormula(SeqExpr):
    """Represents sequence based on a formula

    Elements are generated using a formula

    Examples
    ========

    >>> from sympy import SeqFormula, oo, Symbol
    >>> n = Symbol('n')
    >>> s = SeqFormula(n**2, (n, 0, 5))
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

    sequence starts from negative infinity

    >>> SeqFormula(n**2, (-oo, 0))[0:6]
    [0, 1, 4, 9, 16, 25]

    See Also
    ========

    sympy.series.sequences.SeqPer
    """

    def __new__(cls, formula, limits=None):
        formula = sympify(formula)

        def _find_x(formula):
            free = formula.free_symbols
            if len(formula.free_symbols) == 1:
                return free.pop()
            elif len(formula.free_symbols) == 0:
                return Dummy('k')
            else:
                raise ValueError(
                    " specify dummy variables for %s. If the formula contains"
                    " more than one free symbol, a dummy variable should be"
                    " supplied explicitly e.g., SeqFormula(m*n**2, (n, 0, 5))"
                    % formula)

        x, start, stop = None, None, None
        if limits is None:
            x, start, stop = _find_x(formula), 0, S.Infinity
        if is_sequence(limits, Tuple):
            if len(limits) == 3:
                x, start, stop = limits
            elif len(limits) == 2:
                x = _find_x(formula)
                start, stop = limits

        if not isinstance(x, (Symbol, Idx)) or start is None or stop is None:
            raise ValueError('Invalid limits given: %s' % str(limits))

        if start is S.NegativeInfinity and stop is S.Infinity:
                raise ValueError("Both the start and end value"
                                 "cannot be unbounded")
        limits = sympify((x, start, stop))

        if Interval(limits[1], limits[2]) is S.EmptySet:
            return S.EmptySequence

        return Basic.__new__(cls, formula, limits)

    @property
    def formula(self):
        return self.gen

    def _eval_coeff(self, pt):
        d = self.variables[0]
        return self.formula.subs(d, pt)

    def _add(self, other):
        """See docstring of SeqBase._add"""
        if isinstance(other, SeqFormula):
            form1, v1 = self.formula, self.variables[0]
            form2, v2 = other.formula, other.variables[0]
            formula = form1 + form2.subs(v2, v1)
            start, stop = self._intersect_interval(other)
            return SeqFormula(formula, (v1, start, stop))

    def _mul(self, other):
        """See docstring of SeqBase._mul"""
        if isinstance(other, SeqFormula):
            form1, v1 = self.formula, self.variables[0]
            form2, v2 = other.formula, other.variables[0]
            formula = form1 * form2.subs(v2, v1)
            start, stop = self._intersect_interval(other)
            return SeqFormula(formula, (v1, start, stop))

    def coeff_mul(self, coeff):
        """See docstring of SeqBase.coeff_mul"""
        coeff = sympify(coeff)
        formula = self.formula * coeff
        return SeqFormula(formula, self.args[1])


def sequence(seq, limits=None):
    """Returns appropriate sequence object.
    If seq is a sympy sequence, returns SeqPer object
    otherwise returns SeqFormula object

    Examples
    ========

    >>> from sympy import sequence, SeqPer, SeqFormula
    >>> from sympy.abc import n

    >>> sequence(n**2, (n, 0, 5))
    SeqFormula(n**2, (n, 0, 5))

    >>> sequence((1, 2, 3), (n, 0, 5))
    SeqPer((1, 2, 3), (n, 0, 5))

    See Also
    ========

    sympy.series.sequences.SeqPer
    sympy.series.sequences.SeqFormula
    """
    seq = sympify(seq)

    if is_sequence(seq, Tuple):
        return SeqPer(seq, limits)
    else:
        return SeqFormula(seq, limits)


###############################################################################
#                            OPERATIONS                                       #
###############################################################################


class SeqExprOp(SeqBase):
    """Base class for operations on sequences

    Examples
    ========

    >>> from sympy.series.sequences import SeqExprOp, sequence
    >>> from sympy.abc import n
    >>> s1 = sequence(n**2, (n, 0, 10))
    >>> s2 = sequence((1, 2, 3), (n, 5, 10))
    >>> s = SeqExprOp(s1, s2)
    >>> s.gen
    (n**2, (1, 2, 3))
    >>> s.interval
    [5, 10]
    >>> s.length
    6

    See Also
    ========

    sympy.series.sequences.SeqAdd
    sympy.series.sequences.SeqMul
    """
    @property
    def gen(self):
        """Generator for the sequence
        returns a tuple of generators of all the argument sequences
        """
        return tuple(a.gen for a in self.args)

    @property
    def interval(self):
        """Sequence is defined on the intersection
        of all the intervals of respective sequences
        """
        return Intersection(a.interval for a in self.args)

    @property
    def start(self):
        return self.interval.inf

    @property
    def stop(self):
        return self.interval.sup

    @property
    def variables(self):
        """Cumulative of all the bound variables"""
        return tuple(flatten([a.variables for a in self.args]))

    @property
    def length(self):
        return self.stop - self.start + 1


class SeqAdd(SeqExprOp):
    """
    Represents term-wise addition of sequences

    Rules:
        * The interval on which sequence is defined is the intersection
        of respective intervals of sequences.
        * Anything + EmptySequence, remains unchanged
        * Other rules are defined in _add methods of sequence classes

    Examples
    ========

    >>> from sympy import S, oo, SeqAdd, SeqPer, SeqFormula
    >>> from sympy.abc import n
    >>> SeqAdd(SeqPer((1, 2), (n, 0, oo)), S.EmptySequence)
    SeqPer((1, 2), (n, 0, oo))
    >>> SeqAdd(SeqPer((1, 2), (n, 0, 5)), SeqPer((1, 2), (n, 6, 10)))
    EmptySequence()
    >>> SeqAdd(SeqPer((1, 2), (n, 0, oo)), SeqFormula(n**2, (n, 0, oo)))
    SeqAdd(SeqFormula(n**2, (n, 0, oo)), SeqPer((1, 2), (n, 0, oo)))
    >>> SeqAdd(SeqFormula(n**3), SeqFormula(n**2))
    SeqFormula(n**3 + n**2, (n, 0, oo))

    See Also
    ========

    sympy.series.sequences.SeqMul
    """

    def __new__(cls, *args, **kwargs):
        evaluate = kwargs.get('evaluate', global_evaluate[0])

        # flatten inputs
        args = list(args)

        # adapted from sympy.sets.sets.Union
        def _flatten(arg):
            if isinstance(arg, SeqBase):
                if isinstance(arg, SeqAdd):
                    return sum(map(_flatten, arg.args), [])
                else:
                    return [arg]
            if iterable(arg):
                return sum(map(_flatten, arg), [])
            raise TypeError("Input must be Sequences or "
                            " iterables of Sequences")
        args = _flatten(args)

        args = [a for a in args if a is not S.EmptySequence]

        # Addition of no sequences is EmptySequence
        if not args:
            return S.EmptySequence

        if Intersection(a.interval for a in args) is S.EmptySet:
            return S.EmptySequence

        # reduce using known rules
        if evaluate:
            return SeqAdd.reduce(args)

        args = list(ordered(args, SeqBase._start_key))

        return Basic.__new__(cls, *args)

    @staticmethod
    def reduce(args):
        """
        Simplify a :class:`SeqAdd` using known rules

        Then we iterate through all pairs and ask the constituent
        sequences if they can simplify themselves with any other constituent

        Notes
        =====

        adapted from ``Union.reduce``

        """
        new_args = True
        while(new_args):
            for id1, s in enumerate(args):
                new_args = False
                for id2, t in enumerate(args):
                    if id1 == id2:
                        continue
                    new_seq = s._add(t)
                    # This returns None if s does not know how to add
                    # with t. Returns the newly added sequence otherwise
                    if new_seq is not None:
                        new_args = [a for a in args if a not in (s, t)]
                        new_args.append(new_seq)
                        break
                if new_args:
                    args = new_args
                    break

        if len(args) == 1:
            return args.pop()
        else:
            return SeqAdd(args, evaluate=False)

    def _eval_coeff(self, pt):
        """adds up the coefficients of all the sequences at point pt"""
        return sum(a.coeff(pt) for a in self.args)


class SeqMul(SeqExprOp):
    """
    Represents term-wise multiplication of sequences.
    Handles multiplication of sequences only. For multiplication
    with other objects see ``SeqBase.coeff_mul``.

    Rules:
        * The interval on which sequence is defined is the intersection
        of respective intervals of sequences.
        * Anything * EmptySequence returns EmptySequence
        * Other rules are defined in _mul methods of sequence classes

    Examples
    ========

    >>> from sympy import S, oo, SeqMul, SeqPer, SeqFormula
    >>> from sympy.abc import n
    >>> SeqMul(SeqPer((1, 2), (n, 0, oo)), S.EmptySequence)
    EmptySequence()
    >>> SeqMul(SeqPer((1, 2), (n, 0, 5)), SeqPer((1, 2), (n, 6, 10)))
    EmptySequence()
    >>> SeqMul(SeqPer((1, 2), (n, 0, oo)), SeqFormula(n**2))
    SeqMul(SeqFormula(n**2, (n, 0, oo)), SeqPer((1, 2), (n, 0, oo)))
    >>> SeqMul(SeqFormula(n**3), SeqFormula(n**2))
    SeqFormula(n**5, (n, 0, oo))

    See Also
    ========

    sympy.series.sequences.SeqAdd
    """

    def __new__(cls, *args, **kwargs):
        evaluate = kwargs.get('evaluate', global_evaluate[0])

        # flatten inputs
        args = list(args)

        # adapted from sympy.sets.sets.Union
        def _flatten(arg):
            if isinstance(arg, SeqBase):
                if isinstance(arg, SeqMul):
                    return sum(map(_flatten, arg.args), [])
                else:
                    return [arg]
            elif iterable(arg):
                return sum(map(_flatten, arg), [])
            raise TypeError("Input must be Sequences or "
                            " iterables of Sequences")
        args = _flatten(args)

        # Multiplication of no sequences is EmptySequence
        if not args:
            return S.EmptySequence

        if Intersection(a.interval for a in args) is S.EmptySet:
            return S.EmptySequence

        # reduce using known rules
        if evaluate:
            return SeqMul.reduce(args)

        args = list(ordered(args, SeqBase._start_key))

        return Basic.__new__(cls, *args)

    @staticmethod
    def reduce(args):
        """
        Simplify a :class:`SeqMul` using known rules

        Then we iterate through all pairs and ask the constituent
        sequences if they can simplify themselves with any other constituent

        Notes
        =====
        adapted from ``Union.reduce``

        """
        new_args = True
        while(new_args):
            for id1, s in enumerate(args):
                new_args = False
                for id2, t in enumerate(args):
                    if id1 == id2:
                        continue
                    new_seq = s._mul(t)
                    # This returns None if s does not know how to multiply
                    # with t. Returns the newly multiplied sequence otherwise
                    if new_seq is not None:
                        new_args = [a for a in args if a not in (s, t)]
                        new_args.append(new_seq)
                        break
                if new_args:
                    args = new_args
                    break

        if len(args) == 1:
            return args.pop()
        else:
            return SeqMul(args, evaluate=False)

    def _eval_coeff(self, pt):
        """multiplies the coefficients of all the sequences at point pt"""
        val = 1
        for a in self.args:
            val *= a.coeff(pt)
        return val
