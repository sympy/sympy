from __future__ import print_function, division

from sympy.core.expr import Expr
from sympy.core.singleton import (S, Singleton)
from sympy.core.symbol import Dummy
from sympy.core.function import Lambda
from sympy.core.add import Add
from sympy.core.compatibility import (range, integer_types, with_metaclass\
                                      , is_sequence, iterable, ordered)
from sympy.core.sympify import sympify
from sympy.core.containers import Tuple
from sympy.core.evaluate import global_evaluate
from sympy.functions.elementary.integers import ceiling
from sympy.functions.elementary.miscellaneous import Min, Max
from sympy.polys import lcm
from sympy.sets.sets import Interval, Set, Intersection
from sympy.utilities.iterables import flatten


class SeqBase(Expr):
    """Base class for sequences"""

    is_iterable = True

    is_EmptySequence = False
    is_SeqPer = False
    is_SeqFunc = False
    is_SeqFormula = False
    is_SeqAdd = False

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

    @property
    def gen(self):
        """Returns the generator for the sequence"""
        raise NotImplementedError("(%s)._gen" % self)

    @property
    def interval(self):
        """The interval on which the sequence is defined"""
        raise NotImplementedError("(%s)._interval" % self)

    @property
    def start(self):
        """The starting point of the sequence. This point is included"""
        raise NotImplementedError("(%s)._start" % self)

    @property
    def stop(self):
        """The ending point of the sequence. This point is included"""
        raise NotImplementedError("(%s)._stop" % self)

    @property
    def step(self):
        """Increase points by step"""
        raise NotImplementedError("(%s)._step" % self)

    @property
    def length(self):
        """Length of the sequence"""
        raise NotImplementedError("(%s)._length" % self)

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
        >>> SeqFormula((m*n**2, n), (0, 5)).free_symbols
        set([m])
        """
        fsyms = set().union(*[a.free_symbols for a in self.args])
        for d in self.variables:
            if d in fsyms:
                fsyms.remove(d)
        return fsyms

    def coeff(self, pt):
        """Returns the coefficient at point pt"""
        if pt < self.start or pt > self.stop:
            raise IndexError("Index %s out of bounds %s" %(i, self.interval))
        return self._eval_coeff(pt)

    def _eval_coeff(self, pt):
        raise NotImplementedError("The _eval_coeff method should be added to"
                                  " %s to return coefficient so it is available"
                                  " when coeff calls it."
                                  % self.func)

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

    def add(self, other):
        """
        Returns the term-wise addition of 'self' and 'other'.

        Examples
        ========

        As a shortcut it is possible to use the '+' operator:

        >>> from sympy import S, oo, SeqAdd, SeqFormula
        >>> from sympy.abc import n
        >>> SeqAdd(SeqFormula(n**3), SeqFormula(n**2), evaluate=False)
        SeqFormula((n**2, n), ([0, oo), 1)) + SeqFormula((n**3, n), ([0, oo), 1))

        Similarly it is possible to use
        the '-' operator for term-wise subtraction:

        TODO (__neg__ needs to be implemented, will use SeqMul maybe)
        """
        return SeqAdd(self, other)

    def _add(self, other):
        """
        Should only be used internally

        self._add(other) returns a new, term-wise added sequence if self
        knows how to add with other, otherwise it returns ``None``.

        Used within :class:`SeqAdd` class
        """
        return None

    def __add__(self, other):
        return self.add(other)

    def __iter__(self):
        i = 0
        while(i < self.length):
            pt = self._ith_point(i)
            yield self.coeff(pt)
            i += 1

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

    >>> from sympy import S, SeqPer, oo
    >>> S.EmptySequence
    EmptySequence()

    >>> SeqPer((1, 2)).add(S.EmptySequence)
    SeqPer((1, 2), ([0, oo), 1))
    """

    is_EmptySequence = True

    @property
    def interval(self):
        return S.EmptySet

    @property
    def length(self):
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
            raise ValueError("Both the start and end value"
                             " cannot be unbounded")

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

    See Also
    ========

    sympy.series.sequences.SeqPer
    sympy.series.sequences.SeqFormula
    sympy.series.sequences.SeqFunc
    """

    def __new__(cls, gen, interval=(0, S.Infinity, 1)):
        bounds, step = _parse_interval(interval)
        if bounds is S.EmptySet:
            return S.EmptySequence
        gen = cls._validate(gen)
        interval = Tuple(bounds, step)
        return Expr.__new__(cls, gen, interval)

    @staticmethod
    def _validate(gen):
        """Validates the generator
        Should return unnchanged if no change required
        """
        return sympify(gen)

    @property
    def gen(self):
        return self.args[0]

    @property
    def interval(self):
        return self.args[1][0]

    @property
    def start(self):
        return self.interval.inf

    @property
    def stop(self):
        return self.interval.sup

    @property
    def step(self):
        return self.args[1][1]

    @property
    def length(self):
        return ceiling((self.stop - self.start + 1) / self.step)


class SeqPer(SeqExpr):
    """Represents a periodic sequence

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

    is_SeqPer = True

    @staticmethod
    def _validate(periodical):
        """
        Examples
        ========

        >>> from sympy.series.sequences import SeqPer as seq
        >>> seq._validate((1, 2, 3))
        (1, 2, 3)
        >>> seq._validate([1, 2, 3])
        (1, 2, 3)
        >>> seq._validate([1, (2, 3)])
        (1, 2, 3)
        """
        if is_sequence(periodical, Tuple):
            periodical = tuple(flatten(periodical))
            return sympify(periodical)
        raise ValueError("invalid period %s should be something like e.g (1, 2)"
            % periodical)

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
        return self.periodical[idx]

    def _add(self, other):
        """See docstring of SeqBase._add"""
        if other.is_SeqPer:
            start = Max(self.start, other.start)
            stop = Min(self.stop, other.stop)
            step = Max(self.step, other.step)

            per1, lper1 = self.periodical, self.period
            per2, lper2 = other.periodical, other.period

            per_length = lcm(lper1, lper2)

            new_per = []
            for x in range(per_length):
                ele1 = per1[x % lper1]
                ele2 = per2[x % lper2]
                new_per.append(ele1 + ele2)

            return SeqPer(new_per, (start, stop, step))


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

    >>> SeqFormula(n**2, (0, 5, 2))[:]
    [0, 4, 16]

    sequence starts from negative infinity

    >>> SeqFormula(n**2, (-oo, 0))[0:6]
    [0, 1, 4, 9, 16, 25]

    See Also
    ========

    sympy.series.sequences.SeqPer
    sympy.series.sequences.SeqFunc
    """

    is_SeqFormula = True

    @staticmethod
    def _validate(formula):
        """
        Examples
        ========

        >>> from sympy.series.sequences import SeqFormula as seq
        >>> from sympy.abc import n, m
        >>> seq._validate((n**2, n))
        (n**2, n)
        >>> seq._validate(n**2)
        (n**2, n)
        >>> seq._validate((m*n**2, n))
        (m*n**2, n)
        >>> seq._validate(1)
        (1, _k)
        """
        formula = sympify(formula)
        if not is_sequence(formula, Tuple):
            free = formula.free_symbols
            if len(free) == 0:
                k = Dummy('k')
                formula = Tuple(formula, k)
            elif len(free) == 1:
                formula = Tuple(formula, free.pop())
            else:
                raise ValueError(
                    " specify dummy variables for %s. If the formula contains"
                    " more than one free symbol, a dummy variable should be"
                    " supplied explicitly e.g., SeqFormula((m*n**2, n), (0, 5))"
                    % formula)
        return formula

    @property
    def formula(self):
        return self.gen[0]

    @property
    def variables(self):
        """Bounded variables
        In the case of SeqFormula
        it will be a tuple with only one symbol
        """
        return (self.gen[1],)

    def _eval_coeff(self, pt):
        d = self.variables[0]
        return self.formula.subs(d, pt)

    def _add(self, other):
        """See docstring of SeqBase._add"""
        if other.is_SeqFormula:
            start = Max(self.start, other.start)
            stop = Min(self.stop, other.stop)
            step = Max(self.step, other.step)

            form1, v1 = self.formula, self.variables[0]
            form2, v2 = other.formula, other.variables[0]
            formula = (form1 + form2.subs(v2, v1), v1)

            return SeqFormula(formula, (start, stop, step))


class SeqFunc(SeqExpr):
    """Represents sequence based on a function

    Elements are generated by calling a function.
    Only single argument functions are allowed.
    Only SymPy's :class:`Lambda` functions are permitted.

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
    sympy.core.function.Lambda
    """

    is_SeqFunc = True

    @staticmethod
    def _validate(function):
        if len(function.variables) != 1:
            raise ValueError("Only single argument functions are allowed")
        return sympify(function)

    @property
    def function(self):
        return self.gen

    @property
    def variables(self):
        return self.gen.variables

    def _eval_coeff(self, pt):
        return self.function(pt)

    def _add(self, other):
        """See docstring of SeqBase._add"""
        if other.is_SeqFunc:
            start = Max(self.start, other.start)
            stop = Min(self.stop, other.stop)
            step = Max(self.step, other.step)

            func1, v1 = self.function, self.variables[0]
            func2, v2 = other.function, other.variables[0]
            function = Lambda(v1, func1.expr + func2.expr.subs(v2, v1))

            return SeqFunc(function, (start, stop, step))


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


class SeqExprOp(SeqBase):
    """Base class for operations on sequences

    Examples
    ========

    >>> from sympy.series.sequences import SeqExprOp, sequence
    >>> from sympy.abc import n
    >>> s1 = sequence(formula=(n**2, n), interval=(0, 10))
    >>> s2 = sequence(periodical=(1, 2, 3), interval=(5, 10))
    >>> s = SeqExprOp(s1, s2)
    >>> s.gen
    ((n**2, n), (1, 2, 3))
    >>> s.interval
    [5, 10]
    >>> s.length
    6

    See Also
    ========

    sympy.series.sequences.SeqAdd
    """
    @property
    def gen(self):
        """Generator for the sequence
        returns a tuple of generators of all the argument sequences
        """
        return tuple([a.gen for a in self.args])

    @property
    def interval(self):
        """Sequence is defined on the intersection
        of all the intervals of respective sequences
        """
        return Intersection(*[a.interval for a in self.args])

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
    def step(self):
        """By default maximum step size is taken as the step size"""
        return max(a.step for a in self.args)

    @property
    def length(self):
        return ceiling((self.stop - self.start + 1) / self.step)


class SeqAdd(Add, SeqExprOp):
    """
    Represents addition of sequences

    Rules:
        * The interval on which sequence is defined is the intersection
        of respective intervals of sequences.
        * Anything + EmptySequence, remains unchanged
        * Other rules are defined in _add methods of sequence classes

    Examples
    ========

    >>> from sympy import S, oo, SeqAdd, SeqPer, SeqFormula
    >>> from sympy.abc import n
    >>> SeqAdd(SeqPer((1, 2)), S.EmptySequence)
    SeqPer((1, 2), ([0, oo), 1))
    >>> SeqAdd(SeqPer((1, 2), (0, 5)), SeqPer((1, 2), (6, 10)))
    EmptySequence()
    >>> SeqAdd(SeqPer((1, 2)), SeqFormula(n**2))
    SeqFormula((n**2, n), ([0, oo), 1)) + SeqPer((1, 2), ([0, oo), 1))
    """

    is_SeqAdd = True

    def __new__(cls, *args, **kwargs):
        evaluate = kwargs.get('evaluate', global_evaluate[0])

        # flatten inputs
        args = list(args)

        # adapted from sympy.sets.sets.Union
        def _flatten(arg):
            if isinstance(arg, SeqBase):
                if arg.is_SeqAdd:
                    return sum(map(_flatten, arg.args), [])
                else:
                    return [arg]
            if iterable(arg):
                return sum(map(_flatten, arg), [])
            raise TypeError("Input must be Sequences or "
                            " iterables of Sequences")
        args = _flatten(args)

        args = [a for a in args if not a is S.EmptySequence]

        # Addition of no sequences is EmptySequence
        if len(args) == 0:
            return S.EmptySequence

        if Intersection(*[a.interval for a in args]) is S.EmptySet:
            return S.EmptySequence

        # reduce using known rules
        if evaluate:
            return SeqAdd.reduce(args)

        args = list(ordered(args, SeqBase._start_key))

        return Expr.__new__(cls, *args)

    @staticmethod
    def reduce(args):
        """
        Simplify a :class:`SeqAdd` using known rules

        Then we iterate through all pairs and ask the constituent sets if they
        can simplify themselves with any other constituent

        Notes
        =====
        adapted from ``Union.reduce``

        """
        new_args = True
        while(new_args):
            for s in args:
                new_args = False
                for t in args:
                    if t == s:
                        continue
                    new_seq = s._add(t)
                    # This returns None if s does not know how to add
                    # with t. Returns the newly added sequence otherwise
                    if new_seq is not None:
                        new_args = [a for a in args if not a in (s, t)]
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
        val = 0
        for a in self.args:
            val += a.coeff(pt)
        return val
