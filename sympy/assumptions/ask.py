"""Module for querying SymPy objects about assumptions."""
from __future__ import print_function, division

from sympy.core import sympify

from sympy.logic.boolalg import (to_cnf, And, Not, Or, Implies, Equivalent,
    BooleanFunction, BooleanAtom)
from sympy.logic.inference import satisfiable
from sympy.assumptions.assume import (global_assumptions, Predicate,
        AppliedPredicate)
from sympy.core.decorators import deprecated


# Deprecated predicates should be added to this list
deprecated_predicates = [
    'bounded',
    'infinity'
]


class QClass(object):
    """
    This class contains all the supported keys by ``ask``.
    """

    @property
    def hermitian(self):
        return Predicate('Hermitian')

    @property
    def antihermitian(self):
        r"""
        Antihermitian predicate.

        ``ask(Q.antihermitian(x))`` is true iff x belongs to the field of
        antihermitian operators.

        TODO: Add examples
        """
        return Predicate('antihermitian')

    @property
    def real(self):
        r"""
        Real number predicate.

        ``Q.real(x)`` is true iff ``x`` is a real number, i.e., it is in the
        interval `(-\infty, \infty)`.  Note that, in particular the infinities are
        not real. Use ``Q.extended_real`` if you want to consider those as well.

        A few important facts about reals:

        - Every real number is positive, negative, or zero.  Furthermore, because
          these sets are pairwise disjoint, each real number is exactly one of
          those three.

        - Every real number is also complex.

        - Every real number is either rational or irrational.

        - Every real number is either algebraic or transcendental.

        - The facts ``Q.negative``, ``Q.zero``, ``Q.positive``, ``Q.nonnegative``,
          ``Q.nonpositive``, ``Q.nonzero``, ``Q.integer``, ``Q.rational``, and
          ``Q.irrational`` all imply ``Q.real``, as do all facts that imply those
          facts.

        - The facts ``Q.algebraic``, and ``Q.transcendental`` do not imply
          ``Q.real``; they imply ``Q.complex``. An algebraic or transcendental
          number may or may not be real.

        - The "non" facts (i.e., ``Q.nonnegative``, ``Q.nonzero``, and
          ``Q.nonpositive``) are not equivalent to not the fact, but rather, not
          the fact *and* ``Q.real``.  For example, ``Q.nonnegative`` means
          ``~Q.negative & Q.real``. So for example, ``I`` is not nonnegative,
          nonzero, or nonpositive.


        Examples
        ========

        >>> from sympy import Q, ask, symbols
        >>> x = symbols('x')
        >>> ask(Q.real(x), Q.positive(x))
        True
        >>> ask(Q.real(0))
        True

        """
        return Predicate('real')

    @property
    def extended_real(self):
        return Predicate('extended_real')

    @property
    def imaginary(self):
        r"""
        Imaginary number predicate.

        ``ask(Q.imaginary(x))`` is true iff ``x`` can be written as a real
        number multiplied by the imaginary unit ``I``. Please note that ``0``
        is not considered to be an imaginary number.


        Examples
        ========

        >>> from sympy import Q, ask, I
        >>> ask(Q.imaginary(3*I))
        True
        >>> ask(Q.imaginary(2 + 3*I))
        False
        >>> ask(Q.imaginary(0))
        False

        """
        return Predicate('imaginary')

    @property
    def complex(self):
        r"""
        Complex number predicate.

        ``ask(Q.complex(x))`` is true iff ``x`` belongs to the set of complex
        numbers.


        Examples
        ========

        >>> from sympy import Q, Symbol, ask, I, oo
        >>> x = Symbol('x')
        >>> ask(Q.complex(0))
        True
        >>> ask(Q.complex(2 + 3*I))
        True
        >>> ask(Q.complex(oo))
        False

        """
        return Predicate('complex')

    @property
    def algebraic(self):
        return Predicate('algebraic')

    @property
    def integer(self):
        r"""
        Integer predicate.

        ``ask(Q.integer(x))`` is true iff ``x`` belongs to the set of integer numbers.


        Examples
        ========

        >>> from sympy import Q, ask, S
        >>> ask(Q.integer(5))
        True
        >>> ask(Q.integer(S(1)/2))
        False

        """
        return Predicate('integer')

    @property
    def rational(self):
        return Predicate('rational')

    @property
    def irrational(self):
        return Predicate('irrational')

    @property
    def finite(self):
        r"""
        Finite predicate.

        ``ask(Q.finite(x))`` will return ``True`` if ``x`` is neither an infinity
        nor a ``NaN``. In other words, ``ask(Q.finite(x))`` is true for all ``x``
        having a bounded absolute value.

        Examples
        ========

        >>> from sympy import Q, ask, Symbol, S, oo
        >>> x = Symbol('x')
        >>> ask(Q.finite(S.NaN))
        False
        >>> ask(Q.finite(oo))
        False
        >>> ask(Q.finite(1))
        True

        """
        return Predicate('finite')

    @property
    @deprecated(useinstead="finite", issue=9425, deprecated_since_version="0.7.7")
    def bounded(self):
        return Predicate('finite')

    @property
    def infinite(self):
        return Predicate('infinite')

    @property
    @deprecated(useinstead="infinite", issue=9426, deprecated_since_version="0.7.7")
    def infinity(self):
        return Predicate('infinite')

    @property
    def infinitesimal(self):
        return Predicate('infinitesimal')

    @property
    def positive(self):
        r"""
        Positive real number predicate.

        ``Q.positive(x)`` is true iff ``x`` is real and `x > 0`, that is if ``x``
        is in the interval `(0, \infty)`.  In particular, infinity is not
        positive.

        A few important facts about positive numbers:

        - Note that ``Q.nonpositive`` and ``~Q.positive`` are *not* the same
          thing. ``~Q.positive(x)`` simply means that ``x`` is not positive,
          whereas ``Q.nonpositive(x)`` means that ``x`` is real and not
          positive, i.e., ``Q.positive(x)`` is logically equivalent to
          `Q.negative(x) | Q.zero(x)``.  So for example, ``~Q.positive(I)`` is
          true, whereas ``Q.nonpositive(I)`` is false.

        - See the docstring of ``Q.real`` for more information about related facts.


        Examples
        ========

        >>> from sympy import Q, ask, symbols, I
        >>> x = symbols('x')
        >>> ask(Q.positive(x), Q.real(x) & ~Q.negative(x) & ~Q.zero(x))
        True
        >>> ask(Q.positive(1))
        True
        >>> ask(Q.nonpositive(I))
        False
        >>> ask(~Q.positive(I))
        True

        """
        return Predicate('positive')

    @property
    def negative(self):
        r"""
        Negative number predicate.

        ``Q.negative(x)`` is true iff ``x`` is a real number and `x < 0`, that is,
        it is in the interval `(-oo, 0)`.  Note in particular that negative
        infinity is not negative.

        A few important facts about negative numbers:

        - Note that ``Q.nonnegative`` and ``~Q.negative`` are *not* the same
          thing. ``~Q.negative(x)`` simply means that ``x`` is not negative,
          whereas ``Q.nonnegative(x)`` means that ``x`` is real and not
          negative, i.e., ``Q.nonnegative(x)`` is logically equivalent to
          ``Q.zero(x) | Q.positive(x)``.  So for example, ``~Q.negative(I)`` is
          true, whereas ``Q.nonnegative(I)`` is false.

        - See the docstring of ``Q.real`` for more information about related facts.


        Examples
        ========

        >>> from sympy import Q, ask, symbols, I
        >>> x = symbols('x')
        >>> ask(Q.negative(x), Q.real(x) & ~Q.positive(x) & ~Q.zero(x))
        True
        >>> ask(Q.negative(-1))
        True
        >>> ask(Q.nonnegative(I))
        False
        >>> ask(~Q.negative(I))
        True

        """
        return Predicate('negative')

    @property
    def zero(self):
        return Predicate('zero')

    @property
    def nonzero(self):
        r"""
        Nonzero real number predicate.

        ``Q.nonzero(x)`` is true iff ``x`` is real and ``x`` is not zero.  Note in
        particular that ``Q.nonzero(x)`` is false if ``x`` is not real.  Use
        ``~Q.zero(x)`` if you want the negation of being zero without any real
        assumptions.

        A few important facts about nonzero numbers:

        - ``Q.nonzero`` is logically equivalent to ``Q.positive | Q.negative``.

        - See the docstring of ``Q.real`` for more information about related facts.


        Examples
        ========

        >>> from sympy import Q, ask, symbols, I
        >>> x = symbols('x')
        >>> print(ask(Q.nonzero(x), ~Q.zero(x)))
        None
        >>> ask(Q.nonzero(x), Q.positive(x))
        True
        >>> ask(Q.nonzero(x), Q.zero(x))
        False
        >>> ask(Q.nonzero(0))
        False
        >>> ask(Q.nonzero(I))
        False
        >>> ask(~Q.zero(I))
        True

        """
        return Predicate('nonzero')

    @property
    def nonpositive(self):
        return Predicate('nonpositive')

    @property
    def nonnegative(self):
        return Predicate('nonnegative')

    @property
    def even(self):
        r"""
        Even number predicate.

        ``ask(Q.even(x))`` is true iff ``x`` belongs to the set of even numbers.


        Examples
        ========

        >>> from sympy import Q, ask, pi
        >>> ask(Q.even(0))
        True
        >>> ask(Q.even(2))
        True
        >>> ask(Q.even(3))
        False
        >>> ask(Q.even(pi))
        False

        """
        return Predicate('even')

    @property
    def odd(self):
        return Predicate('odd')

    @property
    def prime(self):
        return Predicate('prime')

    @property
    def composite(self):
        r"""
        Composite number predicate.

        ``ask(Q.composite(x))`` is true iff ``x`` is a positive integer and has
        at least one positive divisor other than ``1`` and the number itself.


        Examples
        ========

        >>> from sympy import Q, ask
        >>> ask(Q.composite(0))
        False
        >>> ask(Q.composite(1))
        True
        >>> ask(Q.composite(2))
        False
        >>> ask(Q.composite(20))
        True

        """
        return Predicate('composite')

    @property
    def commutative(self):
        r"""
        Commutative predicate.

        ``ask(Q.commutative(x))`` is true iff ``x`` commutes with any other
        object with respect to multiplication operation.

        TODO: Add examples
        """
        return Predicate('commutative')

    @property
    def is_true(self):
        r"""
        Generic predicate.

        ``Q.is_true(x)`` is true iff ``x`` is true. This only makes sense if ``x`` is a
        predicate.

        Examples
        ========

        >>> from sympy import ask, Q, symbols
        >>> x = symbols('x')
        >>> ask(Q.is_true(True))
        True

        """
        return Predicate('is_true')

    @property
    def symmetric(self):
        return Predicate('symmetric')

    @property
    def invertible(self):
        return Predicate('invertible')

    @property
    def orthogonal(self):
        return Predicate('orthogonal')

    @property
    def unitary(self):
        return Predicate('unitary')

    @property
    def positive_definite(self):
        return Predicate('positive_definite')

    @property
    def upper_triangular(self):
        return Predicate('upper_triangular')

    @property
    def lower_triangular(self):
        return Predicate('lower_triangular')

    @property
    def diagonal(self):
        return Predicate('diagonal')

    @property
    def fullrank(self):
        return Predicate('fullrank')

    @property
    def square(self):
        return Predicate('square')

    @property
    def integer_elements(self):
        return Predicate('integer_elements')

    @property
    def real_elements(self):
        return Predicate('real_elements')

    @property
    def complex_elements(self):
        return Predicate('complex_elements')

    @property
    def transcendental(self):
        return Predicate('transcendental')

    @property
    def singular(self):
        return Predicate('singular')

    @property
    def normal(self):
        return Predicate('normal')

    @property
    def triangular(self):
        return Predicate('triangular')

    @property
    def unit_triangular(self):
        return Predicate('unit_triangular')


Q = QClass()


def _extract_facts(expr, symbol):
    """
    Helper for ask().

    Extracts the facts relevant to the symbol from an assumption.
    Returns None if there is nothing to extract.
    """
    if isinstance(expr, bool):
        return
    if not expr.has(symbol):
        return
    if isinstance(expr, AppliedPredicate):
        if expr.arg == symbol:
            return expr.func
        else:
            return
    if isinstance(expr, Not) and expr.args[0].func in (And, Or):
        cls = Or if expr.args[0] == And else And
        expr = cls(*[~arg for arg in expr.args[0].args])
    args = [_extract_facts(arg, symbol) for arg in expr.args]
    if isinstance(expr, And):
        args = [x for x in args if x is not None]
        if args:
            return expr.func(*args)
    if args and all(x != None for x in args):
        return expr.func(*args)


def ask(proposition, assumptions=True, context=global_assumptions):
    """
    Method for inferring properties about objects.

    **Syntax**

        * ask(proposition)

        * ask(proposition, assumptions)

            where ``proposition`` is any boolean expression

    Examples
    ========

    >>> from sympy import ask, Q, pi
    >>> from sympy.abc import x, y
    >>> ask(Q.rational(pi))
    False
    >>> ask(Q.even(x*y), Q.even(x) & Q.integer(y))
    True
    >>> ask(Q.prime(x*y), Q.integer(x) &  Q.integer(y))
    False

    **Remarks**
        Relations in assumptions are not implemented (yet), so the following
        will not give a meaningful result.

        >>> ask(Q.positive(x), Q.is_true(x > 0)) # doctest: +SKIP

        It is however a work in progress.

    """
    from sympy.assumptions.satask import satask

    if not isinstance(proposition, (BooleanFunction, AppliedPredicate, bool, BooleanAtom)):
        raise TypeError("proposition must be a valid logical expression")

    if not isinstance(assumptions, (BooleanFunction, AppliedPredicate, bool, BooleanAtom)):
        raise TypeError("assumptions must be a valid logical expression")

    if isinstance(proposition, AppliedPredicate):
        key, expr = proposition.func, sympify(proposition.arg)
    else:
        key, expr = Q.is_true, sympify(proposition)

    assumptions = And(assumptions, And(*context))
    assumptions = to_cnf(assumptions)

    local_facts = _extract_facts(assumptions, expr)

    if local_facts and satisfiable(And(local_facts, known_facts_cnf)) is False:
        raise ValueError("inconsistent assumptions %s" % assumptions)

    # direct resolution method, no logic
    res = key(expr)._eval_ask(assumptions)
    if res is not None:
        return bool(res)

    if assumptions == True:
        return

    if local_facts is None:
        return satask(proposition, assumptions=assumptions, context=context)


    # See if there's a straight-forward conclusion we can make for the inference
    if local_facts.is_Atom:
        if key in known_facts_dict[local_facts]:
            return True
        if Not(key) in known_facts_dict[local_facts]:
            return False
    elif local_facts.func is And and all(k in known_facts_dict for k in local_facts.args):
        for assum in local_facts.args:
            if assum.is_Atom:
                if key in known_facts_dict[assum]:
                    return True
                if Not(key) in known_facts_dict[assum]:
                    return False
            elif assum.func is Not and assum.args[0].is_Atom:
                if key in known_facts_dict[assum]:
                    return False
                if Not(key) in known_facts_dict[assum]:
                    return True
    elif (isinstance(key, Predicate) and
            local_facts.func is Not and local_facts.args[0].is_Atom):
        if local_facts.args[0] in known_facts_dict[key]:
            return False

    # Failing all else, we do a full logical inference
    res = ask_full_inference(key, local_facts, known_facts_cnf)
    if res is None:
        return satask(proposition, assumptions=assumptions, context=context)
    return res


def ask_full_inference(proposition, assumptions, known_facts_cnf):
    """
    Method for inferring properties about objects.

    """
    if not satisfiable(And(known_facts_cnf, assumptions, proposition)):
        return False
    if not satisfiable(And(known_facts_cnf, assumptions, Not(proposition))):
        return True
    return None


def register_handler(key, handler):
    """
    Register a handler in the ask system. key must be a string and handler a
    class inheriting from AskHandler::

        >>> from sympy.assumptions import register_handler, ask, Q
        >>> from sympy.assumptions.handlers import AskHandler
        >>> class MersenneHandler(AskHandler):
        ...     # Mersenne numbers are in the form 2**n + 1, n integer
        ...     @staticmethod
        ...     def Integer(expr, assumptions):
        ...         import math
        ...         return ask(Q.integer(math.log(expr + 1, 2)))
        >>> register_handler('mersenne', MersenneHandler)
        >>> ask(Q.mersenne(7))
        True

    """
    if type(key) is Predicate:
        key = key.name
    try:
        getattr(Q.__class__, key).add_handler(handler)
    except AttributeError:
        setattr(Q.__class__, key, Predicate(key, handlers=[handler]))


def remove_handler(key, handler):
    """Removes a handler from the ask system. Same syntax as register_handler"""
    if type(key) is Predicate:
        key = key.name
    getattr(Q.__class__, key).remove_handler(handler)


def single_fact_lookup(known_facts_keys, known_facts_cnf):
    # Compute the quick lookup for single facts
    mapping = {}
    for key in known_facts_keys:
        mapping[key] = set([key])
        for other_key in known_facts_keys:
            if other_key != key:
                if ask_full_inference(other_key, key, known_facts_cnf):
                    mapping[key].add(other_key)
    return mapping


def compute_known_facts(known_facts, known_facts_keys):
    """Compute the various forms of knowledge compilation used by the
    assumptions system.

    This function is typically applied to the variables
    ``known_facts`` and ``known_facts_keys`` defined at the bottom of
    this file.
    """
    from textwrap import dedent, wrap

    fact_string = dedent('''\
    """
    The contents of this file are the return value of
    ``sympy.assumptions.ask.compute_known_facts``.  Do NOT manually
    edit this file.  Instead, run ./bin/ask_update.py.
    """

    from sympy.logic.boolalg import And, Not, Or
    from sympy.assumptions.ask import Q

    # -{ Known facts in CNF }-
    known_facts_cnf = And(
        %s
    )

    # -{ Known facts in compressed sets }-
    known_facts_dict = {
        %s
    }
    ''')
    # Compute the known facts in CNF form for logical inference
    LINE = ",\n    "
    HANG = ' '*8
    cnf = to_cnf(known_facts)
    c = LINE.join([str(a) for a in cnf.args])
    mapping = single_fact_lookup(known_facts_keys, cnf)
    items = sorted(mapping.items(), key=str)
    keys = [str(i[0]) for i in items]
    values = ['set(%s)' % sorted(i[1], key=str) for i in items]
    m = LINE.join(['\n'.join(
        wrap("%s: %s" % (k, v),
            subsequent_indent=HANG,
            break_long_words=False))
        for k, v in zip(keys, values)]) + ','
    return fact_string % (c, m)

# handlers tells us what ask handler we should use
# for a particular key
_val_template = 'sympy.assumptions.handlers.%s'
_handlers = [
    ("antihermitian",     "sets.AskAntiHermitianHandler"),
    ("finite",           "calculus.AskFiniteHandler"),
    ("commutative",       "AskCommutativeHandler"),
    ("complex",           "sets.AskComplexHandler"),
    ("composite",         "ntheory.AskCompositeHandler"),
    ("even",              "ntheory.AskEvenHandler"),
    ("extended_real",     "sets.AskExtendedRealHandler"),
    ("hermitian",         "sets.AskHermitianHandler"),
    ("imaginary",         "sets.AskImaginaryHandler"),
    ("infinitesimal",     "calculus.AskInfinitesimalHandler"),
    ("integer",           "sets.AskIntegerHandler"),
    ("irrational",        "sets.AskIrrationalHandler"),
    ("rational",          "sets.AskRationalHandler"),
    ("negative",          "order.AskNegativeHandler"),
    ("nonzero",           "order.AskNonZeroHandler"),
    ("nonpositive",       "order.AskNonPositiveHandler"),
    ("nonnegative",       "order.AskNonNegativeHandler"),
    ("zero",              "order.AskZeroHandler"),
    ("positive",          "order.AskPositiveHandler"),
    ("prime",             "ntheory.AskPrimeHandler"),
    ("real",              "sets.AskRealHandler"),
    ("odd",               "ntheory.AskOddHandler"),
    ("algebraic",         "sets.AskAlgebraicHandler"),
    ("is_true",           "common.TautologicalHandler"),
    ("symmetric",         "matrices.AskSymmetricHandler"),
    ("invertible",        "matrices.AskInvertibleHandler"),
    ("orthogonal",        "matrices.AskOrthogonalHandler"),
    ("unitary",           "matrices.AskUnitaryHandler"),
    ("positive_definite", "matrices.AskPositiveDefiniteHandler"),
    ("upper_triangular",  "matrices.AskUpperTriangularHandler"),
    ("lower_triangular",  "matrices.AskLowerTriangularHandler"),
    ("diagonal",          "matrices.AskDiagonalHandler"),
    ("fullrank",          "matrices.AskFullRankHandler"),
    ("square",            "matrices.AskSquareHandler"),
    ("integer_elements",  "matrices.AskIntegerElementsHandler"),
    ("real_elements",     "matrices.AskRealElementsHandler"),
    ("complex_elements",  "matrices.AskComplexElementsHandler"),
]

for name, value in _handlers:
    register_handler(name, _val_template % value)

known_facts_keys = [getattr(Q, attr) for attr in Q.__class__.__dict__
                    if not (attr.startswith('__') or attr in deprecated_predicates)]

known_facts = And(
    Implies(Q.infinite, ~Q.finite),
    Implies(Q.real, Q.complex),
    Implies(Q.real, Q.hermitian),
    Equivalent(Q.even, Q.integer & ~Q.odd),
    Equivalent(Q.extended_real, Q.real | Q.infinite),
    Equivalent(Q.odd, Q.integer & ~Q.even),
    Equivalent(Q.prime, Q.integer & Q.positive & ~Q.composite),
    Implies(Q.integer, Q.rational),
    Implies(Q.rational, Q.algebraic),
    Implies(Q.algebraic, Q.complex),
    Equivalent(Q.transcendental, Q.complex & ~Q.algebraic),
    Implies(Q.imaginary, Q.complex & ~Q.real),
    Implies(Q.imaginary, Q.antihermitian),
    Implies(Q.antihermitian, ~Q.hermitian),
    Equivalent(Q.negative, Q.nonzero & ~Q.positive),
    Equivalent(Q.positive, Q.nonzero & ~Q.negative),
    Equivalent(Q.rational, Q.real & ~Q.irrational),
    Equivalent(Q.real, Q.rational | Q.irrational),
    Implies(Q.nonzero, Q.real),
    Equivalent(Q.nonzero, Q.positive | Q.negative),
    Equivalent(Q.nonpositive, ~Q.positive & Q.real),
    Equivalent(Q.nonnegative, ~Q.negative & Q.real),
    Equivalent(Q.zero, Q.real & ~Q.nonzero),
    Implies(Q.zero, Q.even),

    Implies(Q.orthogonal, Q.positive_definite),
    Implies(Q.orthogonal, Q.unitary),
    Implies(Q.unitary & Q.real, Q.orthogonal),
    Implies(Q.unitary, Q.normal),
    Implies(Q.unitary, Q.invertible),
    Implies(Q.normal, Q.square),
    Implies(Q.diagonal, Q.normal),
    Implies(Q.positive_definite, Q.invertible),
    Implies(Q.diagonal, Q.upper_triangular),
    Implies(Q.diagonal, Q.lower_triangular),
    Implies(Q.lower_triangular, Q.triangular),
    Implies(Q.upper_triangular, Q.triangular),
    Implies(Q.triangular, Q.upper_triangular | Q.lower_triangular),
    Implies(Q.upper_triangular & Q.lower_triangular, Q.diagonal),
    Implies(Q.diagonal, Q.symmetric),
    Implies(Q.unit_triangular, Q.triangular),
    Implies(Q.invertible, Q.fullrank),
    Implies(Q.invertible, Q.square),
    Implies(Q.symmetric, Q.square),
    Implies(Q.fullrank & Q.square, Q.invertible),
    Equivalent(Q.invertible, ~Q.singular),
    Implies(Q.integer_elements, Q.real_elements),
    Implies(Q.real_elements, Q.complex_elements),
)

from sympy.assumptions.ask_generated import known_facts_dict, known_facts_cnf
