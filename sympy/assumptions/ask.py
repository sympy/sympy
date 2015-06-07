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
        """
        Hermitian predicate.

        ``ask(Q.hermitian(x))`` is true iff ``x`` belongs to the field of
        hermitian operators.

        TODO: Add examples

        """
        return Predicate('Hermitian')

    @property
    def antihermitian(self):
        """
        Antihermitian predicate.

        ``Q.antihermitian(x)`` is true iff ``x`` belongs to the field of
        antihermitian operators, i.e., operators in the form ``x*I``, where
        ``x`` is Hermitian.

        TODO: Add examples

        """
        return Predicate('antihermitian')

    @property
    def real(self):
        r"""
        Real number predicate.

        ``Q.real(x)`` is true iff ``x`` is a real number, i.e., it is in the
        interval :math:`(-\infty, \infty)`.  Note that, in particular the infinities
        are not real. Use ``Q.extended_real`` if you want to consider those as well.

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

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Real_number

        """
        return Predicate('real')

    @property
    def extended_real(self):
        r"""
        Extended real predicate.

        ``Q.extended_real(x)`` is true iff ``x`` is a real number including
        :math:`{-\infty, \infty}`.

        See documentation of ``Q.real`` for more information about related facts.

        Examples
        ========

        >>> from sympy import ask, Q, oo, I
        >>> ask(Q.extended_real(1))
        True
        >>> ask(Q.extended_real(I))
        False
        >>> ask(Q.extended_real(oo))
        True

        """
        return Predicate('extended_real')

    @property
    def imaginary(self):
        """
        Imaginary number predicate.

        ``Q.imaginary(x)`` is true iff ``x`` can be written as a real
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

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Imaginary_number

        """
        return Predicate('imaginary')

    @property
    def complex(self):
        """
        Complex number predicate.

        ``Q.complex(x)`` is true iff ``x`` belongs to the set of complex
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

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Complex_number

        """
        return Predicate('complex')

    @property
    def algebraic(self):
        """
        Algebraic number predicate.

        ``Q.algebraic(x)`` is true iff ``x`` belongs to the set of
        algebraic numbers.

        Examples
        ========

        >>> from sympy import ask, Q, sqrt, I, pi
        >>> ask(Q.algebraic(sqrt(2)))
        True
        >>> ask(Q.algebraic(I))
        True
        >>> ask(Q.algebraic(pi))
        False

        References
        ==========

        .. [1] http://en.wikipedia.org/wiki/Algebraic_number
        """
        return Predicate('algebraic')

    @property
    def transcendental(self):
        """
        Transcedental number predicate.

        ``Q.algebraic(x)`` is true iff ``x`` belongs to the set of
        transcendental numbers. A transcendental number is a real
        or complex number that is not algebraic.

        TODO: Add examples

        """
        return Predicate('transcendental')

    @property
    def integer(self):
        """
        Integer predicate.

        ``Q.integer(x)`` is true iff ``x`` belongs to the set of integer numbers.

        Examples
        ========

        >>> from sympy import Q, ask, S
        >>> ask(Q.integer(5))
        True
        >>> ask(Q.integer(S(1)/2))
        False

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Integer

        """
        return Predicate('integer')

    @property
    def rational(self):
        """
        Rational number predicate.

        ``Q.rational(x)`` is true iff ``x`` belongs to the set of
        rational numbers.

        Examples
        ========

        >>> from sympy import ask, Q, pi, S
        >>> ask(Q.rational(0))
        True
        >>> ask(Q.rational(S(1)/2))
        True
        >>> ask(Q.rational(pi))
        False

        References
        ==========

        https://en.wikipedia.org/wiki/Rational_number

        """
        return Predicate('rational')

    @property
    def irrational(self):
        """
        Irrational number predicate.

        ``Q.irrational(x)`` is true iff ``x``  is any real number that
        cannot be expressed as a ratio of integers.

        Examples
        ========

        >>> from sympy import ask, Q, pi, S
        >>> ask(Q.irrational(0))
        False
        >>> ask(Q.irrational(S(1)/2))
        False
        >>> ask(Q.irrational(pi))
        True

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Irrational_number

        """
        return Predicate('irrational')

    @property
    def finite(self):
        """
        Finite predicate.

        ``Q.finite(x)`` is true if ``x`` is neither an infinity
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

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Finite

        """
        return Predicate('finite')

    @property
    @deprecated(useinstead="finite", issue=9425, deprecated_since_version="0.7.7")
    def bounded(self):
        """
        See documentation of ``Q.finite``.
        """
        return Predicate('finite')

    @property
    def infinite(self):
        """
        Infinite number predicate.

        ``Q.infinite(x)`` is true iff the absolute value of ``x`` is
        arbitrarily large.

        """
        return Predicate('infinite')

    @property
    @deprecated(useinstead="infinite", issue=9426, deprecated_since_version="0.7.7")
    def infinity(self):
        """
        See documentation of ``Q.infinite``.
        """
        return Predicate('infinite')

    @property
    def infinitesimal(self):
        """
        See documentation of ``Q.zero``.
        """
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

        - See the documentation of ``Q.real`` for more information about
          related facts.

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

        ``Q.negative(x)`` is true iff ``x`` is a real number and :math:`x < 0`, that is,
        it is in the interval :math:`(-\infty, 0)`.  Note in particular that negative
        infinity is not negative.

        A few important facts about negative numbers:

        - Note that ``Q.nonnegative`` and ``~Q.negative`` are *not* the same
          thing. ``~Q.negative(x)`` simply means that ``x`` is not negative,
          whereas ``Q.nonnegative(x)`` means that ``x`` is real and not
          negative, i.e., ``Q.nonnegative(x)`` is logically equivalent to
          ``Q.zero(x) | Q.positive(x)``.  So for example, ``~Q.negative(I)`` is
          true, whereas ``Q.nonnegative(I)`` is false.

        - See the documentation of ``Q.real`` for more information about
          related facts.

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
        """
        Zero number predicate.

        ``ask(Q.zero(x))`` is true iff the value of ``x`` is zero.

        Examples
        ========

        >>> from sympy import ask, Q, oo
        >>> ask(Q.zero(0))
        True
        >>> ask(Q.zero(1/oo))
        True
        >>> ask(Q.zero(0*oo))
        False
        >>> ask(Q.zero(1))
        False

        """
        return Predicate('zero')

    @property
    def nonzero(self):
        """
        Nonzero real number predicate.

        ``ask(Q.nonzero(x))`` is true iff ``x`` is real and ``x`` is not zero.  Note in
        particular that ``Q.nonzero(x)`` is false if ``x`` is not real.  Use
        ``~Q.zero(x)`` if you want the negation of being zero without any real
        assumptions.

        A few important facts about nonzero numbers:

        - ``Q.nonzero`` is logically equivalent to ``Q.positive | Q.negative``.

        - See the documentation of ``Q.real`` for more information about
          related facts.

        Examples
        ========

        >>> from sympy import Q, ask, symbols, I
        >>> x = symbols('x')
        >>> ask(Q.nonzero(x), ~Q.zero(x))
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
        """
        Nonpositive real number predicate.

        ``ask(Q.nonpositive(x))`` is true iff ``x`` belongs to the set of
        negative numbers including zero.

        Examples
        ========

        >>> from sympy import Q, ask, I
        >>> ask(Q.nonpositive(-1))
        True
        >>> ask(Q.nonpositive(0))
        True
        >>> ask(Q.nonpositive(1))
        False
        >>> ask(Q.nonpositive(I))
        False
        >>> ask(Q.nonpositive(-I))
        False

        """
        return Predicate('nonpositive')

    @property
    def nonnegative(self):
        """
        Nonnegative real number predicate.

        ``ask(Q.nonnegative(x))`` is true iff ``x`` belongs to the set of
        positive numbers including zero.

        Examples
        ========

        >>> from sympy import Q, ask, I
        >>> ask(Q.nonnegative(1))
        True
        >>> ask(Q.nonnegative(0))
        True
        >>> ask(Q.nonnegative(-1))
        False
        >>> ask(Q.nonnegative(I))
        False
        >>> ask(Q.nonnegative(-I))
        False

        """
        return Predicate('nonnegative')

    @property
    def even(self):
        """
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
        """
        Odd number predicate.

        ``ask(Q.odd(x))`` is true iff ``x`` belongs to the set of odd numbers.

        Examples
        ========

        >>> from sympy import Q, ask, pi
        >>> ask(Q.odd(0))
        False
        >>> ask(Q.odd(2))
        False
        >>> ask(Q.odd(3))
        True
        >>> ask(Q.odd(pi))
        False

        """
        return Predicate('odd')

    @property
    def prime(self):
        """
        Prime number predicate.

        ``ask(Q.prime(x))`` is true iff ``x`` is a natural number greater
        than 1 that has no positive divisors other than ``1`` and the
        number itself.

        Examples
        ========

        >>> from sympy import Q, ask
        >>> ask(Q.prime(0))
        False
        >>> ask(Q.prime(1))
        False
        >>> ask(Q.prime(2))
        True
        >>> ask(Q.prime(20))
        False

        """
        return Predicate('prime')

    @property
    def composite(self):
        """
        Composite number predicate.

        ``ask(Q.composite(x))`` is true iff ``x`` is a positive integer and has
        at least one positive divisor other than ``1`` and the number itself.

        Examples
        ========

        >>> from sympy import Q, ask
        >>> ask(Q.composite(0))
        False
        >>> ask(Q.composite(1))
        False
        >>> ask(Q.composite(2))
        False
        >>> ask(Q.composite(20))
        True

        """
        return Predicate('composite')

    @property
    def commutative(self):
        """
        Commutative predicate.

        ``ask(Q.commutative(x))`` is true iff ``x`` commutes with any other
        object with respect to multiplication operation.

        TODO: Add examples

        """
        return Predicate('commutative')

    @property
    def is_true(self):
        """
        Generic predicate.

        ``ask(Q.is_true(x))`` is true iff ``x`` is true. This only makes
        sense if ``x`` is a predicate.

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
        """
        Symmetric matrix predicate.

        ``Q.symmetric(x)`` is true iff the square matrix ``x`` is equal to
        its transpose. Every square diagonal matrix is a symmetric matrix.

        Examples
        ========

        >>> from sympy import Q, ask, MatrixSymbol
        >>> X = MatrixSymbol('X', 2, 2)
        >>> Y = MatrixSymbol('Y', 2, 3)
        >>> Z = MatrixSymbol('Z', 2, 2)
        >>> ask(Q.symmetric(X*Z), Q.symmetric(X) & Q.symmetric(Z))
        True
        >>> ask(Q.symmetric(X + Z), Q.symmetric(X) & Q.symmetric(Z))
        True
        >>> ask(Q.symmetric(Y))
        False

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Symmetric_matrix

        """
        return Predicate('symmetric')

    @property
    def invertible(self):
        """
        Invertible matrix predicate.

        ``Q.invertible(x)`` is true iff ``x`` is an invertible matrix.

        Examples
        ========

        >>> from sympy import Q, ask, MatrixSymbol
        >>> X = MatrixSymbol('X', 2, 2)
        >>> Y = MatrixSymbol('Y', 2, 3)
        >>> Z = MatrixSymbol('Z', 2, 2)
        >>> ask(Q.invertible(X*Y), Q.invertible(X))
        False
        >>> ask(Q.invertible(X*Z), Q.invertible(X) & Q.invertible(Z))
        True
        >>> ask(Q.invertible(X), Q.fullrank(X) & Q.square(X))
        True

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Invertible_matrix

        """
        return Predicate('invertible')

    @property
    def orthogonal(self):
        """
        Orthogonal matrix predicate.

        ``Q.orthogonal(x)`` is true iff ``x`` is an orthogonal matrix.

        Examples
        ========

        >>> from sympy import Q, ask, MatrixSymbol, Identity
        >>> X = MatrixSymbol('X', 2, 2)
        >>> Y = MatrixSymbol('Y', 2, 3)
        >>> Z = MatrixSymbol('Z', 2, 2)
        >>> ask(Q.orthogonal(Y))
        False
        >>> ask(Q.orthogonal(X*Z*X), Q.orthogonal(X) & Q.orthogonal(Z))
        True
        >>> ask(Q.orthogonal(Identity(3)))
        True

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Orthogonal_matrix

        """
        return Predicate('orthogonal')

    @property
    def unitary(self):
        """
        Unitary matrix predicate.

        Unitary matrix is an analogue to orthogonal matrix. A square
        matrix ``M`` with complex elements is unitary if ``M'M = MM' = I``
        where ``M'`` is the conjugate transpose matrix of ``M``.

        Examples
        ========

        >>> from sympy import Q, ask, MatrixSymbol, Identity
        >>> X = MatrixSymbol('X', 2, 2)
        >>> Y = MatrixSymbol('Y', 2, 3)
        >>> Z = MatrixSymbol('Z', 2, 2)
        >>> ask(Q.unitary(Y))
        False
        >>> ask(Q.unitary(X*Z*X), Q.unitary(X) & Q.unitary(Z))
        True
        >>> ask(Q.unitary(Identity(3)))
        True

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Unitary_matrix

        """
        return Predicate('unitary')

    @property
    def positive_definite(self):
        """
        Positive definite matrix predicate.

        If ``M`` is a ``nxn`` symmetric real matrix, it is said
        to be positive definite if :math:`Z^TMZ` is positive for
        every non-zero column vector ``Z`` of ``n`` real numbers.

        Examples
        ========

        >>> from sympy import Q, ask, MatrixSymbol, Identity
        >>> X = MatrixSymbol('X', 2, 2)
        >>> Y = MatrixSymbol('Y', 2, 3)
        >>> Z = MatrixSymbol('Z', 2, 2)
        >>> ask(Q.positive_definite(Y))
        False
        >>> ask(Q.positive_definite(Identity(3)))
        True
        >>> ask(Q.positive_definite(X + Z), Q.positive_definite(X) &
        ...     Q.positive_definite(Z))
        True

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Positive-definite_matrix

        """
        return Predicate('positive_definite')

    @property
    def upper_triangular(self):
        """
        Upper triangular matrix predicate.

        A matrix ``M`` is called upper triangular matrix if :math:`a_{ij}=0`
        for :math:`i>j`.

        Examples
        ========

        >>> from sympy import Q, ask, ZeroMatrix, Identity
        >>> ask(Q.upper_triangular(Identity(3)))
        True
        >>> ask(Q.upper_triangular(ZeroMatrix(3, 3)))
        True

        References
        ==========

        .. [1] http://mathworld.wolfram.com/UpperTriangularMatrix.html

        """
        return Predicate('upper_triangular')

    @property
    def lower_triangular(self):
        """
        Lower triangular matrix predicate.

        A matrix ``M`` is called lower triangular matrix if :math:`a_{ij}=0`
        for :math:`i>j`.

        Examples
        ========

        >>> from sympy import Q, ask, ZeroMatrix, Identity
        >>> ask(Q.lower_triangular(Identity(3)))
        True
        >>> ask(Q.lower_triangular(ZeroMatrix(3, 3)))
        True

        References
        ==========

        .. [1] http://mathworld.wolfram.com/LowerTriangularMatrix.html
        """
        return Predicate('lower_triangular')

    @property
    def diagonal(self):
        """
        Diagonal matrix predicate.

        ``Q.diagonal(x)`` is true iff ``x`` is a diagonal matrix.

        Examples
        ========

        >>> from sympy import Q, ask, MatrixSymbol, ZeroMatrix
        >>> X = MatrixSymbol('X', 2, 2)
        >>> ask(Q.diagonal(ZeroMatrix(3, 3)))
        True
        >>> ask(Q.diagonal(X), Q.lower_triangular(X) &
        ...     Q.upper_triangular(X))
        True

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Diagonal_matrix

        """
        return Predicate('diagonal')

    @property
    def fullrank(self):
        """
        Fullrank matrix predicate.

        ``Q.fullrank(x)`` is true iff ``x`` is a full rank matrix.

        Examples
        ========

        >>> from sympy import Q, ask, MatrixSymbol, ZeroMatrix, Identity
        >>> X = MatrixSymbol('X', 2, 2)
        >>> ask(Q.fullrank(X.T), Q.fullrank(X))
        True
        >>> ask(Q.fullrank(ZeroMatrix(3, 3)))
        False
        >>> ask(Q.fullrank(Identity(3)))
        True

        """
        return Predicate('fullrank')

    @property
    def square(self):
        """
        Square matrix predicate.

        ``Q.square(x)`` is true iff ``x`` is a square matrix.

        Examples
        ========

        >>> from sympy import Q, ask, MatrixSymbol, ZeroMatrix, Identity
        >>> X = MatrixSymbol('X', 2, 2)
        >>> Y = MatrixSymbol('X', 2, 3)
        >>> ask(Q.square(X))
        True
        >>> ask(Q.square(Y))
        False
        >>> ask(Q.square(ZeroMatrix(3, 3)))
        True
        >>> ask(Q.square(Identity(3)))
        True

        References
        ==========

        .. [1] https://en.wikipedia.org/wiki/Square_matrix

        """
        return Predicate('square')

    @property
    def integer_elements(self):
        """
        Integer elements matrix predicate.

        ``Q.integer_elements(x)`` is true iff all the elements of ``x``
        are integers.

        Examples
        ========

        >>> from sympy import Q, ask, MatrixSymbol
        >>> X = MatrixSymbol('X', 4, 4)
        >>> ask(Q.integer(X[1, 2]), Q.integer_elements(X))
        True

        """
        return Predicate('integer_elements')

    @property
    def real_elements(self):
        """
        Real elements matrix predicate.

        ``Q.real_elements(x)`` is true iff all the elements of ``x``
        are real numbers.

        Examples
        ========

        >>> from sympy import Q, ask, MatrixSymbol
        >>> X = MatrixSymbol('X', 4, 4)
        >>> ask(Q.real(X[1, 2]), Q.real_elements(X))
        True

        """
        return Predicate('real_elements')

    @property
    def complex_elements(self):
        """
        Complex elements matrix predicate.

        ``Q.complex_elements(x)`` is true iff all the elements of ``x``
        are complex numbers.

        Examples
        ========

        >>> from sympy import Q, ask, MatrixSymbol
        >>> X = MatrixSymbol('X', 4, 4)
        >>> ask(Q.complex(X[1, 2]), Q.complex_elements(X))
        True
        >>> ask(Q.complex_elements(X), Q.integer_elements(X))
        True

        """
        return Predicate('complex_elements')

    @property
    def singular(self):
        """
        Singular matrix predicate.

        A matrix is singular iff the value of it's determinant is 0.

        Examples
        ========

        >>> from sympy import Q, ask, MatrixSymbol
        >>> X = MatrixSymbol('X', 4, 4)
        >>> ask(Q.singular(X), Q.invertible(X))
        False
        >>> ask(Q.singular(X), ~Q.invertible(X))
        True

        References
        ==========

        .. [1] http://mathworld.wolfram.com/SingularMatrix.html

        """
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
