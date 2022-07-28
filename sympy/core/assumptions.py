"""
This module contains the machinery handling assumptions.
Do also consider the guide :ref:`assumptions-guide`.

All symbolic objects have assumption attributes that can be accessed via
``.is_<assumption name>`` attribute.

Assumptions determine certain properties of symbolic objects and can
have 3 possible values: ``True``, ``False``, ``None``.  ``True`` is returned if the
object has the property and ``False`` is returned if it does not or cannot
(i.e. does not make sense):

    >>> from sympy import I
    >>> I.is_algebraic
    True
    >>> I.is_real
    False
    >>> I.is_prime
    False

When the property cannot be determined (or when a method is not
implemented) ``None`` will be returned. For example,  a generic symbol, ``x``,
may or may not be positive so a value of ``None`` is returned for ``x.is_positive``.

By default, all symbolic values are in the largest set in the given context
without specifying the property. For example, a symbol that has a property
being integer, is also real, complex, etc.

Here follows a list of possible assumption names:

.. glossary::

    commutative
        object commutes with any other object with
        respect to multiplication operation. See [12]_.

    complex
        object can have only values from the set
        of complex numbers. See [13]_.

    imaginary
        object value is a number that can be written as a real
        number multiplied by the imaginary unit ``I``.  See
        [3]_.  Please note that ``0`` is not considered to be an
        imaginary number, see
        `issue #7649 <https://github.com/sympy/sympy/issues/7649>`_.

    real
        object can have only values from the set
        of real numbers.

    extended_real
        object can have only values from the set
        of real numbers, ``oo`` and ``-oo``.

    integer
        object can have only values from the set
        of integers.

    odd
    even
        object can have only values from the set of
        odd (even) integers [2]_.

    prime
        object is a natural number greater than 1 that has
        no positive divisors other than 1 and itself.  See [6]_.

    composite
        object is a positive integer that has at least one positive
        divisor other than 1 or the number itself.  See [4]_.

    zero
        object has the value of 0.

    nonzero
        object is a real number that is not zero.

    rational
        object can have only values from the set
        of rationals.

    algebraic
        object can have only values from the set
        of algebraic numbers [11]_.

    transcendental
        object can have only values from the set
        of transcendental numbers [10]_.

    irrational
        object value cannot be represented exactly by :class:`~.Rational`, see [5]_.

    finite
    infinite
        object absolute value is bounded (arbitrarily large).
        See [7]_, [8]_, [9]_.

    negative
    nonnegative
        object can have only negative (nonnegative)
        values [1]_.

    positive
    nonpositive
        object can have only positive (nonpositive) values.

    extended_negative
    extended_nonnegative
    extended_positive
    extended_nonpositive
    extended_nonzero
        as without the extended part, but also including infinity with
        corresponding sign, e.g., extended_positive includes ``oo``

    hermitian
    antihermitian
        object belongs to the field of Hermitian
        (antihermitian) operators.

Examples
========

    >>> from sympy import Symbol
    >>> x = Symbol('x', real=True); x
    x
    >>> x.is_real
    True
    >>> x.is_complex
    True

See Also
========

.. seealso::

    :py:class:`sympy.core.numbers.ImaginaryUnit`
    :py:class:`sympy.core.numbers.Zero`
    :py:class:`sympy.core.numbers.One`
    :py:class:`sympy.core.numbers.Infinity`
    :py:class:`sympy.core.numbers.NegativeInfinity`
    :py:class:`sympy.core.numbers.ComplexInfinity`

Notes
=====

The fully-resolved assumptions for any SymPy expression
can be obtained as follows:

    >>> from sympy.core.assumptions import assumptions
    >>> x = Symbol('x',positive=True)
    >>> assumptions(x + I)
    {'commutative': True, 'complex': True, 'composite': False, 'even':
    False, 'extended_negative': False, 'extended_nonnegative': False,
    'extended_nonpositive': False, 'extended_nonzero': False,
    'extended_positive': False, 'extended_real': False, 'finite': True,
    'imaginary': False, 'infinite': False, 'integer': False, 'irrational':
    False, 'negative': False, 'noninteger': False, 'nonnegative': False,
    'nonpositive': False, 'nonzero': False, 'odd': False, 'positive':
    False, 'prime': False, 'rational': False, 'real': False, 'zero':
    False}

Developers Notes
================

The current (and possibly incomplete) values are stored
in the ``obj._assumptions dictionary``; queries to getter methods
(with property decorators) or attributes of objects/classes
will return values and update the dictionary.

    >>> eq = x**2 + I
    >>> eq._assumptions
    {}
    >>> eq.is_finite
    True
    >>> eq._assumptions
    {'finite': True, 'infinite': False}

For a :class:`~.Symbol`, there are two locations for assumptions that may
be of interest. The ``assumptions0`` attribute gives the full set of
assumptions derived from a given set of initial assumptions. The
latter assumptions are stored as ``Symbol._assumptions.generator``

    >>> Symbol('x', prime=True, even=True)._assumptions.generator
    {'even': True, 'prime': True}

The ``generator`` is not necessarily canonical nor is it filtered
in any way: it records the assumptions used to instantiate a Symbol
and (for storage purposes) represents a more compact representation
of the assumptions needed to recreate the full set in
``Symbol.assumptions0``.


References
==========

.. [1] https://en.wikipedia.org/wiki/Negative_number
.. [2] https://en.wikipedia.org/wiki/Parity_%28mathematics%29
.. [3] https://en.wikipedia.org/wiki/Imaginary_number
.. [4] https://en.wikipedia.org/wiki/Composite_number
.. [5] https://en.wikipedia.org/wiki/Irrational_number
.. [6] https://en.wikipedia.org/wiki/Prime_number
.. [7] https://en.wikipedia.org/wiki/Finite
.. [8] https://docs.python.org/3/library/math.html#math.isfinite
.. [9] http://docs.scipy.org/doc/numpy/reference/generated/numpy.isfinite.html
.. [10] https://en.wikipedia.org/wiki/Transcendental_number
.. [11] https://en.wikipedia.org/wiki/Algebraic_number
.. [12] https://en.wikipedia.org/wiki/Commutative_property
.. [13] https://en.wikipedia.org/wiki/Complex_number

"""

from .facts import FactRules, FactKB
from .core import BasicMeta
from .sympify import sympify

from sympy.core.random import shuffle

if 0:
    # normal generation of _assume_rules
 _assume_rules = FactRules([

    'integer        ->  rational',
    'rational       ->  real',
    'rational       ->  algebraic',
    'algebraic      ->  complex',
    'transcendental ==  complex & !algebraic',
    'real           ->  hermitian',
    'imaginary      ->  complex',
    'imaginary      ->  antihermitian',
    'extended_real  ->  commutative',
    'complex        ->  commutative',
    'complex        ->  finite',

    'odd            ==  integer & !even',
    'even           ==  integer & !odd',

    'real           ->  complex',
    'extended_real  ->  real | infinite',
    'real           ==  extended_real & finite',

    'extended_real        ==  extended_negative | zero | extended_positive',
    'extended_negative    ==  extended_nonpositive & extended_nonzero',
    'extended_positive    ==  extended_nonnegative & extended_nonzero',

    'extended_nonpositive ==  extended_real & !extended_positive',
    'extended_nonnegative ==  extended_real & !extended_negative',

    'real           ==  negative | zero | positive',
    'negative       ==  nonpositive & nonzero',
    'positive       ==  nonnegative & nonzero',

    'nonpositive    ==  real & !positive',
    'nonnegative    ==  real & !negative',

    'positive       ==  extended_positive & finite',
    'negative       ==  extended_negative & finite',
    'nonpositive    ==  extended_nonpositive & finite',
    'nonnegative    ==  extended_nonnegative & finite',
    'nonzero        ==  extended_nonzero & finite',

    'zero           ->  even & finite',
    'zero           ==  extended_nonnegative & extended_nonpositive',
    'zero           ==  nonnegative & nonpositive',
    'nonzero        ->  real',

    'prime          ->  integer & positive',
    'composite      ->  integer & positive & !prime',
    '!composite     ->  !positive | !even | prime',

    'irrational     ==  real & !rational',

    'imaginary      ->  !extended_real',

    'infinite       ==  !finite',
    'noninteger     ==  extended_real & !integer',
    'extended_nonzero == extended_real & !zero',
])
else:
    # pre-generated data
    _pre_calculated_assumptions = {'full_implications': {('prime', True): {('zero', False), ('extended_nonnegative', True), ('negative', False), ('nonzero', True), ('extended_positive', True), ('extended_real', True), ('extended_nonzero', True), ('composite', False), ('infinite', False), ('positive', True), ('extended_nonpositive', False), ('imaginary', False), ('commutative', True), ('noninteger', False), ('irrational', False), ('hermitian', True), ('real', True), ('finite', True), ('nonnegative', True), ('extended_negative', False), ('nonpositive', False), ('integer', True), ('algebraic', True), ('transcendental', False), ('complex', True), ('rational', True)}, ('zero', True): {('extended_nonnegative', True), ('transcendental', False), ('negative', False), ('extended_real', True), ('odd', False), ('composite', False), ('infinite', False), ('imaginary', False), ('commutative', True), ('nonpositive', True), ('extended_positive', False), ('nonzero', False), ('noninteger', False), ('irrational', False), ('extended_nonzero', False), ('positive', False), ('prime', False), ('hermitian', True), ('real', True), ('finite', True), ('nonnegative', True), ('extended_negative', False), ('integer', True), ('algebraic', True), ('complex', True), ('extended_nonpositive', True), ('even', True), ('rational', True)}, ('positive', True): {('infinite', False), ('zero', False), ('extended_nonnegative', True), ('extended_real', True), ('finite', True), ('extended_negative', False), ('extended_nonpositive', False), ('nonpositive', False), ('nonnegative', True), ('imaginary', False), ('negative', False), ('hermitian', True), ('commutative', True), ('nonzero', True), ('extended_positive', True), ('complex', True), ('extended_nonzero', True), ('real', True)}, ('negative', True): {('zero', False), ('nonnegative', False), ('nonzero', True), ('extended_real', True), ('extended_nonzero', True), ('composite', False), ('infinite', False), ('imaginary', False), ('extended_nonnegative', False), ('commutative', True), ('nonpositive', True), ('extended_positive', False), ('extended_negative', True), ('positive', False), ('prime', False), ('hermitian', True), ('real', True), ('finite', True), ('extended_nonpositive', True), ('complex', True)}, ('extended_real', False): {('zero', False), ('real', False), ('negative', False), ('nonnegative', False), ('integer', False), ('odd', False), ('composite', False), ('even', False), ('extended_nonpositive', False), ('extended_nonnegative', False), ('rational', False), ('extended_positive', False), ('nonzero', False), ('noninteger', False), ('irrational', False), ('extended_nonzero', False), ('positive', False), ('prime', False), ('extended_negative', False), ('nonpositive', False)}, ('complex', False): {('composite', False), ('irrational', False), ('zero', False), ('even', False), ('nonpositive', False), ('algebraic', False), ('imaginary', False), ('positive', False), ('rational', False), ('real', False), ('negative', False), ('nonnegative', False), ('integer', False), ('prime', False), ('transcendental', False), ('odd', False), ('nonzero', False)}, ('odd', True): {('noninteger', False), ('irrational', False), ('zero', False), ('even', False), ('extended_real', True), ('finite', True), ('infinite', False), ('imaginary', False), ('integer', True), ('algebraic', True), ('hermitian', True), ('commutative', True), ('nonzero', True), ('transcendental', False), ('complex', True), ('extended_nonzero', True), ('rational', True), ('real', True)}, ('imaginary', True): {('zero', False), ('antihermitian', True), ('real', False), ('negative', False), ('nonnegative', False), ('integer', False), ('odd', False), ('composite', False), ('even', False), ('infinite', False), ('extended_nonpositive', False), ('extended_nonnegative', False), ('rational', False), ('commutative', True), ('extended_positive', False), ('nonzero', False), ('noninteger', False), ('irrational', False), ('extended_real', False), ('extended_nonzero', False), ('positive', False), ('prime', False), ('finite', True), ('extended_negative', False), ('nonpositive', False), ('complex', True)}, ('extended_negative', True): {('composite', False), ('zero', False), ('extended_real', True), ('imaginary', False), ('extended_nonnegative', False), ('positive', False), ('nonnegative', False), ('prime', False), ('commutative', True), ('extended_nonpositive', True), ('extended_nonzero', True), ('extended_positive', False)}, ('extended_positive', True): {('zero', False), ('extended_nonnegative', True), ('extended_real', True), ('extended_negative', False), ('extended_nonpositive', False), ('nonpositive', False), ('imaginary', False), ('negative', False), ('commutative', True), ('extended_nonzero', True)}, ('rational', False): {('composite', False), ('zero', False), ('even', False), ('integer', False), ('prime', False), ('odd', False)}, ('irrational', True): {('composite', False), ('zero', False), ('even', False), ('finite', True), ('extended_real', True), ('infinite', False), ('imaginary', False), ('rational', False), ('integer', False), ('prime', False), ('hermitian', True), ('noninteger', True), ('commutative', True), ('nonzero', True), ('complex', True), ('extended_nonzero', True), ('odd', False), ('real', True)}, ('commutative', False): {('zero', False), ('real', False), ('negative', False), ('nonnegative', False), ('integer', False), ('odd', False), ('algebraic', False), ('composite', False), ('complex', False), ('even', False), ('extended_nonpositive', False), ('imaginary', False), ('extended_nonnegative', False), ('rational', False), ('extended_positive', False), ('nonzero', False), ('noninteger', False), ('irrational', False), ('extended_real', False), ('extended_nonzero', False), ('positive', False), ('prime', False), ('extended_negative', False), ('nonpositive', False), ('transcendental', False)}, ('extended_real', True): {('imaginary', False), ('commutative', True)}, ('finite', True): {('infinite', False)}, ('infinite', True): {('zero', False), ('real', False), ('negative', False), ('finite', False), ('nonnegative', False), ('integer', False), ('algebraic', False), ('odd', False), ('composite', False), ('complex', False), ('even', False), ('imaginary', False), ('rational', False), ('nonzero', False), ('irrational', False), ('positive', False), ('prime', False), ('nonpositive', False), ('transcendental', False)}, ('algebraic', False): {('composite', False), ('zero', False), ('even', False), ('rational', False), ('integer', False), ('prime', False), ('odd', False)}, ('extended_nonpositive', False): {('nonpositive', False), ('zero', False), ('extended_negative', False), ('negative', False)}, ('finite', False): {('zero', False), ('real', False), ('negative', False), ('nonnegative', False), ('integer', False), ('algebraic', False), ('odd', False), ('composite', False), ('complex', False), ('even', False), ('imaginary', False), ('rational', False), ('nonzero', False), ('irrational', False), ('positive', False), ('prime', False), ('nonpositive', False), ('infinite', True), ('transcendental', False)}, ('extended_nonnegative', False): {('composite', False), ('zero', False), ('positive', False), ('nonnegative', False), ('prime', False), ('extended_positive', False)}, ('noninteger', True): {('composite', False), ('zero', False), ('even', False), ('extended_real', True), ('imaginary', False), ('integer', False), ('prime', False), ('commutative', True), ('extended_nonzero', True), ('odd', False)}, ('real', False): {('composite', False), ('irrational', False), ('zero', False), ('even', False), ('nonpositive', False), ('rational', False), ('positive', False), ('negative', False), ('integer', False), ('prime', False), ('nonnegative', False), ('odd', False), ('nonzero', False)}, ('hermitian', False): {('composite', False), ('irrational', False), ('zero', False), ('even', False), ('nonpositive', False), ('rational', False), ('positive', False), ('real', False), ('negative', False), ('nonnegative', False), ('prime', False), ('integer', False), ('odd', False), ('nonzero', False)}, ('transcendental', True): {('infinite', False), ('composite', False), ('zero', False), ('finite', True), ('even', False), ('rational', False), ('integer', False), ('prime', False), ('commutative', True), ('complex', True), ('algebraic', False), ('odd', False)}, ('integer', True): {('noninteger', False), ('irrational', False), ('infinite', False), ('finite', True), ('imaginary', False), ('algebraic', True), ('hermitian', True), ('commutative', True), ('transcendental', False), ('complex', True), ('extended_real', True), ('rational', True), ('real', True)}, ('composite', True): {('zero', False), ('extended_nonnegative', True), ('negative', False), ('nonzero', True), ('extended_positive', True), ('extended_real', True), ('extended_nonzero', True), ('infinite', False), ('positive', True), ('extended_nonpositive', False), ('imaginary', False), ('commutative', True), ('noninteger', False), ('irrational', False), ('prime', False), ('hermitian', True), ('real', True), ('finite', True), ('nonnegative', True), ('extended_negative', False), ('nonpositive', False), ('integer', True), ('algebraic', True), ('transcendental', False), ('complex', True), ('rational', True)}, ('real', True): {('infinite', False), ('finite', True), ('imaginary', False), ('hermitian', True), ('commutative', True), ('complex', True), ('extended_real', True)}, ('extended_nonzero', True): {('zero', False), ('imaginary', False), ('extended_real', True), ('commutative', True)}, ('nonzero', True): {('infinite', False), ('zero', False), ('finite', True), ('extended_real', True), ('imaginary', False), ('hermitian', True), ('commutative', True), ('complex', True), ('extended_nonzero', True), ('real', True)}, ('algebraic', True): {('infinite', False), ('finite', True), ('commutative', True), ('transcendental', False), ('complex', True)}, ('extended_nonpositive', True): {('composite', False), ('imaginary', False), ('positive', False), ('prime', False), ('commutative', True), ('extended_real', True), ('extended_positive', False)}, ('nonpositive', True): {('infinite', False), ('composite', False), ('finite', True), ('imaginary', False), ('positive', False), ('prime', False), ('hermitian', True), ('commutative', True), ('extended_nonpositive', True), ('complex', True), ('extended_real', True), ('extended_positive', False), ('real', True)}, ('positive', False): {('composite', False), ('prime', False)}, ('rational', True): {('infinite', False), ('irrational', False), ('finite', True), ('imaginary', False), ('algebraic', True), ('hermitian', True), ('commutative', True), ('transcendental', False), ('complex', True), ('extended_real', True), ('real', True)}, ('complex', True): {('infinite', False), ('finite', True), ('commutative', True)}, ('extended_nonnegative', True): {('extended_negative', False), ('imaginary', False), ('negative', False), ('commutative', True), ('extended_real', True)}, ('nonnegative', True): {('infinite', False), ('finite', True), ('extended_nonnegative', True), ('extended_negative', False), ('imaginary', False), ('negative', False), ('hermitian', True), ('commutative', True), ('complex', True), ('extended_real', True), ('real', True)}, ('even', True): {('noninteger', False), ('irrational', False), ('infinite', False), ('finite', True), ('imaginary', False), ('integer', True), ('algebraic', True), ('hermitian', True), ('commutative', True), ('rational', True), ('transcendental', False), ('complex', True), ('extended_real', True), ('odd', False), ('real', True)}, ('nonpositive', False): {('zero', False), ('negative', False)}, ('nonzero', False): {('composite', False), ('prime', False), ('positive', False), ('negative', False)}, ('extended_nonzero', False): {('composite', False), ('extended_negative', False), ('positive', False), ('negative', False), ('extended_positive', False), ('prime', False), ('nonzero', False)}, ('even', False): {('zero', False)}, ('integer', False): {('composite', False), ('zero', False), ('even', False), ('prime', False), ('odd', False)}, ('nonnegative', False): {('composite', False), ('zero', False), ('prime', False), ('positive', False)}, ('antihermitian', False): {('imaginary', False)}, ('infinite', False): {('finite', True)}, ('extended_positive', False): {('composite', False), ('prime', False), ('positive', False)}, ('extended_negative', False): {('negative', False)}, ('odd', False): set(), ('zero', False): set(), ('negative', False): set(), ('prime', False): set(), ('composite', False): set()}, 'beta_triggers': {('prime', True): [1, 2], ('zero', True): [], ('positive', True): [0, 30, 32, 33, 34, 35], ('negative', True): [0, 34, 35], ('extended_real', False): [], ('complex', False): [3, 4, 7, 8, 9, 35, 36], ('odd', True): [8, 9, 11, 12, 13, 14, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26], ('imaginary', True): [0], ('extended_negative', True): [4, 5, 6, 15, 16, 21, 24, 25, 27, 35], ('extended_positive', True): [4, 5, 6, 15, 17, 22, 23, 26, 27, 35], ('rational', False): [7, 8, 9, 15, 16, 17, 34, 35, 36], ('irrational', True): [0, 8, 9, 11, 12, 13, 14, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26], ('commutative', False): [], ('extended_real', True): [4, 5, 6, 8, 9, 10, 13, 14, 35, 36], ('finite', True): [3, 5, 6, 23, 24, 25, 26, 27], ('infinite', True): [7, 8, 9, 35, 36], ('algebraic', False): [0, 7, 8, 9, 15, 16, 17, 34, 35, 36], ('extended_nonpositive', False): [7, 8, 14, 15, 17, 22, 36], ('finite', False): [7, 8, 9, 35, 36], ('extended_nonnegative', False): [7, 9, 13, 15, 16, 21, 36], ('noninteger', True): [4, 5, 6, 8, 9, 11, 12, 13, 14, 15, 16, 17, 27], ('real', False): [3, 4, 7, 8, 9, 35, 36], ('hermitian', False): [3, 4, 7, 8, 9, 35, 36], ('transcendental', True): [3, 5, 6, 7, 8, 9, 15, 16, 17, 23, 24, 25, 26, 27, 34, 35, 36], ('integer', True): [1, 2, 8, 9, 10, 13, 14, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 36], ('composite', True): [1, 2], ('real', True): [0, 8, 9, 10, 13, 14, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 34, 35, 36], ('extended_nonzero', True): [4, 5, 6, 8, 9, 11, 12, 13, 14, 15, 16, 17, 27, 35], ('nonzero', True): [0, 8, 9, 11, 12, 13, 14, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 34, 35], ('algebraic', True): [3, 5, 6, 23, 24, 25, 26, 27], ('extended_nonpositive', True): [4, 5, 6, 9, 10, 11, 14, 15, 16, 18, 21, 25, 28, 35, 36], ('nonpositive', True): [0, 9, 10, 11, 14, 16, 18, 19, 22, 24, 26, 27, 28, 29, 34, 35, 36], ('positive', False): [15, 16, 18, 21], ('rational', True): [8, 9, 10, 13, 14, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 35, 36], ('complex', True): [0, 3, 5, 6, 23, 24, 25, 26, 27], ('extended_nonnegative', True): [4, 5, 6, 8, 10, 12, 13, 15, 17, 18, 22, 26, 28, 35, 36], ('nonnegative', True): [0, 8, 10, 12, 13, 17, 18, 20, 21, 23, 25, 27, 28, 29, 34, 35, 36], ('even', True): [8, 9, 10, 13, 14, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 30, 31, 33, 36], ('nonpositive', False): [7, 8, 9, 15, 17, 22, 36], ('nonzero', False): [15, 18, 21, 22], ('extended_nonzero', False): [7, 10, 13, 14, 15, 18, 21, 22], ('even', False): [1, 7, 8, 9, 15, 16, 17, 36], ('integer', False): [7, 8, 9, 15, 16, 17, 35, 36], ('nonnegative', False): [7, 8, 9, 15, 16, 21, 36], ('antihermitian', False): [], ('infinite', False): [3, 5, 6, 23, 24, 25, 26, 27], ('extended_positive', False): [7, 9, 10, 13, 15, 16, 18, 21], ('extended_negative', False): [7, 8, 10, 14, 15, 17, 18, 22], ('odd', False): [2], ('zero', False): [7, 8, 9, 15, 16, 17, 36], ('negative', False): [15, 17, 18, 22], ('prime', False): [30, 31, 32], ('composite', False): [31, 32, 33]}, 'prereq': {'zero': {'extended_nonzero', 'negative', 'odd', 'nonzero', 'hermitian', 'real', 'rational', 'extended_real', 'extended_nonpositive', 'extended_nonnegative', 'imaginary', 'infinite', 'commutative', 'extended_positive', 'finite', 'irrational', 'even', 'nonpositive', 'algebraic', 'positive', 'complex', 'composite', 'extended_negative', 'prime', 'noninteger', 'nonnegative', 'integer', 'transcendental'}, 'extended_nonnegative': {'nonnegative', 'extended_positive', 'negative', 'positive', 'composite', 'prime', 'extended_real', 'extended_negative', 'imaginary', 'commutative', 'zero'}, 'negative': {'extended_nonzero', 'nonzero', 'hermitian', 'real', 'extended_nonpositive', 'extended_real', 'extended_nonnegative', 'imaginary', 'infinite', 'commutative', 'extended_positive', 'zero', 'finite', 'nonpositive', 'positive', 'complex', 'composite', 'extended_negative', 'prime', 'nonnegative'}, 'nonzero': {'extended_nonzero', 'negative', 'odd', 'hermitian', 'real', 'extended_real', 'imaginary', 'infinite', 'commutative', 'zero', 'finite', 'irrational', 'positive', 'complex', 'composite', 'prime'}, 'extended_positive': {'extended_nonnegative', 'extended_nonzero', 'nonpositive', 'negative', 'positive', 'composite', 'extended_nonpositive', 'prime', 'extended_real', 'extended_negative', 'imaginary', 'commutative', 'zero'}, 'extended_real': {'extended_nonzero', 'negative', 'odd', 'nonzero', 'real', 'rational', 'extended_nonpositive', 'extended_nonnegative', 'imaginary', 'commutative', 'extended_positive', 'zero', 'irrational', 'even', 'nonpositive', 'positive', 'composite', 'extended_negative', 'prime', 'noninteger', 'nonnegative', 'integer'}, 'extended_nonzero': {'irrational', 'extended_positive', 'noninteger', 'negative', 'odd', 'positive', 'composite', 'nonzero', 'prime', 'extended_real', 'extended_negative', 'imaginary', 'commutative', 'zero'}, 'composite': {'extended_nonzero', 'negative', 'nonzero', 'hermitian', 'real', 'rational', 'extended_real', 'extended_nonpositive', 'extended_nonnegative', 'imaginary', 'infinite', 'commutative', 'zero', 'extended_positive', 'finite', 'irrational', 'nonpositive', 'algebraic', 'positive', 'complex', 'extended_negative', 'prime', 'noninteger', 'nonnegative', 'integer', 'transcendental'}, 'infinite': {'negative', 'odd', 'nonzero', 'real', 'rational', 'imaginary', 'zero', 'finite', 'irrational', 'even', 'nonpositive', 'algebraic', 'positive', 'composite', 'complex', 'prime', 'nonnegative', 'integer', 'transcendental'}, 'positive': {'extended_nonzero', 'negative', 'nonzero', 'hermitian', 'real', 'extended_nonpositive', 'extended_real', 'extended_nonnegative', 'imaginary', 'infinite', 'commutative', 'zero', 'extended_positive', 'finite', 'nonpositive', 'composite', 'complex', 'extended_negative', 'prime', 'nonnegative'}, 'extended_nonpositive': {'nonpositive', 'extended_positive', 'negative', 'positive', 'composite', 'prime', 'extended_real', 'extended_negative', 'imaginary', 'commutative', 'zero'}, 'imaginary': {'extended_nonzero', 'negative', 'odd', 'nonzero', 'real', 'rational', 'extended_real', 'extended_nonpositive', 'extended_nonnegative', 'infinite', 'commutative', 'extended_positive', 'zero', 'finite', 'irrational', 'even', 'antihermitian', 'nonpositive', 'positive', 'complex', 'composite', 'extended_negative', 'prime', 'noninteger', 'nonnegative', 'integer'}, 'commutative': {'extended_nonzero', 'negative', 'odd', 'nonzero', 'real', 'rational', 'extended_real', 'extended_nonpositive', 'extended_nonnegative', 'imaginary', 'extended_positive', 'zero', 'irrational', 'even', 'nonpositive', 'algebraic', 'positive', 'composite', 'complex', 'extended_negative', 'prime', 'noninteger', 'nonnegative', 'integer', 'transcendental'}, 'noninteger': {'irrational', 'even', 'odd', 'composite', 'prime', 'extended_real', 'imaginary', 'commutative', 'zero', 'integer'}, 'irrational': {'finite', 'even', 'odd', 'composite', 'complex', 'infinite', 'hermitian', 'real', 'rational', 'prime', 'extended_real', 'imaginary', 'commutative', 'zero', 'integer'}, 'hermitian': {'irrational', 'even', 'nonpositive', 'negative', 'odd', 'positive', 'composite', 'nonzero', 'real', 'rational', 'prime', 'nonnegative', 'zero', 'integer'}, 'real': {'negative', 'odd', 'nonzero', 'hermitian', 'rational', 'extended_real', 'imaginary', 'infinite', 'commutative', 'zero', 'finite', 'irrational', 'even', 'nonpositive', 'positive', 'complex', 'composite', 'prime', 'nonnegative', 'integer'}, 'finite': {'negative', 'odd', 'nonzero', 'real', 'rational', 'imaginary', 'infinite', 'zero', 'irrational', 'even', 'nonpositive', 'algebraic', 'positive', 'composite', 'complex', 'prime', 'nonnegative', 'integer', 'transcendental'}, 'nonnegative': {'finite', 'extended_nonnegative', 'negative', 'positive', 'complex', 'composite', 'infinite', 'hermitian', 'real', 'prime', 'extended_real', 'extended_negative', 'imaginary', 'commutative', 'zero'}, 'extended_negative': {'extended_nonnegative', 'extended_nonzero', 'extended_positive', 'negative', 'positive', 'composite', 'extended_nonpositive', 'prime', 'extended_real', 'nonnegative', 'imaginary', 'commutative', 'zero'}, 'nonpositive': {'finite', 'extended_positive', 'negative', 'positive', 'complex', 'composite', 'infinite', 'hermitian', 'real', 'extended_nonpositive', 'prime', 'extended_real', 'imaginary', 'commutative', 'zero'}, 'integer': {'odd', 'hermitian', 'real', 'rational', 'extended_real', 'imaginary', 'infinite', 'commutative', 'zero', 'finite', 'irrational', 'even', 'algebraic', 'composite', 'complex', 'prime', 'noninteger', 'transcendental'}, 'algebraic': {'finite', 'even', 'odd', 'composite', 'complex', 'infinite', 'rational', 'prime', 'commutative', 'zero', 'integer', 'transcendental'}, 'transcendental': {'finite', 'even', 'odd', 'algebraic', 'composite', 'complex', 'infinite', 'rational', 'prime', 'commutative', 'zero', 'integer'}, 'complex': {'negative', 'odd', 'nonzero', 'real', 'rational', 'imaginary', 'infinite', 'commutative', 'zero', 'finite', 'irrational', 'even', 'nonpositive', 'algebraic', 'positive', 'composite', 'prime', 'nonnegative', 'integer', 'transcendental'}, 'rational': {'odd', 'hermitian', 'real', 'extended_real', 'imaginary', 'infinite', 'commutative', 'zero', 'finite', 'irrational', 'even', 'algebraic', 'composite', 'complex', 'prime', 'integer', 'transcendental'}, 'odd': {'hermitian', 'real', 'rational', 'extended_real', 'imaginary', 'infinite', 'commutative', 'zero', 'finite', 'irrational', 'even', 'algebraic', 'complex', 'noninteger', 'integer', 'transcendental'}, 'prime': {'extended_nonzero', 'negative', 'nonzero', 'hermitian', 'real', 'rational', 'extended_real', 'extended_nonpositive', 'extended_nonnegative', 'imaginary', 'infinite', 'commutative', 'zero', 'extended_positive', 'finite', 'irrational', 'nonpositive', 'algebraic', 'composite', 'complex', 'positive', 'extended_negative', 'noninteger', 'nonnegative', 'integer', 'transcendental'}, 'even': {'odd', 'hermitian', 'real', 'rational', 'extended_real', 'imaginary', 'infinite', 'commutative', 'zero', 'finite', 'irrational', 'algebraic', 'complex', 'noninteger', 'integer', 'transcendental'}, 'antihermitian': {'imaginary'}}, 'beta_rules': [({('complex', True), ('algebraic', False)}, ('transcendental', True)), ({('even', False), ('integer', True)}, ('odd', True)), ({('integer', True), ('odd', False)}, ('even', True)), ({('infinite', False), ('real', False)}, ('extended_real', False)), ({('extended_real', True), ('real', False)}, ('infinite', True)), ({('infinite', False), ('extended_real', True)}, ('real', True)), ({('extended_real', True), ('finite', True)}, ('real', True)), ({('zero', False), ('extended_positive', False), ('extended_negative', False)}, ('extended_real', False)), ({('zero', False), ('extended_real', True), ('extended_negative', False)}, ('extended_positive', True)), ({('zero', False), ('extended_real', True), ('extended_positive', False)}, ('extended_negative', True)), ({('extended_real', True), ('extended_positive', False), ('extended_negative', False)}, ('zero', True)), ({('extended_nonpositive', True), ('extended_nonzero', True)}, ('extended_negative', True)), ({('extended_nonzero', True), ('extended_nonnegative', True)}, ('extended_positive', True)), ({('extended_real', True), ('extended_positive', False)}, ('extended_nonpositive', True)), ({('extended_real', True), ('extended_negative', False)}, ('extended_nonnegative', True)), ({('zero', False), ('positive', False), ('negative', False)}, ('real', False)), ({('zero', False), ('positive', False), ('real', True)}, ('negative', True)), ({('zero', False), ('negative', False), ('real', True)}, ('positive', True)), ({('positive', False), ('negative', False), ('real', True)}, ('zero', True)), ({('nonpositive', True), ('nonzero', True)}, ('negative', True)), ({('nonnegative', True), ('nonzero', True)}, ('positive', True)), ({('positive', False), ('real', True)}, ('nonpositive', True)), ({('negative', False), ('real', True)}, ('nonnegative', True)), ({('finite', True), ('extended_positive', True)}, ('positive', True)), ({('finite', True), ('extended_negative', True)}, ('negative', True)), ({('extended_nonpositive', True), ('finite', True)}, ('nonpositive', True)), ({('extended_nonnegative', True), ('finite', True)}, ('nonnegative', True)), ({('extended_nonzero', True), ('finite', True)}, ('nonzero', True)), ({('extended_nonpositive', True), ('extended_nonnegative', True)}, ('zero', True)), ({('nonpositive', True), ('nonnegative', True)}, ('zero', True)), ({('even', True), ('prime', False), ('positive', True)}, ('composite', True)), ({('composite', False), ('prime', False), ('even', True)}, ('positive', False)), ({('composite', False), ('prime', False), ('positive', True)}, ('even', False)), ({('composite', False), ('even', True), ('positive', True)}, ('prime', True)), ({('rational', False), ('real', True)}, ('irrational', True)), ({('integer', False), ('extended_real', True)}, ('noninteger', True)), ({('zero', False), ('extended_real', True)}, ('extended_nonzero', True))], 'defined_facts': ['extended_nonzero', 'negative', 'hermitian', 'real', 'rational', 'extended_real', 'extended_nonnegative', 'extended_positive', 'zero', 'nonpositive', 'positive', 'complex', 'composite', 'extended_negative', 'nonnegative', 'integer', 'transcendental', 'odd', 'nonzero', 'extended_nonpositive', 'imaginary', 'infinite', 'commutative', 'finite', 'irrational', 'even', 'antihermitian', 'algebraic', 'prime', 'noninteger']}
    _assume_rules=FactRules.from_python(_pre_calculated_assumptions)
    
_assume_defined = _assume_rules.defined_facts.copy()
_assume_defined.add('polar')
_assume_defined = frozenset(_assume_defined)


def assumptions(expr, _check=None):
    """return the T/F assumptions of ``expr``"""
    n = sympify(expr)
    if n.is_Symbol:
        rv = n.assumptions0  # are any important ones missing?
        if _check is not None:
            rv = {k: rv[k] for k in set(rv) & set(_check)}
        return rv
    rv = {}
    for k in _assume_defined if _check is None else _check:
        v = getattr(n, 'is_{}'.format(k))
        if v is not None:
            rv[k] = v
    return rv


def common_assumptions(exprs, check=None):
    """return those assumptions which have the same True or False
    value for all the given expressions.

    Examples
    ========

    >>> from sympy.core import common_assumptions
    >>> from sympy import oo, pi, sqrt
    >>> common_assumptions([-4, 0, sqrt(2), 2, pi, oo])
    {'commutative': True, 'composite': False,
    'extended_real': True, 'imaginary': False, 'odd': False}

    By default, all assumptions are tested; pass an iterable of the
    assumptions to limit those that are reported:

    >>> common_assumptions([0, 1, 2], ['positive', 'integer'])
    {'integer': True}
    """
    check = _assume_defined if check is None else set(check)
    if not check or not exprs:
        return {}

    # get all assumptions for each
    assume = [assumptions(i, _check=check) for i in sympify(exprs)]
    # focus on those of interest that are True
    for i, e in enumerate(assume):
        assume[i] = {k: e[k] for k in set(e) & check}
    # what assumptions are in common?
    common = set.intersection(*[set(i) for i in assume])
    # which ones hold the same value
    a = assume[0]
    return {k: a[k] for k in common if all(a[k] == b[k]
        for b in assume)}


def failing_assumptions(expr, **assumptions):
    """
    Return a dictionary containing assumptions with values not
    matching those of the passed assumptions.

    Examples
    ========

    >>> from sympy import failing_assumptions, Symbol

    >>> x = Symbol('x', positive=True)
    >>> y = Symbol('y')
    >>> failing_assumptions(6*x + y, positive=True)
    {'positive': None}

    >>> failing_assumptions(x**2 - 1, positive=True)
    {'positive': None}

    If *expr* satisfies all of the assumptions, an empty dictionary is returned.

    >>> failing_assumptions(x**2, positive=True)
    {}

    """
    expr = sympify(expr)
    failed = {}
    for k in assumptions:
        test = getattr(expr, 'is_%s' % k, None)
        if test is not assumptions[k]:
            failed[k] = test
    return failed  # {} or {assumption: value != desired}


def check_assumptions(expr, against=None, **assume):
    """
    Checks whether assumptions of ``expr`` match the T/F assumptions
    given (or possessed by ``against``). True is returned if all
    assumptions match; False is returned if there is a mismatch and
    the assumption in ``expr`` is not None; else None is returned.

    Explanation
    ===========

    *assume* is a dict of assumptions with True or False values

    Examples
    ========

    >>> from sympy import Symbol, pi, I, exp, check_assumptions
    >>> check_assumptions(-5, integer=True)
    True
    >>> check_assumptions(pi, real=True, integer=False)
    True
    >>> check_assumptions(pi, negative=True)
    False
    >>> check_assumptions(exp(I*pi/7), real=False)
    True
    >>> x = Symbol('x', positive=True)
    >>> check_assumptions(2*x + 1, positive=True)
    True
    >>> check_assumptions(-2*x - 5, positive=True)
    False

    To check assumptions of *expr* against another variable or expression,
    pass the expression or variable as ``against``.

    >>> check_assumptions(2*x + 1, x)
    True

    To see if a number matches the assumptions of an expression, pass
    the number as the first argument, else its specific assumptions
    may not have a non-None value in the expression:

    >>> check_assumptions(x, 3)
    >>> check_assumptions(3, x)
    True

    ``None`` is returned if ``check_assumptions()`` could not conclude.

    >>> check_assumptions(2*x - 1, x)

    >>> z = Symbol('z')
    >>> check_assumptions(z, real=True)

    See Also
    ========

    failing_assumptions

    """
    expr = sympify(expr)
    if against is not None:
        if assume:
            raise ValueError(
                'Expecting `against` or `assume`, not both.')
        assume = assumptions(against)
    known = True
    for k, v in assume.items():
        if v is None:
            continue
        e = getattr(expr, 'is_' + k, None)
        if e is None:
            known = None
        elif v != e:
            return False
    return known


class StdFactKB(FactKB):
    """A FactKB specialized for the built-in rules

    This is the only kind of FactKB that Basic objects should use.
    """
    def __init__(self, facts=None):
        super().__init__(_assume_rules)
        # save a copy of the facts dict
        if not facts:
            self._generator = {}
        elif not isinstance(facts, FactKB):
            self._generator = facts.copy()
        else:
            self._generator = facts.generator
        if facts:
            self.deduce_all_facts(facts)

    def copy(self):
        return self.__class__(self)

    @property
    def generator(self):
        return self._generator.copy()


def as_property(fact):
    """Convert a fact name to the name of the corresponding property"""
    return 'is_%s' % fact


def make_property(fact):
    """Create the automagic property corresponding to a fact."""

    def getit(self):
        try:
            return self._assumptions[fact]
        except KeyError:
            if self._assumptions is self.default_assumptions:
                self._assumptions = self.default_assumptions.copy()
            return _ask(fact, self)

    getit.func_name = as_property(fact)
    return property(getit)


def _ask(fact, obj):
    """
    Find the truth value for a property of an object.

    This function is called when a request is made to see what a fact
    value is.

    For this we use several techniques:

    First, the fact-evaluation function is tried, if it exists (for
    example _eval_is_integer). Then we try related facts. For example

        rational   -->   integer

    another example is joined rule:

        integer & !odd  --> even

    so in the latter case if we are looking at what 'even' value is,
    'integer' and 'odd' facts will be asked.

    In all cases, when we settle on some fact value, its implications are
    deduced, and the result is cached in ._assumptions.
    """
    # FactKB which is dict-like and maps facts to their known values:
    assumptions = obj._assumptions

    # A dict that maps facts to their handlers:
    handler_map = obj._prop_handler

    # This is our queue of facts to check:
    facts_to_check = [fact]
    facts_queued = {fact}

    # Loop over the queue as it extends
    for fact_i in facts_to_check:

        # If fact_i has already been determined then we don't need to rerun the
        # handler. There is a potential race condition for multithreaded code
        # though because it's possible that fact_i was checked in another
        # thread. The main logic of the loop below would potentially skip
        # checking assumptions[fact] in this case so we check it once after the
        # loop to be sure.
        if fact_i in assumptions:
            continue

        # Now we call the associated handler for fact_i if it exists.
        fact_i_value = None
        handler_i = handler_map.get(fact_i)
        if handler_i is not None:
            fact_i_value = handler_i(obj)

        # If we get a new value for fact_i then we should update our knowledge
        # of fact_i as well as any related facts that can be inferred using the
        # inference rules connecting the fact_i and any other fact values that
        # are already known.
        if fact_i_value is not None:
            assumptions.deduce_all_facts(((fact_i, fact_i_value),))

        # Usually if assumptions[fact] is now not None then that is because of
        # the call to deduce_all_facts above. The handler for fact_i returned
        # True or False and knowing fact_i (which is equal to fact in the first
        # iteration) implies knowing a value for fact. It is also possible
        # though that independent code e.g. called indirectly by the handler or
        # called in another thread in a multithreaded context might have
        # resulted in assumptions[fact] being set. Either way we return it.
        fact_value = assumptions.get(fact)
        if fact_value is not None:
            return fact_value

        # Extend the queue with other facts that might determine fact_i. Here
        # we randomise the order of the facts that are checked. This should not
        # lead to any non-determinism if all handlers are logically consistent
        # with the inference rules for the facts. Non-deterministic assumptions
        # queries can result from bugs in the handlers that are exposed by this
        # call to shuffle. These are pushed to the back of the queue meaning
        # that the inference graph is traversed in breadth-first order.
        new_facts_to_check = list(_assume_rules.prereq[fact_i] - facts_queued)
        shuffle(new_facts_to_check)
        facts_to_check.extend(new_facts_to_check)
        facts_queued.update(new_facts_to_check)

    # The above loop should be able to handle everything fine in a
    # single-threaded context but in multithreaded code it is possible that
    # this thread skipped computing a particular fact that was computed in
    # another thread (due to the continue). In that case it is possible that
    # fact was inferred and is now stored in the assumptions dict but it wasn't
    # checked for in the body of the loop. This is an obscure case but to make
    # sure we catch it we check once here at the end of the loop.
    if fact in assumptions:
        return assumptions[fact]

    # This query can not be answered. It's possible that e.g. another thread
    # has already stored None for fact but assumptions._tell does not mind if
    # we call _tell twice setting the same value. If this raises
    # InconsistentAssumptions then it probably means that another thread
    # attempted to compute this and got a value of True or False rather than
    # None. In that case there must be a bug in at least one of the handlers.
    # If the handlers are all deterministic and are consistent with the
    # inference rules then the same value should be computed for fact in all
    # threads.
    assumptions._tell(fact, None)
    return None


class ManagedProperties(BasicMeta):
    """Metaclass for classes with old-style assumptions"""
    def __init__(cls, *args, **kws):
        BasicMeta.__init__(cls, *args, **kws)

        local_defs = {}
        for k in _assume_defined:
            attrname = as_property(k)
            v = cls.__dict__.get(attrname, '')
            if isinstance(v, (bool, int, type(None))):
                if v is not None:
                    v = bool(v)
                local_defs[k] = v

        defs = {}
        for base in reversed(cls.__bases__):
            assumptions = getattr(base, '_explicit_class_assumptions', None)
            if assumptions is not None:
                defs.update(assumptions)
        defs.update(local_defs)

        cls._explicit_class_assumptions = defs
        cls.default_assumptions = StdFactKB(defs)

        cls._prop_handler = {}
        for k in _assume_defined:
            eval_is_meth = getattr(cls, '_eval_is_%s' % k, None)
            if eval_is_meth is not None:
                cls._prop_handler[k] = eval_is_meth

        # Put definite results directly into the class dict, for speed
        for k, v in cls.default_assumptions.items():
            setattr(cls, as_property(k), v)

        # protection e.g. for Integer.is_even=F <- (Rational.is_integer=F)
        derived_from_bases = set()
        for base in cls.__bases__:
            default_assumptions = getattr(base, 'default_assumptions', None)
            # is an assumption-aware class
            if default_assumptions is not None:
                derived_from_bases.update(default_assumptions)

        for fact in derived_from_bases - set(cls.default_assumptions):
            pname = as_property(fact)
            if pname not in cls.__dict__:
                setattr(cls, pname, make_property(fact))

        # Finally, add any missing automagic property (e.g. for Basic)
        for fact in _assume_defined:
            pname = as_property(fact)
            if not hasattr(cls, pname):
                setattr(cls, pname, make_property(fact))
