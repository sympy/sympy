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


def _load_pre_generated_assumption_rules():
    """ Load the assumption rules from pre-generated data

    To update the pre-generated data, see :method::`_generate_assumption_rules`
    """
    _assume_rules=FactRules._from_python(_pre_calculated_assumptions)
    return _assume_rules

def _generate_assumption_rules():
    """ Generate the default assumption rules

    This method should only be called to update the pre-generated assumption rules.
    These are stored in the variable `_pre_calculated_assumptions`
    To update the pre-generated assumptions run:

    >>> from sympy.core.assumptions import _generate_assumption_rules
    >>> generated_assumptions = _generate_assumption_rules()
    >>> output = str(generated_assumptions._to_python())

    And copy the generated output to `_pre_calculated_assumptions`

    """
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
    return _assume_rules

_pre_calculated_assumptions = {'full_implications': {('extended_real', True): {('commutative', True), ('imaginary', False)}, ('extended_nonzero', True): {('zero', False), ('commutative', True), ('imaginary', False), ('extended_real', True)}, ('rational', True): {('hermitian', True), ('complex', True), ('real', True), ('finite', True), ('algebraic', True), ('irrational', False), ('commutative', True), ('imaginary', False), ('extended_real', True), ('infinite', False), ('transcendental', False)}, ('imaginary', True): {('complex', True), ('nonpositive', False), ('real', False), ('extended_nonpositive', False), ('zero', False), ('negative', False), ('extended_negative', False), ('odd', False), ('commutative', True), ('extended_nonnegative', False), ('prime', False), ('positive', False), ('extended_positive', False), ('antihermitian', True), ('irrational', False), ('infinite', False), ('rational', False), ('even', False), ('finite', True), ('extended_nonzero', False), ('noninteger', False), ('integer', False), ('nonnegative', False), ('extended_real', False), ('composite', False), ('nonzero', False)}, ('extended_nonpositive', False): {('extended_negative', False), ('zero', False), ('nonpositive', False), ('negative', False)}, ('noninteger', True): {('even', False), ('zero', False), ('extended_nonzero', True), ('commutative', True), ('imaginary', False), ('extended_real', True), ('prime', False), ('integer', False), ('odd', False), ('composite', False)}, ('composite', True): {('hermitian', True), ('complex', True), ('extended_nonnegative', True), ('nonpositive', False), ('positive', True), ('extended_nonpositive', False), ('zero', False), ('extended_positive', True), ('negative', False), ('extended_negative', False), ('rational', True), ('commutative', True), ('prime', False), ('transcendental', False), ('nonnegative', True), ('extended_nonzero', True), ('irrational', False), ('integer', True), ('extended_real', True), ('nonzero', True), ('infinite', False), ('real', True), ('finite', True), ('algebraic', True), ('imaginary', False), ('noninteger', False)}, ('odd', True): {('hermitian', True), ('complex', True), ('real', True), ('finite', True), ('even', False), ('extended_nonzero', True), ('algebraic', True), ('irrational', False), ('integer', True), ('imaginary', False), ('rational', True), ('extended_real', True), ('commutative', True), ('nonzero', True), ('noninteger', False), ('zero', False), ('infinite', False), ('transcendental', False)}, ('finite', False): {('nonpositive', False), ('real', False), ('algebraic', False), ('zero', False), ('negative', False), ('complex', False), ('infinite', True), ('prime', False), ('positive', False), ('transcendental', False), ('irrational', False), ('rational', False), ('even', False), ('imaginary', False), ('integer', False), ('nonnegative', False), ('odd', False), ('composite', False), ('nonzero', False)}, ('positive', True): {('hermitian', True), ('complex', True), ('extended_positive', True), ('finite', True), ('real', True), ('negative', False), ('extended_nonnegative', True), ('extended_negative', False), ('nonnegative', True), ('extended_nonzero', True), ('nonpositive', False), ('extended_nonpositive', False), ('commutative', True), ('imaginary', False), ('extended_real', True), ('nonzero', True), ('zero', False), ('infinite', False)}, ('real', True): {('hermitian', True), ('complex', True), ('finite', True), ('commutative', True), ('imaginary', False), ('extended_real', True), ('infinite', False)}, ('even', True): {('hermitian', True), ('complex', True), ('real', True), ('finite', True), ('algebraic', True), ('irrational', False), ('integer', True), ('rational', True), ('commutative', True), ('extended_real', True), ('imaginary', False), ('noninteger', False), ('infinite', False), ('odd', False), ('transcendental', False)}, ('negative', True): {('hermitian', True), ('complex', True), ('nonnegative', False), ('zero', False), ('commutative', True), ('extended_nonnegative', False), ('prime', False), ('positive', False), ('extended_positive', False), ('extended_nonzero', True), ('extended_real', True), ('nonzero', True), ('infinite', False), ('real', True), ('finite', True), ('nonpositive', True), ('extended_nonpositive', True), ('imaginary', False), ('extended_negative', True), ('composite', False)}, ('nonpositive', True): {('hermitian', True), ('complex', True), ('real', True), ('finite', True), ('extended_positive', False), ('extended_nonpositive', True), ('commutative', True), ('imaginary', False), ('extended_real', True), ('prime', False), ('infinite', False), ('positive', False), ('composite', False)}, ('hermitian', False): {('rational', False), ('even', False), ('negative', False), ('zero', False), ('positive', False), ('nonpositive', False), ('real', False), ('irrational', False), ('prime', False), ('integer', False), ('nonnegative', False), ('odd', False), ('composite', False), ('nonzero', False)}, ('integer', True): {('hermitian', True), ('complex', True), ('real', True), ('finite', True), ('algebraic', True), ('irrational', False), ('rational', True), ('commutative', True), ('imaginary', False), ('extended_real', True), ('noninteger', False), ('infinite', False), ('transcendental', False)}, ('zero', True): {('hermitian', True), ('complex', True), ('extended_nonnegative', True), ('negative', False), ('extended_negative', False), ('rational', True), ('commutative', True), ('prime', False), ('even', True), ('positive', False), ('transcendental', False), ('extended_positive', False), ('nonnegative', True), ('irrational', False), ('integer', True), ('extended_real', True), ('infinite', False), ('real', True), ('finite', True), ('nonpositive', True), ('extended_nonpositive', True), ('algebraic', True), ('imaginary', False), ('extended_nonzero', False), ('noninteger', False), ('odd', False), ('composite', False), ('nonzero', False)}, ('rational', False): {('even', False), ('zero', False), ('prime', False), ('integer', False), ('odd', False), ('composite', False)}, ('extended_nonzero', False): {('negative', False), ('extended_positive', False), ('extended_negative', False), ('prime', False), ('positive', False), ('composite', False), ('nonzero', False)}, ('prime', True): {('hermitian', True), ('complex', True), ('extended_nonnegative', True), ('nonpositive', False), ('positive', True), ('extended_nonpositive', False), ('zero', False), ('extended_positive', True), ('negative', False), ('extended_negative', False), ('rational', True), ('commutative', True), ('transcendental', False), ('nonnegative', True), ('extended_nonzero', True), ('irrational', False), ('integer', True), ('extended_real', True), ('nonzero', True), ('infinite', False), ('real', True), ('finite', True), ('algebraic', True), ('imaginary', False), ('noninteger', False), ('composite', False)}, ('real', False): {('rational', False), ('even', False), ('negative', False), ('zero', False), ('positive', False), ('nonpositive', False), ('irrational', False), ('prime', False), ('integer', False), ('nonnegative', False), ('odd', False), ('composite', False), ('nonzero', False)}, ('nonnegative', False): {('zero', False), ('positive', False), ('composite', False), ('prime', False)}, ('infinite', True): {('nonpositive', False), ('real', False), ('finite', False), ('algebraic', False), ('zero', False), ('negative', False), ('complex', False), ('prime', False), ('positive', False), ('transcendental', False), ('irrational', False), ('rational', False), ('even', False), ('imaginary', False), ('integer', False), ('nonnegative', False), ('odd', False), ('composite', False), ('nonzero', False)}, ('irrational', True): {('hermitian', True), ('rational', False), ('complex', True), ('real', True), ('finite', True), ('even', False), ('zero', False), ('noninteger', True), ('extended_nonzero', True), ('commutative', True), ('imaginary', False), ('extended_real', True), ('prime', False), ('nonzero', True), ('integer', False), ('infinite', False), ('odd', False), ('composite', False)}, ('extended_real', False): {('nonpositive', False), ('real', False), ('extended_nonpositive', False), ('zero', False), ('negative', False), ('extended_negative', False), ('extended_nonnegative', False), ('prime', False), ('positive', False), ('extended_positive', False), ('irrational', False), ('rational', False), ('even', False), ('extended_nonzero', False), ('noninteger', False), ('integer', False), ('nonnegative', False), ('odd', False), ('composite', False), ('nonzero', False)}, ('integer', False): {('even', False), ('prime', False), ('zero', False), ('odd', False), ('composite', False)}, ('extended_positive', True): {('negative', False), ('extended_nonnegative', True), ('extended_negative', False), ('nonpositive', False), ('extended_nonzero', True), ('extended_nonpositive', False), ('commutative', True), ('imaginary', False), ('extended_real', True), ('zero', False)}, ('commutative', False): {('nonpositive', False), ('real', False), ('extended_nonpositive', False), ('algebraic', False), ('zero', False), ('negative', False), ('extended_negative', False), ('complex', False), ('odd', False), ('extended_nonnegative', False), ('prime', False), ('positive', False), ('transcendental', False), ('extended_positive', False), ('irrational', False), ('rational', False), ('even', False), ('imaginary', False), ('extended_nonzero', False), ('noninteger', False), ('integer', False), ('nonnegative', False), ('extended_real', False), ('composite', False), ('nonzero', False)}, ('extended_positive', False): {('positive', False), ('composite', False), ('prime', False)}, ('algebraic', False): {('rational', False), ('even', False), ('zero', False), ('prime', False), ('integer', False), ('odd', False), ('composite', False)}, ('finite', True): {('infinite', False)}, ('algebraic', True): {('complex', True), ('finite', True), ('commutative', True), ('infinite', False), ('transcendental', False)}, ('extended_nonnegative', True): {('negative', False), ('extended_negative', False), ('commutative', True), ('imaginary', False), ('extended_real', True)}, ('complex', False): {('rational', False), ('even', False), ('transcendental', False), ('negative', False), ('zero', False), ('positive', False), ('nonpositive', False), ('real', False), ('irrational', False), ('imaginary', False), ('prime', False), ('algebraic', False), ('integer', False), ('nonnegative', False), ('odd', False), ('composite', False), ('nonzero', False)}, ('extended_nonpositive', True): {('extended_positive', False), ('commutative', True), ('imaginary', False), ('extended_real', True), ('prime', False), ('positive', False), ('composite', False)}, ('transcendental', True): {('rational', False), ('complex', True), ('even', False), ('finite', True), ('zero', False), ('commutative', True), ('prime', False), ('algebraic', False), ('integer', False), ('infinite', False), ('odd', False), ('composite', False)}, ('extended_negative', True): {('extended_positive', False), ('extended_nonpositive', True), ('extended_nonzero', True), ('commutative', True), ('imaginary', False), ('extended_real', True), ('prime', False), ('extended_nonnegative', False), ('zero', False), ('nonnegative', False), ('positive', False), ('composite', False)}, ('nonzero', True): {('hermitian', True), ('complex', True), ('real', True), ('finite', True), ('extended_nonzero', True), ('commutative', True), ('imaginary', False), ('extended_real', True), ('zero', False), ('infinite', False)}, ('nonnegative', True): {('hermitian', True), ('complex', True), ('real', True), ('finite', True), ('negative', False), ('extended_nonnegative', True), ('extended_negative', False), ('commutative', True), ('imaginary', False), ('extended_real', True), ('infinite', False)}, ('infinite', False): {('finite', True)}, ('extended_nonnegative', False): {('extended_positive', False), ('prime', False), ('zero', False), ('nonnegative', False), ('positive', False), ('composite', False)}, ('even', False): {('zero', False)}, ('antihermitian', False): {('imaginary', False)}, ('extended_negative', False): {('negative', False)}, ('complex', True): {('infinite', False), ('finite', True), ('commutative', True)}, ('nonzero', False): {('positive', False), ('composite', False), ('prime', False), ('negative', False)}, ('positive', False): {('composite', False), ('prime', False)}, ('nonpositive', False): {('zero', False), ('negative', False)}, ('odd', False): set(), ('zero', False): set(), ('negative', False): set(), ('prime', False): set(), ('composite', False): set()}, 'beta_triggers': {('extended_real', True): [4, 5, 6, 8, 9, 10, 13, 14, 35, 36], ('extended_nonzero', True): [4, 5, 6, 8, 10, 11, 12, 13, 14, 15, 16, 18, 27, 35], ('rational', True): [8, 9, 10, 13, 14, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 35, 36], ('imaginary', True): [0], ('extended_nonpositive', False): [7, 10, 14, 15, 18, 22, 36], ('noninteger', True): [4, 5, 6, 8, 10, 11, 12, 13, 14, 15, 16, 18, 27], ('composite', True): [1, 2], ('odd', True): [8, 10, 11, 12, 13, 14, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26], ('finite', False): [7, 8, 10, 35, 36], ('positive', True): [0, 30, 31, 33, 34, 35], ('real', True): [0, 8, 9, 10, 13, 14, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 34, 35, 36], ('even', True): [8, 9, 10, 13, 14, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 30, 32, 33, 36], ('negative', True): [0, 34, 35], ('nonpositive', True): [0, 8, 9, 11, 14, 16, 17, 19, 22, 24, 26, 27, 28, 29, 34, 35, 36], ('hermitian', False): [3, 5, 7, 8, 10, 35, 36], ('integer', True): [1, 2, 8, 9, 10, 13, 14, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 36], ('zero', True): [], ('rational', False): [7, 8, 10, 15, 16, 18, 34, 35, 36], ('extended_nonzero', False): [7, 9, 13, 14, 15, 17, 21, 22], ('prime', True): [1, 2], ('real', False): [3, 5, 7, 8, 10, 35, 36], ('nonnegative', False): [7, 8, 10, 15, 16, 21, 36], ('infinite', True): [7, 8, 10, 35, 36], ('irrational', True): [0, 8, 10, 11, 12, 13, 14, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26], ('extended_real', False): [], ('integer', False): [7, 8, 10, 15, 16, 18, 35, 36], ('extended_positive', True): [4, 5, 6, 15, 18, 22, 23, 26, 27, 35], ('commutative', False): [], ('extended_positive', False): [7, 8, 9, 13, 15, 16, 17, 21], ('algebraic', False): [0, 7, 8, 10, 15, 16, 18, 34, 35, 36], ('finite', True): [3, 4, 6, 23, 24, 25, 26, 27], ('algebraic', True): [3, 4, 6, 23, 24, 25, 26, 27], ('extended_nonnegative', True): [4, 5, 6, 9, 10, 12, 13, 15, 17, 18, 22, 26, 28, 35, 36], ('complex', False): [3, 5, 7, 8, 10, 35, 36], ('extended_nonpositive', True): [4, 5, 6, 8, 9, 11, 14, 15, 16, 17, 21, 25, 28, 35, 36], ('transcendental', True): [3, 4, 6, 7, 8, 10, 15, 16, 18, 23, 24, 25, 26, 27, 34, 35, 36], ('extended_negative', True): [4, 5, 6, 15, 16, 21, 24, 25, 27, 35], ('nonzero', True): [0, 8, 10, 11, 12, 13, 14, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 34, 35], ('nonnegative', True): [0, 9, 10, 12, 13, 17, 18, 20, 21, 23, 25, 27, 28, 29, 34, 35, 36], ('infinite', False): [3, 4, 6, 23, 24, 25, 26, 27], ('extended_nonnegative', False): [7, 8, 13, 15, 16, 21, 36], ('even', False): [1, 7, 8, 10, 15, 16, 18, 36], ('antihermitian', False): [], ('extended_negative', False): [7, 9, 10, 14, 15, 17, 18, 22], ('complex', True): [0, 3, 4, 6, 23, 24, 25, 26, 27], ('nonzero', False): [15, 17, 21, 22], ('positive', False): [15, 16, 17, 21], ('nonpositive', False): [7, 8, 10, 15, 18, 22, 36], ('odd', False): [2], ('zero', False): [7, 8, 10, 15, 16, 18, 36], ('negative', False): [15, 17, 18, 22], ('prime', False): [30, 31, 32], ('composite', False): [31, 32, 33]}, 'prereq': {'commutative': {'extended_nonzero', 'extended_real', 'nonnegative', 'complex', 'algebraic', 'extended_nonnegative', 'integer', 'extended_positive', 'imaginary', 'real', 'zero', 'negative', 'prime', 'nonpositive', 'composite', 'odd', 'extended_negative', 'nonzero', 'even', 'rational', 'noninteger', 'extended_nonpositive', 'positive', 'transcendental', 'irrational'}, 'imaginary': {'extended_nonzero', 'extended_real', 'complex', 'nonnegative', 'finite', 'antihermitian', 'extended_nonnegative', 'integer', 'extended_positive', 'real', 'negative', 'prime', 'nonpositive', 'composite', 'odd', 'extended_negative', 'infinite', 'nonzero', 'rational', 'commutative', 'even', 'noninteger', 'extended_nonpositive', 'positive', 'zero', 'irrational'}, 'zero': {'extended_nonzero', 'extended_real', 'nonnegative', 'complex', 'finite', 'algebraic', 'extended_nonnegative', 'integer', 'extended_positive', 'imaginary', 'real', 'negative', 'prime', 'nonpositive', 'composite', 'odd', 'extended_negative', 'infinite', 'nonzero', 'rational', 'commutative', 'hermitian', 'even', 'noninteger', 'extended_nonpositive', 'positive', 'transcendental', 'irrational'}, 'extended_real': {'extended_nonzero', 'nonnegative', 'extended_nonnegative', 'integer', 'extended_positive', 'imaginary', 'real', 'negative', 'prime', 'nonpositive', 'composite', 'odd', 'extended_negative', 'nonzero', 'even', 'rational', 'commutative', 'noninteger', 'extended_nonpositive', 'positive', 'zero', 'irrational'}, 'hermitian': {'negative', 'prime', 'nonpositive', 'nonnegative', 'odd', 'composite', 'nonzero', 'even', 'rational', 'integer', 'positive', 'zero', 'real', 'irrational'}, 'complex': {'nonnegative', 'transcendental', 'finite', 'algebraic', 'integer', 'imaginary', 'real', 'negative', 'prime', 'nonpositive', 'composite', 'odd', 'nonzero', 'infinite', 'even', 'rational', 'commutative', 'positive', 'zero', 'irrational'}, 'real': {'extended_real', 'complex', 'nonnegative', 'finite', 'integer', 'imaginary', 'negative', 'prime', 'nonpositive', 'composite', 'odd', 'nonzero', 'infinite', 'even', 'rational', 'commutative', 'hermitian', 'positive', 'zero', 'irrational'}, 'finite': {'nonnegative', 'complex', 'transcendental', 'algebraic', 'integer', 'imaginary', 'real', 'negative', 'prime', 'nonpositive', 'composite', 'odd', 'nonzero', 'infinite', 'even', 'rational', 'positive', 'zero', 'irrational'}, 'algebraic': {'prime', 'complex', 'transcendental', 'finite', 'odd', 'composite', 'infinite', 'even', 'rational', 'integer', 'commutative', 'zero'}, 'irrational': {'prime', 'extended_real', 'complex', 'finite', 'odd', 'composite', 'infinite', 'even', 'rational', 'integer', 'hermitian', 'commutative', 'imaginary', 'zero', 'real'}, 'infinite': {'nonnegative', 'complex', 'transcendental', 'finite', 'algebraic', 'integer', 'imaginary', 'real', 'negative', 'prime', 'nonpositive', 'composite', 'odd', 'nonzero', 'even', 'rational', 'positive', 'zero', 'irrational'}, 'transcendental': {'prime', 'complex', 'finite', 'algebraic', 'odd', 'composite', 'infinite', 'even', 'rational', 'integer', 'commutative', 'zero'}, 'nonpositive': {'negative', 'prime', 'extended_real', 'complex', 'finite', 'positive', 'composite', 'infinite', 'commutative', 'hermitian', 'extended_positive', 'extended_nonpositive', 'imaginary', 'zero', 'real'}, 'extended_nonpositive': {'negative', 'prime', 'nonpositive', 'extended_real', 'positive', 'composite', 'extended_negative', 'commutative', 'extended_positive', 'imaginary', 'zero'}, 'negative': {'extended_nonzero', 'extended_real', 'complex', 'nonnegative', 'finite', 'extended_nonnegative', 'extended_positive', 'imaginary', 'real', 'prime', 'nonpositive', 'composite', 'extended_negative', 'nonzero', 'infinite', 'commutative', 'hermitian', 'extended_nonpositive', 'positive', 'zero'}, 'extended_negative': {'extended_nonzero', 'negative', 'prime', 'extended_real', 'nonnegative', 'positive', 'extended_nonnegative', 'composite', 'commutative', 'extended_positive', 'extended_nonpositive', 'imaginary', 'zero'}, 'odd': {'extended_real', 'complex', 'finite', 'algebraic', 'integer', 'imaginary', 'real', 'zero', 'infinite', 'even', 'rational', 'commutative', 'hermitian', 'noninteger', 'transcendental', 'irrational'}, 'extended_nonnegative': {'negative', 'prime', 'extended_real', 'nonnegative', 'positive', 'composite', 'extended_negative', 'commutative', 'extended_positive', 'imaginary', 'zero'}, 'prime': {'extended_nonzero', 'extended_real', 'nonnegative', 'complex', 'finite', 'algebraic', 'extended_nonnegative', 'integer', 'extended_positive', 'imaginary', 'real', 'zero', 'negative', 'nonpositive', 'composite', 'extended_negative', 'nonzero', 'infinite', 'rational', 'commutative', 'hermitian', 'noninteger', 'extended_nonpositive', 'positive', 'transcendental', 'irrational'}, 'positive': {'extended_nonzero', 'extended_real', 'nonnegative', 'complex', 'finite', 'extended_nonnegative', 'extended_positive', 'imaginary', 'real', 'negative', 'prime', 'nonpositive', 'composite', 'extended_negative', 'nonzero', 'infinite', 'commutative', 'hermitian', 'extended_nonpositive', 'zero'}, 'extended_positive': {'extended_nonzero', 'negative', 'prime', 'nonpositive', 'extended_real', 'positive', 'extended_nonnegative', 'composite', 'extended_negative', 'commutative', 'extended_nonpositive', 'imaginary', 'zero'}, 'antihermitian': {'imaginary'}, 'rational': {'extended_real', 'complex', 'finite', 'algebraic', 'integer', 'imaginary', 'real', 'zero', 'prime', 'odd', 'composite', 'infinite', 'even', 'commutative', 'hermitian', 'transcendental', 'irrational'}, 'even': {'extended_real', 'complex', 'finite', 'algebraic', 'integer', 'imaginary', 'real', 'zero', 'odd', 'infinite', 'rational', 'commutative', 'hermitian', 'noninteger', 'transcendental', 'irrational'}, 'extended_nonzero': {'negative', 'prime', 'extended_real', 'positive', 'composite', 'odd', 'extended_negative', 'nonzero', 'commutative', 'extended_positive', 'noninteger', 'imaginary', 'zero', 'irrational'}, 'noninteger': {'prime', 'extended_real', 'odd', 'composite', 'integer', 'even', 'commutative', 'imaginary', 'zero', 'irrational'}, 'integer': {'extended_real', 'complex', 'finite', 'algebraic', 'imaginary', 'real', 'zero', 'prime', 'composite', 'odd', 'infinite', 'even', 'rational', 'commutative', 'hermitian', 'noninteger', 'transcendental', 'irrational'}, 'nonnegative': {'negative', 'prime', 'extended_real', 'complex', 'finite', 'positive', 'extended_nonnegative', 'composite', 'extended_negative', 'infinite', 'commutative', 'hermitian', 'imaginary', 'zero', 'real'}, 'composite': {'extended_nonzero', 'extended_real', 'nonnegative', 'complex', 'finite', 'algebraic', 'extended_nonnegative', 'integer', 'extended_positive', 'imaginary', 'real', 'zero', 'negative', 'prime', 'nonpositive', 'extended_negative', 'nonzero', 'infinite', 'rational', 'commutative', 'hermitian', 'noninteger', 'extended_nonpositive', 'positive', 'transcendental', 'irrational'}, 'nonzero': {'extended_nonzero', 'extended_real', 'complex', 'finite', 'imaginary', 'real', 'negative', 'prime', 'odd', 'composite', 'infinite', 'commutative', 'hermitian', 'positive', 'zero', 'irrational'}}, 'beta_rules': [({('algebraic', False), ('complex', True)}, ('transcendental', True)), ({('integer', True), ('even', False)}, ('odd', True)), ({('integer', True), ('odd', False)}, ('even', True)), ({('real', False), ('infinite', False)}, ('extended_real', False)), ({('infinite', False), ('extended_real', True)}, ('real', True)), ({('real', False), ('extended_real', True)}, ('infinite', True)), ({('finite', True), ('extended_real', True)}, ('real', True)), ({('extended_negative', False), ('zero', False), ('extended_positive', False)}, ('extended_real', False)), ({('extended_positive', False), ('zero', False), ('extended_real', True)}, ('extended_negative', True)), ({('extended_negative', False), ('extended_positive', False), ('extended_real', True)}, ('zero', True)), ({('extended_negative', False), ('zero', False), ('extended_real', True)}, ('extended_positive', True)), ({('extended_nonpositive', True), ('extended_nonzero', True)}, ('extended_negative', True)), ({('extended_nonnegative', True), ('extended_nonzero', True)}, ('extended_positive', True)), ({('extended_positive', False), ('extended_real', True)}, ('extended_nonpositive', True)), ({('extended_negative', False), ('extended_real', True)}, ('extended_nonnegative', True)), ({('zero', False), ('positive', False), ('negative', False)}, ('real', False)), ({('real', True), ('positive', False), ('zero', False)}, ('negative', True)), ({('real', True), ('positive', False), ('negative', False)}, ('zero', True)), ({('real', True), ('zero', False), ('negative', False)}, ('positive', True)), ({('nonzero', True), ('nonpositive', True)}, ('negative', True)), ({('nonzero', True), ('nonnegative', True)}, ('positive', True)), ({('real', True), ('positive', False)}, ('nonpositive', True)), ({('real', True), ('negative', False)}, ('nonnegative', True)), ({('extended_positive', True), ('finite', True)}, ('positive', True)), ({('extended_negative', True), ('finite', True)}, ('negative', True)), ({('extended_nonpositive', True), ('finite', True)}, ('nonpositive', True)), ({('extended_nonnegative', True), ('finite', True)}, ('nonnegative', True)), ({('finite', True), ('extended_nonzero', True)}, ('nonzero', True)), ({('extended_nonnegative', True), ('extended_nonpositive', True)}, ('zero', True)), ({('nonnegative', True), ('nonpositive', True)}, ('zero', True)), ({('even', True), ('prime', False), ('positive', True)}, ('composite', True)), ({('prime', False), ('composite', False), ('positive', True)}, ('even', False)), ({('even', True), ('prime', False), ('composite', False)}, ('positive', False)), ({('even', True), ('composite', False), ('positive', True)}, ('prime', True)), ({('rational', False), ('real', True)}, ('irrational', True)), ({('integer', False), ('extended_real', True)}, ('noninteger', True)), ({('zero', False), ('extended_real', True)}, ('extended_nonzero', True))], 'defined_facts': ['extended_real', 'finite', 'algebraic', 'extended_nonnegative', 'integer', 'extended_positive', 'imaginary', 'real', 'nonpositive', 'composite', 'rational', 'commutative', 'hermitian', 'noninteger', 'extended_nonpositive', 'extended_nonzero', 'nonnegative', 'complex', 'transcendental', 'antihermitian', 'negative', 'prime', 'odd', 'extended_negative', 'infinite', 'nonzero', 'even', 'positive', 'zero', 'irrational']}
_assume_rules = _load_pre_generated_assumption_rules()

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
