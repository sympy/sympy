"""
This module contains the machinery handling assumptions.

All symbolic objects have assumption attributes that can be accessed via
.is_<assumption name> attribute.

Assumptions determine certain properties of symbolic objects. Assumptions
can have 3 possible values: True, False, None.  None is returned when it is
impossible to say something about the property. For example, a generic Symbol
is not known beforehand to be positive.

By default, all symbolic values are in the largest set in the given context
without specifying the property. For example, a symbol that has a property
being integer, is also real, complex, etc.

Here follows a list of possible assumption names:

    - commutative   - object commutes with any other object with
                        respect to multiplication operation.
    - real          - object can have only values from the set
                        of real numbers
    - integer       - object can have only values from the set
                        of integers
    - bounded       - object absolute value is bounded
    - positive      - object can have only positive values
    - negative      - object can have only negative values
    - nonpositive      - object can have only nonpositive values
    - nonnegative      - object can have only nonnegative values
    - irrational    - object value cannot be represented exactly by Rational
    - unbounded     - object value is arbitrarily large
    - infinitesimal - object value is infinitesimal


Implementation note: assumption values are stored in
._assumptions dictionary or are returned by getter methods (with
property decorators) or are attributes of objects/classes.

Examples
========

    >>> from sympy import Symbol
    >>> Symbol('x', real = True)
    x

"""
from __future__ import print_function, division

from sympy.core.facts import FactRules, FactKB
from sympy.core.core import BasicMeta
from sympy.core.compatibility import integer_types, with_metaclass

# This are the rules under which our assumptions function
#
# References
# ----------
#
# negative,     -- http://en.wikipedia.org/wiki/Negative_number
# nonnegative
#
# even, odd     -- http://en.wikipedia.org/wiki/Parity_(mathematics)
# imaginary     -- http://en.wikipedia.org/wiki/Imaginary_number
# composite     -- http://en.wikipedia.org/wiki/Composite_number
# finite        -- http://en.wikipedia.org/wiki/Finite
# infinitesimal -- http://en.wikipedia.org/wiki/Infinitesimal
# irrational    -- http://en.wikipedia.org/wiki/Irrational_number
# ...

_assume_rules = FactRules([

    'integer        ->  rational',
    'rational       ->  real',
    'real           ->  complex',
    'real           ->  hermitian',
    'imaginary      ->  complex',
    'imaginary      ->  antihermitian',
    'complex        ->  commutative',

    'odd            ==  integer & !even',
    'even           ==  integer & !odd',

    'real           ==  negative | zero | positive',

    'negative       ==  nonpositive & nonzero',
    'positive       ==  nonnegative & nonzero',
    'zero           ==  nonnegative & nonpositive',

    'nonpositive    ==  real & !positive',
    'nonnegative    ==  real & !negative',

    'zero           ->  infinitesimal & even',

    'prime          ->  integer & positive',
    'composite      ==  integer & positive & !prime',

    'irrational     ==  real & !rational',

    'imaginary      ->  !real',


    '!bounded     ==  unbounded',
    'noninteger     ==  real & !integer',
    '!zero        ==  nonzero',

    # XXX do we need this ?
    'finite     ->  bounded',       # XXX do we need this?
    'finite     ->  !zero',         # XXX wrong?
    'infinitesimal ->  !finite',    # XXX is this ok?
])

_assume_defined = _assume_rules.defined_facts.copy()
_assume_defined.add('polar')
_assume_defined = frozenset(_assume_defined)


class StdFactKB(FactKB):
    """A FactKB specialised for the built-in rules

    This is the only kind of FactKB that Basic objects should use.
    """
    rules = _assume_rules

    def __init__(self, facts=None):
        if facts:
            self.deduce_all_facts(facts)

    def copy(self):
        return self.__class__(self)


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
    assumptions = obj._assumptions
    handler_map = obj._prop_handler

    # Store None into the assumptions so that recursive attempts at
    # evaluating the same fact don't trigger infinite recursion.
    assumptions._tell(fact, None)

    # First try the assumption evaluation function if it exists
    try:
        evaluate = handler_map[fact]
    except KeyError:
        pass
    else:
        a = evaluate(obj)
        if a is not None:
            assumptions.deduce_all_facts(((fact, a),))
            return a

    # Try assumption's prerequisites
    for pk in _assume_rules.prereq[fact]:
        if pk in assumptions:
            continue
        if pk in handler_map:
            _ask(pk, obj)

            # we might have found the value of fact
            ret_val = assumptions.get(fact)
            if ret_val is not None:
                return ret_val

    # Note: the result has already been cached
    return None


class ManagedProperties(with_metaclass(BasicMeta, BasicMeta)):
    """Metaclass for classes with old-style assumptions"""
    def __init__(cls, *args, **kws):
        BasicMeta.__init__(cls, *args, **kws)

        local_defs = {}
        for k in _assume_defined:
            attrname = as_property(k)
            v = cls.__dict__.get(attrname, '')
            if isinstance(v, (bool, integer_types, type(None))):
                if v is not None:
                    v = bool(v)
                local_defs[k] = v

        defs = {}
        for base in reversed(cls.__bases__):
            try:
                defs.update(base._explicit_class_assumptions)
            except AttributeError:
                pass
        defs.update(local_defs)

        cls._explicit_class_assumptions = defs
        cls.default_assumptions = StdFactKB(defs)

        cls._prop_handler = {}
        for k in _assume_defined:
            try:
                cls._prop_handler[k] = getattr(cls, '_eval_is_%s' % k)
            except AttributeError:
                pass

        # Put definite results directly into the class dict, for speed
        for k, v in cls.default_assumptions.items():
            setattr(cls, as_property(k), v)

        # protection e.g. for Integer.is_even=F <- (Rational.is_integer=F)
        derived_from_bases = set()
        for base in cls.__bases__:
            try:
                derived_from_bases |= set(base.default_assumptions)
            except AttributeError:
                continue  # not an assumption-aware class
        for fact in derived_from_bases - set(cls.default_assumptions):
            pname = as_property(fact)
            if pname not in cls.__dict__:
                setattr(cls, pname, make_property(fact))

        # Finally, add any missing automagic property (e.g. for Basic)
        for fact in _assume_defined:
            pname = as_property(fact)
            if not hasattr(cls, pname):
                setattr(cls, pname, make_property(fact))
