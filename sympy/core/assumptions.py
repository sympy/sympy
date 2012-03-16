from sympy.core.facts import FactRules, FactKB
from sympy.core.core import BasicMeta

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
    'imaginary      ->  complex',
    'complex        ->  commutative',

    'odd            ==  integer & !even',
    'even           ==  integer & !odd',

    'real           ==  negative | zero | positive',

    'positive       ->  real & !negative & !zero',
    'negative       ->  real & !positive & !zero',

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

class PropertyManager(StdFactKB):
    """This object is responsible for ensuring the consistency of class
    properties obeying a set of logic rules."""
    def __init__(self, cls):
        local_defs = {}
        for k in _assume_defined:
            attrname = as_property(k)
            v = cls.__dict__.get(attrname, '')
            if isinstance(v, (bool, int, long, type(None))):
                if v is not None:
                    v = bool(v)
                local_defs[k] = v

        defs = {}
        for base in reversed(cls.__bases__):
            try:
                defs.update(base.default_assumptions.definitions)
            except AttributeError:
                pass
        defs.update(local_defs)

        self.definitions = defs
        self.deduce_all_facts(defs)

    def copy(self):
        return StdFactKB(self)

    def make_property(manager, fact):
        """Create the automagic property corresponding to a fact."""

        def getit(self):
            try:
                return self._assumptions[fact]
            except KeyError:
                if self._assumptions is self.default_assumptions:
                    self._assumptions = self.default_assumptions.copy()
                return self._what_known_about(fact)

        getit.func_name = as_property(fact)
        return property(getit)


class WithAssumptions(BasicMeta):
    """Metaclass for classes with old-style assumptions"""
    __metaclass__ = BasicMeta

    def __new__(mcl, name, bases, attrdict):
        if not any(issubclass(base, AssumeMixin) for base in bases):
            bases = (AssumeMixin,) + bases
            if '__slots__' in attrdict:
                attrdict['__slots__'] += AssumeMixin._assume_slots
        return super(WithAssumptions, mcl).__new__(mcl, name, bases, attrdict)

    def __init__(cls, *args, **kws):
        BasicMeta.__init__(cls, *args, **kws)

        pmanager = PropertyManager(cls)
        cls.default_assumptions = pmanager

        # Put definite results directly into the class dict, for speed
        for k, v in pmanager.iteritems():
            setattr(cls, as_property(k), v)

        # protection e.g. for Integer.is_even=F <- (Rational.is_integer=F)
        derived_from_bases = set()
        for base in cls.__bases__:
            try:
                derived_from_bases |= set(base.default_assumptions)
            except AttributeError:
                continue        #not an assumption-aware class
        for fact in derived_from_bases - set(pmanager):
            pname = as_property(fact)
            if pname not in cls.__dict__:
                setattr(cls, pname, pmanager.make_property(fact))

        # Finally, add any missing automagic property (e.g. for Basic)
        for fact in _assume_defined:
            pname = as_property(fact)
            if not hasattr(cls, pname):
                setattr(cls, pname, pmanager.make_property(fact))


class AssumeMixin(object):
    """ Define default assumption methods.

    AssumeMeths should be used to derive Basic class only.

    All symbolic objects have assumption attributes that can be accessed via
    .is_<assumption name> attribute.

    Assumptions determine certain properties of symbolic objects. Assumptions
    can have 3 possible values: True, False, None.  None is returned when it is
    impossible to say something about the property. For example, a Symbol is
    not know beforehand to be positive.

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
        - comparable    - object.evalf() returns Number object.
        - irrational    - object value cannot be represented exactly by Rational
        - unbounded     - object value is arbitrarily large
        - infinitesimal - object value is infinitesimal


    Example rules:

      positive=T            ->  nonpositive=F, real=T
      real=T & positive=F   ->  nonpositive=T

      unbounded=F|T         ->  bounded=not unbounded   XXX ok?
      irrational=T          ->  real=T


    Implementation note: assumption values are stored in
    ._assumption dictionary or are returned by getter methods (with
    property decorators) or are attributes of objects/classes.

    Examples
    ========

        - True, when we are sure about a property. For example, when we are
        working only with real numbers:
        >>> from sympy import Symbol
        >>> Symbol('x', real = True)
        x

        - False

        - None (if you don't know if the property is True or false)
    """
    _assume_slots = ['_assumptions']
    try:
        # This particular __slots__ definition breaks SymPy in Jython.
        # See issue 1233.
        import java
    except ImportError:
        __slots__ = []

    def  _init_assumptions(self):
        self._assumptions  = self.default_assumptions

    # XXX better name?
    @property
    def assumptions0(self):
        """
        Return object `type` assumptions.

        For example:

          Symbol('x', real=True)
          Symbol('x', integer=True)

        are different objects. In other words, besides Python type (Symbol in
        this case), the initial assumptions are also forming their typeinfo.

        Examples
        ========

        >>> from sympy import Symbol
        >>> from sympy.abc import x
        >>> x.assumptions0
        {'commutative': True}
        >>> x = Symbol("x", positive=True)
        >>> x.assumptions0
        {'commutative': True, 'complex': True, 'imaginary': False,
        'negative': False, 'nonnegative': True, 'nonpositive': False,
        'nonzero': True, 'positive': True, 'real': True, 'zero': False}

        """
        return {}

    def _what_known_about(self, k):
        """tries hard to give an answer to: what is known about fact `k`

           NOTE: You should not use this directly -- see make__get_assumption
                 instead

           This function is called when a request is made to see what a fact
           value is.

           If we are here, it means that the asked-for fact is not known, and
           we should try to find a way to find its value.

           For this we use several techniques:

           1. _eval_is_<fact>
           ------------------

           first fact-evaluation function is tried,  for example
           _eval_is_integer


           2. relations
           ------------

           if the first step did not succeeded (no such function, or its return
           is None) then we try related facts. For example

                       means
             rational   -->   integer

           another example is joined rule:

             integer & !odd  --> even

           so in the latter case if we are looking at what 'even' value is,
           'integer' and 'odd' facts will be asked.

           In all cases when we settle on some fact value, it is given to
           _learn_new_facts to deduce all its implications, and also the result
           is cached in ._assumptions for later quick access.
        """
        assumptions = self._assumptions

        # Store None into the assumptions so that recursive attempts at
        # evaluating the same fact don't trigger infinite recursion.
        assumptions.deduce_all_facts(((k, None),))

        # First try the assumption evaluation function if it exists
        try:
            a = getattr(self, '_eval_is_' + k)()
        except AttributeError:
            pass
        else:
            if a is not None:
                assumptions.deduce_all_facts(((k, a),))
                return a

        # Try assumption's prerequisites
        for pk in _assume_rules.prereq[k]:
            if pk in assumptions:
                continue
            if hasattr(self, '_eval_is_' + pk):
                self._what_known_about(pk)

                # we might have found the value of k
                ret_val = assumptions.get(k)
                if ret_val is not None:
                    return ret_val

        # Note: the result has already been cached
        return None
