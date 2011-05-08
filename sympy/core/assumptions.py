from sympy.core.compatibility import cmp
from sympy.core.facts import FactRules
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
_assume_defined.add('comparable')
_assume_defined = frozenset(_assume_defined)


###################################
# positive/negative from .evalf() #
###################################

# properties that indicate ordering on real axis
_real_ordering = set(['negative', 'nonnegative', 'positive', 'nonpositive'])

# what can be said from cmp(x.evalf(),0)
# if x.evalf() is zero we can say nothing so nonpositive is the same as negative
_real_cmp0_table= {
        'positive': {1: True,  -1: False, 0: None},
        'negative': {1: False, -1: True,  0: None},
        }
_real_cmp0_table['nonpositive'] = _real_cmp0_table['negative']
_real_cmp0_table['nonnegative'] = _real_cmp0_table['positive']

class CycleDetected(Exception):
    """(internal) used to detect cycles when evaluating assumptions
       through prerequisites
    """
    pass

def as_property(fact):
    """Convert a fact name to the name of the corresponding property"""
    return 'is_%s' % fact

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

        cls_premises = {}
        for k in _assume_defined:
            attrname = as_property(k)
            v = cls.__dict__.get(attrname, '')
            if isinstance(v, (bool, int, long, type(None))):
                if v is not None:
                    v = bool(v)
                cls_premises[k] = v

        cls._default_premises = {}
        for base in reversed(cls.__bases__):
            try:
                cls._default_premises.update(base._default_premises)
            except AttributeError:
                pass
        cls._default_premises.update(cls_premises)

        # deduce all consequences from default assumptions -- make it complete
        default_assumptions = _assume_rules.deduce_all_facts(cls._default_premises)

        # and store completed set into cls -- this way we'll avoid rededucing
        # extensions of class default assumptions each time on instance
        # creation -- we keep it prededuced already.
        for k, v in default_assumptions.iteritems():
            setattr(cls, as_property(k), v)
        cls.default_assumptions = default_assumptions

        # protection e.g. for Initeger.is_even=F <- (Rational.is_integer=F)
        base_derived_premises = set()
        for base in cls.__bases__:
            try:
                base_derived_premises |= (set(base.default_assumptions) -
                                                set(base._default_premises))
            except AttributeError:
                continue        #not an assumption-aware class

        for fact in base_derived_premises:
            if as_property(fact) not in cls.__dict__:
                cls.add_property(fact)

        for fact in _assume_defined:
            if not hasattr(cls, as_property(fact)):
                cls.add_property(fact)

    def add_property(cls, fact):
        """Add to the class the automagic property corresponding to a fact."""

        def getit(self):
            try:
                return self._assumptions[fact]
            except KeyError:
                return self._what_known_about(fact)

        getit.func_name = '%s__is_%s' % (cls.__name__, fact)
        setattr(cls, as_property(fact), property(getit))


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

    Examples:

        - True, when we are sure about a property. For example, when we are
        working only with real numbers:
        >>> from sympy import Symbol
        >>> Symbol('x', real = True)
        x

        - False

        - None (if you don't know if the property is True or false)
    """
    _assume_slots = ['_assumptions',    # assumptions
                 '_a_inprogress',   # already-seen requests (when deducing
                                    # through prerequisites -- see CycleDetected)
                 '_assume_type_keys', # assumptions typeinfo keys
                ]
    __slots__ = []

    def  _init_assumptions(self, assumptions):
        # initially assumptions are shared between instances and class
        self._assumptions  = self.default_assumptions
        self._a_inprogress = []

        # NOTE this could be made lazy -- probably not all instances will need
        # fully derived assumptions?
        if assumptions:
            self._learn_new_facts(assumptions)
            #                      ^
            # FIXME this is slow   |    another NOTE: speeding this up is *not*
            #        |             |    important. say for %timeit x+y most of
            # .------'             |    the time is spent elsewhere
            # |                    |
            # |  XXX _learn_new_facts  could be asked about what *new* facts have
            # v  XXX been learned -- we'll need this to append to _hashable_content
            basek = set(self.default_assumptions.keys())
            k2    = set(self._assumptions.keys())
            newk  = k2.difference(basek)

            self._assume_type_keys = frozenset(newk)
        else:
            self._assume_type_keys = None

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

        Example:

        >>> from sympy import Symbol
        >>> from sympy.abc import x
        >>> x.assumptions0
        {}
        >>> x = Symbol("x", positive=True)
        >>> x.assumptions0
        {'commutative': True, 'complex': True, 'imaginary': False,
        'negative': False, 'nonnegative': True, 'nonpositive': False,
        'nonzero': True, 'positive': True, 'real': True, 'zero': False}

        """
        cls = type(self)
        A   = self._assumptions

        # assumptions shared:
        if A is cls.default_assumptions or (self._assume_type_keys is None):
            assumptions0 = {}
        else:
            assumptions0 = dict( (k, A[k]) for k in self._assume_type_keys )

        return assumptions0

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


           3. evalf() for comparable
           -------------------------

           as a last resort for comparable objects we get their numerical value
           -- this helps to determine facts like 'positive' and 'negative'



           In all cases when we settle on some fact value, it is given to
           _learn_new_facts to deduce all its implications, and also the result
           is cached in ._assumptions for later quick access.
        """
        if k not in _assume_defined:
            raise AttributeError('undefined assumption %r' % (k))

        seen = self._a_inprogress
        if k in seen:
            raise CycleDetected
        seen.append(k)

        try:
            # First try the assumption evaluation function if it exists
            if hasattr(self, '_eval_is_' + k):
                try:
                    a = getattr(self, '_eval_is_' + k)()
                except CycleDetected:
                    pass
                else:
                    if a is not None:
                        self._learn_new_facts( ((k, a),) )
                        return a

            # Try assumption's prerequisites
            for pk in _assume_rules.prereq.get(k, ()):
                if hasattr(self, '_eval_is_' + pk):
                    # cycle
                    if pk in seen:
                        continue
                    a = getattr(self, 'is_' + pk)
                    if a is not None:
                        self._learn_new_facts( ((pk,a),) )
                        # it is possible that we either know or don't know k at
                        # this point
                        try:
                            return self._assumptions[k]
                        except KeyError:
                            pass
        finally:
            seen.pop()

        # For positive/negative try to ask evalf
        if k in _real_ordering and self.is_comparable:
            #FIXME-py3k: this fails for complex numbers, when we define cmp
            #FIXME-py3k: as (a>b) - (a<b)
            a = _real_cmp0_table[k][cmp(self.evalf(), 0)]
            if a is not None:
                self._learn_new_facts( ((k,a),) )
                return a

        # No result -- unknown
        # cache it  (NB ._learn_new_facts(k, None) to learn other properties,
        # and because assumptions may not be detached)
        self._learn_new_facts( ((k,None),) )
        return None


    def _learn_new_facts(self, facts):
        """Learn new facts about self.

           *******************************************************************
           * internal routine designed to be used only from assumptions core *
           *******************************************************************

           Given new facts and already present knowledge (._assumptions) we ask
           inference engine to derive full set of new facts which follow from
           this combination.

           The result is stored back into ._assumptions
        """
        if not facts:
            return

        base = self._assumptions
        if base is self.default_assumptions:
            base = base.copy()
            self._assumptions = base
        # NOTE it modifies base inplace
        _assume_rules.deduce_all_facts(facts, base)

