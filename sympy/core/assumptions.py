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


def make__get_assumption(classname, name):
    """Cooks function which will get named assumption

       e.g.

       class C:

           is_xxx = make__get_assumption('C', 'xxx')
           is_yyy = property( make__get_assumption('C', 'yyy'))


       then

       c = C()

       c.is_xxx()   # note braces -- it's a function call
       c.is_yyy     # no braces   -- it's a property
    """
    def getit(self):
        try:
            return self._assumptions[name]
        except KeyError:
            return self._what_known_about(name)

    getit.func_name = '%s__is_%s' % (classname, name)
    return getit


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

        # initialize default_assumptions dictionary
        default_assumptions = {}

        for k,v in cls.__dict__.iteritems():
            if not k.startswith('is_'):
                continue

            # this is not an assumption (e.g. is_Integer)
            if k[3:] not in _assume_defined:
                continue

            k = k[3:]
            if isinstance(v,(bool,int,long,type(None))):
                if v is not None:
                    v = bool(v)
                default_assumptions[k] = v

        # XXX maybe we should try to keep ._default_premises out of class ?
        # XXX __slots__ in class ?
        cls._default_premises = default_assumptions

        for base in cls.__bases__:
            try:
                base_premises = base._default_premises
            except AttributeError:
                continue    # no ._default_premises is ok

            for k,v in base_premises.iteritems():

                # if an assumption is already present in child, we should ignore base
                # e.g. Integer.is_integer=T, but Rational.is_integer=F (for speed)
                if k in default_assumptions:
                    continue

                default_assumptions[k] = v

        # deduce all consequences from default assumptions -- make it complete
        xass = _assume_rules.deduce_all_facts(default_assumptions)

        # and store completed set into cls -- this way we'll avoid rededucing
        # extensions of class default assumptions each time on instance
        # creation -- we keep it prededuced already.
        cls.default_assumptions = xass

        # first we need to collect derived premises
        derived_premises = {}
        for k,v in xass.iteritems():
            if k not in default_assumptions:
                derived_premises[k] = v
        cls._derived_premises = derived_premises

        for k,v in xass.iteritems():
            assert v == cls.__dict__.get('is_'+k, v),  (cls,k,v)
            # NOTE: this way Integer.is_even = False (inherited from Rational)
            # NOTE: the next code blocks add 'protection-properties' to overcome this
            setattr(cls, 'is_'+k, v)

        # protection e.g. for Initeger.is_even=F <- (Rational.is_integer=F)
        for base in cls.__bases__:
            try:
                base_derived_premises = base._derived_premises
            except AttributeError:
                continue    # no ._derived_premises is ok

            for k,v in base_derived_premises.iteritems():
                if ('is_'+k) not in cls.__dict__:
                    is_k = make__get_assumption(cls.__name__, k)
                    setattr(cls, 'is_'+k, property(is_k))


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
                ]
    __slots__ = []


    def __getstate__(self, cls=None):
        if cls is None:
            # This is the case for the instance that gets pickled
            cls = self.__class__

        d = {}
        # Get all data that should be stored from super classes
        for c in cls.__bases__:
            if hasattr(c, "__getstate__"):
                d.update(c.__getstate__(self, c))

        # Get all information that should be stored from cls and return the dic
        for name in cls.__slots__:
            if hasattr(self, name):
                d[name] = getattr(self, name)
        return d

    def __setstate__(self, d):
        # All values that were pickled are now assigned to a fresh instance
        for name, value in d.iteritems():
            try:
                setattr(self, name, value)
            except:
                pass


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

