"""Base class for all objects in sympy"""

import sympy.mpmath as mpmath

from assumptions import AssumeMeths, make__get_assumption
from sympify import _sympify, _sympifyit, sympify, SympifyError
from cache import cacheit, Memoizer, MemoizerArg

# from numbers  import Number, Integer, Rational, Real /cyclic/
# from interval import Interval /cyclic/
# from symbol   import Symbol, Wild, Temporary /cyclic/
# from add      import Add  /cyclic/
# from mul      import Mul  /cyclic/
# from power    import Pow  /cyclic/
# from function import Derivative, FunctionClass   /cyclic/
# from relational import Equality, Unequality, Inequality, StrictInequality /cyclic/
# from sympy.functions.elementary.complexes import abs as abs_   /cyclic/
# from sympy.printing import StrPrinter

# used for canonical ordering of symbolic sequences
# via __cmp__ method:
# FIXME this is *so* irrelevant and outdated!
ordering_of_classes = [
    # singleton numbers
    'Zero', 'One','Half','Infinity','NaN','NegativeOne','NegativeInfinity',
    # numbers
    'Integer','Rational','Real',
    # singleton symbols
    'Exp1','Pi','ImaginaryUnit',
    # symbols
    'Symbol','Wild','Temporary',
    # Functions that should come before Pow/Add/Mul
    'ApplyConjugate', 'ApplyAbs',
    # arithmetic operations
    'Pow', 'Mul', 'Add',
    # function values
    'Apply',
    'ApplyExp','ApplyLog',
    'ApplySin','ApplyCos','ApplyTan','ApplyCot',
    'ApplyASin','ApplyACos','ApplyATan','ApplyACot',
    'ApplySinh','ApplyCosh','ApplyTanh','ApplyCoth',
    'ApplyASinh','ApplyACosh','ApplyATanh','ApplyACoth',
    'ApplyRisingFactorial','ApplyFallingFactorial',
    'ApplyFactorial','ApplyBinomial',
    'ApplyFloor', 'ApplyCeiling',
    'ApplyRe','ApplyIm', 'ApplyArg',
    'ApplySqrt','ApplySign',
    'ApplyMrvLog',
    'ApplyGamma','ApplyLowerGamma','ApplyUpperGamma','ApplyPolyGamma',
    'ApplyErf',
    'ApplyChebyshev','ApplyChebyshev2',
    'Derivative','Integral',
    # defined singleton functions
    'Abs','Sign','Sqrt',
    'Floor', 'Ceiling',
    'Re', 'Im', 'Arg',
    'Conjugate',
    'Exp','Log','MrvLog',
    'Sin','Cos','Tan','Cot','ASin','ACos','ATan','ACot',
    'Sinh','Cosh','Tanh','Coth','ASinh','ACosh','ATanh','ACoth',
    'RisingFactorial','FallingFactorial',
    'Factorial','Binomial',
    'Gamma','LowerGamma','UpperGamma','PolyGamma',
    'Erf',
    # special polynomials
    'Chebyshev','Chebyshev2',
    # undefined functions
    'Function','WildFunction',
    # anonymous functions
    'Lambda',
    # operators
    'FDerivative','FApply',
    # composition of functions
    'FPow', 'Composition',
    # Landau O symbol
    'Order',
    # relational operations
    'Equality', 'Unequality', 'StrictInequality', 'Inequality',
    ]


class BasicType(type):
    pass

class BasicMeta(BasicType):

    classnamespace = {}
    singleton = {}

    def __init__(cls, *args, **kws):
        n = cls.__name__
        c = BasicMeta.classnamespace.get(n)
        BasicMeta.classnamespace[n] = cls
        super(BasicMeta, cls).__init__(cls)

        # --- assumptions ---

        # initialize default_assumptions dictionary
        default_assumptions = {}

        for k,v in cls.__dict__.iteritems():
            if not k.startswith('is_'):
                continue

            # this is not an assumption (e.g. is_Integer)
            if k[3:] not in AssumeMeths._assume_defined:
                continue

            k = k[3:]
            if isinstance(v,(bool,int,long,type(None))):
                if v is not None:
                    v = bool(v)
                default_assumptions[k] = v
                #print '  %r <-- %s' % (k,v)


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
        xass = AssumeMeths._assume_rules.deduce_all_facts(default_assumptions)

        # and store completed set into cls -- this way we'll avoid rededucing
        # extensions of class default assumptions each time on instance
        # creation -- we keep it prededuced already.
        cls.default_assumptions = xass

        #print '\t(%2i)  %s' % (len(default_assumptions), default_assumptions)
        #print '\t(%2i)  %s' % (len(xass), xass)



        # let's store new derived assumptions back into class.
        # this will result in faster access to this attributes.
        #
        # Timings
        # -------
        #
        # a = Integer(5)
        # %timeit a.is_zero     -> 20 us (without this optimization)
        # %timeit a.is_zero     ->  2 us (with    this optimization)
        #
        #
        # BTW: it is very important to study the lessons learned here --
        #      we could drop Basic.__getattr__ completely (!)
        #
        # %timeit x.is_Add      -> 2090 ns  (Basic.__getattr__  present)
        # %timeit x.is_Add      ->  825 ns  (Basic.__getattr__  absent)
        #
        # so we may want to override all assumptions is_<xxx> methods and
        # remove Basic.__getattr__


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
                if not cls.__dict__.has_key('is_'+k):
                    #print '%s -- overriding: %s' % (cls.__name__, k)
                    is_k = make__get_assumption(cls.__name__, k)
                    setattr(cls, 'is_'+k, property(is_k))



    def __cmp__(cls, other):
        # If the other object is not a Basic subclass, then we are not equal to
        # it.
        if not isinstance(other, BasicType):
            return -1
        n1 = cls.__name__
        n2 = other.__name__
        c = cmp(n1,n2)
        if not c: return 0

        UNKNOWN = len(ordering_of_classes)+1
        try:
            i1 = ordering_of_classes.index(n1)
        except ValueError:
            #print 'Add',n1,'to basic.ordering_of_classes list'
            #return c
            i1 = UNKNOWN
        try:
            i2 = ordering_of_classes.index(n2)
        except ValueError:
            #print 'Add',n2,'to basic.ordering_of_classes list'
            #return c
            i2 = UNKNOWN
        if i1 == UNKNOWN and i2 == UNKNOWN:
            return c
        return cmp(i1,i2)

    def __lt__(cls, other):
        if cls.__cmp__(other)==-1:
            return True
        return False

    def __gt__(cls, other):
        if cls.__cmp__(other)==1:
            return True
        return False




class Basic(AssumeMeths):
    """
    Base class for all objects in sympy.

    Conventions:

    1)
    When you want to access parameters of some instance, always use .args:
    Example:

    >>> from sympy import symbols, cot
    >>> x, y = symbols('xy')

    >>> cot(x).args
    (x,)

    >>> cot(x).args[0]
    x

    >>> (x*y).args
    (x, y)

    >>> (x*y).args[1]
    y


    2) Never use internal methods or variables (the ones prefixed with "_").
    Example:

    >>> cot(x)._args    #don't use this, use cot(x).args instead
    (x,)


    """

    __metaclass__ = BasicMeta

    __slots__ = ['_mhash',              # hash value
                 '_args',               # arguments
                 '_assume_type_keys',   # assumptions typeinfo keys
                ]

    # To be overridden with True in the appropriate subclasses
    is_Atom = False
    is_Symbol = False
    is_Function = False
    is_Add = False
    is_Mul = False
    is_Pow = False
    is_Number = False
    is_Real = False
    is_Rational = False
    is_Integer = False
    is_NumberSymbol = False
    is_Order = False
    is_Derivative   = False

    def __new__(cls, *args, **assumptions):
        obj = object.__new__(cls)

        # FIXME we are slowed a *lot* by Add/Mul passing is_commutative as the
        # only assumption.
        #
        # .is_commutative is not an assumption -- it's like typeinfo!!!
        # we should remove it.

        # initially assumptions are shared between instances and class
        obj._assumptions  = cls.default_assumptions
        obj._a_inprogress = []

        # NOTE this could be made lazy -- probably not all instances will need
        # fully derived assumptions?
        if assumptions:
            obj._learn_new_facts(assumptions)
            #                      ^
            # FIXME this is slow   |    another NOTE: speeding this up is *not*
            #        |             |    important. say for %timeit x+y most of
            # .------'             |    the time is spent elsewhere
            # |                    |
            # |  XXX _learn_new_facts  could be asked about what *new* facts have
            # v  XXX been learned -- we'll need this to append to _hashable_content
            basek = set(cls.default_assumptions.keys())
            k2    = set(obj._assumptions.keys())
            newk  = k2.difference(basek)

            obj._assume_type_keys = frozenset(newk)
        else:
            obj._assume_type_keys = None

        obj._mhash = None # will be set by __hash__ method.
        obj._args = args  # all items in args must be Basic objects
        return obj


    # XXX better name?
    @property
    def assumptions0(self):
        """
        Return object `type` assumptions.

        For example:

          Symbol('x', real=True)
          Symbol('x', integer=True)

        are different objects, and besides Python type (Symbol), initial
        assumptions, are too forming their typeinfo.

        Example:

        >>> from sympy import Symbol
        >>> x = Symbol("x")
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


    def new(self, *args):
        """
        Create new 'similar' object.

        this is conceptually equivalent to:

          type(self) (*args)

        but takes type assumptions into account. See also: assumptions0

        Example:

        >>> x = Symbol("x")
        >>> x.new("x")
        x

        """
        obj = type(self) (*args, **self.assumptions0)
        return obj


    # NOTE NOTE NOTE
    # --------------
    #
    # new-style classes + __getattr__ is *very* slow!

    # def __getattr__(self, name):
    #     raise Warning('no way, *all* attribute access will be 2.5x slower')

    # here is what we do instead:
    for k in AssumeMeths._assume_defined:
        exec "is_%s  = property(make__get_assumption('Basic', '%s'))" % (k,k)
    del k

    # NB: there is no need in protective __setattr__

    def __getnewargs__(self):
        """ Pickling support.
        """
        return tuple(self.args)

    def __hash__(self):
        # hash cannot be cached using cache_it because infinite recurrence
        # occurs as hash is needed for setting cache dictionary keys
        h = self._mhash
        if h is None:
            h = (type(self).__name__,) + self._hashable_content()

            if self._assume_type_keys is not None:
                a = []
                kv= self._assumptions
                for k in sorted(self._assume_type_keys):
                    a.append( (k, kv[k]) )

                h = hash( h + tuple(a) )

            else:
                h = hash( h )


            self._mhash = h
            return h

        else:
            return h

    def _hashable_content(self):
        # If class defines additional attributes, like name in Symbol,
        # then this method should be updated accordingly to return
        # relevant attributes as tuple.
        return self._args

    def __nonzero__(self):
        """Tests if 'self' is an instance of Zero class.

           This should be understand as an idiom:

               [1] bool(x) <=> bool(x is not S.Zero)

               [2] bool(not x) <=> bool(x is S.Zero)

           Allowing definition of __nonzero__ method is important in
           algorithms where uniform handling of int, long values and
           and sympy expressions is required.

           >>> from sympy import *
           >>> x,y = symbols('xy')

           >>> bool(0)
           False
           >>> bool(1)
           True

           >>> bool(S.Zero)
           False
           >>> bool(S.One)
           True

           >>> bool(x*y)
           True
           >>> bool(x + y)
           True

        """
        return self is not S.Zero

    def compare(self, other):
        """
        Return -1,0,1 if the object is smaller, equal, or greater than other.

        Not always in mathematical sense. If the object is of different type
        from other then their classes are ordered according to sorted_classes
        list.

        Example:

        >>> from sympy import *
        >>> x, y = symbols("x y")
        >>> x.compare(y)
        -1
        >>> x.compare(x)
        0
        >>> y.compare(x)
        1

        """
        # all redefinitions of __cmp__ method should start with the
        # following three lines:
        if self is other: return 0
        c = cmp(self.__class__, other.__class__)
        if c: return c
        #
        st = self._hashable_content()
        ot = other._hashable_content()
        c = cmp(len(st),len(ot))
        if c: return c
        for l,r in zip(st,ot):
            if isinstance(l, Basic):
                c = l.compare(r)
            else:
                c = cmp(l, r)
            if c: return c
        return 0

    @staticmethod
    def _compare_pretty(a, b):
        from sympy.series.order import Order
        if isinstance(a, Order) and not isinstance(b, Order):
            return 1
        if not isinstance(a, Order) and isinstance(b, Order):
            return -1

        # FIXME this produces wrong ordering for 1 and 0
        # e.g. the ordering will be 1 0 2 3 4 ...
        # because 1 = x^0, but 0 2 3 4 ... = x^1
        p1, p2, p3 = Wild("p1"), Wild("p2"), Wild("p3")
        r_a = a.match(p1 * p2**p3)
        r_b = b.match(p1 * p2**p3)
        if r_a is not None and r_b is not None:
            c = Basic.compare(r_a[p3], r_b[p3])
            if c!=0:
                return c

        return Basic.compare(a,b)

    @staticmethod
    def compare_pretty(a, b):
        """
        Is a>b in the sense of ordering in printing?

        yes ..... return 1
        no ...... return -1
        equal ... return 0

        Strategy:

        It uses Basic.compare as a fallback, but improves it in many cases,
        like x**3, x**4, O(x**3) etc. In those simple cases, it just parses the
        expression and returns the "sane" ordering such as:

        1 < x < x**2 < x**3 < O(x**4) etc.

        Example:

        >>> x = Symbol("x")
        >>> Basic.compare_pretty(x, x**2)
        -1
        >>> Basic.compare_pretty(x**2, x**2)
        0
        >>> Basic.compare_pretty(x**3, x**2)
        1

        """
        try:
            a = _sympify(a)
        except SympifyError:
            pass

        try:
            b = _sympify(b)
        except SympifyError:
            pass

        # both objects are non-SymPy
        if (not isinstance(a, Basic)) and (not isinstance(b, Basic)):
            return cmp(a,b)

        if not isinstance(a, Basic):
            return -1   # other < sympy

        if not isinstance(b, Basic):
            return +1   # sympy > other

        # now both objects are from SymPy, so we can proceed to usual comparison
        return Basic._compare_pretty(a, b)


    def __eq__(self, other):
        """a == b  -> Compare two symbolic trees and see whether they are equal

           this is the same as:

             a.compare(b) == 0

           but faster
        """

        if type(self) is not type(other):
            try:
                other = _sympify(other)
            except SympifyError:
                return False    # sympy != other

            if type(self) is not type(other):
                return False

        # type(self) == type(other)
        st = self._hashable_content()
        ot = other._hashable_content()

        return (st == ot)

    def __ne__(self, other):
        """a != b  -> Compare two symbolic trees and see whether they are different

           this is the same as:

             a.compare(b) != 0

           but faster
        """

        if type(self) is not type(other):
            try:
                other = _sympify(other)
            except SympifyError:
                return True     # sympy != other

            if type(self) is not type(other):
                return True

        # type(self) == type(other)
        st = self._hashable_content()
        ot = other._hashable_content()

        return (st != ot)

    # TODO all comparison methods should return True/False directly (?)
    # see #153
    #
    # OTOH Py3k says
    #
    #    Comparisons other than == and != between disparate types will raise an
    #    exception unless explicitly supported by the type
    #
    # references:
    #
    # http://www.python.org/dev/peps/pep-0207/
    # http://www.python.org/dev/peps/pep-3100/#id18
    # http://mail.python.org/pipermail/python-dev/2004-June/045111.html

    @_sympifyit('other', False) # sympy >  other
    def __lt__(self, other):
        #return sympify(other) > self
        return StrictInequality(self, other)

    @_sympifyit('other', True)  # sympy >  other
    def __gt__(self, other):
        return StrictInequality(other, self)
        #return sympify(other) < self

    @_sympifyit('other', False) # sympy >  other
    def __le__(self, other):
        return Inequality(self, other)

    @_sympifyit('other', True)  # sympy >  other
    def __ge__(self, other):
        return sympify(other) <= self


    # ***************
    # * Arithmetics *
    # ***************

    def __pos__(self):
        return self
    def __neg__(self):
        return Mul(S.NegativeOne, self)
    def __abs__(self):
        return abs_(self)

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        return Add(self, other)
    @_sympifyit('other', NotImplemented)
    def __radd__(self, other):
        return Add(other, self)

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        return Add(self, -other)
    @_sympifyit('other', NotImplemented)
    def __rsub__(self, other):
        return Add(other, -self)

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        return Mul(self, other)
    @_sympifyit('other', NotImplemented)
    def __rmul__(self, other):
        return Mul(other, self)

    @_sympifyit('other', NotImplemented)
    def __pow__(self, other):
        return Pow(self, other)
    @_sympifyit('other', NotImplemented)
    def __rpow__(self, other):
        return Pow(other, self)

    @_sympifyit('other', NotImplemented)
    def __div__(self, other):
        return Mul(self, Pow(other, S.NegativeOne))
    @_sympifyit('other', NotImplemented)
    def __rdiv__(self, other):
        return Mul(other, Pow(self, S.NegativeOne))

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    # *******************
    # * Logic operators *
    # *******************

    def __and__(self, other):
        """Overloading for & operator"""
        from sympy.logic.boolalg import And
        return And(self, other)
    def __or__(self, other):
        """Overloading for |"""
        from sympy.logic.boolalg import Or
        return Or(self, other)
    def __invert__(self):
        """Overloading for ~"""
        from sympy.logic.boolalg import Not
        return Not(self)
    def __rshift__(self, other):
        """Overloading for >>"""
        from sympy.logic.boolalg import Implies
        return Implies(self, other)
    def __lshift__(self, other):
        """Overloading for <<"""
        from sympy.logic.boolalg import Implies
        return Implies(other, self)


    def __repr__(self):
        return StrPrinter.doprint(self)

    def __str__(self):
        return StrPrinter.doprint(self)

    def atoms(self, *types):
        """Returns the atoms that form the current object.
           An atom is the smallest piece in which we can divide an
           expression.


           Examples:

           >>> from sympy import *
           >>> x,y = symbols('xy')

           >>> sorted((x+y**2 + 2*x*y).atoms())
           [2, x, y]

           You can also filter the results by a given type(s) of object:

           >>> sorted((x+y+2+y**2*sin(x)).atoms(Symbol))
           [x, y]

           >>> sorted((x+y+2+y**3*sin(x)).atoms(Number))
           [2, 3]

           >>> sorted((x+y+2+y**2*sin(x)).atoms(Symbol, Number))
           [2, x, y]

           Or by a type of on object in an impliciy way:

           >>> sorted((x+y+2+y**2*sin(x)).atoms(x))
           [x, y]

        """


        def _atoms(expr, typ):
            """Helper function for recursively denesting atoms"""
            if isinstance(expr, Basic):
                if expr.is_Atom and len(typ) == 0: # if we haven't specified types
                    return [expr]
                else:
                    try:
                        if isinstance(expr, typ): return [expr]
                    except TypeError:
                        #if type is in implicit form
                        if isinstance(expr, tuple(map(type, typ))): return [expr]
            result = []
            #search for a suitable iterator
            if isinstance(expr, Basic):
                iter = expr.iter_basic_args()
            elif isinstance(expr, (tuple, list)):
                iter = expr.__iter__()
            else:
                iter = []

            for obj in iter:
                result.extend(_atoms(obj, typ))
            return result

        return set(_atoms(self, typ=types))

    def is_hypergeometric(self, k):
        from sympy.simplify import hypersimp
        return hypersimp(self, k) is not None

    @property
    def is_number(self):
        """Returns True if 'self' is a number.

           >>> from sympy import *
           >>> x,y = symbols('xy')

           >>> x.is_number
           False
           >>> (2*x).is_number
           False
           >>> (2 + log(2)).is_number
           True

        """
        for obj in self.iter_basic_args():
            if not obj.is_number:
                return False
        else:
            return True

    @property
    def func(self):
        """
        The top-level function in an expression.

        The following should hold for all objects::

            >> x == x.func(*x.args)

        Example:

        >>> x = Symbol("x")
        >>> a = 2*x
        >>> a.func
        <class 'sympy.core.mul.Mul'>
        >>> a.args
        (2, x)
        >>> a.func(*a.args)
        2*x
        >>> a == a.func(*a.args)
        True

        """
        return self.__class__

    @property
    def args(self):
        """Returns a tuple of arguments of 'self'.

        Example:

        >>> from sympy import symbols, cot
        >>> x, y = symbols('xy')

        >>> cot(x).args
        (x,)

        >>> cot(x).args[0]
        x

        >>> (x*y).args
        (x, y)

        >>> (x*y).args[1]
        y

        Note for developers: Never use self._args, always use self.args.
        Only when you are creating your own new function, use _args
        in the __new__. Don't override .args() from Basic (so that it's
        easy to change the interface in the future if needed).
        """
        return self._args[:]

    def iter_basic_args(self):
        """
        Iterates arguments of 'self'.

        Example:

        >>> x = Symbol("x")
        >>> a = 2*x
        >>> a.iter_basic_args()
        <tupleiterator object at 0x...>
        >>> list(a.iter_basic_args())
        [2, x]

        """
        return iter(self.args)

    def is_rational_function(self, *syms):
        """
        Test whether function is rational function - ratio of two polynomials.
        When no arguments are present, Basic.atoms(Symbol) is used instead.

        Example:

        >>> from sympy import symbols, sin
        >>> x, y = symbols('xy')

        >>> (x/y).is_rational_function()
        True

        >>> (x**2).is_rational_function()
        True

        >>> (x/sin(y)).is_rational_function(y)
        False

        """
        p, q = self.as_numer_denom()

        if p.is_polynomial(*syms):
            if q.is_polynomial(*syms):
                return True

        return False

    def _eval_is_polynomial(self, syms):
        return

    def is_polynomial(self, *syms):
        if syms:
            syms = map(sympify, syms)
        else:
            syms = list(self.atoms(Symbol))

        if not syms: # constant polynomial
            return True
        else:
            return self._eval_is_polynomial(syms)

    def as_poly(self, *symbols, **flags):
        """Converts 'self' to a polynomial or returns None.

           When constructing a polynomial an exception will be raised in
           case the input expression is not convertible to a polynomial.
           There are situations when it is easier (simpler or prettier)
           to receive None on failure.

           If no symbols were given and 'self' isn't already a polynomial
           then all  available symbols will be collected and used to form
           a new polynomial.

           >>> from sympy import *
           >>> x,y = symbols('xy')

           >>> print (x**2 + x*y).as_poly()
           Poly(x**2 + x*y, x, y)

           >>> print (x**2 + x*y).as_poly(x, y)
           Poly(x**2 + x*y, x, y)

           >>> print (x**2 + sin(y)).as_poly(x, y)
           None

        """
        from sympy.polys import Poly, PolynomialError

        try:
            if not symbols:
                if isinstance(self, Poly):
                    return self
                else:
                    symbols = sorted(self.atoms(Symbol))

            return Poly(self, *symbols, **flags)
        except PolynomialError:
            return None

    def as_basic(self):
        """Converts polynomial to a valid sympy expression.

           >>> from sympy import *
           >>> x,y = symbols('xy')

           >>> p = (x**2 + x*y).as_poly(x, y)

           >>> p.as_basic()
           x*y + x**2

           >>> f = sin(x)

           >>> f.as_basic()
           sin(x)

        """
        return self

    def subs(self, *args):
        """
        Substitutes an expression.

        Calls either _subs_old_new, _subs_dict or _subs_list depending
        if you give it two arguments (old, new), a dictionary or a list.

        Examples:

        >>> from sympy import *
        >>> x,y = symbols('xy')
        >>> (1+x*y).subs(x, pi)
        1 + pi*y
        >>> (1+x*y).subs({x:pi, y:2})
        1 + 2*pi
        >>> (1+x*y).subs([(x,pi), (y,2)])
        1 + 2*pi

        """
        if len(args) == 1:
            sequence = args[0]
            if isinstance(sequence, dict):
                return self._subs_dict(sequence)
            elif isinstance(sequence, (list, tuple)):
                return self._subs_list(sequence)
            else:
                raise TypeError("Not an iterable container")
        elif len(args) == 2:
            old, new = args
            return self._subs_old_new(old, new)
        else:
            raise TypeError("subs accepts either 1 or 2 arguments")

    @cacheit
    def _subs_old_new(self, old, new):
        """Substitutes an expression old -> new."""
        old = sympify(old)
        new = sympify(new)
        return self._eval_subs(old, new)

    def _eval_subs(self, old, new):
        if self==old:
            return new
        return self

    def _subs_list(self, sequence):
        """
        Performs an order sensitive substitution from the
        input sequence list.

        Examples:

        >>> from sympy import *
        >>> x, y = symbols('xy')
        >>> (x+y)._subs_list( [(x, 3),     (y, x**2)] )
        3 + x**2
        >>> (x+y)._subs_list( [(y, x**2),  (x, 3)   ] )
        12

        """
        if not isinstance(sequence, (list, tuple)):
            raise TypeError("Not an iterable container")
        result = self
        for old, new in sequence:
            if hasattr(result, 'subs'):
                result = result.subs(old, new)
        return result

    def _subs_dict(self, sequence):
        """Performs sequential substitution.

           Given a collection of key, value pairs, which correspond to
           old and new expressions respectively,  substitute all given
           pairs handling properly all overlapping keys  (according to
           'in' relation).

           We have to use naive O(n**2) sorting algorithm, as 'in'
           gives only partial order and all asymptotically faster
           fail (depending on the initial order).

           >>> from sympy import *
           >>> x, y = symbols('xy')

           >>> a,b,c,d,e = symbols('abcde')

           >>> A = (sqrt(sin(2*x)), a)
           >>> B = (sin(2*x), b)
           >>> C = (cos(2*x), c)
           >>> D = (x, d)
           >>> E = (exp(x), e)

           >>> expr = sqrt(sin(2*x))*sin(exp(x)*x)*cos(2*x) + sin(2*x)

           >>> expr._subs_dict([A,B,C,D,E])
           b + a*c*sin(d*e)

        """
        if isinstance(sequence, dict):
            sequence = sequence.items()
        elif not isinstance(sequence, (list, tuple)):
            raise TypeError("Not an iterable container")

        subst = []

        for pattern in sequence:
            for i, (expr, _) in enumerate(subst):
                if pattern[0] in expr:
                    subst.insert(i, pattern)
                    break
            else:
                subst.append(pattern)
        subst.reverse()
        return self._subs_list(subst)

    def _seq_subs(self, old, new):
        if self==old:
            return new
        #new functions are initialized differently, than old functions
        if isinstance(self.func, FunctionClass):
            args = self.args[:]
        else:
            args = (self.func,)+self[:]
        return self.__class__(*[s.subs(old, new) for s in args])

    def __contains__(self, what):
        if self == what: return True
        for x in self._args:
            if what in x:
                return True
        return False

    @cacheit
    def has_any_symbols(self, *syms):
        """Return True if 'self' has any of the symbols.

           >>> from sympy import *
           >>> x,y,z = symbols('xyz')

           >>> (x**2 + sin(x*y)).has_any_symbols(z)
           False

           >>> (x**2 + sin(x*y)).has_any_symbols(x, y)
           True

           >>> (x**2 + sin(x*y)).has_any_symbols(x, y, z)
           True

        """
        syms = set(syms)

        if not syms:
            return True
        else:
            def search(expr):
                if expr.is_Atom:
                    if expr.is_Symbol:
                        return expr in syms
                    else:
                        return False
                else:
                    for term in expr.iter_basic_args():
                        if search(term):
                            return True
                    else:
                        return False

            return search(self)

    @cacheit
    def has_all_symbols(self, *syms):
        """Return True if 'self' has all of the symbols.

           >>> from sympy import *
           >>> x,y,z = symbols('xyz')

           >>> (x**2 + sin(x*y)).has_all_symbols(x, y)
           True

           >>> (x**2 + sin(x*y)).has_all_symbols(x, y, z)
           False

        """
        syms = set(syms)

        if not syms:
            return True
        else:
            def search(expr):
                if expr.is_Atom:
                    if expr.is_Symbol and expr in syms:
                        syms.remove(expr)
                else:
                    for term in expr.iter_basic_args():
                        if not syms:
                            break
                        else:
                            search(term)

            search(self)

            return not syms

    def has(self, *patterns):
        """
        Return True if self has any of the patterns.

        Example:
        >>> x = Symbol("x")
        >>> (2*x).has(x)
        True
        >>> (2*x/x).has(x)
        False

        """
        if len(patterns)>1:
            for p in patterns:
                if self.has(p):
                    return True
            return False
        elif not patterns:
            raise TypeError("has() requires at least 1 argument (got none)")
        p = sympify(patterns[0])
        if p.is_Symbol and not isinstance(p, Wild): # speeds up
            return p in self.atoms(p.__class__)
        if isinstance(p, BasicType):
            #XXX this is very fragile:
            if str(self).find(str(p.__name__)) == -1:
                #didn't find p in self, let's try corner cases
                if p is Derivative:
                    if str(self).find("D(") != -1:
                        return True
                return False
            else:
                return True
        if p.matches(self) is not None:
            return True
        if not False:
            args = self.args[:]
        else:
            args = (self.func,)+self.args[:]
        for e in args:
            if e.has(p):
                return True
        return False

    def _eval_interval(self, x, a, b):
        """
        Returns evaluation over an interval.  For most funtions this is:

        self.subs(x, b) - self.subs(x, a),

        possibly using limit() if NaN is returned from subs.

        """
        from sympy.series import limit
        A = self.subs(x, a)

        if A is S.NaN:
            A = limit(self, x, a)
            if A is S.NaN:
                return self

        B = self.subs(x, b)

        if B is S.NaN:
            B = limit(self, x, b)
        if B is S.NaN:
            return self
        return B - A

    def _eval_power(self, other):
        return None

    def _eval_derivative(self, s):
        return

    def _eval_fapply(self, *args, **assumptions):
        return

    def _eval_fpower(b, e):
        return

    def _eval_apply_evalf(self,*args):
        return

    def _eval_eq_nonzero(self, other):
        return

    @classmethod
    def _eval_apply_subs(cls, *args):
        return

    def _eval_conjugate(self):
        if self.is_real:
            return self

    def conjugate(self):
        from sympy.functions.elementary.complexes import conjugate as c
        return c(self)

    def removeO(self):
        "Removes the O(..) symbol if there is one"
        if self.is_Order:
            return Integer(0)
        for i,x in enumerate(self.args):
            if x.is_Order:
                return Add(*(self.args[:i]+self.args[i+1:]))
        return self

    def getO(e):
        "Returns the O(..) symbol, or None if there is none."
        if e.is_Order:
            return e
        if e.is_Add:
            for x in e.args:
                if x.is_Order:
                    return x

    #@classmethod
    def matches(pattern, expr, repl_dict={}, evaluate=False):
        """
        Helper method for match() - switches the pattern and expr.

        Can be used to solve linear equations:
          >>> from sympy import Symbol, Wild, Integer
          >>> a,b = map(Symbol, 'ab')
          >>> x = Wild('x')
          >>> (a+b*x).matches(Integer(0))
          {x_: -a/b}

        """

        # weed out negative one prefixes
        sign = 1
        if pattern.is_Mul and pattern.args[0] == -1:
          pattern = -pattern; sign = -sign
        if expr.is_Mul and expr.args[0] == -1:
          expr = -expr; sign = -sign

        if evaluate:
            pat = pattern
            for old,new in repl_dict.items():
                pat = pat.subs(old, new)
            if pat!=pattern:
                return pat.matches(expr, repl_dict)
        expr = sympify(expr)
        if not isinstance(expr, pattern.__class__):
            # if we can omit the first factor, we can match it to sign * one
            if pattern.is_Mul and Mul(*pattern.args[1:]) == expr:
               return pattern.args[0].matches(Rational(sign), repl_dict, evaluate)
            # two-factor product: if the 2nd factor matches, the first part must be sign * one
            if pattern.is_Mul and len(pattern.args[:]) == 2:
               dd = pattern.args[1].matches(expr, repl_dict, evaluate)
               if dd == None: return None
               dd = pattern.args[0].matches(Rational(sign), dd, evaluate)
               return dd
            return None

        if len(pattern.args[:])==0:
            if pattern==expr:
                return repl_dict
            return None
        d = repl_dict.copy()

        # weed out identical terms
        pp = list(pattern.args)
        ee = list(expr.args)
        for p in pattern.args:
          for e in expr.args:
            if e == p:
              if e in ee: ee.remove(e)
              if p in pp: pp.remove(p)

        # only one symbol left in pattern -> match the remaining expression
        if len(pp) == 1 and isinstance(pp[0], Wild):
          if len(ee) == 1: d[pp[0]] = sign * ee[0]
          else: d[pp[0]] = sign * (type(expr)(*ee))
          return d

        if len(ee) != len(pp):
            return None

        i = 0
        for p,e in zip(pp, ee):
            if i == 0 and sign != 1:
              try: e = sign * e
              except TypeError: return None
            d = p.matches(e, d, evaluate=not i)
            i += 1
            if d is None:
                return None
        return d

    def match(self, pattern):
        """
        Pattern matching.

        Wild symbols match all.

        Return None when expression (self) does not match
        with pattern. Otherwise return a dictionary such that

          pattern.subs(self.match(pattern)) == self

        Example:

        >>> from sympy import symbols
        >>> x, y = symbols("x y")
        >>> p = Wild("p")
        >>> q = Wild("q")
        >>> r = Wild("r")
        >>> e = (x+y)**(x+y)
        >>> e.match(p**p)
        {p_: x + y}
        >>> e.match(p**q)
        {p_: x + y, q_: x + y}
        >>> e = (2*x)**2
        >>> e.match(p*q**r)
        {p_: 4, q_: x, r_: 2}
        >>> (p*q**r).subs(e.match(p*q**r))
        4*x**2

        """
        pattern = sympify(pattern)
        return pattern.matches(self, {})

    def solve4linearsymbol(eqn, rhs, symbols = None):
        """
        Solve equation "eqn == rhs" with respect to some linear symbol in eqn.

        Returns (symbol, solution). If eqn is nonlinear with respect to all
        symbols, then return trivial solution (eqn, rhs).
        """
        if eqn.is_Symbol:
            return (eqn, rhs)
        if symbols is None:
            symbols = eqn.atoms(Symbol)
        if symbols:
            # find  symbol
            for s in symbols:
                deqn = eqn.diff(s)
                if deqn.diff(s) is S.Zero:
                    # eqn = a + b*c, a=eqn(c=0),b=deqn(c=0)
                    return s, (rhs - eqn.subs(s,0))/deqn.subs(s,0)
        # no linear symbol, return trivial solution
        return eqn, rhs

    @cacheit
    def count_ops(self, symbolic=True):
        """ Return the number of operations in expressions.

        Examples:
        >>> (1+a+b**2).count_ops()
        POW + 2 * ADD
        >>> (sin(x)*x+sin(x)**2).count_ops()
        ADD + MUL + POW + 2 * SIN
        """
        return Integer(len(self[:])-1) + sum([t.count_ops(symbolic=symbolic) for t in self])

    def doit(self, **hints):
        """Evaluate objects that are not evaluated by default like limits,
           integrals, sums and products. All objects of this kind will be
           evaluated unless some species were excluded via 'hints'.

           >>> from sympy import *
           >>> x, y = symbols('xy')

           >>> 2*Integral(x, x)
           2*Integral(x, x)

           >>> (2*Integral(x, x)).doit()
           x**2

        """
        terms = [ term.doit(**hints) for term in self.args ]
        return self.new(*terms)

    ###########################################################################
    ###################### EXPRESSION EXPANSION METHODS #######################
    ###########################################################################

    # These should be overridden in subclasses

    def _eval_expand_basic(self, deep=True, **hints):
        return self

    def _eval_expand_power_exp(self, deep=True, **hints):
        return self

    def _eval_expand_power_base(self, deep=True, **hints):
        return self

    def _eval_expand_mul(self, deep=True, **hints):
        return self

    def _eval_expand_multinomial(self, deep=True, **hints):
        return self

    def _eval_expand_log(self, deep=True, **hints):
        return self

    def _eval_expand_complex(self, deep=True, **hints):
        return self

    def _eval_expand_trig(self, deep=True, **hints):
        return self

    def _eval_expand_func(self, deep=True, **hints):
        return self

    def expand(self, deep=True, power_base=True, power_exp=True, mul=True, \
           log=True, multinomial=True, basic=True, **hints):
        """
        Expand an expression using hints.

        See the docstring in function.expand for more information.
        """
        hints['power_base'] = power_base
        hints['power_exp'] = power_exp
        hints['mul'] = mul
        hints['log'] = log
        hints['multinomial'] = multinomial
        hints['basic'] = basic

        expr = self

        for hint in hints:
            if hints[hint] == True:
                func = getattr(expr, '_eval_expand_'+hint, None)
                if func is not None:
                    expr = func(deep=deep, **hints)

        if not 'basic' in hints:
            if not expr.is_Atom:
                result = expr._eval_expand_basic()

                if result is not None:
                    expr = result

        return expr

    def _eval_rewrite(self, pattern, rule, **hints):
        if self.is_Atom:
            return self
        sargs = self.args
        terms = [ t._eval_rewrite(pattern, rule, **hints) for t in sargs ]
        return self.new(*terms)

    def rewrite(self, *args, **hints):
        """Rewrites expression containing applications of functions
           of one kind in terms of functions of different kind. For
           example you can rewrite trigonometric functions as complex
           exponentials or combinatorial functions as gamma function.

           As a pattern this function accepts a list of functions to
           to rewrite (instances of DefinedFunction class). As rule
           you can use string or a destination function instance (in
           this case rewrite() will use the str() function).

           There is also possibility to pass hints on how to rewrite
           the given expressions. For now there is only one such hint
           defined called 'deep'. When 'deep' is set to False it will
           forbid functions to rewrite their contents.

           >>> from sympy import *
           >>> x, y = symbols('xy')

           >>> sin(x).rewrite(sin, exp)
           -I*(exp(I*x) - exp(-I*x))/2

        """
        if self.is_Atom or not args:
            return self
        else:
            pattern, rule = args[:-1], args[-1]

            if not isinstance(rule, str):

                if rule == C.tan:
                    rule = "tan"
                elif rule == C.exp:
                    rule = "exp"
                elif isinstance(rule, FunctionClass):   # new-style functions
                    #print rule
                    rule = rule.__name__  # XXX proper attribute for name?
                    #print rule
                else:
                    rule = str(rule)

            rule = '_eval_rewrite_as_' + rule

            if not pattern:
                return self._eval_rewrite(None, rule, **hints)
            else:
                if isinstance(pattern[0], (tuple, list)):
                    pattern = pattern[0]

                pattern = [ p.__class__ for p in pattern if self.has(p) ]

                if pattern:
                    return self._eval_rewrite(tuple(pattern), rule, **hints)
                else:
                    return self

    def coeff(self, x):
        """
        Returns the coefficient of the term "x" or None if there is no "x".

        Example:

        >>> x = Symbol('x')

        >>> (3+2*x+4*x**2).coeff(1)
        >>> (3+2*x+4*x**2).coeff(x)
        2
        >>> (3+2*x+4*x**2).coeff(x**2)
        4
        >>> (3+2*x+4*x**2).coeff(x**3)
        >>>
        """
        from sympy import collect
        x = sympify(x)
        const = x.as_coeff_terms()[0] # constant multiplying x
        if const != S.One: # get rid of constants
            result = self.coeff(x/const)
            if result is not None:
                return(result/const)
            else:
                return None
        if x.is_Integer:
            return
        self = self.expand() # collect expects it's arguments in expanded form
        result = collect(self, x, evaluate=False)
        if x in result:
            return result[x]
        else:
            return None

    def as_coefficient(self, expr):
        """Extracts symbolic coefficient at the given expression. In
           other words, this functions separates 'self' into product
           of 'expr' and 'expr'-free coefficient. If such separation
           is not possible it will return None.

           >>> from sympy import *
           >>> x, y = symbols('xy')

           >>> E.as_coefficient(E)
           1
           >>> (2*E).as_coefficient(E)
           2

           >>> (2*E + x).as_coefficient(E)
           >>> (2*sin(E)*E).as_coefficient(E)

           >>> (2*pi*I).as_coefficient(pi*I)
           2

           >>> (2*I).as_coefficient(pi*I)

        """
        if expr.is_Add:
            return None
        else:
            w = Wild('w')

            coeff = self.match(w * expr)

            if coeff is not None:
                if expr.is_Mul:
                    expr = expr.args
                else:
                    expr = [expr]

                if coeff[w].has(*expr):
                    return None
                else:
                    return coeff[w]
            else:
                return None

    def as_independent(self, *deps):
        """Returns a pair with separated parts of a given expression
           independent of specified symbols in the first place and
           dependent on them in the other. Both parts are valid
           SymPy expressions.

           >>> from sympy import *
           >>> x, y = symbols('xy')

           >>> (2*x*sin(x)+y+x).as_independent(x)
           (y, x + 2*x*sin(x))

           >>> (x*sin(x)*cos(y)).as_independent(x)
           (cos(y), x*sin(x))

           All other expressions are multiplicative:

           >>> (sin(x)).as_independent(x)
           (1, sin(x))

           >>> (sin(x)).as_independent(y)
           (sin(x), 1)

        """
        indeps, depend = [], []

        if self.is_Add or self.is_Mul:
            for term in self.args[:]:
                if term.has(*deps):
                    depend.append(term)
                else:
                    indeps.append(term)

            return (self.__class__(*indeps),
                    self.__class__(*depend))
        else:
            if self.has(*deps):
                return (S.One, self)
            else:
                return (self, S.One)

    def as_real_imag(self):
        """Performs complex expansion on 'self' and returns a tuple
           containing collected both real and imaginary parts. This
           method can't be confused with re() and im() functions,
           which does not perform complex expansion at evaluation.

           However it is possible to expand both re() and im()
           functions and get exactly the same results as with
           a single call to this function.

           >>> from sympy import *

           >>> x, y = symbols('xy', real=True)

           >>> (x + y*I).as_real_imag()
           (x, y)

           >>> z, w = symbols('zw')

           >>> (z + w*I).as_real_imag()
           (-im(w) + re(z), im(z) + re(w))

        """
        expr = self.expand(complex=True, deep=False)

        if not expr.is_Add:
            expr = [expr]

        re_part, im_part = [], []
        if isinstance(expr, Basic):
            expr = expr.args

        for term in expr:
            coeff = term.as_coefficient(S.ImaginaryUnit)

            if coeff is None:
                re_part.append(term)
            else:
                im_part.append(coeff)

        return (Add(*re_part), Add(*im_part))

    def as_powers_dict(self):
        return { self : S.One }

    def as_base_exp(self):
        # a -> b ** e
        return self, S.One

    def as_coeff_terms(self, x=None):
        # a -> c * t
        if x is not None:
            if not self.has(x):
                return self, tuple()
        return S.One, (self,)

    def as_coeff_factors(self, x=None):
        # a -> c + f
        if x is not None:
            if not self.has(x):
                return self, tuple()
        return S.Zero, (self,)

    def as_numer_denom(self):
        # a/b -> a,b
        base, exp = self.as_base_exp()
        coeff, terms = exp.as_coeff_terms()
        if coeff.is_negative:
            # b**-e -> 1, b**e
            return S.One, base ** (-exp)
        return self, S.One

    def normal(self):
        n, d = self.as_numer_denom()
        if d is S.One:
            return n
        return n/d

    def extract_multiplicatively(self, c):
        """Return None if it's not possible to make self in the form
           c * something in a nice way, i.e. preserving the properties
           of arguments of self.

           >>> from sympy import *

           >>> x, y = symbols('xy', real=True)

           >>> ((x*y)**3).extract_multiplicatively(x**2 * y)
           x*y**2

           >>> ((x*y)**3).extract_multiplicatively(x**4 * y)

           >>> (2*x).extract_multiplicatively(2)
           x

           >>> (2*x).extract_multiplicatively(3)

           >>> (Rational(1,2)*x).extract_multiplicatively(3)
           x/6

        """
        c = sympify(c)
        if c is S.One:
            return self
        elif c == self:
            return S.One
        elif c.is_Mul:
            x = self.extract_multiplicatively(c.as_two_terms()[0])
            if x != None:
                return x.extract_multiplicatively(c.as_two_terms()[1])
        quotient = self / c
        if self.is_Number:
            if self is S.Infinity:
                if c.is_positive:
                    return S.Infinity
            elif self is S.NegativeInfinity:
                if c.is_negative:
                    return S.Infinity
                elif c.is_positive:
                    return S.NegativeInfinity
            elif self is S.ComplexInfinity:
                if not c.is_zero:
                    return S.ComplexInfinity
            elif self is S.NaN:
                return S.NaN
            elif self.is_Integer:
                if not quotient.is_Integer:
                    return None
                elif self.is_positive and quotient.is_negative:
                    return None
                else:
                    return quotient
            elif self.is_Rational:
                if not quotient.is_Rational:
                    return None
                elif self.is_positive and quotient.is_negative:
                    return None
                else:
                    return quotient
            elif self.is_Real:
                if not quotient.is_Real:
                    return None
                elif self.is_positive and quotient.is_negative:
                    return None
                else:
                    return quotient
        elif self.is_NumberSymbol or self.is_Symbol or self is S.ImaginaryUnit:
            if quotient.is_Mul and len(quotient.args) == 2:
                if quotient.args[0].is_Integer and quotient.args[0].is_positive and quotient.args[1] == self:
                    return quotient
            elif quotient.is_Integer:
                return quotient
        elif self.is_Add:
            newargs = []
            for arg in self.args:
                newarg = arg.extract_multiplicatively(c)
                if newarg != None:
                    newargs.append(newarg)
                else:
                    return None
            return C.Add(*newargs)
        elif self.is_Mul:
            for i in xrange(len(self.args)):
                newargs = list(self.args)
                del(newargs[i])
                tmp = C.Mul(*newargs).extract_multiplicatively(c)
                if tmp != None:
                    return tmp * self.args[i]
        elif self.is_Pow:
            if c.is_Pow and c.base == self.base:
                new_exp = self.exp.extract_additively(c.exp)
                if new_exp != None:
                    return self.base ** (new_exp)
            elif c == self.base:
                new_exp = self.exp.extract_additively(1)
                if new_exp != None:
                    return self.base ** (new_exp)

    def extract_additively(self, c):
        """Return None if it's not possible to make self in the form
           something + c in a nice way, i.e. preserving the properties
           of arguments of self.

           >>> from sympy import *

           >>> x, y = symbols('xy', real=True)

           >>> ((x*y)**3).extract_additively(1)

           >>> (x+1).extract_additively(x)
           1

           >>> (x+1).extract_additively(2*x)

           >>> (x+1).extract_additively(-x)
           1 + 2*x

           >>> (-x+1).extract_additively(2*x)
           1 - 3*x

        """
        c = sympify(c)
        if c is S.Zero:
            return self
        elif c == self:
            return S.Zero
        elif self is S.Zero:
            return None
        elif c.is_Add:
            x = self.extract_additively(c.as_two_terms()[0])
            if x != None:
                return x.extract_additively(c.as_two_terms()[1])
        sub = self - c
        if self.is_Number:
            if self.is_Integer:
                if not sub.is_Integer:
                    return None
                elif self.is_positive and sub.is_negative:
                    return None
                else:
                    return sub
            elif self.is_Rational:
                if not sub.is_Rational:
                    return None
                elif self.is_positive and sub.is_negative:
                    return None
                else:
                    return sub
            elif self.is_Real:
                if not sub.is_Real:
                    return None
                elif self.is_positive and sub.is_negative:
                    return None
                else:
                    return sub
        elif self.is_NumberSymbol or self.is_Symbol or self is S.ImaginaryUnit:
            if sub.is_Mul and len(sub.args) == 2:
                if sub.args[0].is_Integer and sub.args[0].is_positive and sub.args[1] == self:
                    return sub
            elif sub.is_Integer:
                return sub
        elif self.is_Add:
            terms = self.as_two_terms()
            subs0 = terms[0].extract_additively(c)
            if subs0 != None:
                return subs0 + terms[1]
            else:
                subs1 = terms[1].extract_additively(c)
                if subs1 != None:
                    return subs1 + terms[0]
        elif self.is_Mul:
            self_coeff, self_terms = self.as_coeff_terms()
            if c.is_Mul:
                c_coeff, c_terms = c.as_coeff_terms()
                if c_terms == self_terms:
                    new_coeff = self_coeff.extract_additively(c_coeff)
                    if new_coeff != None:
                        return new_coeff * C.Mul(*self_terms)
            elif c == self_terms:
                new_coeff = self_coeff.extract_additively(1)
                if new_coeff != None:
                    return new_coeff * C.Mul(*self_terms)

    def could_extract_minus_sign(self):
        """Canonical way to choose an element in the set {e, -e} where
           e is any expression. If the canonical element is e, we have
           e.could_extract_minus_sign() == True, else
           e.could_extract_minus_sign() == False.

           For any expression, the set {e.could_extract_minus_sign(),
           (-e).could_extract_minus_sign()} must be {True, False}.

           >>> from sympy import *

           >>> x, y = symbols("xy")

           >>> (x-y).could_extract_minus_sign() != (y-x).could_extract_minus_sign()
           True

        """
        negative_self = -self
        self_has_minus = (self.extract_multiplicatively(-1) != None)
        negative_self_has_minus = ((negative_self).extract_multiplicatively(-1) != None)
        if self_has_minus != negative_self_has_minus:
            return self_has_minus
        else:
            if self.is_Add:
                # We choose the one with less arguments with minus signs
                arg_signs = [arg.could_extract_minus_sign() for arg in self.args]
                positive_args = arg_signs.count(False)
                negative_args = arg_signs.count(True)
                if positive_args > negative_args:
                    return False
                elif positive_args < negative_args:
                    return True
            elif self.is_Mul:
                num, den = self.as_numer_denom()
                if den != 0:
                    return num.could_extract_minus_sign()

            # As a last resort, we choose the one with greater hash
            return hash(self) < hash(negative_self)

    ###################################################################################
    ##################### DERIVATIVE, INTEGRAL, FUNCTIONAL METHODS ####################
    ###################################################################################

    def diff(self, *symbols, **assumptions):
        new_symbols = map(sympify, symbols)
        if not "evaluate" in assumptions:
            assumptions["evaluate"] = True
        ret = Derivative(self, *new_symbols, **assumptions)
        return ret

    def fdiff(self, *indices):
        # FIXME FApply -> ?
        return C.FApply(C.FDerivative(*indices), self)

    def integrate(self, *args, **kwargs):
        from sympy.integrals import integrate
        return integrate(self, *args, **kwargs)

    def __call__(self, subsdict):
        """Use call as a shortcut for subs, but only support the dictionary version"""
        if not isinstance(subsdict, dict):
            raise TypeError("argument must be a dictionary")
        return self.subs(subsdict)

    def __float__(self):
        result = self.evalf()
        if result.is_Number:
            return float(result)
        else:
            raise ValueError("Symbolic value, can't compute")

    def __complex__(self):
        result = self.evalf()
        re, im = result.as_real_imag()
        return complex(float(re), float(im))

    def _evalf(self, prec):
        """Helper for evalf. Does the same thing but takes binary precision"""
        r = self._eval_evalf(prec)
        if r is None:
            r = self
        return r

    def _eval_evalf(self, prec):
        return

    def _seq_eval_evalf(self, prec):
        return self.__class__(*[s._evalf(prec) for s in self.args])

    def _to_mpmath(self, prec, allow_ints=True):
        # mpmath functions accept ints as input
        errmsg = "cannot convert to mpmath number"
        if allow_ints and self.is_Integer:
            return self.p
        v = self._eval_evalf(prec)
        if v is None:
            raise ValueError(errmsg)
        if v.is_Real:
            return mpmath.make_mpf(v._mpf_)
        # Number + Number*I is also fine
        re, im = v.as_real_imag()
        if allow_ints and re.is_Integer:
            re = mpmath.libmpf.from_int(re.p)
        elif re.is_Real:
            re = re._mpf_
        else:
            raise ValueError(errmsg)
        if allow_ints and im.is_Integer:
            im = mpmath.libmpf.from_int(im.p)
        elif im.is_Real:
            im = im._mpf_
        else:
            raise ValueError(errmsg)
        return mpmath.make_mpc((re, im))

    @staticmethod
    def _from_mpmath(x, prec):
        if hasattr(x, "_mpf_"):
            return C.Real._new(x._mpf_, prec)
        elif hasattr(x, "_mpc_"):
            re, im = x._mpc_
            re = C.Real._new(re, prec)
            im = C.Real._new(im, prec)*S.ImaginaryUnit
            return re+im
        else:
            raise TypeError("expected mpmath number (mpf or mpc)")

    ###################################################################################
    ##################### SERIES, LEADING TERM, LIMIT, ORDER METHODS ##################
    ###################################################################################

    def series(self, x, point=0, n=6, dir="+"):
        """
        Series expansion of "self" around "point".

        Usage:
            Returns the Taylor (Laurent or generalized) series of "self" around
            the point "point" (default 0) with respect to "x" until the n-th
            term (default n is 6).

            For dir="+" (default) it calculates the series from the right
            and for dir="-" the series from the left.
            For smooth functions this argument doesn't matter.

        Notes:
            This method is the most high level method and it returns the
            series including the O(x**n) term.

            Internally, it executes a method nseries(), see nseries() docstring
            for more information.
        """
        x = sympify(x)
        point = sympify(point)
        if dir == "+":
            return self.nseries(x, point, n)
        elif dir == "-":
            return self.subs(x, -x).nseries(x, -point, n).subs(x, -x)
        else:
            raise ValueError("Dir has to be '+' or '-'")


    def lseries(self, x, x0):
        """
        lseries is a generator yielding terms in the series.

        Example: if you do:

        for term in sin(x).lseries(x, 0):
            print term

        It will print all terms of the sin(x) series (i.e. it never
        terminates).

        The advantage of lseries() over nseries() is that many times you are
        just interested in the next term in the series (i.e. the first term for
        example), but you don't know how many you should ask for in nseries()
        using the "n" parameter.

        See also nseries().
        """
        return self._eval_lseries(x, x0)

    def _eval_lseries(self, x, x0):
        # default implementation of lseries is using nseries(), and adaptively
        # increasing the "n". As you can see, it is not very efficient, because
        # we are calculating the series over and over again. Subclasses should
        # override this method and implement much more efficient yielding of
        # terms.
        n = 0
        e = self.nseries(x, x0, n)
        while e.is_Order:
            n += 1
            e = self.nseries(x, x0, n)
        series = e.removeO()
        yield series
        while 1:
            n += 1
            e = self.nseries(x, x0, n).removeO()
            while series == e:
                n += 1
                e = self.nseries(x, x0, n).removeO()
            term = e - series
            series = e
            yield term

    def nseries(self, x, x0, n):
        """
        Calculates a generalized series expansion.

        nseries calculates "n" terms in the innermost expressions and then
        builds up the final series just by "cross-mutliplying" everything out.

        Advantage -- it's fast, because we don't have to determine how many
        terms we need to calculate in advance.

        Disadvantage -- you may endup with less terms than you may have
        expected, but the O(x**n) term appended will always be correct, so the
        result is correct, but maybe shorter.

        See also lseries().
        """
        return self._eval_nseries(x, x0, n)

    def _eval_nseries(self, x, x0, n):
        """
        This is a method that should be overriden in subclasses. Users should
        never call this method directly (use .nseries() instead), so you don't
        have to write docstrings for _eval_nseries().
        """
        raise NotImplementedError("(%s).nseries(%s, %s, %s)" % (self, x, x0, n))

    def limit(self, x, xlim, direction='+'):
        """ Compute limit x->xlim.
        """
        from sympy.series.limits import limit
        return limit(self, x, xlim, direction)

    @cacheit
    def as_leading_term(self, *symbols):
        """
        Returns the leading term.

        Example:

        >>> x = Symbol("x")
        >>> (1+x+x**2).as_leading_term(x)
        1
        >>> (1/x**2+x+x**2).as_leading_term(x)
        1/x**2

        Note:

        self is assumed to be the result returned by Basic.series().
        """
        from sympy import powsimp
        if len(symbols)>1:
            c = self
            for x in symbols:
                c = c.as_leading_term(x)
            return c
        elif not symbols:
            return self
        x = sympify(symbols[0])
        assert x.is_Symbol, `x`
        if not self.has(x):
            return self
        obj = self._eval_as_leading_term(x)
        if obj is not None:
            return powsimp(obj, deep=True, combine='exp')
        raise NotImplementedError('as_leading_term(%s, %s)' % (self, x))

    def _eval_as_leading_term(self, x):
        return self

    def as_coeff_exponent(self, x):
        """ c*x**e -> c,e where x can be any symbolic expression.
        """
        x = sympify(x)
        wc = Wild('wc')
        we = Wild('we')
        p  = wc*x**we
        from sympy import collect
        self = collect(self, x)
        d = self.match(p)
        if d is not None and we in d:
            return d[wc], d[we]
        return self, S.Zero

    def leadterm(self, x):
        """
        Returns the leading term a*x**b as a tuple (a, b).

        Example:

        >>> x = Symbol("x")
        >>> (1+x+x**2).leadterm(x)
        (1, 0)
        >>> (1/x**2+x+x**2).leadterm(x)
        (1, -2)

        Note:

        self is assumed to be the result returned by Basic.series().
        """
        from sympy import powsimp
        x = sympify(x)
        c,e = self.as_leading_term(x).as_coeff_exponent(x)
        c = powsimp(c, deep=True, combine='exp')
        if not c.has(x):
            return c,e
        raise ValueError("cannot compute leadterm(%s, %s), got c=%s" % (self, x, c))


    ##########################################################################
    ##################### END OF BASIC CLASS #################################
    ##########################################################################

class Atom(Basic):
    """
    A parent class for atomic things. An atom is an expression with no subexpressions.

    Examples: Symbol, Number, Rational, Integer, ...
    But not: Add, Mul, Pow, ...
    """

    is_Atom = True

    __slots__ = []

    def _eval_derivative(self, s):
        if self==s: return S.One
        return S.Zero

    def pattern_match(pattern, expr, repl_dict):
        if pattern==expr:
            return repl_dict
        return None

    def as_numer_denom(self):
        return self, S.One

    def count_ops(self, symbolic=True):
        return S.Zero

    def doit(self, **hints):
        return self

    def _eval_is_polynomial(self, syms):
        return True

    @property
    def is_number(self):
        return True

    def _eval_nseries(self, x, x0, n):
        return self


class SingletonMeta(BasicMeta):
    """Metaclass for all singletons

       All singleton classes should put this into their __metaclass__, and
       _not_ to define __new__

       example:

       class Zero(Integer):
           __metaclass__ = SingletonMeta

           p = 0
           q = 1
    """

    def __init__(cls, *args, **kw):
        BasicMeta.__init__(cls, *args, **kw)

        # we are going to inject singletonic __new__, here it is:
        def cls_new(cls):
            try:
                obj = getattr(SingletonFactory, cls.__name__)

            except AttributeError:
                obj = Basic.__new__(cls, *(), **{})
                setattr(SingletonFactory, cls.__name__, obj)

            return obj

        cls_new.__name__ = '%s.__new__' % (cls.__name__)

        assert not cls.__dict__.has_key('__new__'), \
                'Singleton classes are not allowed to redefine __new__'

        # inject singletonic __new__
        cls.__new__      = staticmethod(cls_new)

        # Inject pickling support.
        def cls_getnewargs(self):
            return ()
        cls_getnewargs.__name__ = '%s.__getnewargs__' % cls.__name__

        assert not cls.__dict__.has_key('__getnewargs__'), \
                'Singleton classes are not allowed to redefine __getnewargs__'
        cls.__getnewargs__ = cls_getnewargs


        # tag the class appropriately (so we could verify it later when doing
        # S.<something>
        cls.is_Singleton = True

class SingletonFactory:
    """
    A map between singleton classes and the corresponding instances.
    E.g. S.Exp == C.Exp()
    """

    def __getattr__(self, clsname):
        if clsname == "__repr__":
            return lambda: "S"

        cls = getattr(C, clsname)
        assert cls.is_Singleton
        obj = cls()

        # store found object in own __dict__, so the next lookups will be
        # serviced without entering __getattr__, and so will be fast
        setattr(self, clsname, obj)
        return obj

S = SingletonFactory()

# S(...) = sympify(...)
S.__call__ = sympify

class ClassesRegistry:
    """Namespace for SymPy classes

       This is needed to avoid problems with cyclic imports.
       To get a SymPy class you do this:

         C.<class_name>

       e.g.

         C.Rational
         C.Add
    """

    def __getattr__(self, name):
        try:
            cls = BasicMeta.classnamespace[name]
        except KeyError:
            raise AttributeError("No SymPy class '%s'" % name)

        setattr(self, name, cls)
        return cls

C = ClassesRegistry()

# XXX this is ugly, but needed for Memoizer('str', ...) to work
import cache
cache.C = C
del cache

# /cyclic/
import sympify as _
_.Basic     = Basic
_.BasicType = BasicType
_.S         = S
del _
