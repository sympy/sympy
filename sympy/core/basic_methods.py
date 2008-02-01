""" Implementation of Basic low-level methods.
"""

import decimal
from assumptions import AssumeMeths

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

#

def repr_level(flag=None, _cache=[1]):
    if flag is None:
        return _cache[0]
    old_flag = _cache[0]
    _cache[0] = max(0, min(2, int(flag))) # restrict to 0,1,2
    return old_flag

class BasicType(type):
    pass

class MetaBasicMeths(BasicType):

    classnamespace = {}
    repr_level = 0        # defines the output of repr()
    singleton = {}

    def __init__(cls,*args,**kws):
        n = cls.__name__
        c = MetaBasicMeths.classnamespace.get(n)
        if c is None:
            MetaBasicMeths.classnamespace[n] = cls
        else:
            print 'Ignoring redefinition of %s: %s defined earlier than %s' % (n, c, cls)
        type.__init__(cls, *args, **kws)

        # initialize default_assumptions dictionary
        default_assumptions = {}
        for k in dir(cls):
            if not k.startswith('is_'):
                continue
            v = getattr(cls, k)
            k = k[3:]
            if isinstance(v,(bool,int,long)):
                default_assumptions[k] = bool(v)
        cls.default_assumptions = default_assumptions

    def __cmp__(cls, other):
        try:
            other = cls.sympify(other)
        except ValueError:
            #if we cannot sympify it, other is definitely not equal to cls
            return -1
        n1 = cls.__name__
        n2 = other.__name__
        c = cmp(n1,n2)
        if not c: return 0

        UNKNOWN = len(ordering_of_classes)+1
        try:
            i1 = ordering_of_classes.index(n1)
        except ValueError:
            #print 'Add',n1,'to basic_methods.ordering_of_classes list'
            #return c
            i1 = UNKNOWN
        try:
            i2 = ordering_of_classes.index(n2)
        except ValueError:
            #print 'Add',n2,'to basic_methods.ordering_of_classes list'
            #return c
            i2 = UNKNOWN
        if i1 == UNKNOWN and i2 == UNKNOWN:
            return c
        return cmp(i1,i2)


class BasicMeths(AssumeMeths):

    __metaclass__ = MetaBasicMeths

    Lambda_precedence = 1
    Add_precedence = 40
    Mul_precedence = 50
    Pow_precedence = 60
    Apply_precedence = 70
    Item_precedence = 75
    Atom_precedence = 1000

    def __getattr__(self, name):
        # if it's not an assumption -- we don't have it
        if not name.startswith('is_'):
            # it is important to return shortly for speed reasons:
            # we have *lots* of non-'is_' attribute access, e.g.
            # '_eval_<smth>', and a lot of them does *not* exits.
            #
            # if we are here -- it surely does not exist,
            # so let's get out of here as fast as possible.
            raise AttributeError(name)

        return self._get_assumption(name)

    def __setattr__(self, name, val):
        if name.startswith('is_'):
            raise AttributeError("Modification of assumptions is not allowed")
        else:
            AssumeMeths.__setattr__(self, name, val)

    def __hash__(self):
        # hash cannot be cached using cache_it because infinite recurrence
        # occurs as hash is needed for setting cache dictionary keys
        h = self._mhash
        if h is None:
            a = self._assume_hashable_content()
            self._mhash = h = hash((self.__class__.__name__,) + self._hashable_content() + a)
        return h

    def _hashable_content(self):
        # If class defines additional attributes, like name in Symbol,
        # then this method should be updated accordingly to return
        # relevant attributes as tuple.
        return self._args

    @property
    def precedence(self):
        return 0

    def tostr(self, level=0):
        return self.torepr()

    def torepr(self):
        l = []
        for o in self.args:
            try:
                l.append(o.torepr())
            except AttributeError:
                l.append(repr(o))
        return self.__class__.__name__ + '(' + ', '.join(l) + ')'

    def __str__(self):
        return self.tostr()

    @staticmethod
    def set_repr_level(flag = None):
        """
        Set the representation level used for repr() printing,
        returning the current level. The available levels are:
            0: Lowest level printing. Expressions printing should be
               be able to be evaluated through Python's eval()
               function
            1: Higher level printing. Expressions are printed in a
               one-dimensional fashion, are easier to read than
               level 1, but cannot be parsed through eval()
            2: Highest level printing. Expressions are simply
               two-dimensional, "pretty" versions of the expressions
               that are only useful for readability purposes.

        Notes:
        ======
            - Level 2 printing is done through the printing module in
              smpy.printing.pretty.
        """
        return repr_level(flag)

    def __repr__(self):
        plevel = repr_level()
        if plevel == 1:
            return self.tostr()
        elif plevel == 2:
            from sympy.printing.pretty import pretty
            # in fact, we should just return pretty(self) -- it would be right,
            # in the real world, the situation is somewhat complicated:
            # - there is a bug in python2.4 -- unicode result from __repr__ is
            #   wrongly handled: http://bugs.python.org/issue1459029
            # - interactive interpreter will try to encode unicode strings with
            #   sys.getdefaultencoding() encoding. site.py just deletes
            #   sys.setdefaultencoding and thus, we are out of chance to change
            #   it from 'ascii' to something unicode-aware.
            #
            #   So, by default, python is unable to handle unicode repr's in
            #   interactive sessions.
            #
            #   we could change default site.py to set default encoding based
            #   on locale, but it is not convenient to force users to change
            #   system-wide python setup.
            #
            #   It's ugly, but we are going to workaround this.
            #   See #425 for motivation.
            pstr = pretty(self)
            if isinstance(pstr, unicode):
                import sys
                try:
                    pstr = pstr.encode(sys.stdout.encoding)
                except UnicodeEncodeError:
                    print 'W: unicode problem in __repr__, will use ascii as fallback'
                    pstr = pretty(self, use_unicode=False)

            return pstr

        return self.torepr()

    def __nonzero__(self):
        # prevent using constructs like:
        #   a = Symbol('a')
        #   if a: ..
        raise AssertionError("only Equality|Unequality can define __nonzero__ method, %r" % (self.__class__))

    def compare(self, other):
        """
        Return -1,0,1 if the object is smaller, equal, or greater than other
        (not always in mathematical sense).
        If the object is of different type from other then their classes
        are ordered according to sorted_classes list.
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
