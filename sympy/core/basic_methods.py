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
