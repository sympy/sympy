""" Implementation of Basic low-level methods.
"""

import decimal
from assumptions import AssumeMeths

# used for canonical ordering of symbolic sequences
# via __cmp__ method:
ordering_of_classes = [
    # singleton numbers
    'Zero', 'One','Half','Infinity','NaN','NegativeOne','NegativeInfinity',
    # numbers
    'Integer','Rational','Real',
    # singleton symbols
    'Exp1','Pi','ImaginaryUnit',
    # symbols
    'Symbol','Wild','Temporary',
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
    'ApplyFloor', 'ApplyCeiling',
    'ApplyRe','ApplyIm', 'ApplyArg',
    'ApplySqrt','ApplyAbs','ApplySign',
    'ApplyConjugate',
    'ApplyMrvLog',
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
    # special polynomials
    'Chebyshev','Chebyshev2','UnevaluatedFactorial',
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

def interactive(flag = None, _cache=[False]):
    if flag is None:
        return _cache[0]
    old_flag = _cache[0]
    _cache[0] = bool(flag)
    return old_flag

def mycopy(obj, level=0):
    if isinstance(obj, (list, tuple)):
        return obj.__class__(map(mycopy, obj))
    elif isinstance(obj, dict):
        d = obj.__class__()
        for k,v in obj.items():
            d[mycopy(k)] = mycopy(v)
        return d
    return obj


def cache_it_fast(func):
    func._cache_it_cache = func_cache_it_cache = {}
    def wrapper(*args, **kw_args):
        if kw_args:
            keys = kw_args.keys()
            keys.sort()
            items = [(k+'=',kw_args[k]) for k in keys]
            k = args + tuple(items)
        else:
            k = args
        cache_flag = False
        try:
            r = func_cache_it_cache[k]
        except KeyError:
            r = func(*args, **kw_args)
            cache_flag = True
        if cache_flag:
            func_cache_it_cache[k] = r
        return mycopy(r)
    return wrapper

def cache_it_immutable(func):
    func._cache_it_cache = func_cache_it_cache = {}
    def wrapper(*args, **kw_args):
        if kw_args:
            keys = kw_args.keys()
            keys.sort()
            items = [(k+'=',kw_args[k]) for k in keys]
            k = args + tuple(items)
        else:
            k = args
        try:
            return func_cache_it_cache[k]
        except KeyError:
            pass
        func_cache_it_cache[k] = r = func(*args, **kw_args)
        return r
    return wrapper

def cache_it_debug(func):
    func._cache_it_cache = func_cache_it_cache = {}
    func._cache_it_cache_repr = func_cache_it_cache_repr = {}
    def wrapper(*args, **kw_args):
        if kw_args:
            keys = kw_args.keys()
            keys.sort()
            items = [(k+'=',kw_args[k]) for k in keys]
            k = args + tuple(items)
        else:
            k = args
        cache_flag = False
        try:
            r = func_cache_it_cache[k]
        except KeyError:
            r = func(*args, **kw_args)
            cache_flag = True
        if cache_flag:
            func_cache_it_cache[k] = r
            f = interactive(False)
            func_cache_it_cache_repr[k] = repr(r)
            interactive(f)
        else:
            s = func_cache_it_cache_repr[k]
            f = interactive(False)
            new_s = repr(r)
            interactive(f)
            # check that cache values have not changed
            assert new_s==s,`func,s,r, args[0].__class__`
        return mycopy(r)
    return wrapper

cache_it = cache_it_fast
#cache_it = cache_it_debug # twice slower

def cache_it_nondummy(func):
    func._cache_it_cache = func_cache_it_cache = {}
    def wrapper(*args, **kw_args):
        if kw_args:
            try:
                dummy = kw_args['dummy']
            except KeyError:
                dummy = None
            if dummy:
                return func(*args, **kw_args)
            keys = kw_args.keys()
            keys.sort()
            items = [(k+'=',kw_args[k]) for k in keys]
            k = args + tuple(items)
        else:
            k = args
        try:
            return func_cache_it_cache[k]
        except KeyError:
            pass
        func_cache_it_cache[k] = r = func(*args, **kw_args)
        return r
    return wrapper

class MetaBasicMeths(type):

    classnamespace = {}
    interactive = False        # defines the output of repr()
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

    def __getattr__(cls, name):
        try: return MetaBasicMeths.classnamespace[name]
        except KeyError: pass
        raise AttributeError("'%s' object has no attribute '%s'"%
                             (cls.__name__, name))

    def __cmp__(cls, other):
        n1 = cls.__name__
        n2 = other.__name__
        c = cmp(n1,n2)
        if not c: return 0
        try:
            i1 = ordering_of_classes.index(n1)
        except ValueError:
            print 'Add',n1,'to basic_methods.ordering_of_classes list'
            return c
        try:
            i2 = ordering_of_classes.index(n2)
        except ValueError:
            print 'Add',n2,'to basic_methods.ordering_of_classes list'
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
        try:
            return self._get_assumption(name)
        except AttributeError:
            pass

        if BasicMeths.classnamespace.has_key(name):
            return BasicMeths.classnamespace[name]
        else:
            raise AttributeError("'%s' object has no attribute '%s'"%
                                 (self.__class__.__name__, name))

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
        for o in self:
            try:
                l.append(o.torepr())
            except AttributeError:
                l.append(repr(o))
        return self.__class__.__name__ + '(' + ', '.join(l) + ')'

    def __str__(self):
        return self.tostr()

    @staticmethod
    def set_interactive(flag = None):
        return interactive(flag)

    def __repr__(self):
        if interactive():
            return self.tostr()
        return self.torepr()

    def __getitem__(self, iter):
        return self._args[iter]

    def __contains__(self, what):
        if self == what: return True
        for x in self._args:
            if what in x:
                return True
        return False

    @staticmethod
    def set_precision(prec = None):
        """
        Set precision for Decimal number operations and return previous precision value.
        """
        context = decimal.getcontext()
        oldprec = context.prec
        if prec is not None:
            context.prec = prec
        return oldprec

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
        Basic = self.__class__.Basic
        for l,r in zip(st,ot):
            if isinstance(l, Basic):
                c = l.compare(r)
            else:
                c = cmp(l, r)
            if c: return c
        return 0
