"""Base class for all objects in sympy"""

type_class = type

import decimal
from assumptions import AssumeMeths
from sympify import _sympify, sympify, SympifyError
from cache import cacheit, Memoizer, MemoizerArg

# from numbers  import Number, Integer, Rational, Real /cyclic/
# from interval import Interval /cyclic/
# from symbol   import Symbol, Wild, Temporary /cyclic/
# from add      import Add  /cyclic/
# from mul      import Mul  /cyclic/
# from power    import Pow  /cyclic/
# from function import Derivative, FunctionClass   /cyclic/

def repr_level(flag=None, _cache=[1]):
    if flag is None:
        return _cache[0]
    old_flag = _cache[0]
    _cache[0] = max(0, min(2, int(flag))) # restrict to 0,1,2
    return old_flag


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

class BasicMeta(BasicType):

    classnamespace = {}
    repr_level = 0        # defines the output of repr()
    singleton = {}

    def __init__(cls,*args,**kws):
        n = cls.__name__
        c = BasicMeta.classnamespace.get(n)
        if c is None:
            BasicMeta.classnamespace[n] = cls
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
            other = sympify(other)
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



class Basic(AssumeMeths):
    """
    Base class for all objects in sympy.

    Conventions
    ===========

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

    __slots__ = ['_mhash', '_args']

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

    def __new__(cls, *args, **assumptions):
        obj = object.__new__(cls)
        obj.assume(**assumptions)
        obj._mhash = None # will be set by __hash__ method.
        obj._args = args  # all items in args must be Basic objects
        return obj

    def __getattr__(self, name):
        # if it's not an assumption -- we don't have it
        if name[:3] != 'is_':
            # it is important to return shortly for speed reasons:
            # we have *lots* of non-'is_' attribute access, e.g.
            # '_eval_<smth>', and a lot of them does *not* exits.
            #
            # if we are here -- it surely does not exist,
            # so let's get out of here as fast as possible.
            raise AttributeError(name)

        else:
            return self._get_assumption(name)

    # NB: there is no need in protective __setattr__

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



    ##############
    # STR / REPR #
    ##############

    Lambda_precedence = 1
    Add_precedence = 40
    Mul_precedence = 50
    Pow_precedence = 60
    Apply_precedence = 70
    Item_precedence = 75
    Atom_precedence = 1000

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



    @Memoizer('Basic', MemoizerArg((type, type(None), tuple), name='type'), return_value_converter = lambda obj: obj.copy())
    def atoms(self, type=None):
        """Returns the atoms that form current object.

        Example:
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> y = Symbol('y')
        >>> (x+y**2+ 2*x*y).atoms()
        set([2, x, y])

        You can also filter the results by a given type(s) of object
        >>> (x+y+2+y**2*sin(x)).atoms(type=Symbol)
        set([x, y])

        >>> (x+y+2+y**2*sin(x)).atoms(type=Number)
        set([2])

        >>> (x+y+2+y**2*sin(x)).atoms(type=(Symbol, Number))
        set([2, x, y])
        """
        result = set()
        if type is not None and not isinstance(type, (type_class, tuple)):
            type = sympify(type).__class__
        if self.is_Atom:
            if type is None or isinstance(self, type):
                result.add(self)
        else:
            for obj in self.args:
                result = result.union(obj.atoms(type=type))
        return result

    def is_hypergeometric(self, arg):
        from sympy.simplify import hypersimp
        return hypersimp(self, arg, simplify=False) is not None

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
        for obj in self.args:
            if isinstance(obj, Basic):
                if not obj.is_number:
                    return False
            else:
                raise TypeError
        else:
            return True

    @property
    def func(self):
        """
        The top-level function in an expression.

        The following should hold for all objects:

            x = x.func(*x.args)

        """
        return self.__class__

    @property
    def args(self):
        """Returns a tuple of arguments of "self".

        Example
        -------

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

    def is_fraction(self, *syms):
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
            syms = list(self.atoms(type=Symbol))

        if not syms: # constant polynomial
            return True
        else:
            return self._eval_is_polynomial(syms)

    def as_polynomial(self, *syms, **kwargs): # remove
        return C.Polynomial(self, var=(syms or None), **kwargs)

    def as_poly(self, *symbols, **flags):
        """Converts 'self' to a polynomial or returns None.

           When constructing a polynomial an exception will be raised in
           case the input expression is not convertible to a polynomial.
           There are situations when it is easier (simpler or prettier)
           to receive None on failure.

           >>> from sympy import *
           >>> x,y = symbols('xy')

           >>> print (x**2 + x*y).as_poly(x, y)
           Poly((1, 1), ((2, 0), (1, 1)), (x, y), 'grlex')

           >>> print (x**2 + sin(y)).as_poly(x, y)
           None

        """
        from sympy.polys import Poly, PolynomialError

        try:
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
            raise Exception("subs accept either 1 or 2 arguments")

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
                if expr.is_Symbol:
                    return expr in syms
                else:
                    for term in expr.args:
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
                if expr.is_Symbol:
                    if expr in syms:
                        syms.remove(expr)
                else:
                    for term in expr.args:
                        if not syms:
                            break
                        else:
                            search(term)

            search(self)

            return not syms

    def has(self, *patterns):
        """
        Return True if self has any of the patterns.
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
            #XXX hack, this is very fragile:
            if str(self).find(str(p.__name__)) == -1:
                #didn't find p in self
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

    def _eval_power(self, other):
        return None

    def _eval_derivative(self, s):
        return

    def _eval_integral(self, s):
        return

    def _eval_defined_integral(self, s, a, b):
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

    def _calc_apply_positive(self, *args):
        return

    def _calc_apply_real(self, *args):
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

    #@classmethod
    def matches(pattern, expr, repl_dict={}, evaluate=False):
        """
        Helper method for match() - switches the pattern and expr.

        Can be used to solve linear equations:
          >>> from sympy import Symbol, Wild
          >>> a,b = map(Symbol, 'ab')
          >>> x = Wild('x')
          >>> (a+b*x).matches(0)
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

        """
        pattern = sympify(pattern)
        return pattern.matches(self, {})

    def solve4linearsymbol(eqn, rhs, symbols = None):
        """ Solve equation
          eqn == rhs
        with respect to some linear symbol in eqn.
        Returns (symbol, solution). If eqn is nonlinear
        with respect to all symbols, then return
        trivial solution (eqn, rhs).
        """
        if eqn.is_Symbol:
            return (eqn, rhs)
        if symbols is None:
            symbols = eqn.atoms(type=Symbol)
        if symbols:
            # find  symbol
            for s in symbols:
                deqn = eqn.diff(s)
                if deqn.diff(s) is S.Zero:
                    # eqn = a + b*c, a=eqn(c=0),b=deqn(c=0)
                    return s, (rhs - eqn.subs(s,0))/deqn.subs(s,0)
        # no linear symbol, return trivial solution
        return eqn, rhs

    def _calc_splitter(self, d):
        if d.has_key(self):
            return d[self]
        r = self.__class__(*[t._calc_splitter(d) for t in self])
        if d.has_key(r):
            return d[r]
        s = d[r] = Temporary()
        return s

    def splitter(self):
        d = {}
        r = self._calc_splitter(d)
        l = [(s.dummy_index,s,e) for e,s in d.items()]
        l.sort()
        return [(s,e) for i,s,e in l]

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
        return self.__class__(*terms, **self._assumptions)

    ###########################################################################
    ################# EXPRESSION REPRESENTATION METHODS #######################
    ###########################################################################

    def _eval_expand_basic(self):
        terms, rewrite = [], False

        for term in self.args:
            if not isinstance(term, Basic) or \
                   isinstance(term, Atom):
                terms.append(term)
            else:
                T = term._eval_expand_basic()

                if T is None:
                    terms.append(term)
                else:
                    terms.append(T)
                    rewrite = True

        if rewrite:
            return self.__class__(*terms, **self._assumptions)
        else:
            return None

    def _eval_expand_power(self, *args):
        if self.is_Atom:
            return self
        if not isinstance(self, C.Apply):   # FIXME Apply -> Function
            sargs = self[:]
        else:
            sargs = (self.func,)+self[:]
        terms = [ term._eval_expand_power(*args) for term in sargs ]
        return self.__class__(*terms, **self._assumptions)

    def _eval_expand_complex(self, *args):
        if self.is_Atom:
            return self
        sargs = self.args[:]
        terms = [ term._eval_expand_complex(*args) for term in sargs ]
        return self.__class__(*terms, **self._assumptions)

    def _eval_expand_trig(self, *args):
        if self.is_Atom:
            return self
        sargs = self.args[:]
        terms = [ term._eval_expand_trig(*args) for term in sargs ]
        return self.__class__(*terms, **self._assumptions)

    def _eval_expand_func(self, *args):
        if self.is_Atom:
            return self
        sargs = self.args
        terms = [ term._eval_expand_func(*args) for term in sargs ]
        return self.__class__(*terms, **self._assumptions)

    def expand(self, **hints):
        """Expand an expression using hints.

           Currently supported hints are basic, power, complex, trig
           and func.  Hints are applied with arbitrary order so your
           code shouldn't depend on the way hints are passed to this
           method. Expand 'basic' is the default and run always,
           provided that it isn't turned off by the user.

           >>> from sympy import *
           >>> x,y = symbols('xy')

           >>> (y*(x + y)**2).expand()
           y*x**2 + 2*x*y**2 + y**3

           >>> (x+y).expand(complex=True)
           I*im(x) + I*im(y) + re(x) + re(y)

        """
        expr = self

        for hint in hints:
            if hints[hint] == True:
                func = getattr(expr, '_eval_expand_'+hint, None)

                if func is not None:
                    expr = func()

        if not hints.has_key('basic'):
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
        return self.__class__(*terms, **self._assumptions)

    def rewrite(self, *args, **hints):
        """Rewrites expression containing applications of functions
           of one kind in terms of functions of different kind. For
           example you can rewrite trigonometric functions as complex
           exponentials or combinatorial functions as gamma function.

           As a pattern this function accepts a list of functions to
           to rewrite (instances of DefinedFunction class). As rule
           you can use string or a destination function instance (in
           this cas rewrite() will use tostr() method).

           There is also possibility to pass hints on how to rewrite
           the given expressions. For now there is only one such hint
           defined called 'deep'. When 'deep' is set to False it will
           forbid functions to rewrite their contents.

           >>> from sympy import *
           >>> x, y = symbols('xy')

           >>> sin(x).rewrite(sin, exp)
           -1/2*I*(exp(I*x) - exp(-I*x))

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
                    rule = rule.tostr()

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
        expr = self.expand(complex=True)

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

    def as_indep_terms(self, x):
        coeff, terms = self.as_coeff_terms()
        indeps = [coeff]
        new_terms = []
        for t in terms:
            if t.has(x):
                new_terms.append(x)
            else:
                indeps.append(x)
        return Mul(*indeps), Mul(*new_terms)

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

    def as_expr_orders(self):
        """ Split expr + Order(..) to (expr, Order(..)).
        """
        l1 = []
        l2 = []
        if self.is_Add:
            for f in self:
                if isinstance(f, C.Order):
                    l2.append(f)
                else:
                    l1.append(f)
        elif self.is_Order:
            l2.append(self)
        else:
            l1.append(self)
        return Add(*l1), Add(*l2)

    def normal(self):
        n, d = self.as_numer_denom()
        if d is S.One:
            return n
        return n/d

    ###################################################################################
    ##################### DERIVATIVE, INTEGRAL, FUNCTIONAL METHODS ####################
    ###################################################################################

    def diff(self, *symbols, **assumptions):
        new_symbols = map(sympify, symbols)
        if not assumptions.has_key("evaluate"):
            assumptions["evaluate"] = True
        ret = Derivative(self, *new_symbols, **assumptions)
        return ret

    def fdiff(self, *indices):
        # FIXME FApply -> ?
        return C.FApply(C.FDerivative(*indices), self)

    def integral(self, *symbols, **assumptions):
        new_symbols = []
        for s in symbols:
            s = sympify(s)
            if s.is_Integer and new_symbols:
                last_s = new_symbols[-1]
                i = int(s)
                new_symbols += [last_s] * (i-1)
            elif s.is_Symbol:
                new_symbols.append(s)
            else:
                raise TypeError(".integral() argument must be Symbol|Integer|Equality instance (got %s)" % (s.__class__.__name__))
        return C.Integral(self, *new_symbols, **assumptions)

    #XXX fix the removeme
    def __call__(self, *args, **removeme):
        return Function(self[0])(*args)

    def _eval_evalf(self):
        return

    def _seq_eval_evalf(self):
        return self.__class__(*[s.evalf() for s in self.args])

    def __float__(self):
        result = self.evalf(precision=16)

        if result.is_Number:
            return float(result)
        else:
            raise ValueError("Symbolic value, can't compute")

    def evalf(self, precision=None):
        if precision is None:
            r = self._eval_evalf()
        else:
            old_precision = Basic.set_precision(precision)
            r = self._eval_evalf()
            Basic.set_precision(old_precision)
        if r is None:
            r = self
        return r

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


    ###################################################################################
    ##################### SERIES, LEADING TERM, LIMIT, ORDER METHODS ##################
    ###################################################################################

    def series(self, x, point=0, n=6):
        """
        Usage
        =====
            Returns the Taylor (Laurent or generalized) series of "self" around
            the point "point" (default 0) with respect to "x" until the n-th
            term (default n is 6).

        Notes
        =====
            This method is the most high level method and it returns the
            series including the O(x**n) term.

            Internally, it executes a method oseries(), which takes an
            O instance as the only parameter and it is responsible for
            returning a series (without the O term) up to the given order.
        """
        x = sympify(x)
        point = sympify(point)
        if point != 0:
            raise NotImplementedError("series expansion around arbitrary point")
            #self = self.subs(x, x + point)
        o = C.Order(x**n,x)
        r = self.oseries(o)
        if r==self:
            return self
        return r + o

    @cacheit
    def oseries(self, order):
        """
        Return the series of an expression up to given Order symbol (without the
        actual O term).

        The general philosophy is this: simply start with the most simple
        Taylor (Laurent) term and calculate one be one and use
        order.contains(term) method to determine if your term is still
        significant and should be added to the series, or we should stop.
        """
        order = C.Order(order)
        if order is S.Zero:
            return self
        if order.contains(self):
            return S.Zero
        if len(order.symbols)>1:
            r = self
            for s in order.symbols:
                o = C.Order(order.expr, s)
                r = r.oseries(o)
            return r
        x = order.symbols[0]
        if not self.has(x):
            return self
        obj = self._eval_oseries(order)
        if obj is not None:
            #obj2 = obj.expand(trig=True)
            obj2 = obj.expand()
            if obj2 != obj:
                r = obj2.oseries(order)
                return r
            return obj
        raise NotImplementedError('(%s).oseries(%s)' % (self, order))

    def nseries(self, x, x0, n):
        """
        Calculates a generalized series expansion.

        The difference between oseries and nseries is that nseries calculates
        "n" terms in the innermost expressions and then builds up the final
        series just by "cross-mutliplying" everything out.

        Advantage -- it's fast, because we don't have to determine how many
        terms we need to calculate in advance.

        Disadvantage -- you may endup with less terms than you may have
        expected, but the O(x**n) term appended will always be correct, so the
        result is correct, but maybe shorter.
        """
        raise NotImplementedError("(%s).nseries(%s, %s, %s)" % (self, x, x0, n))

    def _eval_oseries(self, order):
        return

    def _compute_oseries(self, arg, order, taylor_term, unevaluated_func, correction = 0):
        """
        compute series sum(taylor_term(i, arg), i=0..n-1) such
        that order.contains(taylor_term(n, arg)). Assumes that arg->0 as x->0.
        """
        x = order.symbols[0]
        ln = C.log
        o = C.Order(arg, x)
        if o is S.Zero:
            return unevaluated_func(arg)
        if o.expr==1:
            e = ln(order.expr*x)/ln(x)
        else:
            e = ln(order.expr)/ln(o.expr)
        n = e.limit(x,0) + 1 + correction
        if n.is_unbounded:
            # requested accuracy gives infinite series,
            # order is probably nonpolynomial e.g. O(exp(-1/x), x).
            return unevaluated_func(arg)
        try:
            n = int(n)
        except TypeError:
            #well, the n is something more complicated (like 1+log(2))
            n = int(n.evalf()) + 1
        assert n>=0,`n`
        l = []
        g = None
        for i in xrange(n+2):
            g = taylor_term(i, arg, g)
            g = g.oseries(order)
            l.append(g)
        return Add(*l)

    def limit(self, x, xlim, direction='+'):
        """ Compute limit x->xlim.
        """
        from sympy.series.limits import limit
        return limit(self, x, xlim, direction)

    @cacheit
    def as_leading_term(self, *symbols):
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
        expr = self.expand(trig=True)
        obj = expr._eval_as_leading_term(x)
        if obj is not None:
            return obj
        raise NotImplementedError('as_leading_term(%s, %s)' % (self, x))

    def as_coeff_exponent(self, x):
        """ c*x**e -> c,e where x can be any symbolic expression.
        """
        x = sympify(x)
        wc = Wild('wc')
        we = Wild('we')
        c, terms = self.as_coeff_terms()
        p  = wc*x**we
        d = self.match(p)
        if d is not None:
            return d[wc], d[we]
        return self, S.Zero

    def ldegree(self, x):
        x = sympify(x)
        c,e = self.as_leading_term(x).as_coeff_exponent(x)
        if not c.has(x):
            return e
        raise ValueError("cannot compute ldegree(%s, %s), got c=%s" % (self, x, c))

    def leadterm(self, x):
        x = sympify(x)
        c,e = self.as_leading_term(x).as_coeff_exponent(x)
        if not c.has(x):
            return c,e
        raise ValueError("cannot compute ldegree(%s, %s), got c=%s" % (self, x, c))


    ##########################################################################
    ##################### END OF BASIC CLASS #################################
    ##########################################################################

class Atom(Basic):
    """
    A parent class for atomic things.

    Examples: Symbol, Number, Rational, Integer, ...
    But not: Add, Mul, Pow, ...
    """

    is_Atom = True

    precedence = Basic.Atom_precedence

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

    def _calc_splitter(self, d):
        return self

    def count_ops(self, symbolic=True):
        return S.Zero

    def doit(self, **hints):
        return self

    def _eval_integral(self, s):
        if s==self:
            return self**2/2
        return self*s

    def _eval_defined_integral(self, s, a, b):
        if s==self:
            return (b**2-a**2)/2
        return self*(b-a)

    def _eval_is_polynomial(self, syms):
        return True

    def _eval_oseries(self, order):
        # .oseries() method checks for order.contains(self)
        return self

    def _eval_as_leading_term(self, x):
        return self

    @property
    def is_number(self):
        return True

    def nseries(self, x, x0, n):
        return self


class Singleton(Basic):
    """ Singleton object.
    """

    __slots__ = []

    def __new__(cls, *args, **assumptions):
        # if you need to overload __new__, then
        # use the same code as below to ensure
        # that only one instance of Singleton
        # class is created.
        obj = Singleton.__dict__.get(cls.__name__)
        if obj is None:
            obj = Basic.__new__(cls,*args,**assumptions)
            setattr(Singleton, cls.__name__, obj)
        return obj

class SingletonFactory:
    """
    A map between singleton classes and the corresponding instances.
    E.g. S.Exp == C.Exp()
    """

    def __getattr__(self, clsname):
        if clsname == "__repr__":
            return lambda: "S"
        obj = Singleton.__dict__.get(clsname)
        if obj is None:
            cls = getattr(C, clsname)
            assert issubclass(cls, Singleton),`cls`
            obj = cls()

        # store found object in own __dict__, so the next lookups will be
        # serviced without entering __getattr__, and so will be fast
        setattr(self, clsname, obj)
        return obj

S = SingletonFactory()

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
