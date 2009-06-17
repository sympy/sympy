"""
There are two types of functions:
1) defined function like exp or sin that has a name and body
   (in the sense that function can be evaluated).
    e = exp
2) undefined function with a name but no body. Undefined
  functions can be defined using a Function class as follows:
    f = Function('f')
  (the result will be Function instance)
3) this isn't implemented yet: anonymous function or lambda function that has
no name but has body with dummy variables. An anonymous function
   object creation examples:
    f = Lambda(x, exp(x)*x)
    f = Lambda(exp(x)*x)  # free symbols in the expression define the number of arguments
    f = exp * Lambda(x,x)
4) isn't implemented yet: composition of functions, like (sin+cos)(x), this
works in sympycore, but needs to be ported back to SymPy.


Example:
    >>> from sympy import *
    >>> f = Function("f")
    >>> x = Symbol("x")
    >>> f(x)
    f(x)
    >>> print srepr(f(x).func)
    Function('f')
    >>> f(x).args
    (x,)
"""

from basic import Basic, Atom, S, C, sympify
from basic import BasicType, BasicMeta
from operations import AssocOp
from cache import cacheit
from itertools import repeat
from numbers import Rational, Integer
from symbol import Symbol
from add    import Add
from multidimensional import vectorize

from sympy import mpmath
from sympy.utilities.decorator import deprecated

class PoleError(Exception):
    pass

class FunctionClass(BasicMeta):
    """
    Base class for function classes. FunctionClass is a subclass of type.

    Use Function('<function name>' [ , signature ]) to create
    undefined function classes.
    """

    _new = type.__new__

    def __new__(cls, arg1, arg2, arg3=None, **options):
        assert not options,`options`
        if isinstance(arg1, type):
            # the following code gets executed when one types
            # FunctionClass(Function, "f")
            # i.e. cls = FunctionClass, arg1 = Function, arg2 = "f"
            # and we simply do an equivalent of:
            # class f(Function):
            #     ...
            # return f
            ftype, name, signature = arg1, arg2, arg3
            #XXX this probably needs some fixing:
            assert ftype.__name__.endswith('Function'),`ftype`
            attrdict = ftype.__dict__.copy()
            attrdict['undefined_Function'] = True
            if signature is not None:
                attrdict['signature'] = signature
            bases = (ftype,)
            return type.__new__(cls, name, bases, attrdict)
        else:
            name, bases, attrdict = arg1, arg2, arg3
            return type.__new__(cls, name, bases, attrdict)

    def __repr__(cls):
        return cls.__name__

class Function(Basic):
    """
    Base class for applied functions.
    Constructor of undefined classes.

    """

    __metaclass__ = FunctionClass

    is_Function = True

    nargs = None

    @vectorize(1)
    @cacheit
    def __new__(cls, *args, **options):
        # NOTE: this __new__ is twofold:
        #
        # 1 -- it can create another *class*, which can then be instantiated by
        #      itself e.g. Function('f') creates a new class f(Function)
        #
        # 2 -- on the other hand, we instantiate -- that is we create an
        #      *instance* of a class created earlier in 1.
        #
        # So please keep, both (1) and (2) in mind.

        # (1) create new function class
        #     UC: Function('f')
        if cls is Function:
            #when user writes Function("f"), do an equivalent of:
            #taking the whole class Function(...):
            #and rename the Function to "f" and return f, thus:
            #In [13]: isinstance(f, Function)
            #Out[13]: False
            #In [14]: isinstance(f, FunctionClass)
            #Out[14]: True

            if len(args) == 1 and isinstance(args[0], str):
                #always create Function
                return FunctionClass(Function, *args)
                return FunctionClass(Function, *args, **options)
            else:
                print args
                print type(args[0])
                raise TypeError("You need to specify exactly one string")

        # (2) create new instance of a class created in (1)
        #     UC: Function('f')(x)
        #     UC: sin(x)
        args = map(sympify, args)
        # these lines should be refactored
        for opt in ["nargs", "dummy", "comparable", "noncommutative", "commutative"]:
            if opt in options:
                del options[opt]
        # up to here.
        if options.get('evaluate') is False:
            return Basic.__new__(cls, *args, **options)
        evaluated = cls.eval(*args)
        if evaluated is not None: return evaluated
        # Just undefined functions have nargs == None
        if not cls.nargs and hasattr(cls, 'undefined_Function'):
            r = Basic.__new__(cls, *args, **options)
            r.nargs = len(args)
            return r
        return Basic.__new__(cls, *args, **options)

    @property
    def is_commutative(self):
        return True

    @classmethod
    @deprecated
    def canonize(cls, *args):
        return cls.eval(*args)

    @classmethod
    def eval(cls, *args):
        """
        Returns a canonical form of cls applied to arguments args.

        The eval() method is called when the class cls is about to be
        instantiated and it should return either some simplified instance
        (possible of some other class), or if the class cls should be
        unmodified, return None.

        Example of eval() for the function "sign"
        ---------------------------------------------

        @classmethod
        def eval(cls, arg):
            if arg is S.NaN:
                return S.NaN
            if arg is S.Zero: return S.Zero
            if arg.is_positive: return S.One
            if arg.is_negative: return S.NegativeOne
            if isinstance(arg, C.Mul):
                coeff, terms = arg.as_coeff_terms()
                if coeff is not S.One:
                    return cls(coeff) * cls(C.Mul(*terms))

        """
        return

    @property
    def func(self):
        return self.__class__

    def _eval_subs(self, old, new):
        if self == old:
            return new
        elif old.is_Function and new.is_Function:
            if old == self.func:
                if self.nargs is new.nargs or not new.nargs:
                    return new(*self.args[:])
                # Written down as an elif to avoid a super-long line
                elif isinstance(new.nargs,tuple) and self.nargs in new.nargs:
                    return new(*self.args[:])
        obj = self.func._eval_apply_subs(*(self.args[:] + (old,) + (new,)))
        if obj is not None:
            return obj
        return Basic._seq_subs(self, old, new)

    def _eval_evalf(self, prec):
        # Lookup mpmath function based on name
        fname = self.func.__name__
        try:
            if not hasattr(mpmath, fname):
                from sympy.utilities.lambdify import MPMATH_TRANSLATIONS
                fname = MPMATH_TRANSLATIONS[fname]
            func = getattr(mpmath, fname)
        except (AttributeError, KeyError):
            return

        # Convert all args to mpf or mpc
        try:
            args = [arg._to_mpmath(prec) for arg in self.args]
        except ValueError:
            return

        # Set mpmath precision and apply. Make sure precision is restored
        # afterwards
        orig = mpmath.mp.prec
        try:
            mpmath.mp.prec = prec
            v = func(*args)
        finally:
            mpmath.mp.prec = orig

        return Basic._from_mpmath(v, prec)

    def _eval_is_comparable(self):
        if self.is_Function:
            r = True
            for s in self.args:
                c = s.is_comparable
                if c is None: return
                if not c: r = False
            return r
        return

    def _eval_derivative(self, s):
        # f(x).diff(s) -> x.diff(s) * f.fdiff(1)(s)
        i = 0
        l = []
        r = S.Zero
        for a in self.args:
            i += 1
            da = a.diff(s)
            if da is S.Zero:
                continue
            if isinstance(self.func, FunctionClass):
                df = self.fdiff(i)
                l.append(df * da)
        return Add(*l)

    def _eval_is_commutative(self):
        r = True
        for a in self._args:
            c = a.is_commutative
            if c is None: return None
            if not c: r = False
        return r

    def _eval_eq_nonzero(self, other):
        if isinstance(other.func, self.func.__class__) and len(self)==len(other):
            for a1,a2 in zip(self,other):
                if not (a1==a2):
                    return False
            return True

    def as_base_exp(self):
        return self, S.One

    def count_ops(self, symbolic=True):
        #      f()             args
        return 1   + Add(*[t.count_ops(symbolic) for t in self.args])

    def _eval_nseries(self, x, x0, n):
        assert len(self.args) == 1
        arg = self.args[0]
        arg0 = arg.limit(x, 0)
        from sympy import oo
        if arg0 in [-oo, oo]:
            raise PoleError("Cannot expand around %s" % (arg))
        if arg0 is not S.Zero:
            e = self.func(arg)
            e1 = e.expand()
            if e==e1:
                #for example when e = sin(x+1) or e = sin(cos(x))
                #let's try the general algorithm
                term = e.subs(x, S.Zero)
                series = term
                fact = S.One
                for i in range(n-1):
                    i += 1
                    fact *= Rational(i)
                    e = e.diff(x)
                    term = e.subs(x, S.Zero)*(x**i)/fact
                    term = term.expand()
                    series += term
                return series + C.Order(x**n, x)
            return self.nseries(x, x0, n)
        l = []
        g = None
        for i in xrange(n+2):
            g = self.taylor_term(i, arg, g)
            g = g.nseries(x, x0, n)
            l.append(g)
        return Add(*l) + C.Order(x**n, x)

    def _eval_is_polynomial(self, syms):
        for arg in self.args:
            if arg.has(*syms):
                return False
        return True

    def _eval_expand_basic(self, deep=True, **hints):
        if not deep:
            return self
        sargs, terms = self.args[:], []
        for term in sargs:
            if hasattr(term, '_eval_expand_basic'):
                newterm = term._eval_expand_basic(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_power_exp(self, deep=True, **hints):
        if not deep:
            return self
        sargs, terms = self.args[:], []
        for term in sargs:
            if hasattr(term, '_eval_expand_power_exp'):
                newterm = term._eval_expand_power_exp(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_power_base(self, deep=True, **hints):
        if not deep:
            return self
        sargs, terms = self.args[:], []
        for term in sargs:
            if hasattr(term, '_eval_expand_power_base'):
                newterm = term._eval_expand_power_base(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_mul(self, deep=True, **hints):
        if not deep:
            return self
        sargs, terms = self.args[:], []
        for term in sargs:
            if hasattr(term, '_eval_expand_mul'):
                newterm = term._eval_expand_mul(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_multinomial(self, deep=True, **hints):
        if not deep:
            return self
        sargs, terms = self.args[:], []
        for term in sargs:
            if hasattr(term, '_eval_expand_multinomail'):
                newterm = term._eval_expand_multinomial(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_log(self, deep=True, **hints):
        if not deep:
            return self
        sargs, terms = self.args[:], []
        for term in sargs:
            if hasattr(term, '_eval_expand_log'):
                newterm = term._eval_expand_log(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_complex(self, deep=True, **hints):
        if deep:
            func = self.func(*[ a.expand(deep, **hints) for a in self.args ])
        else:
            func = self.func(*self.args)
        return C.re(func) + S.ImaginaryUnit * C.im(func)

    def _eval_expand_trig(self, deep=True, **hints):
        sargs, terms = self.args[:], []
        for term in sargs:
            if hasattr(term, '_eval_expand_trig'):
                newterm = term._eval_expand_trig(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_func(self, deep=True, **hints):
        sargs, terms = self.args[:], []
        for term in sargs:
            if hasattr(term, '_eval_expand_func'):
                newterm = term._eval_expand_func(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_rewrite(self, pattern, rule, **hints):
        if hints.get('deep', False):
            args = [ a._eval_rewrite(pattern, rule, **hints) for a in self ]
        else:
            args = self.args[:]

        if pattern is None or isinstance(self.func, pattern):
            if hasattr(self, rule):
                rewritten = getattr(self, rule)(*args)

                if rewritten is not None:
                    return rewritten

        return self.func(*args, **self._assumptions)

    def fdiff(self, argindex=1):
        if self.nargs is not None:
            if isinstance(self.nargs, tuple):
                nargs = self.nargs[-1]
            else:
                nargs = self.nargs
            if not (1<=argindex<=nargs):
                raise TypeError("argument index %r is out of range [1,%s]" % (argindex,nargs))
        return Derivative(self,self.args[argindex-1],evaluate=False)

    @classmethod
    def _eval_apply_evalf(cls, arg):
        arg = arg.evalf(prec)

        #if cls.nargs == 1:
        # common case for functions with 1 argument
        #if arg.is_Number:
        if arg.is_number:
            func_evalf = getattr(arg, cls.__name__)
            return func_evalf()

    def _eval_as_leading_term(self, x):
        """General method for the leading term"""
        arg = self.args[0].as_leading_term(x)

        if C.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

    @classmethod
    def taylor_term(cls, n, x, *previous_terms):
        """General method for the taylor term.

        This method is slow, because it differentiates n-times.  Subclasses can
        redefine it to make it faster by using the "previous_terms".
        """
        x = sympify(x)
        return cls(x).diff(x, n).subs(x, 0) * x**n / C.Factorial(n)

class WildFunction(Function, Atom):
    """
    WildFunction() matches any expression but another WildFunction()
    XXX is this as intended, does it work ?
    """

    nargs = 1

    def __new__(cls, name=None, **assumptions):
        if name is None:
            name = 'Wf%s' % (Symbol.dummycount + 1) # XXX refactor dummy counting
            Symbol.dummycount += 1
        obj = Function.__new__(cls, name, **assumptions)
        obj.name = name
        return obj

    def matches(pattern, expr, repl_dict={}, evaluate=False):
        for p,v in repl_dict.items():
            if p==pattern:
                if v==expr: return repl_dict
                return None
        if pattern.nargs is not None:
            if not hasattr(expr,'nargs') or pattern.nargs != expr.nargs:
                return None
        repl_dict = repl_dict.copy()
        repl_dict[pattern] = expr
        return repl_dict

    @classmethod
    def _eval_apply_evalf(cls, arg):
        return

    @property
    def is_number(self):
        return False

class Derivative(Basic):
    """
    Carries out differentation of the given expression with respect to symbols.

    expr must define ._eval_derivative(symbol) method that returns
    the differentation result or None.

    Examples:

    Derivative(Derivative(expr, x), y) -> Derivative(expr, x, y)
    Derivative(expr, x, 3)  -> Derivative(expr, x, x, x)
    """

    is_Derivative   = True

    @staticmethod
    def _symbolgen(*symbols):
        """
        Generator of all symbols in the argument of the Derivative.

        Example:
        >> ._symbolgen(x, 3, y)
        (x, x, x, y)
        >> ._symbolgen(x, 10**6)
        (x, x, x, x, x, x, x, ...)

        The second example shows why we don't return a list, but a generator,
        so that the code that calls _symbolgen can return earlier for special
        cases, like x.diff(x, 10**6).
        """
        last_s = sympify(symbols[len(symbols)-1])
        for i in xrange(len(symbols)):
            s = sympify(symbols[i])
            next_s = None
            if s != last_s:
                next_s = sympify(symbols[i+1])

            if isinstance(s, Integer):
                continue
            elif isinstance(s, Symbol):
                # handle cases like (x, 3)
                if isinstance(next_s, Integer):
                    # yield (x, x, x)
                    for copy_s in repeat(s,int(next_s)):
                        yield copy_s
                else:
                    yield s
            else:
                yield s

    def __new__(cls, expr, *symbols, **assumptions):
        expr = sympify(expr)
        if not symbols: return expr
        symbols = Derivative._symbolgen(*symbols)
        if expr.is_commutative:
            assumptions["commutative"] = True
        if "evaluate" in assumptions:
            evaluate = assumptions["evaluate"]
            del assumptions["evaluate"]
        else:
            evaluate = False
        if not evaluate and not isinstance(expr, Derivative):
            obj = Basic.__new__(cls, expr, *symbols, **assumptions)
            return obj
        unevaluated_symbols = []
        for s in symbols:
            s = sympify(s)
            if not isinstance(s, Symbol):
                raise ValueError('Invalid literal: %s is not a valid variable' % s)
            if not expr.has(s):
                return S.Zero
            obj = expr._eval_derivative(s)
            if obj is None:
                unevaluated_symbols.append(s)
            elif obj is S.Zero:
                return S.Zero
            else:
                expr = obj

        if not unevaluated_symbols:
            return expr
        return Basic.__new__(cls, expr, *unevaluated_symbols, **assumptions)

    def _eval_derivative(self, s):
        #print
        #print self
        #print s
        #stop
        if s not in self.symbols:
            obj = self.expr.diff(s)
            if isinstance(obj, Derivative):
                return Derivative(obj.expr, *(obj.symbols+self.symbols))
            return Derivative(obj, *self.symbols)
        return Derivative(self.expr, *((s,)+self.symbols), **{'evaluate': False})

    def doit(self):
        return Derivative(self.expr, *self.symbols,**{'evaluate': True})

    @property
    def expr(self):
        return self._args[0]

    @property
    def symbols(self):
        return self._args[1:]

    def _eval_subs(self, old, new):
        if self==old:
            return new
        return Derivative(*map(lambda x: x._eval_subs(old, new), self.args), **{'evaluate': True})

    def matches(pattern, expr, repl_dict={}, evaluate=False):
        # this method needs a cleanup.

        #print "?   :",pattern, expr, repl_dict, evaluate
        #if repl_dict:
        #    return repl_dict
        for p,v in repl_dict.items():
            if p==pattern:
                if v==expr: return repl_dict
                return None
        assert isinstance(pattern, Derivative)
        if isinstance(expr, Derivative):
            if len(expr.symbols) == len(pattern.symbols):
                    #print "MAYBE:",pattern, expr, repl_dict, evaluate
                    return Basic.matches(pattern, expr, repl_dict, evaluate)
        #print "NONE:",pattern, expr, repl_dict, evaluate
        return None
        #print pattern, expr, repl_dict, evaluate
        stop
        if pattern.nargs is not None:
            if pattern.nargs != expr.nargs:
                return None
        repl_dict = repl_dict.copy()
        repl_dict[pattern] = expr
        return repl_dict

    def _eval_lseries(self, x, x0):
        stop
        arg = self.args[0]
        dx = self.args[1]
        for term in arg.lseries(x, x0):
            yield term.diff(dx)

    def _eval_nseries(self, x, x0, n):
        arg = self.args[0]
        arg = arg.nseries(x, x0, n)
        o = arg.getO()
        dx = self.args[1]
        if o:
            return arg.removeO().diff(dx) + arg.getO()/dx
        else:
            return arg.removeO().diff(dx)

class Lambda(Function):
    """
    Lambda(x, expr) represents a lambda function similar to Python's
    'lambda x: expr'. A function of several variables is written as
    Lambda((x, y, ...), expr).

    A simple example:
        >>> from sympy import Symbol
        >>> x = Symbol('x')
        >>> f = Lambda(x, x**2)
        >>> f(4)
        16

    For multivariate functions, use:
        >>> x = Symbol('x')
        >>> y = Symbol('y')
        >>> z = Symbol('z')
        >>> t = Symbol('t')
        >>> f2 = Lambda(x,y,z,t,x+y**z+t**z)
        >>> f2(1,2,3,4)
        73

    Multivariate functions can be curries for partial applications:
        >>> sum2numbers = Lambda(x,y,x+y)
        >>> sum2numbers(1,2)
        3
        >>> plus1 = sum2numbers(1)
        >>> plus1(3)
        4

    A handy shortcut for lots of arguments:
        >>> from sympy import *
        >>> var('x y z')
        (x, y, z)
        >>> p = x, y, z
        >>> f = Lambda(p, x + y*z)
        >>> f(*p)
        x + y*z

    """

    # a minimum of 2 arguments (parameter, expression) are needed
    nargs = 2
    def __new__(cls,*args):
        assert len(args) >= 2,"Must have at least one parameter and an expression"
        if len(args) == 2 and isinstance(args[0], (list, tuple)):
            args = tuple(args[0])+(args[1],)
        obj = Function.__new__(cls,*args)
        obj.nargs = len(args)-1
        return obj

    @classmethod
    @deprecated
    def canonize(cls, *args):
        return cls.eval(*args)

    @classmethod
    def eval(cls,*args):
        obj = Basic.__new__(cls, *args)
        #use dummy variables internally, just to be sure
        nargs = len(args)-1

        expression = args[nargs]
        funargs = [Symbol(arg.name, dummy=True) for arg in args[:nargs]]
        #probably could use something like foldl here
        for arg,funarg in zip(args[:nargs],funargs):
            expression = expression.subs(arg,funarg)
        funargs.append(expression)
        obj._args = tuple(funargs)

        return obj

    def apply(self, *args):
        """Applies the Lambda function "self" to the arguments given.
        This supports partial application.

        Example:
            >>> from sympy import Symbol
            >>> x = Symbol('x')
            >>> y = Symbol('y')
            >>> f = Lambda(x, x**2)
            >>> f.apply(4)
            16
            >>> sum2numbers = Lambda(x,y,x+y)
            >>> sum2numbers(1,2)
            3
            >>> plus1 = sum2numbers(1)
            >>> plus1(3)
            4

        """

        nparams = self.nargs
        assert nparams >= len(args),"Cannot call function with more parameters than function variables: %s (%d variables) called with %d arguments" % (str(self),nparams,len(args))


        #replace arguments
        expression = self.args[self.nargs]
        for arg,funarg in zip(args,self.args[:nparams]):
            expression = expression.subs(funarg,arg)

        #curry the rest
        if nparams != len(args):
            unused_args = list(self.args[len(args):nparams])
            unused_args.append(expression)
            return Lambda(*tuple(unused_args))
        return expression

    def __call__(self, *args):
        return self.apply(*args)

    def __eq__(self, other):
        if isinstance(other, Lambda):
            if not len(self.args) == len(other.args):
                return False

            selfexpr = self.args[self.nargs]
            otherexpr = other.args[other.nargs]
            for selfarg,otherarg in zip(self.args[:self.nargs],other.args[:other.nargs]):
                otherexpr = otherexpr.subs(otherarg,selfarg)
            if selfexpr == otherexpr:
                return True
           # if self.args[1] == other.args[1].subs(other.args[0], self.args[0]):
           #     return True
        return False

@vectorize(0,1,2)
def diff(f, x, times = 1, evaluate=True):
    """Differentiate f with respect to x

    It's just a wrapper to unify .diff() and the Derivative class,
    it's interface is similar to that of integrate()

    see http://documents.wolfram.com/v5/Built-inFunctions/AlgebraicComputation/Calculus/D.html
    """

    return Derivative(f,x,times, **{'evaluate':evaluate})

@vectorize(0)
def expand(e, deep=True, power_base=True, power_exp=True, mul=True, \
           log=True, multinomial=True, basic=True, **hints):
    """
    Expand an expression using hints.

    Currently supported hints are log, power_exp, power_base,
    multinomial, log, mul, complex, trig, func, and basic.  log, power_exp,
    power_base, multinomial, mul, and basic are evaluated by default.
    If you want to not evaluate those expansions, you need to set them
    all to false.

    basic is a generic keyword for methods that want to be expanded
    automatically.  Currently, no class that is part of SymPy uses
    basic.  If you want your class expand methods to run automatically and
    they don't fit one of the already automatic methods, wrap it around
    _eval_expand_basic.

    Also see expand_log, expand_mul, expand_complex,
    expand_trig, and expand_func, which are wrappers around those
    expansion methods.  Hints are applied with arbitrary order so your
    code shouldn't depend on the way hints are passed to this
    method.

    If deep is set to True, things like arguments of functions are
    recursively expanded.  User deep=False to only expand on the top
    level.

    >>> from sympy import *
    >>> x,y = symbols('xy')

    mul - Distributes multiplication over addition.
    >>> (y*(x + z)).expand(mul=True)
    x*y + y*z

    complex - Split and expression into real and imaginary parts
    >>> (x+y).expand(complex=True)
    I*im(x) + I*im(y) + re(x) + re(y)
    >>> cos(x).expand(complex=True)
    cos(re(x))*cosh(im(x)) - I*sin(re(x))*sinh(im(x))

    power_exp - Expand addition in exponents into multiplied bases.
    >>> exp(x+y).expand(power_exp=True)
    exp(x)*exp(y)
    >>> (2**(x+y)).expand(power_exp=True
    2**x*2**y

    power_base - Split powers of multiplied bases.
    >>> ((x*y)**z)).expand(power_base=True)
    x**z*y**z

    log - Pull out powers of logs as exponents and split multiplied arguments
    of logs into addition of logs.  Note that this only works if the
    arguments of the log function have the proper assumptions such that this
    is a correct expansion.  Namely, the arguments must be positive and the
    exponents must be real.
    >>> log(x**2*y).expand(log=True)
    log(x**2*y)
    >>> x, y = symbols('xy', positive=True)
    >>> log(x**2*y).expand(log=True)
    2*log(x)+log(y)

    trig - Do trigonometric expansion.
    >>> cos(x+y).expand(trig=True)
    cos(x)*cos(y) - sin(x)*sin(y)

    func - Expand other functions
    >>> gamma(x+1).expand(func=True)
    x*gamma(x)

    multinomial - Expand (x + y + ...)**n where n is a positive integer
    >>> ((x+y+z)**2).expand(multinomial=True)
    2*x*y + 2*x*z + 2*y*z + x**2 + y**2 + z**2

    You can pick and choose the methods.
    >>> (exp(x+y)*(x+y)).expand()
    x*exp(x)*exp(y) + y*exp(x)*exp(y)
    >>> (exp(x+y)*(x+y)).expand(power_exp=False)
    x*exp(x + y) + y*exp(x + y)
    >>> (exp(x+y)*(x+y)).expand(mul=False)
    (x + y)*exp(x)*exp(y)

    Use deep=False to only expand on the top level.
    >>> exp(x+exp(x+y)).expand()
    exp(x)*exp(exp(x)*exp(y))
    >>> exp(x+exp(x+y)).expand(deep=False)
    exp(x)*exp(exp(x + y))

    Note that because hints are applied in arbitrary order, some hints may
    prevent expansion by other hints if they are applied first.  In
    particular, mul may distribute multiplications and prevent log and
    power_base from expanding them.  Also, if mul is applied before multinomial,
    the expression might not be fully distributed.  The solution is to expand
    with mul=False first, then run expand_mul if you need further expansion.
    Examples:
    >>> x, y, z = symbols('xyz', positive=True)
    #>>> expand(log(x*(y+z))) # could be either one below
    log(x*y + x*z)
    log(x) + log(y + z)
    >>> expand_log(log(x*y+x*z))
    log(x*y+x*z)
    >>> expand(log(x*(y+z), mul=False)
    log(x) + log(y + z)

    #>>> expand((x*(y+z))**x) # could be either one below
    (x*y + x*z)**x
    x**x*(y + z)**x
    >>> expand((x*(y+z))**x, mul=False)
    x**x*(y + z)**x

    #>>> expand(x*(y+z)**2) # could be either one below
    2*x*y*z + x*y**2 + x*z**2
    x*(y + z)**2
    >>> expand(x*(y+z)**2, mul=False
    x*(2*y*z + y**2 + z**2)
    >>> expand_mul(_)
    2*x*y*z + x*y**2 + x*z**2
       """
    hints['power_base'] = power_base
    hints['power_exp'] = power_exp
    hints['mul'] = mul
    hints['log'] = log
    hints['multinomial'] = multinomial
    hints['basic'] = basic
    return sympify(e).expand(deep=deep, **hints)

# These are simple wrappers around single hints.  Fell free to add ones for
# power_exp, power_base, multinomial, or basic if you need them.
def expand_mul(expr, deep=True):
    """
    Wrapper around expand that only uses the mul hint.  See the expand
    docstring for more information.

    Example:
    >>> from sympy import *
    >>> x, y = symbols('xy', positive=True)
    >>> expand_mul(exp(x+y)*(x+y)*log(x*y**2))
    x*exp(x + y)*log(x*y**2) + y*exp(x + y)*log(x*y**2)
    """
    return sympify(expr).expand(deep=deep, mul=True, power_exp=False,\
    power_base=False, basic=False, multinomial=False, log=False)

def expand_log(expr, deep=True):
    """
    Wrapper around expand that only uses the log hint.  See the expand
    docstring for more information.

    Example:
    >>> from sympy import *
    >>> x, y = symbols('xy', positive=True)
    >>> expand_log(exp(x+y)*(x+y)*log(x*y**2))
    (x + y)*(2*log(y) + log(x))*exp(x + y)
    """
    return sympify(expr).expand(deep=deep, log=True, mul=False,\
    power_exp=False, power_base=False, multinomial=False, basic=False)

def expand_func(expr, deep=True):
    """
    Wrapper around expand that only uses the func hint.  See the expand
    docstring for more information.

    Example:
    >>> from sympy import *
    >>> x = Symbol('x')
    >>> expand_func(gamma(x + 2))
    x*(1 + x)*gamma(x)
    """
    return sympify(expr).expand(deep=deep, func=True, basic=False,\
    log=False, mul=False, power_exp=False, power_base=False, multinomial=False)

def expand_trig(expr, deep=True):
    """
    Wrapper around expand that only uses the trig hint.  See the expand
    docstring for more information.

    Example:
    >>> from sympy import *
    >>> x, y = symbols('xy')
    >>> expand_trig(sin(x+y)*(x+y))
    (x + y)*(cos(x)*sin(y) + cos(y)*sin(x))
    """
    return sympify(expr).expand(deep=deep, trig=True, basic=False,\
    log=False, mul=False, power_exp=False, power_base=False, multinomial=False)

def expand_complex(expr, deep=True):
    """
    Wrapper around expand that only uses the complex hint.  See the expand
    docstring for more information.

    Example:
    >>> from sympy import *
    >>> z = Symbol('z')
    >>> expand_complex(z**(2*I))
    I*im(z**(2*I)) + re(z**(2*I))
    """
    return sympify(expr).expand(deep=deep, complex=True, basic=False,\
    log=False, mul=False, power_exp=False, power_base=False, multinomial=False)



# /cyclic/
import basic as _
_.Derivative    = Derivative
_.FunctionClass = FunctionClass
del _

import add as _
_.FunctionClass = FunctionClass
del _

import mul as _
_.FunctionClass = FunctionClass
_.WildFunction  = WildFunction
del _

import operations as _
_.Lambda        = Lambda
_.WildFunction  = WildFunction
del _

import symbol as _
_.Function      = Function
_.WildFunction  = WildFunction
del _

import numbers as _
_.FunctionClass = FunctionClass
del _
