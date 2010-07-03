"""
There are two types of functions:
1) defined function like exp or sin that has a name and body
   (in the sense that function can be evaluated).
    e = exp
2) undefined function with a name but no body. Undefined
   functions can be defined using a Function class as follows:
       f = Function('f')
   (the result will be a Function instance)
3) this isn't implemented yet: anonymous function or lambda function that has
   no name but has body with dummy variables. Examples of anonymous function
   creation:
       f = Lambda(x, exp(x)*x)
       f = Lambda(exp(x)*x)  # free symbols in the expression define the number of arguments
       f = exp * Lambda(x,x)
4) isn't implemented yet: composition of functions, like (sin+cos)(x), this
   works in sympy core, but needs to be ported back to SymPy.


Example:
    >>> import sympy
    >>> f = sympy.Function("f")
    >>> from sympy.abc import x
    >>> f(x)
    f(x)
    >>> print sympy.srepr(f(x).func)
    Function('f')
    >>> f(x).args
    (x,)

"""

from core import BasicMeta, C
from basic import Basic
from singleton import S
from expr import Expr, AtomicExpr

from cache import cacheit
#from numbers import Rational, Integer
#from symbol import Symbol, Dummy
from sympy.utilities.decorator import deprecated
from sympy.utilities import all, any

from sympy import mpmath

from itertools import repeat

class PoleError(Exception):
    pass

class FunctionClass(BasicMeta):
    """
    Base class for function classes. FunctionClass is a subclass of type.

    Use Function('<function name>' [ , signature ]) to create
    undefined function classes.
    """
    __metaclass__ = BasicMeta

    _new = type.__new__

    def __repr__(cls):
        return cls.__name__

class UndefinedFunction(FunctionClass):
    """
    The (meta)class of undefined functions.
    """
    def __new__(mcl, name, signature=None):
        attrdict = {'undefined_Function': True}
        if signature is not None:
            attrdict['signature'] = signature
        bases = (Function,)
        return BasicMeta.__new__(mcl, name, bases, attrdict)

class Application(Basic):
    """
    Base class for applied functions.

    Instances of Application represent the result of applying an application of
    any type to any object.
    """
    __metaclass__ = FunctionClass
    __slots__ = []

    is_Function = True

    nargs = None

    @cacheit
    def __new__(cls, *args, **options):
        args = map(sympify, args)
        # these lines should be refactored
        for opt in ["nargs", "dummy", "comparable", "noncommutative", "commutative"]:
            if opt in options:
                del options[opt]
        # up to here.
        if not options.pop('evaluate', True):
            return super(Application, cls).__new__(cls, *args, **options)
        evaluated = cls.eval(*args)
        if evaluated is not None:
            return evaluated
        # Just undefined functions have nargs == None
        if not cls.nargs and hasattr(cls, 'undefined_Function'):
            r = super(Application, cls).__new__(cls, *args, **options)
            r.nargs = len(args)
            return r
        return super(Application, cls).__new__(cls, *args, **options)

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
                coeff, terms = arg.as_coeff_mul()
                if coeff is not S.One:
                    return cls(coeff) * cls(arg._new_rawargs(*terms))

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
                if self.nargs == new.nargs or not new.nargs:
                    return new(*self.args)
                # Written down as an elif to avoid a super-long line
                elif isinstance(new.nargs, tuple) and self.nargs in new.nargs:
                    return new(*self.args)
        return self.func(*[s.subs(old, new) for s in self.args])


class Function(Application, Expr):
    """
    Base class for applied numeric functions.
    Constructor of undefined function classes.

    """

    @cacheit
    def __new__(cls, *args, **options):
        # NOTE: this __new__ is twofold:
        #
        # 1 -- it can create another *class*, which can then be instantiated by
        #      itself e.g. Function('f') creates a new class f(Function)
        #
        # 2 -- on the other hand, we instantiate -- that is we create an
        #      *instance* of a class created earlier in 1.

        # (1) create new function class
        #     UC: Function('f')
        if cls is Function:
            return UndefinedFunction(*args)

        # (2) create new instance of a class created in (1)
        #     UC: Function('f')(x)
        #     UC: sin(x)
        return Application.__new__(cls, *args, **options)


    @property
    def is_commutative(self):
        if all(getattr(t, 'is_commutative') for t in self.args):
            return True
        else:
            return False

    @classmethod
    @deprecated
    def canonize(cls, *args):
        return cls.eval(*args)

    def _eval_evalf(self, prec):
        # Lookup mpmath function based on name
        fname = self.func.__name__
        try:
            if not hasattr(mpmath, fname):
                from sympy.utilities.lambdify import MPMATH_TRANSLATIONS
                fname = MPMATH_TRANSLATIONS[fname]
            func = getattr(mpmath, fname)
        except (AttributeError, KeyError):
            try:
                return C.Real(self._imp_(*self.args), prec)
            except (AttributeError, TypeError):
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

        return Expr._from_mpmath(v, prec)

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

    def as_base_exp(self):
        return self, S.One

    def _eval_nseries(self, x, n):
        if self.func.nargs != 1:
            raise NotImplementedError('series for user-defined and \
multi-arg functions are not supported.')
        arg = self.args[0]
        arg0 = arg.limit(x, 0)
        from sympy import oo
        if arg0 == S.NaN or arg0.is_bounded == False:
            raise PoleError("Cannot expand %s around 0" % (arg))
        if arg0:
            e = self
            e1 = e.expand()
            if e == e1:
                #for example when e = sin(x+1) or e = sin(cos(x))
                #let's try the general algorithm
                term = e.subs(x, S.Zero)
                if arg0 == S.NaN or term.is_bounded == False:
                    raise PoleError("Cannot expand %s around 0" % (self))
                series = term
                fact = S.One
                for i in range(n-1):
                    i += 1
                    fact *= Rational(i)
                    e = e.diff(x)
                    subs = e.subs(x, S.Zero)
                    if subs == S.NaN:
                        # try to evaluate a limit if we have to
                        subs = e.limit(x, S.Zero)
                        if subs.is_bounded == False:
                            raise PoleError("Cannot expand %s around 0" % (self))
                    term = subs*(x**i)/fact
                    term = term.expand()
                    series += term
                return series + C.Order(x**n, x)
            return e1.nseries(x, n=n)
        l = []
        g = None
        for i in xrange(n+2):
            g = self.taylor_term(i, arg, g)
            g = g.nseries(x, n=n)
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
        sargs, terms = self.args, []
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
        sargs, terms = self.args, []
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
        sargs, terms = self.args, []
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
        sargs, terms = self.args, []
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
        sargs, terms = self.args, []
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
        sargs, terms = self.args, []
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
        sargs, terms = self.args, []
        for term in sargs:
            if hasattr(term, '_eval_expand_trig'):
                newterm = term._eval_expand_trig(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_func(self, deep=True, **hints):
        sargs, terms = self.args, []
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
            args = self.args

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
        u = self.args[argindex - 1]
        if u.is_Symbol:
            uself = self
        else:
            u = C.Dummy('u')
            uself = self.func(u)
        return Derivative(uself, u, evaluate=False)

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

class WildFunction(Function, AtomicExpr):
    """
    WildFunction() matches any expression but another WildFunction()
    XXX is this as intended, does it work ?
    """

    nargs = 1

    def __new__(cls, name, **assumptions):
        obj = Function.__new__(cls, name, **assumptions)
        obj.name = name
        return obj

    def matches(self, expr, repl_dict={}, evaluate=False):
        if self in repl_dict:
            if repl_dict[self] == expr:
                return repl_dict
            else:
                return None
        if self.nargs is not None:
            if not hasattr(expr,'nargs') or self.nargs != expr.nargs:
                return None
        repl_dict = repl_dict.copy()
        repl_dict[self] = expr
        return repl_dict

    @property
    def is_number(self):
        return False

class Derivative(Expr):
    """
    Carries out differentiation of the given expression with respect to symbols.

    expr must define ._eval_derivative(symbol) method that returns
    the differentiation result. This function only needs to consider the
    non-trivial case where expr contains symbol and it should call the diff()
    method interally (not _eval_derivative); Derivative should be the only
    one to call _eval_derivative.

    Examples:

    Derivative(Derivative(expr, x), y) -> Derivative(expr, x, y)
    Derivative(expr, x, 3)  -> Derivative(expr, x, x, x)

    """

    is_Derivative   = True

    def __new__(cls, expr, *symbols, **assumptions):
        expr = sympify(expr)
        if not symbols:
            return expr

        # standardize symbols
        symbols = list(sympify(symbols))
        if not symbols[-1].is_Integer:
            symbols.append(S.One)
        symbol_count = []
        all_zero = True
        i = 0
        while i < len(symbols) - 1: # process up to final Integer
            s, count = symbols[i: i+2]
            iwas = i
            if s.is_Symbol:
                if count.is_Symbol:
                    count = 1
                    i += 1
                elif count.is_Integer:
                    count = int(count)
                    i += 2

            if i == iwas: # didn't get an update because of bad input
                raise ValueError('Derivative expects Symbol [, Integer] args but got %s, %s' % (s, count))

            symbol_count.append((s, count))
            if all_zero and not count == 0:
                all_zero = False

        # We make a special case for 0th derivative, because there
        # is no good way to unambiguously print this.
        if all_zero:
            return expr

        evaluate = assumptions.pop('evaluate', False)

        # look for a quick exit if there are symbols that are not in the free symbols
        if evaluate:
            if set(sc[0] for sc in symbol_count
                  ).difference(expr.free_symbols):
                return S.Zero

        # We make a generator so as to only generate a symbol when necessary.
        # If a high order of derivative is requested and the expr becomes 0
        # after a few differentiations, then we won't need the other symbols
        symbolgen = (s for s, count in symbol_count for i in xrange(count))

        if expr.is_commutative:
            assumptions['commutative'] = True

        if (not (hasattr(expr, '_eval_derivative') and
                 evaluate) and
            not isinstance(expr, Derivative)):
            symbols = list(symbolgen)
            obj = Expr.__new__(cls, expr, *symbols, **assumptions)
            return obj

        # compute the derivative now
        unevaluated_symbols = []
        for s in symbolgen:
            obj = expr._eval_derivative(s)
            if obj is None:
                unevaluated_symbols.append(s)
            elif obj is S.Zero:
                return S.Zero
            else:
                expr = obj

        if not unevaluated_symbols:
            return expr

        return Expr.__new__(cls, expr, *unevaluated_symbols, **assumptions)

    def _eval_derivative(self, s):
        if s not in self.variables:
            obj = self.expr.diff(s)
            if not obj:
                return obj
            if isinstance(obj, Derivative):
                return Derivative(obj.expr, *(self.variables + obj.variables))
            return Derivative(obj, *self.variables)
        return Derivative(self.expr, *(self.variables + (s, )), **{'evaluate': False})

    def doit(self, **hints):
        expr = self.expr
        if hints.get('deep', True):
            expr = expr.doit(**hints)
        hints['evaluate'] = True
        return Derivative(expr, *self.variables, **hints)

    @property
    def expr(self):
        return self._args[0]

    @property
    def variables(self):
        return self._args[1:]

    @property
    def free_symbols(self):
        return self.expr.free_symbols

    def _eval_subs(self, old, new):
        if self==old:
            return new
        return Derivative(*map(lambda x: x._eval_subs(old, new), self.args))

    def matches(self, expr, repl_dict={}, evaluate=False):
        if self in repl_dict:
            if repl_dict[self] == expr:
                return repl_dict
        elif isinstance(expr, Derivative):
            if len(expr.variables) == len(self.variables):
                return Expr.matches(self, expr, repl_dict, evaluate)

    def _eval_lseries(self, x):
        dx = self.args[1:]
        for term in self.args[0].lseries(x):
            yield Derivative(term, *dx)

    def _eval_nseries(self, x, n):
        arg = self.args[0].nseries(x, n=n)
        o = arg.getO()
        dx = self.args[1:]
        rv = [Derivative(a, *dx) for a in Add.make_args(arg.removeO())]
        if o:
            rv.append(o/x)
        return Add(*rv)

    def _eval_as_leading_term(self, x):
        return self.args[0].as_leading_term(x)

class Lambda(Function):
    """
    Lambda(x, expr) represents a lambda function similar to Python's
    'lambda x: expr'. A function of several variables is written as
    Lambda((x, y, ...), expr).

    A simple example:
        >>> from sympy import Lambda
        >>> from sympy.abc import x
        >>> f = Lambda(x, x**2)
        >>> f(4)
        16

    For multivariate functions, use:
        >>> from sympy.abc import y, z, t
        >>> f2 = Lambda(x, y, z, t, x + y**z + t**z)
        >>> f2(1, 2, 3, 4)
        73

    Multivariate functions can be curries for partial applications:
        >>> sum2numbers = Lambda(x, y, x+y)
        >>> sum2numbers(1,2)
        3
        >>> plus1 = sum2numbers(1)
        >>> plus1(3)
        4

    A handy shortcut for lots of arguments:
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
    def eval(cls,*args):
        obj = Expr.__new__(cls, *args)
        #use dummy variables internally, just to be sure
        nargs = len(args)-1

        expression = args[nargs]
        funargs = [C.Dummy(arg.name) for arg in args[:nargs]]
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
            >>> from sympy import Lambda
            >>> from sympy.abc import x, y
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

    def __hash__(self):
        return super(Lambda, self).__hash__()

    @property
    def is_identity(self):
        """Return ``True`` if this ``Lambda`` is an identity function. """
        if len(self.args) == 2:
            return self.args[0] == self.args[1]
        else:
            return None

def diff(f, *symbols, **kwargs):
    """
    Differentiate f with respect to symbols.

    This is just a wrapper to unify .diff() and the Derivative class; its
    interface is similar to that of integrate().  You can use the same
    shortcuts for multiple variables as with Derivative.  For example,
    diff(f(x), x, x, x) and diff(f(x), x, 3) both return the third derivative
    of f(x).

    You can pass evaluate=False to get an unevaluated Derivative class.  Note
    that if there are 0 symbols (such as diff(f(x), x, 0), then the result will
    be the function (the zeroth derivative), even if evaluate=False.

    Examples:
    >>> from sympy import sin, cos, Function, diff
    >>> from sympy.abc import x, y
    >>> f = Function('f')

    >>> diff(sin(x), x)
    cos(x)
    >>> diff(f(x), x, x, x)
    D(f(x), x, x, x)
    >>> diff(f(x), x, 3)
    D(f(x), x, x, x)
    >>> diff(sin(x)*cos(y), x, 2, y, 2)
    cos(y)*sin(x)

    >>> type(diff(sin(x), x))
    cos
    >>> type(diff(sin(x), x, evaluate=False))
    <class 'sympy.core.function.Derivative'>
    >>> type(diff(sin(x), x, 0))
    sin
    >>> type(diff(sin(x), x, 0, evaluate=False))
    sin

    See Also
    http://documents.wolfram.com/v5/Built-inFunctions/AlgebraicComputation/Calculus/D.html

    """
    kwargs.setdefault('evaluate', True)
    return Derivative(f, *symbols, **kwargs)

def expand(e, deep=True, power_base=True, power_exp=True, mul=True, \
           log=True, multinomial=True, basic=True, **hints):
    """
    Expand an expression using methods given as hints.

    Hints are applied with arbitrary order so your code shouldn't
    depend on the way hints are passed to this method.

    Hints evaluated unless explicitly set to False are:
      basic, log, multinomial, mul, power_base, and power_exp
    The following hints are supported but not applied unless set to True:
      complex, func, and trig.

    basic is a generic keyword for methods that want to be expanded
    automatically.  For example, Integral uses expand_basic to expand the
    integrand.  If you want your class expand methods to run automatically and
    they don't fit one of the already automatic methods, wrap it around
    _eval_expand_basic.

    If deep is set to True, things like arguments of functions are
    recursively expanded.  Use deep=False to only expand on the top
    level.

    Also see expand_log, expand_mul, expand_complex, expand_trig,
    and expand_func, which are wrappers around those expansion methods.

    >>> from sympy import cos, exp
    >>> from sympy.abc import x, y, z

    mul - Distributes multiplication over addition.
    >>> (y*(x + z)).expand(mul=True)
    x*y + y*z

    complex - Split an expression into real and imaginary parts.
    >>> (x+y).expand(complex=True)
    I*im(x) + I*im(y) + re(x) + re(y)
    >>> cos(x).expand(complex=True)
    cos(re(x))*cosh(im(x)) - I*sin(re(x))*sinh(im(x))

    power_exp - Expand addition in exponents into multiplied bases.
    >>> exp(x+y).expand(power_exp=True)
    exp(x)*exp(y)
    >>> (2**(x+y)).expand(power_exp=True)
    2**x*2**y

    power_base - Split powers of multiplied bases.
    >>> ((x*y)**z).expand(power_base=True)
    x**z*y**z

    log - Pull out power of an argument as a coefficient and split logs products
    into sums of logs.  Note that these only work if the arguments of the log
    function have the proper assumptions: the arguments must be positive and the
    exponents must be real.
    >>> from sympy import log, symbols
    >>> log(x**2*y).expand(log=True)
    log(y*x**2)
    >>> x, y = symbols('x,y', positive=True)
    >>> log(x**2*y).expand(log=True)
    2*log(x) + log(y)

    trig - Do trigonometric expansions.
    >>> cos(x+y).expand(trig=True)
    cos(x)*cos(y) - sin(x)*sin(y)

    func - Expand other functions.
    >>> from sympy import gamma
    >>> gamma(x+1).expand(func=True)
    x*gamma(x)

    multinomial - Expand (x + y + ...)**n where n is a positive integer.
    >>> ((x+y+z)**2).expand(multinomial=True)
    2*x*y + 2*x*z + 2*y*z + x**2 + y**2 + z**2

    You can shut off methods that you don't want.
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

    Note: because hints are applied in arbitrary order, some hints may
    prevent expansion by other hints if they are applied first.  In
    particular, mul may distribute multiplications and prevent log and
    power_base from expanding them.  Also, if mul is applied before multinomial,
    the expression might not be fully distributed.  The solution is to expand
    with mul=False first, then run expand_mul if you need further expansion.

    Examples:
    >>> from sympy import expand_log, expand, expand_mul
    >>> x, y, z = symbols('x,y,z', positive=True)

    >> expand(log(x*(y+z))) # could be either one below
    log(x*y + x*z)
    log(x) + log(y + z)

    >>> expand_log(log(x*y+x*z))
    log(x*y + x*z)

    >> expand(log(x*(y+z)), mul=False)
    log(x) + log(y + z)


    >> expand((x*(y+z))**x) # could be either one below
    (x*y + x*z)**x
    x**x*(y + z)**x

    >>> expand((x*(y+z))**x, mul=False)
    x**x*(y + z)**x


    >> expand(x*(y+z)**2) # could be either one below
    2*x*y*z + x*y**2 + x*z**2
    x*(y + z)**2

    >>> expand(x*(y+z)**2, mul=False)
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

# These are simple wrappers around single hints.  Feel free to add ones for
# power_exp, power_base, multinomial, or basic if you need them.
def expand_mul(expr, deep=True):
    """
    Wrapper around expand that only uses the mul hint.  See the expand
    docstring for more information.

    Example:
    >>> from sympy import symbols, expand_mul, exp, log
    >>> x, y = symbols('x,y', positive=True)
    >>> expand_mul(exp(x+y)*(x+y)*log(x*y**2))
    x*exp(x + y)*log(x*y**2) + y*exp(x + y)*log(x*y**2)

    """
    return sympify(expr).expand(deep=deep, mul=True, power_exp=False,\
    power_base=False, basic=False, multinomial=False, log=False)

def expand_multinomial(expr, deep=True):
    """
    Wrapper around expand that only uses the multinomial hint.  See the expand
    docstring for more information.

    Example:
    >>> from sympy import symbols, expand_multinomial, exp
    >>> x, y = symbols('xy', positive=True)
    >>> expand_multinomial((x + exp(x + 1))**2)
    x**2 + 2*x*exp(1 + x) + exp(2 + 2*x)

    """
    return sympify(expr).expand(deep=deep, mul=False, power_exp=False,\
    power_base=False, basic=False, multinomial=True, log=False)


def expand_log(expr, deep=True):
    """
    Wrapper around expand that only uses the log hint.  See the expand
    docstring for more information.

    Example:
    >>> from sympy import symbols, expand_log, exp, log
    >>> x, y = symbols('x,y', positive=True)
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
    >>> from sympy import expand_func, gamma
    >>> from sympy.abc import x
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
    >>> from sympy import expand_trig, sin, cos
    >>> from sympy.abc import x, y
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
    >>> from sympy import expand_complex, I, im, re
    >>> from sympy.abc import z
    >>> expand_complex(z**(2*I))
    I*im(z**(2*I)) + re(z**(2*I))

    """
    return sympify(expr).expand(deep=deep, complex=True, basic=False,\
    log=False, mul=False, power_exp=False, power_base=False, multinomial=False)

def count_ops(expr, visual=False):
    """
    Return a representation (integer or expression) of the operations in expr.

    If `visual` is False (default) then the sum of the coefficients of the
    visual expression will be returned.

    If `visual` is True then the number of each type of operation is shown
    with the core class types (or their virtual equivalent) multiplied by the
    number of times they occur.

    If expr is an iterable, the sum of the op counts of the
    items will be returned.

    Examples:
        >>> from sympy.abc import a, b, x, y
        >>> from sympy import sin, count_ops

    Although there isn't a SUB object, minus signs are interpreted as
    either negations or subtractions:
        >>> (x - y).count_ops(visual=True)
        SUB
        >>> (-x).count_ops(visual=True)
        NEG

    Here, there are two Adds and a Pow:
        >>> (1 + a + b**2).count_ops(visual=True)
        POW + 2*ADD

    In the following, an Add, Mul, Pow and two functions:
        >>> (sin(x)*x + sin(x)**2).count_ops(visual=True)
        ADD + MUL + POW + 2*SIN

    for a total of 5:
        >>> (sin(x)*x + sin(x)**2).count_ops(visual=False)
        5

    Note that "what you type" is not always what you get. The expression
    1/x/y is translated by sympy into 1/(x*y) so it gives a DIV and MUL rather
    than two DIVs:
        >>> (1/x/y).count_ops(visual=True)
        DIV + MUL

    The visual option can be used to demonstrate the difference in
    operations for expressions in different forms. Here, the Horner
    representation is compared with the expanded form of a polynomial:
        >>> eq=x*(1 + x*(2 + x*(3 + x)))
        >>> count_ops(eq.expand(), visual=True) - count_ops(eq, visual=True)
        -MUL + 3*POW

    The count_ops function also handles iterables:
        >>> count_ops([x, sin(x), None, True, x + 2], visual=False)
        2
        >>> count_ops([x, sin(x), None, True, x + 2], visual=True)
        ADD + SIN
        >>> count_ops({x: sin(x), x + 2: y + 1}, visual=True)
        SIN + 2*ADD

    """
    from sympy.simplify.simplify import fraction

    expr = sympify(expr)
    if isinstance(expr, Expr):

        ops = []
        args = [expr]
        NEG = C.Symbol('NEG')
        DIV = C.Symbol('DIV')
        SUB = C.Symbol('SUB')
        ADD = C.Symbol('ADD')
        def isneg(a):
            c = a.as_coeff_mul()[0]
            return c.is_Number and c.is_negative
        while args:
            a = args.pop()
            if a.is_Rational:
                #-1/3 = NEG + DIV
                if a is not S.One:
                    if a.p < 0:
                        ops.append(NEG)
                    if a.q != 1:
                        ops.append(DIV)
                    continue
            elif a.is_Mul:
                if isneg(a):
                    ops.append(NEG)
                    if a.args[0] is S.NegativeOne:
                        a = a.as_two_terms()[1]
                    else:
                        a = -a
                n, d = fraction(a)
                if n.is_Integer:
                    ops.append(DIV)
                    if n < 0:
                        ops.append(NEG)
                    args.append(d)
                    continue # won't be -Mul but could be Add
                elif d is not S.One:
                    if not d.is_Integer:
                        args.append(d)
                    ops.append(DIV)
                    args.append(n)
                    continue # could be -Mul
            elif a.is_Add:
                aargs = list(a.args)
                negs = 0
                for i, ai in enumerate(aargs):
                    if isneg(ai):
                        negs += 1
                        args.append(-ai)
                        if i > 0:
                            ops.append(SUB)
                    else:
                        args.append(ai)
                        if i > 0:
                            ops.append(ADD)
                if negs == len(aargs): # -x - y = NEG + SUB
                    ops.append(NEG)
                elif isneg(aargs[0]): # -x + y = SUB, but we already recorded an ADD
                    ops.append(SUB - ADD)
                continue
            if a.is_Pow and a.exp is S.NegativeOne:
                ops.append(DIV)
                args.append(a.base) # won't be -Mul but could be Add
                continue
            if (a.is_Mul or
                a.is_Pow or
                a.is_Function or
                isinstance(a, Derivative) or
                isinstance(a, C.Integral)):

                o = C.Symbol(a.func.__name__.upper())
                # count the args
                if (a.is_Mul or
                    isinstance(a, C.LatticeOp)):
                   ops.append(o*(len(a.args) - 1))
                else:
                    ops.append(o)
            args.extend(a.args)

    elif type(expr) is dict:
        ops = [count_ops(k, visual=visual) +
               count_ops(v, visual=visual) for k, v in expr.iteritems()]
    elif hasattr(expr, '__iter__'):
        ops = [count_ops(i, visual=visual) for i in expr]
    elif not isinstance(expr, Basic):
        ops = []
    else: # it's Basic not isinstance(expr, Expr):
        assert isinstance(expr, Basic)
        ops = [count_ops(a, visual=visual) for a in expr.args]

    if not ops:
        if visual:
            return S.Zero
        return 0

    ops = Add(*ops)

    if visual:
        return ops

    if ops.is_Number:
        return int(ops)

    return sum(int((a.args or [1])[0]) for a in Add.make_args(ops))

from numbers import Rational, Integer
from sympify import sympify
from add    import Add
