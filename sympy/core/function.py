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
    >>> f(x).func
    <class 'sympy.core.function.f'>
    >>> f(x).args
    (x,)
"""

from basic import Basic, Singleton, Atom, S, C
from basic_methods import BasicType, MetaBasicMeths
from methods import ArithMeths, NoRelMeths, RelMeths
from operations import AssocOp
from cache import cache_it

class FunctionClass(MetaBasicMeths):
    """
    Base class for function classes. FunctionClass is a subclass of type.

    Use Function('<function name>' [ , signature ]) to create
    undefined function classes.
    """

    _new = type.__new__

    def __new__(cls, arg1, arg2, arg3=None, **options):
        assert not options,`options`
        if isinstance(arg1, type):
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

    def torepr(cls):
        return cls.__name__

class Function(Basic, ArithMeths, RelMeths):
    """
    Base class for applied functions.
    Constructor of undefined classes.

    """

    __metaclass__ = FunctionClass

    precedence = Basic.Apply_precedence

    nargs = None

    @cache_it
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
                raise Exception("You need to specify exactly one string")

        # (2) create new instance of a class created in (1)
        #     UC: Function('f')(x)
        #     UC: sin(x)
        args = map(Basic.sympify, args)
        # these lines should be refactored
        for opt in ["nargs", "dummy", "comparable", "noncommutative", "commutative"]:
            if opt in options:
                del options[opt]
        # up to here.
        r = cls.canonize(*args, **options)
        if isinstance(r, Basic):
            return r
        elif r is None:
            pass
        elif not isinstance(r, tuple):
            args = (r,)
        return Basic.__new__(cls, *args, **options)

    @property
    def is_comparable(self):
        return True

    @property
    def is_commutative(self):
        return True

    @classmethod
    def canonize(cls, *args, **options):
        """
        Returns a canonical form of cls applied to arguments args.

        The canonize() method is called when the class cls is about to be
        instantiated and it should return either some simplified instance
        (possible of some other class), or if the class cls should be
        unmodified, return None.

        Example of canonize() for the function "sign"
        ---------------------------------------------

        @classmethod
        def canonize(cls, arg):
            if arg is S.NaN:
                return S.NaN
            if arg is S.Zero: return S.One
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
        elif isinstance(old, FunctionClass) and isinstance(new, FunctionClass):
            if old == self.func and old.nargs == new.nargs:
                return new(*self.args[:])
        obj = self.func._eval_apply_subs(*(self.args[:] + (old,) + (new,)))
        if obj is not None:
            return obj
        return Basic._seq_subs(self, old, new)

    def _eval_expand_basic(self, *args):
        return self

    def _eval_evalf(self):
        obj = self.func._eval_apply_evalf(*self.args[:])
        if obj is None:
            return self
        return obj

    def _eval_is_comparable(self):
        if isinstance(self.func, DefinedFunction):
            r = True
            for s in self:
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
        return C.Add(*l)

    def _eval_is_commutative(self):
        r = True
        for a in self._args:
            c = a.is_commutative
            if c is None: return None
            if not c: r = False
        return r

    def _calc_positive(self):
        return self.func._calc_apply_positive(*self[:])

    def _calc_real(self):
        return self.func._calc_apply_real(*self[:])

    def _calc_unbounded(self):
        return self.func._calc_apply_unbounded(*self[:])

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
        return 1   + C.Add(*[t.count_ops(symbolic) for t in self.args])

    def _eval_oseries(self, order):
        assert self.func.nargs==1,`self.func`
        arg = self.args[0]
        x = order.symbols[0]
        if not C.Order(1,x).contains(arg):
            return self.func(arg)
        arg0 = arg.limit(x, 0)
        if arg0 is not S.Zero:
            e = self.func(arg)
            e1 = e.expand()
            if e==e1:
                #for example when e = sin(x+1) or e = sin(cos(x))
                #let's try the general algorithm
                term = e.subs(x, S.Zero)
                series = S.Zero
                fact = S.One
                i = 0
                while not order.contains(term) or term == 0:
                    series += term
                    i += 1
                    fact *= C.Rational(i)
                    e = e.diff(x)
                    term = e.subs(x, S.Zero)*(x**i)/fact
                return series

                #print '%s(%s).oseries(%s) is unevaluated' % (self.func,arg,order)
            return e1.oseries(order)
        return self._compute_oseries(arg, order, self.func.taylor_term, self.func)

    def _eval_is_polynomial(self, syms):
        for arg in self.args:
            if arg.has(*syms):
                return False
        return True

    def _eval_expand_complex(self, *args):
        func = self.func(*[ a._eval_expand_complex(*args) for a in self.args ])
        return C.re(func) + S.ImaginaryUnit * C.im(func)

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
                raise TypeError("argument index %r is out of range [1,%s]" % (i,nargs))
        return Derivative(self,self.args[argindex-1],evaluate=False)

    def torepr(self):
        r = '%s(%r)' % (self.func.__base__.__name__, self.func.__name__)
        r+= '(%s)' % ', '.join([a.torepr() for a in self.args])
        return r

    def tostr(self, level=0):
        p = self.precedence
        r = '%s(%s)' % (self.func.__name__, ', '.join([a.tostr() for a in
            self.args]))
        if p <= level:
            return '(%s)' % (r)
        return r

    @classmethod
    def _eval_apply_evalf(cls, arg):
        arg = arg.evalf()

        #if cls.nargs == 1:
        # common case for functions with 1 argument
        #if isinstance(arg, C.Number):
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
        x = Basic.sympify(x)
        return cls(x).diff(x, n).subs(x, 0) * x**n / C.Factorial(n)

class WildFunction(Function, Atom):
    """
    WildFunction() matches any expression but another WildFunction()
    XXX is this as intended, does it work ?
    """

    nargs = 1

    def __new__(cls, name=None, **assumptions):
        if name is None:
            name = 'Wf%s' % (C.Symbol.dummycount + 1) # XXX refactor dummy counting
            C.Symbol.dummycount += 1
        obj = Function.__new__(cls, name, **assumptions)
        obj.name = name
        return obj

    def matches(pattern, expr, repl_dict={}, evaluate=False):
        for p,v in repl_dict.items():
            if p==pattern:
                if v==expr: return repl_dict
                return None
        if pattern.nargs is not None:
            if pattern.nargs != expr.nargs:
                return None
        repl_dict = repl_dict.copy()
        repl_dict[pattern] = expr
        return repl_dict

    def torepr(self):
        return '%s(%r)' % (self.__class__.__name__, self.name)

    def tostr(self, level=0):
        return self.name + '_'

    @classmethod
    def _eval_apply_evalf(cls, arg):
        return

class Derivative(Basic, ArithMeths, RelMeths):
    """
    Carries out differentation of the given expression with respect to symbols.

    expr must define ._eval_derivative(symbol) method that returns differentation result or None.

    Derivative(Derivative(expr, x), y) -> Derivative(expr, x, y)
    """

    precedence = Basic.Apply_precedence

    def __new__(cls, expr, *symbols, **assumptions):
        expr = Basic.sympify(expr)
        if not symbols: return expr
        symbols = map(Basic.sympify, symbols)

        if not assumptions.get("evaluate", False):
            obj = Basic.__new__(cls, expr, *symbols)
            return obj

        for s in symbols:
            assert isinstance(s, C.Symbol),`s`
            if not expr.has(s):
                return S.Zero

        unevaluated_symbols = []
        for s in symbols:
            obj = expr._eval_derivative(s)
            if obj is None:
                unevaluated_symbols.append(s)
            else:
                expr = obj

        if not unevaluated_symbols:
            return expr
        return Basic.__new__(cls, expr, *unevaluated_symbols)

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

    def tostr(self, level=0):
        r = 'D' + `tuple(self.args)`
        if self.precedence <= level:
            r = '(%s)' % (r)
        return r

    def _eval_subs(self, old, new):
        return Derivative(self.args[0].subs(old, new), *self.args[1:], **{'evaluate': True})

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
    """

    #XXX currently only one argument Lambda is supported

    nargs = 2

    @classmethod
    def canonize(cls, x, expr):
        obj = Basic.__new__(cls, x, expr)
        #use dummy variables internally, just to be sure
        tmp = C.Symbol("x", dummy=True)
        obj._args = (tmp, expr.subs(x, tmp))
        return obj

    def apply(self, x):
        """Applies the Lambda function "self" to the "x"

        Example:
            >>> from sympy import Symbol
            >>> x = Symbol('x')
            >>> f = Lambda(x, x**2)
            >>> f.apply(4)
            16
        
        """
        return self.args[1].subs(self.args[0], x)

    def __call__(self, *args):
        return self.apply(*args)

    def __eq__(self, other):
        if isinstance(other, Lambda):
            if self.args[1] == other.args[1].subs(other.args[0], self.args[0]):
                return True
        return False
        


def diff(f, x, times = 1, evaluate=True):
    """Differentiate f with respect to x

    It's just a wrapper to unify .diff() and the Derivative class,
    it's interface is similar to that of integrate()

    see http://documents.wolfram.com/v5/Built-inFunctions/AlgebraicComputation/Calculus/D.html
    """
    f = Basic.sympify(f)
    if evaluate == True:
        for i in range(0,times):
            f = f.diff(x)
        return f
    else:
        return Derivative(f, x, evaluate=evaluate)
