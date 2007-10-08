"""
NOTE: this doctest is for now obsoleted interface.

There are different types of functions:
1) defined function like exp or sin that has a name and body
   (in the sense that function can be evaluated).
    e = exp
2) undefined function with a name but no body. Undefined
  functions can be defined using Symbol class as follows:
    f = Symbol('f', function=True)
  (the result will be Function instance)
  or
    f = Function('f')
3) anonymous function or lambda function that has no name
   but has body with dummy variables. An anonymous function
   object creation examples:
    f = Lambda(x, exp(x)*x)
    f = Lambda(exp(x)*x)  # free symbols in the expression define the number of arguments
    f = exp * Lambda(x,x)
4) composition of functions, inverse functions

One can perform certain operations with functions, the
result will be a lambda function. Allowed operations are
addition and multiplication. Function multiplication is
elementise ie (f*g)(x) is equivalent to f(x) * g(x).
Multiplication by a number is equivalent to multiplication
by a constant function.
Composition of functions is achived via Composition class@
Eg

  f+Composition(2,exp,sin) -> lambda _x: f(x)+2*exp(sin(x))

In the above functions did not have arguments, then
it is said that functions are in unevaluated form.
When calling a function with arguments, then we get
an instance of the function value. Eg

  exp(1) -> 1
  (2*f)(x)  -> 2 * f(x)
  Lambda(x, exp(x)*x)(y) -> exp(y)*y

One can construct undefined function values from Symbol
object:
  f = Symbol('f')
  fx = f(x)
  fx.func -> Function('f')
  fx[:] -> (Symbol('x'),)
As seen above, function values are Apply instances and
have attributes .func and [:].
"""

from basic import Basic, Singleton, Atom, cache_it, S
from basic_methods import BasicType, MetaBasicMeths
from methods import ArithMeths, NoRelMeths, RelMeths
from operations import AssocOp

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

class Function(Basic, RelMeths):
    """
    Base class for applied functions.
    Constructor of undefined classes.

    """

    __metaclass__ = FunctionClass

    precedence = Basic.Apply_precedence

    nofargs = None

    @cache_it
    def __new__(cls, *args, **options):
        if cls is SingleValuedFunction or cls is Function:
            #when user writes SingleValuedFunction("f"), do an equivalent of:
            #taking the whole class SingleValuedFunction(...):
            #and rename the SingleValuedFunction to "f" and return f, thus:
            #In [13]: isinstance(f, SingleValuedFunction)
            #Out[13]: False
            #In [14]: isinstance(f, FunctionClass)
            #Out[14]: True

            if len(args) == 1 and isinstance(args[0], str):
                #always create SingleValuedFunction
                return FunctionClass(SingleValuedFunction, *args)
                return FunctionClass(SingleValuedFunction, *args, **options)
            else:
                print args
                print type(args[0])
                raise Exception("You need to specify exactly one string")
        args = map(Basic.sympify, args)
        # these lines should be refactored
        if "nofargs" in options:
            del options["nofargs"]
        if "dummy" in options:
            del options["dummy"]
        if "comparable" in options:
            del options["comparable"]
        if "noncommutative" in options:
            del options["noncommutative"]
        if "commutative" in options:
            del options["commutative"]
        # up to here.
        r = cls._eval_apply(*args, **options)
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
    def canonize(cls, args, **options):
        return

    @classmethod
    def _eval_apply(self, *args):
        return

    @property
    def func(self):
        return self.__class__

    def _eval_subs(self, old, new):
        if self == old:
            return new
        elif isinstance(old, FunctionClass) and isinstance(new, FunctionClass):
            if old == self.func and old.nofargs == new.nofargs:
                return new(*self[:])
        obj = self.func._eval_apply_subs(*(self[:] + (old,) + (new,)))
        if obj is not None:
            return obj
        return Basic._seq_subs(self, old, new)

    def _eval_expand_basic(self, *args):
        return self

    def _eval_evalf(self):
        obj = self.func._eval_apply_evalf(*self[:])
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
        # Apply(f(x), x).diff(s) -> x.diff(s) * f.fdiff(1)(s)
        i = 0
        l = []
        r = Basic.Zero()
        for a in self:
            i += 1
            da = a.diff(s)
            if isinstance(da, Basic.Zero):
                continue
            if isinstance(self.func, FunctionClass):
                df = self.fdiff(i)
                l.append(df * da)
            #else:
            #    df = self.func.fdiff(i)
            #    l.append(Apply(df,*self[:]) * da)
        return Basic.Add(*l)

    def _eval_power(b, e):
        if len(b[:])==1:
            return b.func._eval_apply_power(b[0], e)
        return

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
        return self, Basic.One()

    def count_ops(self, symbolic=True):
        #      f()             args
        return 1   + Basic.Add(*[t.count_ops(symbolic) for t in self])

    def _eval_oseries(self, order):
        assert self.func.nofargs==1,`self.func`
        arg = self[0]
        x = order.symbols[0]
        if not Basic.Order(1,x).contains(arg):
            return self.func(arg)
        arg0 = arg.limit(x, 0)
        if not isinstance(arg0, Basic.Zero):
            e = self.func(arg)
            e1 = e.expand()
            if e==e1:
                #one example is e = sin(x+1)
                #let's try the general algorithm
                term = e.subs(x, Basic.Rational(0))
                series = Basic.Rational(0)
                fact = Basic.Rational(1)
                i = 0
                while not order.contains(term):
                    series += term
                    i += 1
                    fact *= Basic.Rational(i)
                    e = e.diff(x)
                    term = e.subs(x, Basic.Rational(0))*(x**i)/fact
                return series

                #print '%s(%s).oseries(%s) is unevaluated' % (self.func,arg,order)
            return e1.oseries(order)
        return self._compute_oseries(arg, order, self.func.taylor_term, self.func)

    def _eval_is_polynomial(self, syms):
        for arg in self:
            if arg.has(*syms):
                return False
        return True

    def _eval_expand_complex(self, *args):
        func = self.func(*[ a._eval_expand_complex(*args) for a in self ])
        return Basic.Re()(func) + S.ImaginaryUnit * Basic.Im()(func)

    def _eval_rewrite(self, pattern, rule, **hints):
        if hints.get('deep', False):
            args = [ a._eval_rewrite(pattern, rule, **hints) for a in self ]
        else:
            args = self[:]

        if pattern is None or isinstance(self.func, pattern):
            if hasattr(self, rule):
                rewritten = getattr(self, rule)(*args)

                if rewritten is not None:
                    return rewritten

        return self.func(*args, **self._assumptions)

    def fdiff(self, argindex=1):
        if self.nofargs is not None:
            if isinstance(self.nofargs, tuple):
                nofargs = self.nofargs[-1]
            else:
                nofargs = self.nofargs
            if not (1<=argindex<=nofargs):
                raise TypeError("argument index %r is out of range [1,%s]" % (i,nofargs))
        return Derivative(self,self[argindex-1],evaluate=False)

    def tostr(self, level=0):
        p = self.precedence
        r = '%s(%s)' % (self.func.__name__, ', '.join([a.tostr() for a in self]))
        if p <= level:
            return '(%s)' % (r)
        return r

    @classmethod
    def _eval_apply_evalf(cls, arg):
        return

class WildFunction(Function, Atom):

    nofargs = 1

    def matches(pattern, expr, repl_dict={}, evaluate=False):
        for p,v in repl_dict.items():
            if p==pattern:
                if v==expr: return repl_dict
                return None
        if pattern.nofargs is not None:
            if pattern.nofargs != expr.nofargs:
                return None
        repl_dict = repl_dict.copy()
        repl_dict[pattern] = expr
        return repl_dict

    def tostr(self, level=0):
        return self.name + '_'

    @classmethod
    def _eval_apply_evalf(cls, arg):
        return

class Lambda(Function):
    """
    Lambda(expr, arg1, arg2, ...) -> lambda arg1, arg2,... : expr

    Lambda instance has the same assumptions as its body.

    """
    precedence = Basic.Lambda_precedence
    name = None
    has_derivative = True

    def __new__(cls, expr, *args):
        expr = Basic.sympify(expr)
        args = tuple(map(Basic.sympify, args))
        #if isinstance(expr, Apply):
        #    if expr[:]==args:
        #        return expr.func
        dummy_args = []
        for a in args:
            if not isinstance(a, Basic.Symbol):
                raise TypeError("%s %s-th argument must be Symbol instance (got %r)" \
                                % (cls.__name__, len(dummy_args)+1,a))
            d = a.as_dummy()
            expr = expr.subs(a, d)
            dummy_args.append(d)
        obj = Basic.__new__(cls, expr, *dummy_args, **expr._assumptions)
        return obj

    def _hashable_content(self):
        return self._args

    @property
    def nofargs(self):
        return len(self._args)-1

    def __getitem__(self, iter):
        return self._args[1:][iter]

    def __len__(self):
        return len(self[:])

    @property
    def body(self):
        return self._args[0]

    def tostr(self, level=0):
        precedence = self.precedence
        r = 'lambda %s: %s' % (', '.join([a.tostr() for a in self]),
                               self.body.tostr(precedence))
        if precedence <= level:
            return '(%s)' % r
        return r

    def torepr(self):
        return '%s(%s)' % (self.__class__.__name__, ', '.join([a.torepr() for a in self]))

    def as_coeff_terms(self, x=None):
        c,t = self.body.as_coeff_terms(x)
        return c, [Lambda(Basic.Mul(*t),*self[:])]

    def _eval_power(b, e):
        """
        (lambda x:f(x))**e -> (lambda x:f(x)**e)
        """
        return Lambda(b.body**e, *b[:])

    def _eval_fpower(b, e):
        """
        FPow(lambda x:f(x), 2) -> lambda x:f(f(x)))
        """
        if isinstance(e, Basic.Integer) and e.is_positive and e.p < 10 and len(b)==1:
            r = b.body
            for i in xrange(e.p-1):
                r = b(r)
            return Lambda(r, *b[:])

    def with_dummy_arguments(self, args = None):
        if args is None:
            args = tuple([a.as_dummy() for a in self])
        if len(args) != len(self):
            raise TypeError("different number of arguments in Lambda functions: %s, %s" % (len(args), len(self)))
        expr = self.body
        for a,na in zip(self, args):
            expr = expr.subs(a, na)
        return expr, args

    def _eval_expand_basic(self, *args):
        return Lambda(self.body._eval_expand_basic(*args), *self[:])

    def diff(self, *symbols):
        return Lambda(self.body.diff(*symbols), *self[:])

    def fdiff(self, argindex=1):
        if not (1<=argindex<=len(self)):
            raise TypeError("%s.fderivative() argindex %r not in the range [1,%s]"\
                            % (self.__class__, argindex, len(self)))
        s = self[argindex-1]
        expr = self.body.diff(s)
        return Lambda(expr, *self[:])

    _eval_subs = Basic._seq_subs

    def _eval_apply(self, *args):
        n = self.nofargs
        if n!=len(args):
            raise TypeError('%s takes exactly %s arguments (got %s)'\
                            % (self, n, len(args)))
        expr = self.body
        for da,a in zip(self, args):
            expr = expr.subs(da,a)
        return expr

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

        if not assumptions.get("evaluate", True):
            obj = Basic.__new__(cls, expr, *symbols)
            return obj

        for s in symbols:
            assert isinstance(s, Basic.Symbol),`s`
            if not expr.has(s):
                return Basic.Zero()

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

    def xas_apply(self):
        # Derivative(f(x),x) -> Apply(Lambda(f(_x),_x), x)
        symbols = []
        indices = []
        for s in self.symbols:
            if s not in symbols:
                symbols.append(s)
                indices.append(len(symbols))
            else:
                indices.append(symbols.index(s)+1)
        stop
        return Apply(FApply(FDerivative(*indices), Lambda(self.expr, *symbols)), *symbols)

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
        return Derivative(self.expr, *self.symbols)

    @property
    def expr(self):
        return self._args[0]

    @property
    def symbols(self):
        return self._args[1:]

    def tostr(self, level=0):
        r = 'D' + `tuple(self)`
        if self.precedence <= level:
            r = '(%s)' % (r)
        return r

    def _eval_subs(self, old, new):
        return Derivative(self[0].subs(old, new), *self[1:])

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
        if pattern.nofargs is not None:
            if pattern.nofargs != expr.nofargs:
                return None
        repl_dict = repl_dict.copy()
        repl_dict[pattern] = expr
        return repl_dict


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



# TODO rename me to something more appropriate? e.g. ArithFunction (or just
# Function?)
class SingleValuedFunction(ArithMeths, Function):
    """
    Single-valued functions.
    """

    @classmethod
    def _eval_apply_evalf(cls, arg):
        arg = arg.evalf()

        #if cls.nofargs == 1:
        # common case for functions with 1 argument
        #if isinstance(arg, Basic.Number):
        if arg.is_number:
            func_evalf = getattr(arg, cls.__name__)
            return func_evalf()
