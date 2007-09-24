"""
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

class Apply(Basic, ArithMeths, RelMeths):
    """ Represents unevaluated function value.

    Apply(func, arg1, arg2, ..., **assumptions) <-> func(arg1, arg2, .., **assumptions)
    """

    precedence = Basic.Apply_precedence

    @cache_it
    def __new__(cls, *args, **kwargs):
        args = map(Basic.sympify, args)
        func = args[0]
        func_args = args[1:]

        # f(g+O(h)) -> f(g) + O(f'(g)*h) if f'!=0, otherwise f(g) + O(f''(g)*h**2), etc.
        i = -1
        for a in func_args:
            i += 1
            a0,o0 = a.as_expr_orders()
            if not isinstance(o0, Basic.Zero):
                new_args = func_args[:i] + [a0,] + func_args[i+1:]
                f = func
                df = f
                dfa = Basic.Zero()
                n = 0
                while isinstance(dfa, Basic.Zero) and not (isinstance(df, Lambda) and isinstance(df.body, Basic.Zero)):
                    df = df.fdiff(i+1)
                    dfa = df(*new_args)
                    n += 1
                return f(*new_args) + dfa * o0**n

        obj = func._eval_apply(*func_args)
        if obj is None:
            assert isinstance(func, Function),`args`
            cls = getattr(Basic,'Apply'+func.__class__.__name__, cls)
            obj = Basic.__new__(cls, *args, **kwargs)
            obj._func = func
        return obj

    @property
    def func(self):
        return self._func

    def __getitem__(self, iter):
        return self._args[1:][iter]

    def __len__(self):
        return len(self[:])

    def tostr(self, level=0):
        p = self.precedence
        r = '%s(%s)' % (self.func.tostr(p), ', '.join([a.tostr() for a in self]))
        if p <= level:
            return '(%s)' % (r)
        return r

    def _eval_subs(self, old, new):
        if self == old:
            return new
        elif isinstance(old, Apply) and old[:] == self[:]:
            try:
                newfunc = Lambda(new, *old[:])
                func = self.func.subs(old.func, newfunc)

                if func != self.func:
                    return func(*self[:])
            except TypeError:
                pass
        elif isinstance(old, Function) and isinstance(new, Function):
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
            elif isinstance(self.func, Basic.Sin):
                #Just a temporary workaround for Sin
                df = Basic.cos(self[0])
                l.append(df * da)
            else:
                df = self.func.fdiff(i)
                l.append(Apply(df,*self[:]) * da)
        return Basic.Add(*l)

    def _eval_power(b, e):
        if len(b)==1:
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
        if isinstance(other.func, self.func.__class__) and len(self[:])==len(other[:]):
            for a1,a2 in zip(self,other):
                if not (a1==a2):
                    return False
            return True

    def as_base_exp(self):
        return self, Basic.One()

    def count_ops(self, symbolic=True):
        return Basic.Add(*[t.count_ops(symbolic) for t in self])

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
                print '%s(%s).oseries(%s) is unevaluated' % (self.func,arg,order)
                return
            return e1.oseries(order)
        return self._compute_oseries(arg, order, self.func.taylor_term, self.func)

    def _eval_is_polynomial(self, syms):
        for arg in self:
            if arg.has(*syms):
                return False
        return True

    def _eval_expand_complex(self, *args):
        func = self.func(*[ a._eval_expand_complex(*args) for a in self])
        return Basic.Re()(func) + S.ImaginaryUnit * Basic.Im()(func)

    def _eval_rewrite(self, pattern, rule, **hints):
        if hints.get('deep', False):
            args = [ a._eval_rewrite(pattern, rule, **hints) for a in self]
        else:
            args = self[:]

        if pattern is None or isinstance(self.func, pattern):
            if hasattr(self, rule):
                rewritten = getattr(self, rule)(*args)

                if rewritten is not None:
                    return rewritten

        return self.func(*args, **self._assumptions)

class Function(Basic, ArithMeths, NoRelMeths):
    """ Base class for function objects, represents also undefined functions.

    Undefined functions:
      f = Function('f', nofargs=1)
      f.fdiff(1) -> FApply(DF(1), f)
      f(x).diff(x) -> Apply(f,x).diff(x) -> Apply(FApply(DF(1),f),x)

    Defined functions:
      exp = Exp()
      f.fdiff(1) -> exp
      log = Log()
      log.fdiff(1) -> Lambda(1/_x, _x)

    Anonymous/lambda functions:
      f = Lambda(_x**2 + 1, _x)
      f.fdiff(1) -> Lambda((_x**2+1).diff(_x), _x) -> Lambda(2*_x, _x)
      f = Lambda(g(_x), _x) -> Lambda(Apply(g,_x),_x) -> g
      f = Lambda(g(2*_x), _x) -> Lambda(Apply(g,2*_x),_x)
      f.fdiff(1) -> Lambda(Apply(FApply(DF(1), g), 2*_x) * (2*_x).diff(_x), _x)

    Function derivatives (applying operators to functions):
      df = DF(1)(f) == FApply(DF(1), f)
      df.fdiff(1) -> FApply(DF(1), FApply(DF(1), f)) -> FApply(DF(1,1), f)

    Composition of functions:
      e2 = Composition(exp, sin)
      e2(x) -> exp(sin(x)) -> Apply(exp, Apply(sin, x))
      e2.fdiff(1) -> Lambda(exp(sin(_x)).diff(_x), _x) -> Lambda(Apply(exp, Apply(sin, _x)).diff(_x), _x)
                  -> Lambda(Apply(Compostion(exp, sin), _x).diff(_x), _x)
                  -> Lambda(Apply(FApply(DF(1),exp), Apply(sin, _x)) * Apply(FApply(DF(1), sin),_x), _x)

    """

    nofargs = None
    args = None
    body = None
    has_derivative = False

    def __new__(cls, name, nofargs = None, commutative=True, **assumptions):
        obj = Basic.singleton.get(name)
        if obj is not None:
            return obj()
        obj = Basic.__new__(cls, commutative=commutative, **assumptions)
        obj.name = name
        obj.nofargs = nofargs
        return obj

    def _hashable_content(self):
        return (self.name, self.nofargs)

    def tostr(self, level=0):
        return self.name

    def torepr(self):
        return '%s(%r, nofargs = %r)' % (self.__class__.__name__, self.name, self.nofargs)

    def _eval_apply(self, *args):
        return

    def __call__(self, *args, **assumptions):
        n = self.nofargs
        if n is not None:
            if isinstance(n, int):
                n1,n2 = n,n
            else:
                n1,n2 = n
            if not (n1<=len(args)<=n2):
                raise TypeError('%s takes %s to %s arguments (got %s)'\
                            % (self, n1,n2, len(args)))
        return Apply(self, *args, **assumptions)

    def as_lambda(self):
        body, args = self.with_dummy_arguments()
        return Lambda(body, *args)

    def with_dummy_arguments(self, args = None):
        """ Return a function value with dummy arguments.
        """
        if args is None:
            if self.nofargs is None:
                raise TypeError('unknown number of arguments for undefined function %r' % (self))
            else:
                args = tuple([Basic.Symbol('x%s'%(i+1),dummy=True) for i in range(self.nofargs)])
        else:
            if self.nofargs is None:
                self.nofargs = len(args)
        return self(*args), args

    def fdiff(self, argindex=1):
        if self.nofargs is not None:
            if isinstance(self.nofargs, tuple):
                nofargs = self.nofargs[-1]
            else:
                nofargs = self.nofargs
            if not (1<=argindex<=nofargs):
                raise TypeError("argument index %r is out of range [1,%s]" % (i,nofargs))
        return FDerivative(argindex)(self)

    def diff(self, *args, **assumptions):
        raise TypeError('%s.diff not usable for functions, use .fdiff() method instead.' % (self.__class__.__name__))

    def _eval_power(b, e):
        if b.nofargs is None:
            return
        body, args = b.with_dummy_arguments()
        return Lambda(body**e, *args)

    def inverse(self):
        return FPow(self, -1)

    def _calc_splitter(self, d):
        return self

    def count_ops(self, symbolic):
        if symbolic:
            return Basic.Symbol(self.__class__.__name__.upper())

        if isinstance(self.nofargs, (tuple, list)):
            return Basic.Integer(self.nofargs[0])
        else:
            return Basic.Integer(self.nofargs)

    def taylor_term(self, n, x, *previous_terms):
        raise NotImplementedError('%s.taylor_term(..)' % (self))

class WildFunction(Function, Atom):

    def matches(pattern, expr, repl_dict={}, evaluate=False):
        for p,v in repl_dict.items():
            if p==pattern:
                if v==expr: return repl_dict
                return None
        if pattern.nofargs is not None:
            if isinstance(expr, Apply):
                if pattern.nofargs != expr.func.nofargs:
                    return None
            elif pattern.nofargs != expr.nofargs:
                return None
        repl_dict = repl_dict.copy()
        repl_dict[pattern] = expr
        return repl_dict

    def tostr(self, level=0):
        return self.name + '_'

class FApply(Function):
    """
    Defines n-ary operator that acts on symbolic functions.

    DF(1)(f) -> FApply(DF(1), f)
    DF(2)(FApply(DF(1), f)) -> FApply(DF(2), FApply(DF(1),f)) -> FApply(DF(1,2), f)
    """
    def __new__(cls, operator, *funcs, **assumptions):
        operator = Basic.sympify(operator)
        funcs = map(Basic.sympify, funcs)
        obj = operator._eval_fapply(*funcs, **assumptions)
        if obj is not None:
            return obj
        return Basic.__new__(cls, operator, *funcs, **assumptions)

    def _hashable_content(self):
        return self._args

    @property
    def operator(self):
        return self._args[0]

    @property
    def func(self):
        return self._args[0]

    def __getitem__(self, iter):
        return self._args[1:][iter]

    @property
    def funcs(self):
        return self._args[1:]

    def tostr(self, level=0):
        p = self.precedence
        r = '%s(%s)' % (self.operator.tostr(p), ', '.join([a.tostr() for a in self.funcs]))
        if p <= level:
            return '(%s)' % (r)
        return r

    _eval_subs = Basic._seq_subs


class Lambda(Function):
    """
    Lambda(expr, arg1, arg2, ...) -> lambda arg1, arg2,... : expr

    Lambda instance has the same assumptions as its body.

    Lambda(Apply(g,x),x) -> g
    """
    precedence = Basic.Lambda_precedence
    name = None
    has_derivative = True

    def __new__(cls, expr, *args):
        expr = Basic.sympify(expr)
        args = tuple(map(Basic.sympify, args))
        if isinstance(expr, Apply):
            if expr[:]==args:
                return expr.func
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


class FPow(Function):
    precedence = Basic.Apply_precedence

    def __new__(cls, a, b, **assumptions):
        a = Basic.sympify(a)
        b = Basic.sympify(b)
        if isinstance(b, Basic.Zero):
            return Basic.One()
        if isinstance(b, Basic.One):
            return a
        obj = a._eval_fpower(b)
        if obj is None:
            obj = Basic.__new__(cls, a, b, **assumptions)
        return obj

    def _hashable_content(self):
        return self._args

    @property
    def base(self):
        return self._args[0]

    @property
    def exp(self):
        return self._args[1]

    def _eval_fpower(self, other):
        if isinstance(other, Basic.Number):
            if isinstance(self.exp, Basic.Number):
                # (a ** 2) ** 3 -> a ** (2 * 3)
                return FPow(self.base, self.exp * other)
        return

    def _eval_apply(self, *args):
        if isinstance(self.exp, Basic.Integer) and self.exp.is_positive and self.exp.p < 10:
            r = self.base(*args)
            for i in xrange(self.exp.p-1):
                r = self.base(r)
            return r

    def tostr(self, level=0):
        r = '%s(%s, %s)' % (self.__class__.__name__, self.base.tostr(), self.exp.tostr())
        if self.precedence<=level:
            return '(%s)' % (r)
        return r

    def _eval_subs(self, old, new):
        if self==old: return new
        #elif exp(self.exp * log(self.base)) == old: return new
        return self.base.subs(old, new) ** self.exp.subs(old, new)

    def as_base_exp(self):
        return self.base, self.exp


class Composition(AssocOp, Function):
    """ Composition of functions.

    Composition(f1,f2,f3)(x) -> f1(f2(f3(x)))
    >>> from sympy import exp,log
    >>> l1 = Lambda('x**2','x')
    >>> l2 = Lambda('x+1','x')
    >>> Composition(l1,l2)('y')
    (1 + y)**2
    >>> Composition(l2,l1)('y')
    1 + y**2
    >>> Composition(exp,log,exp,exp)('x')
    exp(exp(x))
    """

    @classmethod
    def flatten(cls, seq):
        c_part = []
        nc_part = []
        while seq:
            o = seq.pop(0)
            if isinstance(o, Basic.One):
                continue
            if o.__class__ is cls:
                # associativity
                seq = list(o._args) + seq
                continue
            if not nc_part:
                nc_part.append(o)
                continue
            # try to combine last terms: a**b * a ** c -> a ** (b+c)
            o1 = nc_part.pop()
            b1,e1 = o1.as_base_exp()
            b2,e2 = o.as_base_exp()
            if b1==b2:
                seq.insert(0, Basic.FPow(b1,e1 + e2))
            elif o1.is_homogeneous and isinstance(o, Basic.Number):
                seq.insert(0,o1)
                seq.insert(0,o)
            elif isinstance(o1, Basic.Number) and isinstance(o, Basic.Number):
                nc_part.append(o1 * o)
            elif isinstance(o1, Basic.Lambda) and isinstance(o, Basic.Lambda):
                nc_part.append(Lambda(o1(o.body),*o[:]))
            else:
                nc_part.append(o1)
                nc_part.append(o)
        return [], nc_part, None, None

    def tostr(self, level=0):
        return '%s(%s)' % (self.__class__.__name__,', '.join([f.tostr() for f in self]))

    def _eval_apply(self, *args):
        l = list(self)
        r = l.pop()(*args)
        while l:
            r = l.pop()(r)
        return r

    def _hashable_content(self):
        return self._args

    def _matches_simple(pattern, expr, repl_dict):
        return


class FDerivative(Function):

    is_homogeneous = True

    def __new__(cls, *indices, **assumptions):
        indices = map(Basic.sympify, indices)
        indices.sort(Basic.compare)
        return Basic.__new__(cls, *indices, **assumptions)

    def _hashable_content(self):
        return self._args

    @property
    def func(self):
        return self._args[0]

    def __getitem__(self, iter):
        return self._args[1:][iter]

    @property
    def indices(self):
        return self._args

    def tostr(self, level=0):
        return 'FD%s' % (repr(tuple(self.indices)))

    def _eval_fapply(self, *funcs):
        if len(funcs)==1:
            func = funcs[0]
            if not self.indices:
                return func
            if isinstance(func, FApply) and isinstance(func.operator, FDerivative):
                # FApply(FD(1),FApply(FD(2),f)) -> FApply(FD(1,2), f)
                return FApply(FDerivative(*(self.indices+func.operator.indices)), *func.funcs)
            if isinstance(func, Lambda):
                # FApply(FD(1),Lambda _x: f(_x)) -> FApply(Lambda _x: f(_x).diff(_x))
                expr = func.body
                unevaluated_indices = []
                for i in self.indices:
                    if isinstance(i, Basic.Integer):
                        a = func[int(i)-1]
                        obj = expr.diff(a)
                        if isinstance(obj, Derivative):
                            unevaluated_indices.append(i)
                        else:
                            expr = obj
                if len(self.indices)==len(unevaluated_indices):
                    return
                if not unevaluated_indices:
                    return Lambda(expr, *func[:])
                return FApply(FDerivative(*unevaluated_indices), Lambda(expr, *func))
            if isinstance(func, DefinedFunction):
                for i in range(len(self.indices)):
                    di = self.indices[i]
                    if isinstance(di, Basic.Integer):
                        dfunc = func.fdiff(int(di))
                        new_indices = self.indices[:i] + self.indices[i+1:]
                        return FDerivative(*new_indices)(dfunc)
        return

    def _eval_power(b, e):
        if isinstance(e, Basic.Integer) and e.is_positive:
            return FDerivative(*(int(e)*b.indices))

    _eval_subs = Basic._seq_subs

    def __call__(self, func):
        return FApply(self, func)


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

    def as_apply(self):
        # Derivative(f(x),x) -> Apply(Lambda(f(_x),_x), x)
        symbols = []
        indices = []
        for s in self.symbols:
            if s not in symbols:
                symbols.append(s)
                indices.append(len(symbols))
            else:
                indices.append(symbols.index(s)+1)
        return Apply(FApply(FDerivative(*indices), Lambda(self.expr, *symbols)), *symbols)

    def _eval_derivative(self, s):
        if s not in self.symbols:
            obj = self.expr.diff(s)
            if isinstance(obj, Derivative):
                return Derivative(obj.expr, *(obj.symbols+self.symbols))
            return Derivative(obj, *self.symbols)
        return Derivative(self.expr, *((s,)+self.symbols))

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
        # XXX Is this required?
        #for a in list(old.atoms(Basic.Symbol)):
        #    if a in self.symbols:
        #        return self.as_apply().subs(old, new)
        return self.__class__(*[s.subs(old, new) for s in self], **{'evaluate': False})


class DefinedFunction(Function, Singleton, Atom):
    """ Base class for defined functions.
    """

    is_commutative = True # the values of functions are commutative

    def __new__(cls, **assumptions):
        obj = Basic.__new__(cls,**assumptions)
        obj.name = obj.__class__.__name__.lower()
        return obj

    def torepr(self):
        return '%s()' % (self.__class__.__name__)

Basic.singleton['D'] = lambda : Derivative
Basic.singleton['FD'] = lambda : FDerivative
#Basic.singleton['Int'] = lambda : Integral

def diff(f, x, times = 1, evaluate=True):
    """Derivate f with respect to x

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

class Function2(Basic, RelMeths):
    """
    Base class for applied functions.
    Constructor of undefined classes.

    We are moving to this more simple scheme. When all functions are moved, we
    simply delete Function and rename Function2 -> Function
    """

    __metaclass__ = FunctionClass

    precedence = Basic.Apply_precedence

    def __new__(cls, *args, **options):
        args = map(Basic.sympify, args)
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
        
    @classmethod
    def canonize(cls, args, **options):
        return

    #@cache_it
    def x__new__(cls, *args, **kwargs):
        args = map(Basic.sympify, args)
        func = args[0]
        func_args = args[1:]

        # f(g+O(h)) -> f(g) + O(f'(g)*h) if f'!=0, otherwise f(g) + O(f''(g)*h**2), etc.
        i = -1
        for a in func_args:
            i += 1
            a0,o0 = a.as_expr_orders()
            if not isinstance(o0, Basic.Zero):
                new_args = func_args[:i] + [a0,] + func_args[i+1:]
                f = func
                df = f
                dfa = Basic.Zero()
                n = 0
                while isinstance(dfa, Basic.Zero) and not (isinstance(df, Lambda) and isinstance(df.body, Basic.Zero)):
                    df = df.fdiff(i+1)
                    dfa = df(*new_args)
                    n += 1
                return f(*new_args) + dfa * o0**n

        obj = func._eval_apply(*func_args)
        if obj is None:
            assert isinstance(func, Function),`args`
            cls = getattr(Basic,'Apply'+func.__class__.__name__, cls)
            obj = Basic.__new__(cls, *args, **kwargs)
        return obj

    @property
    def func(self):
        return self.__class__

    def _eval_subs(self, old, new):
        if self == old:
            return new
        elif isinstance(old, Apply) and old[:] == self[:]:
            try:
                newfunc = Lambda(new, *old[:])
                func = self.func.subs(old.func, newfunc)

                if func != self.func:
                    return func(*self[:])
            except TypeError:
                pass
        elif isinstance(old, Function) and isinstance(new, Function):
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
            else:
                df = self.func.fdiff(i)
                l.append(Apply(df,*self[:]) * da)
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
        return Basic.Add(*[t.count_ops(symbolic) for t in self])

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
                print '%s(%s).oseries(%s) is unevaluated' % (self.func,arg,order)
                return
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
        return FDerivative(argindex)(self)

class SingleValuedFunction(ArithMeths, Function2):
    """
    Single-valued functions.
    """

    def series(self, x, n):
        s = Basic.Rational(0)
        for i in range(n+1):
            s += diff(self, x, i).subs(x, 0) * x**i / Basic.Factorial()(i)
        return s
