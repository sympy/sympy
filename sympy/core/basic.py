"""Base class for all objects in sympy"""

type_class = type

import decimal
from basic_methods import BasicMeths, cache_it, cache_it_immutable, BasicType

class MemoizerArg:
    """ See Memoizer.
    """

    def __init__(self, allowed_types, converter = None, name = None):
        self._allowed_types = allowed_types
        self.converter = converter
        self.name = name

    def fix_allowed_types(self, have_been_here={}):
        i = id(self)
        if have_been_here.get(i): return
        allowed_types = self._allowed_types
        if isinstance(allowed_types, str):
            self.allowed_types = getattr(Basic, allowed_types)
        elif isinstance(allowed_types, (tuple, list)):
            new_allowed_types = []
            for t in allowed_types:
                if isinstance(t, str):
                    t = getattr(Basic, t)
                new_allowed_types.append(t)
            self.allowed_types = tuple(new_allowed_types)
        else:
            self.allowed_types = allowed_types
        have_been_here[i] = True
        return

    def process(self, obj, func, index = None):
        if isinstance(obj, self.allowed_types):
            if self.converter is not None:
                obj = self.converter(obj)
            return obj
        func_src = '%s:%s:function %s' % (func.func_code.co_filename, func.func_code.co_firstlineno, func.func_name)
        if index is None:
            raise ValueError('%s return value must be of type %r but got %r' % (func_src, self.allowed_types, obj))
        if isinstance(index, (int,long)):
            raise ValueError('%s %s-th argument must be of type %r but got %r' % (func_src, index, self.allowed_types, obj))
        if isinstance(index, str):
            raise ValueError('%s %r keyword argument must be of type %r but got %r' % (func_src, index, self.allowed_types, obj))
        raise NotImplementedError(`index,type(index)`)

class Memoizer:
    """ Memoizer function decorator generator.

    Features:
      - checks that function arguments have allowed types
      - optionally apply converters to arguments
      - cache the results of function calls
      - optionally apply converter to function values

    Usage:

      @Memoizer(<allowed types for argument 0>,
                MemoizerArg(<allowed types for argument 1>),
                MemoizerArg(<allowed types for argument 2>, <convert argument before function call>),
                MemoizerArg(<allowed types for argument 3>, <convert argument before function call>, name=<kw argument name>),
                ...
                return_value_converter = <None or converter function, usually makes a copy>
                )
      def function(<arguments>, <kw_argumnets>):
          ...

    Details:
      - if allowed type is string object then there Basic must have attribute
        with the string name that is used as the allowed type --- this is needed
        for applying Memoizer decorator to Basic methods when Basic definition
        is not defined.

    Restrictions:
      - arguments must be immutable
      - when function values are mutable then one must use return_value_converter to
        deep copy the returned values

    Ref: http://en.wikipedia.org/wiki/Memoization
    """

    def __init__(self, *arg_templates, **kw_arg_templates):
        new_arg_templates = []
        for t in arg_templates:
            if not isinstance(t, MemoizerArg):
                t = MemoizerArg(t)
            new_arg_templates.append(t)
        self.arg_templates = tuple(new_arg_templates)
        return_value_converter = kw_arg_templates.pop('return_value_converter', None)
        self.kw_arg_templates = kw_arg_templates.copy()
        for template in self.arg_templates:
            if template.name is not None:
                self.kw_arg_templates[template.name] = template
        if return_value_converter is None:
            self.return_value_converter = lambda obj: obj
        else:
            self.return_value_converter = return_value_converter

    def fix_allowed_types(self, have_been_here={}):
        i = id(self)
        if have_been_here.get(i): return
        for t in self.arg_templates:
            t.fix_allowed_types()
        for k,t in self.kw_arg_templates.items():
            t.fix_allowed_types()
        have_been_here[i] = True

    def __call__(self, func):
        cache = {}
        value_cache = {}
        def wrapper(*args, **kw_args):
            kw_items = tuple(kw_args.items())
            try:
                return self.return_value_converter(cache[args,kw_items])
            except KeyError:
                pass
            self.fix_allowed_types()
            new_args = tuple([template.process(a,func,i) for (a, template, i) in zip(args, self.arg_templates, range(len(args)))])
            assert len(args)==len(new_args)
            new_kw_args = {}
            for k, v in kw_items:
                template = self.kw_arg_templates[k]
                v = template.process(v, func, k)
                new_kw_args[k] = v
            new_kw_items = tuple(new_kw_args.items())
            try:
                return self.return_value_converter(cache[new_args, new_kw_items])
            except KeyError:
                r = func(*new_args, **new_kw_args)
                try:
                    try:
                        r = value_cache[r]
                    except KeyError:
                        value_cache[r] = r
                except TypeError:
                    pass
                cache[new_args, new_kw_items] = cache[args, kw_items] = r
                return self.return_value_converter(r)
        return wrapper

#####

class Basic(BasicMeths):
    """
    Base class for all objects in sympy.

    Conventions
    ===========

    1)
    When you want to access parameters of some instance, always use [].
    Example:

    In [2]: cot(x)[:]
    Out[2]: (x,)

    In [3]: cot(x)[0]
    Out[3]: x

    In [4]: (x*y)[:]
    Out[4]: (x, y)

    In [5]: (x*y)[1]
    Out[5]: y


    2) Never use internal methods or variables (the ones prefixed with "_").
    Example:

    In [6]: cot(x)._args    #don't use this, use cot(x)[:] instead
    Out[6]: (x,)


    """

    def __new__(cls, *args, **assumptions):
        obj = object.__new__(cls)
        obj.assume(**assumptions)
        obj._mhash = None # will be set by BasicMeths.__hash__ method.
        obj._args = args  # all items in args must be Basic objects
        return obj

    @staticmethod
    def sympify(a, sympify_lists=False):
        """Converts an arbitrary expression to a type that can be used
           inside sympy. For example, it will convert python int's into
           instance of sympy.Rational, floats into intances of sympy.Real,
           etc. It is also able to coerce symbolic expressions which does
           inherit after Basic. This can be useful in cooperation with SAGE.

           It currently accepts as arguments:
               - any object defined in sympy (except maybe matrices [TODO])
               - standard numeric python types: int, long, float, Decimal
               - strings (like "0.09" or "2e-19")

           If sympify_lists is set to True then sympify will also accept
           lists, tuples and sets. It will return the same type but with
           all of the entries sympified.

           If the argument is already a type that sympy understands, it will do
           nothing but return that value. This can be used at the begining of a
           function to ensure you are working with the correct type.

           >>> from sympy import *

           >>> sympify(2).is_integer
           True
           >>> sympify(2).is_real
           True

           >>> sympify(2.0).is_real
           True
           >>> sympify("2.0").is_real
           True
           >>> sympify("2e-45").is_real
           True

        """
        if isinstance(a, BasicType):
            return a
        if isinstance(a, Basic):
            return a
        elif isinstance(a, bool):
            raise NotImplementedError("bool support")
        elif isinstance(a, (int, long)):
            return Basic.Integer(a)
        elif isinstance(a, (float, decimal.Decimal)):
            return Basic.Real(a)
        elif isinstance(a, complex):
            real, imag = map(Basic.sympify, (a.real, a.imag))
            ireal, iimag = int(real), int(imag)

            if ireal + iimag*1j == a:
                return ireal + iimag*Basic.ImaginaryUnit()
            return real + Basic.ImaginaryUnit() * imag
        elif isinstance(a, (list,tuple)) and len(a) == 2:
            return Basic.Interval(*a)
        elif isinstance(a, (list,tuple,set)) and sympify_lists:
            return type(a)([Basic.sympify(x, True) for x in a])
        else:
            # XXX this is here because of cyclic-import issues
            from sympy.matrices import Matrix

            if isinstance(a, Matrix):
                raise NotImplementedError('matrix support')

            if not isinstance(a, str):
                # At this point we were given an arbitrary expression
                # which does not inherit after Basic. This may be
                # SAGE's expression (or something alike) so take
                # its normal form via str() and try to parse it.
                a = str(a)

            try:
                return parser.Expr(a).tosymbolic()
            except:
                pass
        raise ValueError("%r is NOT a valid SymPy expression" % a)

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
            type = Basic.sympify(type).__class__
        if isinstance(self, Atom):
            if type is None or isinstance(self, type):
                result.add(self)
        else:
            for obj in self:
                result = result.union(obj.atoms(type=type))
        return result

    def is_hypergeometric(self, arg):
        from sympy.simplify import hypersimp
        return hypersimp(self, arg, simplify=False) is not None

    def is_fraction(self, syms):
        p, q = self.as_numer_denom()

        if p.is_polynomial(*syms):
            if q.is_polynomial(*syms):
                return True

        return False

    def _eval_is_polynomial(self, syms):
        return

    def is_polynomial(self, *syms):
        if syms:
            syms = map(Basic.sympify, syms)
        else:
            syms = list(self.atoms(type=Basic.Symbol))

        if not syms: # constant polynomial
            return True
        else:
            return self._eval_is_polynomial(syms)

    def as_polynomial(self, *syms, **kwargs):
        return Basic.Polynomial(self, var=(syms or None), **kwargs)

    def _eval_subs(self, old, new):
        if self==old:
            return new
        return self

    @cache_it_immutable
    def subs(self, old, new):
        """Substitutes an expression old -> new."""
        old = Basic.sympify(old)
        new = Basic.sympify(new)

        # TODO This code is a start for issue 264. Currently, uncommenting
        #      this code will break A LOT of tests!
        #
        #if not old.is_dummy:
        #    exclude = ['dummy', 'comparable']
        #    for x in self._assume_defined:
        #        if x in exclude: continue
        #        old_val = getattr(old, 'is_' + x)
        #        if old_val is not None and old_val != getattr(new, 'is_' + x):
        #            raise ValueError("Cannot substitute '%s' for '%s' because assumptions do not match" % (str(new), str(old)))

        return self._eval_subs(old, new)

    def _seq_subs(self, old, new):
        if self==old:
            return new
        #new functions are initialized differently, than old functions
        from sympy.core.function import FunctionClass
        if isinstance(self.func, FunctionClass):
            args = self[:]
        else:
            args = (self.func,)+self[:]
        return self.__class__(*[s.subs(old, new) for s in args])

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
        p = Basic.sympify(patterns[0])
        if isinstance(p, Basic.Symbol) and not isinstance(p, Basic.Wild): # speeds up
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
        if not isinstance(self, Basic.Apply):
            args = self[:]
        else:
            args = (self.func,)+self[:]
        for e in args:
            if e.has(p):
                return True
        return False

    def _eval_derivative(self, s):
        return

    def _eval_integral(self, s):
        return

    def _eval_defined_integral(self, s, a, b):
        return

    def _eval_apply(self, *args, **assumptions):
        return

    def _eval_fapply(self, *args, **assumptions):
        return

    def _eval_fpower(b, e):
        return

    def _eval_apply_power(self,b,e):
        return

    def _eval_apply_evalf(self,*args):
        return

    def _eval_eq_nonzero(self, other):
        return

    def _eval_apply_subs(self, *args):
        return

    def _calc_apply_positive(self, *args):
        return

    def _calc_apply_real(self, *args):
        return

    def _eval_conjugate(self):
        if self.is_real:
            return self

    def conjugate(self):
        return S.Conjugate(self)

    def subs_dict(self, old_new_dict):
        r = self
        for old,new in old_new_dict.items():
            r = r.subs(old,new)
        return r

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
        from sympy.core.mul import Mul
        from sympy.core.power import Pow

        # weed out negative one prefixes
        sign = 1
        if isinstance(pattern,Mul) and pattern[0] == -1:
          pattern = -pattern; sign = -sign
        if isinstance(expr, Mul) and expr[0] == -1:
          expr = -expr; sign = -sign

        if evaluate:
            pat = pattern
            for old,new in repl_dict.items():
                pat = pat.subs(old, new)
            if pat!=pattern:
                return pat.matches(expr, repl_dict)
        expr = Basic.sympify(expr)
        if not isinstance(expr, pattern.__class__):
            from sympy.core.numbers import Rational
            # if we can omit the first factor, we can match it to sign * one
            if isinstance(pattern, Mul) and Mul(*pattern[1:]) == expr:
               return pattern[0].matches(Rational(sign), repl_dict, evaluate)
            # two-factor product: if the 2nd factor matches, the first part must be sign * one
            if isinstance(pattern, Mul) and len(pattern[:]) == 2:
               dd = pattern[1].matches(expr, repl_dict, evaluate)
               if dd == None: return None
               dd = pattern[0].matches(Rational(sign), dd, evaluate)
               return dd
            return None

        if len(pattern[:])==0:
            if pattern==expr:
                return repl_dict
            return None
        d = repl_dict.copy()

        # weed out identical terms
        if isinstance(expr, (Basic.Apply, Basic.FApply)):
            pp = list((pattern.func,)+pattern[:])
        else:
            pp = list(pattern)
        if isinstance(expr, (Basic.Apply, Basic.FApply)):
            ee = list((expr.func,)+expr[:])
        else:
            ee = list(expr)
        for p in pattern:
          for e in expr:
            if e == p:
              if e in ee: ee.remove(e)
              if p in pp: pp.remove(p)

        # only one symbol left in pattern -> match the remaining expression
        from sympy.core.symbol import Wild
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

          pattern.subs_dict(self.match(pattern)) == self

        """
        pattern = Basic.sympify(pattern)
        return pattern.matches(self, {})

    def solve4linearsymbol(eqn, rhs, symbols = None):
        """ Solve equation
          eqn == rhs
        with respect to some linear symbol in eqn.
        Returns (symbol, solution). If eqn is nonlinear
        with respect to all symbols, then return
        trivial solution (eqn, rhs).
        """
        if isinstance(eqn, Basic.Symbol):
            return (eqn, rhs)
        if symbols is None:
            symbols = eqn.atoms(type=Basic.Symbol)
        if symbols:
            # find  symbol
            for s in symbols:
                deqn = eqn.diff(s)
                if isinstance(deqn.diff(s), Basic.Zero):
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
        s = d[r] = Basic.Temporary()
        return s

    def splitter(self):
        d = {}
        r = self._calc_splitter(d)
        l = [(s.dummy_index,s,e) for e,s in d.items()]
        l.sort()
        return [(s,e) for i,s,e in l]

    @cache_it_immutable
    def count_ops(self, symbolic=True):
        """ Return the number of operations in expressions.

        Examples:
        >>> (1+a+b**2).count_ops()
        POW + 2 * ADD
        >>> (sin(x)*x+sin(x)**2).count_ops()
        ADD + MUL + POW + 2 * SIN
        """
        return Basic.Integer(len(self[:])-1) + sum([t.count_ops(symbolic=symbolic) for t in self])

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
        terms = [ term.doit(**hints) for term in self ]
        return self.__class__(*terms, **self._assumptions)

    ###########################################################################
    ################# EXPRESSION REPRESENTATION METHODS #######################
    ###########################################################################

    def _eval_expand_basic(self, *args):
        if isinstance(self, Atom):
            return self
        if not isinstance(self, Basic.Apply):
            sargs = self[:]
        else:
            sargs = (self.func,)+self[:]
        terms = [ term._eval_expand_basic(*args) for term in sargs ]
        return self.__class__(*terms, **self._assumptions)

    def _eval_expand_power(self, *args):
        if isinstance(self, Atom):
            return self
        terms = [ term._eval_expand_power(*args) for term in self._args ]
        return self.__class__(*terms, **self._assumptions)

    def _eval_expand_complex(self, *args):
        if isinstance(self, Atom):
            return self
        terms = [ term._eval_expand_complex(*args) for term in self._args ]
        return self.__class__(*terms, **self._assumptions)

    def _eval_expand_trig(self, *args):
        if isinstance(self, Atom):
            return self
        terms = [ term._eval_expand_trig(*args) for term in self._args ]
        return self.__class__(*terms, **self._assumptions)

    def _eval_expand_func(self, *args):
        if isinstance(self, Atom):
            return self
        terms = [ term._eval_expand_func(*args) for term in self._args ]
        return self.__class__(*terms, **self._assumptions)

    def expand(self, *args, **hints):
        """Expand an expression based on different hints. Currently
           supported hints are basic, power, complex, trig and func.
        """
        obj = self

        for hint in hints:
            if hints[hint] == True:
                func = getattr(obj, '_eval_expand_'+hint, None)

                if func is not None:
                    obj = func(*args)

        if hints.get('basic', True):
            obj = obj._eval_expand_basic()

        return obj

    def _eval_rewrite(self, pattern, rule, **hints):
        if isinstance(self, Atom):
            return self
        terms = [ t._eval_rewrite(pattern, rule, **hints) for t in self._args ]
        return self.__class__(*terms, **self._assumptions)

    def rewrite(self, *args, **hints):
        """Rewrites expression containing applications of functions
           of one kind in terms of functions of different kind. For
           example you can rewrite trigonometric functions as complex
           exponentials or combinatorial functions as gamma function.

           As a pattern this function accepts a list of functions to
           to rewrite (instances of DefinedFunction class). As rule
           you can use string or a destinaton function instance (in
           this cas rewrite() will use tostr() method).

           There is also possibility to pass hints on how to rewrite
           the given expressions. For now there is only one such hint
           defined called 'deep'. When 'deep' is set to False it will
           forbid functions to rewrite their contents.

           >>> from sympy import *
           >>> x, y = symbols('xy')

           >>> sin(x).rewrite(sin, exp)
           -1/2*I*(-exp(-I*x) + exp(I*x))

        """
        if isinstance(self, Atom) or not args:
            return self
        else:
            pattern, rule = args[:-1], args[-1]

            if not isinstance(rule, str):
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
        if isinstance(expr, Basic.Add):
            return None
        else:
            w = Basic.Wild('w')

            coeff = self.match(w * expr)

            if coeff is not None:
                if isinstance(expr, Basic.Mul):
                    expr = expr[:]
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
           dependend on them in the other. Both parts are valid
           SymPy expressions.

           >>> from sympy import *
           >>> x, y = symbols('xy')

           >>> (x*sin(x)*cos(y)).as_independent(x)
           (cos(y), x*sin(x))

        """
        indeps, depend = [], []

        if isinstance(self, (Basic.Add, Basic.Mul)):
            terms = self[:]
        else:
            terms = [ self ]

        for term in terms:
            if term.has(*deps):
                depend.append(term)
            else:
                indeps.append(term)

        return Basic.Mul(*indeps), Basic.Mul(*depend)

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
           (-im(w) + re(z), re(w) + im(z))

        """
        expr = self.expand(complex=True)

        if not isinstance(expr, Basic.Add):
            expr = [expr]

        re_part, im_part = [], []

        for term in expr:
            coeff = term.as_coefficient(S.ImaginaryUnit)

            if coeff is None:
                re_part.append(term)
            else:
                im_part.append(coeff)

        return (Basic.Add(*re_part), Basic.Add(*im_part))

    def as_powers_dict(self):
        return { self : S.One }

    def as_base_exp(self):
        # a -> b ** e
        return self, Basic.One()

    def as_coeff_terms(self, x=None):
        # a -> c * t
        if x is not None:
            if not self.has(x):
                return self, []
        return Basic.One(), [self]

    def as_indep_terms(self, x):
        coeff, terms = self.as_coeff_terms()
        indeps = [coeff]
        new_terms = []
        for t in terms:
            if t.has(x):
                new_terms.append(x)
            else:
                indeps.append(x)
        return Basic.Mul(*indeps), Basic.Mul(*new_terms)

    def as_coeff_factors(self, x=None):
        # a -> c + f
        if x is not None:
            if not self.has(x):
                return self, []
        return Basic.Zero(), [self]

    def as_numer_denom(self):
        # a/b -> a,b
        base, exp = self.as_base_exp()
        coeff, terms = exp.as_coeff_terms()
        if coeff.is_negative:
            # b**-e -> 1, b**e
            return Basic.One(), base ** (-exp)
        return self, Basic.One()

    def as_expr_orders(self):
        """ Split expr + Order(..) to (expr, Order(..)).
        """
        l1 = []
        l2 = []
        if isinstance(self, Basic.Add):
            for f in self:
                if isinstance(f, Basic.Order):
                    l2.append(f)
                else:
                    l1.append(f)
        elif isinstance(self, Basic.Order):
            l2.append(self)
        else:
            l1.append(self)
        return Basic.Add(*l1), Basic.Add(*l2)

    def normal(self):
        n, d = self.as_numer_denom()
        if isinstance(d, Basic.One):
            return n
        return n/d

    ###################################################################################
    ##################### DERIVATIVE, INTEGRAL, FUNCTIONAL METHODS ####################
    ###################################################################################

    def diff(self, *symbols, **assumptions):
        new_symbols = []
        for s in symbols:
            s = Basic.sympify(s)
            if isinstance(s, Basic.Integer) and new_symbols:
                last_s = new_symbols.pop()
                i = int(s)
                new_symbols += [last_s] * i
            elif isinstance(s, Basic.Symbol):
                new_symbols.append(s)
            else:
                raise TypeError(".diff() argument must be Symbol|Integer instance (got %s)" % (s.__class__.__name__))
        ret = Basic.Derivative(self, *new_symbols, **assumptions)
        return ret

    def fdiff(self, *indices):
        return Basic.FApply(Basic.FDerivative(*indices), self)

    def integral(self, *symbols, **assumptions):
        new_symbols = []
        for s in symbols:
            s = Basic.sympify(s)
            if isinstance(s, Basic.Integer) and new_symbols:
                last_s = new_symbols[-1]
                i = int(s)
                new_symbols += [last_s] * (i-1)
            elif isinstance(s, (Basic.Symbol, Basic.Equality)):
                new_symbols.append(s)
            else:
                raise TypeError(".integral() argument must be Symbol|Integer|Equality instance (got %s)" % (s.__class__.__name__))
        return Basic.Integral(self, *new_symbols, **assumptions)

    def __call__(self, *args):
        return Basic.Apply(self, *args)

    def _eval_evalf(self):
        return

    def _seq_eval_evalf(self):
        return self.__class__(*[s.evalf() for s in self])

    def __float__(self):
        result = self.evalf(precision=16)

        if isinstance(result, Basic.Number):
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

    ###################################################################################
    ##################### SERIES, LEADING TERM, LIMIT, ORDER METHODS ##################
    ###################################################################################

    def series(self, x, n = 6):
        """
        Usage
        =====
            Return the Taylor series around 0 of self with respect to x until
            the n-th term (default n is 6).

        Notes
        =====
            For computing power series, use oseries() method.
        """
        x = Basic.sympify(x)
        o = Basic.Order(x**n,x)
        r = self.oseries(o)
        if r==self:
            return self
        return r + o

    @cache_it_immutable
    def oseries(self, order, _cache={}):
        """
        Return the series of an expression upto given Order symbol.
        """
        if _cache.has_key((self, order)):
            raise RuntimeError('Detected recursion while computing oseries(%s, %s)' % (self, order))
        order = Basic.Order(order)
        _cache[(self, order)] = 1
        if isinstance(order, Basic.Zero):
            del _cache[(self, order)]
            return self
        o = self.is_order
        if o is not None:
            if o.contains(order):
                del _cache[(self, order)]
                return self
        if order.contains(self):
            del _cache[(self, order)]
            return Basic.Zero()
        if len(order.symbols)>1:
            r = self
            for s in order.symbols:
                o = Basic.Order(order.expr, s)
                r = r.oseries(o)
            del _cache[(self, order)]
            return r
        x = order.symbols[0]
        if not self.has(x):
            del _cache[(self, order)]
            return self
        obj = self._eval_oseries(order)
        if obj is not None:
            obj2 = obj.expand(trig=True)
            if obj2 != obj:
                r = obj2.oseries(order)
                del _cache[(self, order)]
                return r
            del _cache[(self, order)]
            return obj2
        del _cache[(self, order)]
        raise NotImplementedError('(%s).oseries(%s)' % (self, order))

    def _eval_oseries(self, order):
        return

    def _compute_oseries(self, arg, order, taylor_term, unevaluated_func, correction = 0):
        """
        compute series sum(taylor_term(i, arg), i=0..n-1) such
        that order.contains(taylor_term(n, arg)). Assumes that arg->0 as x->0.
        """
        x = order.symbols[0]
        ln = Basic.Log()
        o = Basic.Order(arg, x)
        if isinstance(o, Basic.Zero):
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
        n = int(n)
        assert n>=0,`n`
        l = []
        g = None
        for i in xrange(n+2):
            g = taylor_term(i, arg, g)
            g = g.oseries(order)
            l.append(g)
        return Basic.Add(*l)

    def limit(self, x, xlim, direction='<'):
        """ Compute limit x->xlim.
        """
        return Basic.Limit(self, x, xlim, direction)

    def inflimit(self, x): # inflimit has its own cache
        x = Basic.sympify(x)
        return Basic.InfLimit(self, x)

    @cache_it_immutable
    def as_leading_term(self, *symbols):
        if len(symbols)>1:
            c = self
            for x in symbols:
                c = c.as_leading_term(x)
            return c
        elif not symbols:
            return self
        x = Basic.sympify(symbols[0])
        assert isinstance(x, Basic.Symbol),`x`
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
        x = Basic.sympify(x)
        wc = Basic.Wild()
        we = Basic.Wild()
        c, terms = self.as_coeff_terms()
        p  = wc*x**we
        d = self.match(p)
        if d is not None:
            return d[wc], d[we]
        return self, Basic.Zero()

    def ldegree(self, x):
        x = Basic.sympify(x)
        c,e = self.as_leading_term(x).as_coeff_exponent(x)
        if not c.has(x):
            return e
        raise ValueError("cannot compute ldegree(%s, %s), got c=%s" % (self, x, c))

    def leadterm(self, x):
        x = Basic.sympify(x)
        c,e = self.as_leading_term(x).as_coeff_exponent(x)
        if not c.has(x):
            return c,e
        raise ValueError("cannot compute ldegree(%s, %s), got c=%s" % (self, x, c))

    ##########################################################################
    ##################### END OF BASIC CLASS #################################
    ##########################################################################

class Atom(Basic):

    precedence = Basic.Atom_precedence

    def _eval_derivative(self, s):
        if self==s: return Basic.One()
        return Basic.Zero()

    def pattern_match(pattern, expr, repl_dict):
        if pattern==expr:
            return repl_dict
        return None

    def as_numer_denom(self):
        return self, Basic.One()

    def _calc_splitter(self, d):
        return self

    def count_ops(self, symbolic=True):
        return Basic.Zero()

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


class Singleton(Basic):
    """ Singleton object.
    """

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
    E.g. S.Exp == Basic.Exp()
    """

    def __getattr__(self, clsname):
        if clsname == "__repr__":
            return lambda: "S"
        obj = Singleton.__dict__.get(clsname)
        if obj is None:
            cls = getattr(Basic, clsname)
            assert issubclass(cls, Singleton),`cls`
            obj = cls()
            setattr(self, clsname, obj)
        return obj

S = SingletonFactory()

import parser
