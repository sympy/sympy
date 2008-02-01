"""Base class for all objects in sympy"""

type_class = type

import decimal
from basic_methods import BasicMeths, BasicType, MetaBasicMeths
from assumptions import AssumeMeths
from cache import cache_it, cache_it_immutable, Memoizer, MemoizerArg

class SympifyError(ValueError):
    def __init__(self, expr, base_exc):
        self.expr = expr
        self.base_exc = base_exc
    def __str__(self):
        return "Sympify of expression '%s' failed, because of exception being raised:\n%s: %s" % (self.expr, self.base_exc.__class__.__name__, str(self.base_exc))


def repr_level(flag=None, _cache=[1]):
    if flag is None:
        return _cache[0]
    old_flag = _cache[0]
    _cache[0] = max(0, min(2, int(flag))) # restrict to 0,1,2
    return old_flag



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
        obj._mhash = None # will be set by __hash__ method.
        obj._args = args  # all items in args must be Basic objects
        return obj

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


    @staticmethod
    def sympify(a, sympify_lists=False, locals= {}):
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
        if isinstance(a, Basic):
            return a
        if isinstance(a, BasicType):
            return a
        elif isinstance(a, bool):
            raise NotImplementedError("bool support")
        elif isinstance(a, (int, long)):
            return C.Integer(a)
        elif isinstance(a, (float, decimal.Decimal)):
            return C.Real(a)
        elif isinstance(a, complex):
            real, imag = map(Basic.sympify, (a.real, a.imag))
            ireal, iimag = int(real), int(imag)

            if ireal + iimag*1j == a:
                return ireal + iimag*S.ImaginaryUnit
            return real + S.ImaginaryUnit * imag
        elif (a.__class__ in [list,tuple]) and len(a) == 2:
            # isinstance causes problems in the issue #432, so we use .__class__
            return C.Interval(*a)
        elif isinstance(a, (list,tuple,set)) and sympify_lists:
            return type(a)([Basic.sympify(x, True) for x in a])
        elif hasattr(a, "_sympy_"):
            # the "a" implements _sympy_() method, that returns a SymPy
            # expression (by definition), so we just use it
            return a._sympy_()
        else:
            # XXX this is here because of cyclic-import issues
            from sympy.matrices import Matrix
            from sympy.polynomials import Polynomial

            if isinstance(a, Polynomial):
                return a
            if isinstance(a, Matrix):
                raise NotImplementedError('matrix support')

            if not isinstance(a, str):
                # At this point we were given an arbitrary expression
                # which does not inherit from Basic and doesn't implement
                # _sympy_ (which is a canonical and robust way to convert
                # anything to SymPy expression). 
                # 
                # As a last chance, we try to take "a"'s  normal form via str()
                # and try to parse it. If it fails, then we have no luck and
                # return an exception
                a = str(a)

            try:
                return ast_parser.SymPyParser(local_dict=locals).parse_expr(a)
            except Exception, exc:
                raise SympifyError(a, exc)
        raise ValueError("%r is NOT a valid SymPy expression" % a)


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
            type = Basic.sympify(type).__class__
        if isinstance(self, Atom):
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
        """Returns True if self is a number (like 1, or 1+log(2)), and False
        otherwise (e.g. 1+x).""" 
        return len(self.atoms(C.Symbol)) == 0

    @property
    def args(self):
        """Returns a tuple of arguments of "self".

        Example
        -------
        In [2]: cot(x).args[:]
        Out[2]: (x,)

        In [3]: cot(x).args[0]
        Out[3]: x

        In [4]: (x*y).args[:]
        Out[4]: (x, y)

        In [5]: (x*y).args[1]
        Out[5]: y

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
            syms = map(Basic.sympify, syms)
        else:
            syms = list(self.atoms(type=C.Symbol))

        if not syms: # constant polynomial
            return True
        else:
            return self._eval_is_polynomial(syms)

    def as_polynomial(self, *syms, **kwargs):
        return C.Polynomial(self, var=(syms or None), **kwargs)

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

        #print self, old, new
        return self._eval_subs(old, new)

    def _seq_subs(self, old, new):
        if self==old:
            return new
        #new functions are initialized differently, than old functions
        from sympy.core.function import FunctionClass
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
        if isinstance(p, C.Symbol) and not isinstance(p, C.Wild): # speeds up
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

    def subs_dict(self, old_new_dict):
        """Substitutes "self" using the keys and values from the "old_new_dict".

           This correctly handles "overlapping" keys,
           e.g. when doing substituion like:

             e.subs_dict(x: y, exp(x): y)

           exp(x) is always substituted before x
        """
        r = self

        oldnew = [(o,n) for (o,n) in old_new_dict.items()]


        # Care needs to be taken for cases like these:
        # 
        # consider
        #     e = x*exp(x)
        # 
        # when we do
        #     e.subs_dict(x: y, exp(x): y)
        # 
        # the result depends in which order elementary substitutions are done.
        # 
        # So we *have* to proceess 'exp(x)' first. This is achieved by topologically
        # sorting substitutes by 'in' criteria - if "exp(x)" contains "x", it
        # gets substituted first.


        # let's topologically sort oldnew, so we have more general terms in the
        # beginning     (e.g. [exp(x), x] is not the same as [x, exp(x)] when
        #                doing substitution)
        def in_cmp(a,b):
            if a[0] in b[0]:
                return 1
            else:
                return -1

        oldnew.sort(cmp=in_cmp)

        for old,new in oldnew:
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
        if isinstance(pattern,Mul) and pattern.args[0] == -1:
          pattern = -pattern; sign = -sign
        if isinstance(expr, Mul) and expr.args[0] == -1:
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
            if isinstance(pattern, Mul) and Mul(*pattern.args[1:]) == expr:
               return pattern.args[0].matches(Rational(sign), repl_dict, evaluate)
            # two-factor product: if the 2nd factor matches, the first part must be sign * one
            if isinstance(pattern, Mul) and len(pattern.args[:]) == 2:
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
        if isinstance(eqn, C.Symbol):
            return (eqn, rhs)
        if symbols is None:
            symbols = eqn.atoms(type=C.Symbol)
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
        s = d[r] = C.Temporary()
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
        return C.Integer(len(self[:])-1) + sum([t.count_ops(symbolic=symbolic) for t in self])

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

    def _eval_expand_basic(self, *args):
        if isinstance(self, Atom):
            return self
        sargs = self.args[:]
        terms = [ term._eval_expand_basic(*args) for term in sargs ]
        return self.__class__(*terms, **self._assumptions)

    def _eval_expand_power(self, *args):
        if isinstance(self, Atom):
            return self
        if not isinstance(self, C.Apply):   # FIXME Apply -> Function
            sargs = self[:]
        else:
            sargs = (self.func,)+self[:]
        terms = [ term._eval_expand_power(*args) for term in sargs ]
        return self.__class__(*terms, **self._assumptions)

    def _eval_expand_complex(self, *args):
        if isinstance(self, Atom):
            return self
        sargs = self.args[:]
        terms = [ term._eval_expand_complex(*args) for term in sargs ]
        return self.__class__(*terms, **self._assumptions)

    def _eval_expand_trig(self, *args):
        if isinstance(self, Atom):
            return self
        sargs = self.args[:]
        terms = [ term._eval_expand_trig(*args) for term in sargs ]
        return self.__class__(*terms, **self._assumptions)

    def _eval_expand_func(self, *args):
        if isinstance(self, Atom):
            return self
        sargs = self.args
        terms = [ term._eval_expand_func(*args) for term in sargs ]
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
                # XXX move me out of here (cyclic imports)
                from function import FunctionClass

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
        if isinstance(expr, C.Add):
            return None
        else:
            w = C.Wild('w')

            coeff = self.match(w * expr)

            if coeff is not None:
                if isinstance(expr, C.Mul):
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
           dependend on them in the other. Both parts are valid
           SymPy expressions.

           >>> from sympy import *
           >>> x, y = symbols('xy')

           >>> (x*sin(x)*cos(y)).as_independent(x)
           (cos(y), x*sin(x))

        """
        indeps, depend = [], []

        if isinstance(self, (C.Add, C.Mul)):
            terms = self.args[:]
        else:
            terms = [ self ]

        for term in terms:
            if term.has(*deps):
                depend.append(term)
            else:
                indeps.append(term)

        return C.Mul(*indeps), C.Mul(*depend)

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

        if not isinstance(expr, C.Add):
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

        return (C.Add(*re_part), C.Add(*im_part))

    def as_powers_dict(self):
        return { self : S.One }

    def as_base_exp(self):
        # a -> b ** e
        return self, S.One

    def as_coeff_terms(self, x=None):
        # a -> c * t
        if x is not None:
            if not self.has(x):
                return self, []
        return S.One, [self]

    def as_indep_terms(self, x):
        coeff, terms = self.as_coeff_terms()
        indeps = [coeff]
        new_terms = []
        for t in terms:
            if t.has(x):
                new_terms.append(x)
            else:
                indeps.append(x)
        return C.Mul(*indeps), C.Mul(*new_terms)

    def as_coeff_factors(self, x=None):
        # a -> c + f
        if x is not None:
            if not self.has(x):
                return self, []
        return S.Zero, [self]

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
        if isinstance(self, C.Add):
            for f in self:
                if isinstance(f, C.Order):
                    l2.append(f)
                else:
                    l1.append(f)
        elif isinstance(self, C.Order):
            l2.append(self)
        else:
            l1.append(self)
        return C.Add(*l1), C.Add(*l2)

    def normal(self):
        n, d = self.as_numer_denom()
        if d is S.One:
            return n
        return n/d

    ###################################################################################
    ##################### DERIVATIVE, INTEGRAL, FUNCTIONAL METHODS ####################
    ###################################################################################

    def diff(self, *symbols, **assumptions):
        new_symbols = []
        for s in symbols:
            s = Basic.sympify(s)
            if isinstance(s, C.Integer) and new_symbols:
                last_s = new_symbols.pop()
                i = int(s)
                new_symbols += [last_s] * i
            elif isinstance(s, C.Symbol):
                new_symbols.append(s)
            else:
                raise TypeError(".diff() argument must be Symbol|Integer instance (got %s)" % (s.__class__.__name__))
        if not assumptions.has_key("evaluate"):
            assumptions["evaluate"] = True
        ret = C.Derivative(self, *new_symbols, **assumptions)
        return ret

    def fdiff(self, *indices):
        # FIXME FApply -> ?
        return C.FApply(C.FDerivative(*indices), self)

    def integral(self, *symbols, **assumptions):
        new_symbols = []
        for s in symbols:
            s = Basic.sympify(s)
            if isinstance(s, C.Integer) and new_symbols:
                last_s = new_symbols[-1]
                i = int(s)
                new_symbols += [last_s] * (i-1)
            elif isinstance(s, (C.Symbol, C.Equality)):
                new_symbols.append(s)
            else:
                raise TypeError(".integral() argument must be Symbol|Integer|Equality instance (got %s)" % (s.__class__.__name__))
        return C.Integral(self, *new_symbols, **assumptions)

    #XXX fix the removeme
    def __call__(self, *args, **removeme):
        return C.Function(self[0])(*args)

    def _eval_evalf(self):
        return

    def _seq_eval_evalf(self):
        return self.__class__(*[s.evalf() for s in self.args])

    def __float__(self):
        result = self.evalf(precision=16)

        if isinstance(result, C.Number):
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
        x = Basic.sympify(x)
        point = Basic.sympify(point)
        if point != 0:
            raise NotImplementedError("series expansion around arbitrary point")
            #self = self.subs(x, x + point)
        o = C.Order(x**n,x)
        r = self.oseries(o)
        if r==self:
            return self
        return r + o

    @cache_it_immutable
    def oseries(self, order):
        """
        Return the series of an expression upto given Order symbol (without the
        actual O term).

        The general philosophy is this: simply start with the most simple
        taylor (laurent) term and calculate one be one and use
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
        return C.Add(*l)

    def limit(self, x, xlim, direction='<'):
        """ Compute limit x->xlim.
        """
        from sympy.series.limits_series import Limit
        return Limit(self, x, xlim, direction)

    def inflimit(self, x): # inflimit has its own cache
        x = Basic.sympify(x)
        from sympy.series.limits_series import InfLimit
        return InfLimit(self, x)

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
        assert isinstance(x, C.Symbol),`x`
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
        wc = C.Wild()
        we = C.Wild()
        c, terms = self.as_coeff_terms()
        p  = wc*x**we
        d = self.match(p)
        if d is not None:
            return d[wc], d[we]
        return self, S.Zero

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


# XXX BasicMeths.compare needs Basic
import basic_methods
basic_methods.Basic = Basic
del basic_methods

    ##########################################################################
    ##################### END OF BASIC CLASS #################################
    ##########################################################################

class Atom(Basic):

    precedence = Basic.Atom_precedence

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
            cls = MetaBasicMeths.classnamespace[name]
        except KeyError:
            raise AttributeError("No SymPy class '%s'" % name)

        setattr(self, name, cls)
        return cls

C = ClassesRegistry()

# XXX this is ugly, but needed for Memoizer('str', ...) to work
import cache
cache.C = C
del cache


import ast_parser
