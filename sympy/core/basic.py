"""Base class for all objects in sympy"""

from decorators import _sympifyit
from assumptions import AssumeMeths, make__get_assumption
from cache import cacheit
from core import BasicMeta, BasicType

class Basic(AssumeMeths):
    """
    Base class for all objects in sympy.

    Conventions:

    1)
    When you want to access parameters of some instance, always use .args:
    Example:

    >>> from sympy import symbols, cot
    >>> from sympy.abc import x, y

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

    __slots__ = ['_mhash',              # hash value
                 '_args',               # arguments
                 '_assume_type_keys',   # assumptions typeinfo keys
                ]

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
    is_Derivative = False
    is_Piecewise = False
    is_Poly = False
    is_AlgebraicNumber = False

    def __new__(cls, *args, **assumptions):
        obj = object.__new__(cls)

        # FIXME we are slowed a *lot* by Add/Mul passing is_commutative as the
        # only assumption.
        #
        # .is_commutative is not an assumption -- it's like typeinfo!!!
        # we should remove it.

        # initially assumptions are shared between instances and class
        obj._assumptions  = cls.default_assumptions
        obj._a_inprogress = []

        # NOTE this could be made lazy -- probably not all instances will need
        # fully derived assumptions?
        if assumptions:
            obj._learn_new_facts(assumptions)
            #                      ^
            # FIXME this is slow   |    another NOTE: speeding this up is *not*
            #        |             |    important. say for %timeit x+y most of
            # .------'             |    the time is spent elsewhere
            # |                    |
            # |  XXX _learn_new_facts  could be asked about what *new* facts have
            # v  XXX been learned -- we'll need this to append to _hashable_content
            basek = set(cls.default_assumptions.keys())
            k2    = set(obj._assumptions.keys())
            newk  = k2.difference(basek)

            obj._assume_type_keys = frozenset(newk)
        else:
            obj._assume_type_keys = None

        obj._mhash = None # will be set by __hash__ method.
        obj._args = args  # all items in args must be Basic objects
        return obj


    # XXX better name?
    @property
    def assumptions0(self):
        """
        Return object `type` assumptions.

        For example:

          Symbol('x', real=True)
          Symbol('x', integer=True)

        are different objects. In other words, besides Python type (Symbol in
        this case), the initial assumptions are also forming their typeinfo.

        Example:

        >>> from sympy import Symbol
        >>> from sympy.abc import x
        >>> x.assumptions0
        {}
        >>> x = Symbol("x", positive=True)
        >>> x.assumptions0
        {'commutative': True, 'complex': True, 'imaginary': False,
        'negative': False, 'nonnegative': True, 'nonpositive': False,
        'nonzero': True, 'positive': True, 'real': True, 'zero': False}

        """

        cls = type(self)
        A   = self._assumptions

        # assumptions shared:
        if A is cls.default_assumptions or (self._assume_type_keys is None):
            assumptions0 = {}
        else:
            assumptions0 = dict( (k, A[k]) for k in self._assume_type_keys )

        return assumptions0


    def new(self, *args):
        """
        Create new 'similar' object.

        this is conceptually equivalent to:

          type(self) (*args)

        but takes type assumptions into account. See also: assumptions0

        Example:

        >>> from sympy.abc import x
        >>> x.new("x")
        x

        """
        obj = self.func(*args, **self.assumptions0)
        return obj


    # NOTE NOTE NOTE
    # --------------
    #
    # new-style classes + __getattr__ is *very* slow!

    # def __getattr__(self, name):
    #     raise Warning('no way, *all* attribute access will be 2.5x slower')

    # here is what we do instead:
    for k in AssumeMeths._assume_defined:
        exec "is_%s  = property(make__get_assumption('Basic', '%s'))" % (k,k)
    del k

    # NB: there is no need in protective __setattr__

    def __getnewargs__(self):
        """ Pickling support.
        """
        return tuple(self.args)

    def __hash__(self):
        # hash cannot be cached using cache_it because infinite recurrence
        # occurs as hash is needed for setting cache dictionary keys
        h = self._mhash
        if h is None:
            h = (type(self).__name__,) + self._hashable_content()

            if self._assume_type_keys is not None:
                a = []
                kv= self._assumptions
                for k in sorted(self._assume_type_keys):
                    a.append( (k, kv[k]) )

                h = hash( h + tuple(a) )

            else:
                h = hash( h )


            self._mhash = h
            return h

        else:
            return h

    def _hashable_content(self):
        # If class defines additional attributes, like name in Symbol,
        # then this method should be updated accordingly to return
        # relevant attributes as tuple.
        return self._args

    def compare(self, other):
        """
        Return -1,0,1 if the object is smaller, equal, or greater than other.

        Not in the mathematical sense. If the object is of a different type
        from the "other" then their classes are ordered according to
        the sorted_classes list.

        Example:

        >>> from sympy.abc import x, y
        >>> x.compare(y)
        -1
        >>> x.compare(x)
        0
        >>> y.compare(x)
        1

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
    def _compare_pretty(a, b):
        from sympy.series.order import Order
        if isinstance(a, Order) and not isinstance(b, Order):
            return 1
        if not isinstance(a, Order) and isinstance(b, Order):
            return -1

        # FIXME this produces wrong ordering for 1 and 0
        # e.g. the ordering will be 1 0 2 3 4 ...
        # because 1 = x^0, but 0 2 3 4 ... = x^1
        from sympy.core.symbol import Wild
        p1, p2, p3 = Wild("p1"), Wild("p2"), Wild("p3")
        r_a = a.match(p1 * p2**p3)
        r_b = b.match(p1 * p2**p3)
        if r_a is not None and r_b is not None:
            c = Basic.compare(r_a[p3], r_b[p3])
            if c!=0:
                return c

        return Basic.compare(a,b)

    @staticmethod
    def compare_pretty(a, b):
        """
        Is a > b in the sense of ordering in printing?

        yes ..... return 1
        no ...... return -1
        equal ... return 0

        Strategy:

        It uses Basic.compare as a fallback, but improves it in many cases,
        like x**3, x**4, O(x**3) etc. In those simple cases, it just parses the
        expression and returns the "sane" ordering such as:

        1 < x < x**2 < x**3 < O(x**4) etc.

        Example:

        >>> from sympy.abc import x
        >>> from sympy import Basic
        >>> Basic._compare_pretty(x, x**2)
        -1
        >>> Basic._compare_pretty(x**2, x**2)
        0
        >>> Basic._compare_pretty(x**3, x**2)
        1

        """
        try:
            a = _sympify(a)
        except SympifyError:
            pass

        try:
            b = _sympify(b)
        except SympifyError:
            pass

        # both objects are non-SymPy
        if (not isinstance(a, Basic)) and (not isinstance(b, Basic)):
            return cmp(a,b)

        if not isinstance(a, Basic):
            return -1   # other < sympy

        if not isinstance(b, Basic):
            return +1   # sympy > other

        # now both objects are from SymPy, so we can proceed to usual comparison
        return Basic._compare_pretty(a, b)


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

    def __repr__(self):
        from sympy.printing import sstr
        return sstr(self)

    def __str__(self):
        from sympy.printing import sstr
        return sstr(self)

    def atoms(self, *types):
        """Returns the atoms that form the current object.

           By default, only objects that are truly atomic and can't
           be divided into smaller pieces are returned: symbols, numbers,
           and number symbols like I and pi. It is possible to request
           atoms of any type, however, as demonstrated below.

           Examples::

           >>> from sympy import I, pi, sin
           >>> from sympy.abc import x, y
           >>> list((1+x+2*sin(y+I*pi)).atoms())
           [y, I, 2, x, 1, pi]

           If one or more types are given, the results will contain only
           those types of atoms::

           Examples::

           >>> from sympy import Number, NumberSymbol, Symbol
           >>> sorted((1+x+2*sin(y+I*pi)).atoms(Symbol))
           [x, y]

           >>> sorted((1+x+2*sin(y+I*pi)).atoms(Number))
           [1, 2]

           >>> sorted((1+x+2*sin(y+I*pi)).atoms(Number, NumberSymbol))
           [1, 2, pi]

           >>> sorted((1+x+2*sin(y+I*pi)).atoms(Number, NumberSymbol, I))
           [1, 2, pi, I]

           Note that I (imaginary unit) and zoo (complex infinity) are special
           types of number symbols and are not part of the NumberSymbol class.

           The type can be given implicitly, too::

           >>> sorted((1+x+2*sin(y+I*pi)).atoms(x)) # x is a Symbol
           [x, y]

           Be careful to check your assumptions when using the implicit option
           since S(1).is_Integer = True but type(S(1)) is One, a special type
           of sympy atom, while type(S(2)) is type Integer and will find all
           integers in an expression::

           >>> from sympy import S
           >>> sorted((1+x+2*sin(y+I*pi)).atoms(S(1)))
           [1]

           >>> sorted((1+x+2*sin(y+I*pi)).atoms(S(2)))
           [1, 2]

           Finally, arguments to atoms() can select more than atomic atoms: any
           sympy type (loaded in core/__init__.py) can be listed as an argument
           and those types of "atoms" as found in scanning the arguments of the
           expression nonrecursively::

           >>> from sympy import Function, Mul
           >>> sorted((1+x+2*sin(y+I*pi)).atoms(Function))
           [sin(y + pi*I)]

           >>> sorted((1+x+2*sin(y+I*pi)).atoms(Mul))
           [2*sin(y + pi*I)]


        """


        def _atoms(expr, typ):
            """Helper function for recursively denesting atoms"""
            if isinstance(expr, Basic):
                if expr.is_Atom and len(typ) == 0: # if we haven't specified types
                    return [expr]
                else:
                    try:
                        if isinstance(expr, typ): return [expr]
                    except TypeError:
                        #one or more types is in implicit form
                        for t in typ:
                            if type(type(t)) is type:
                                if isinstance(expr, t): return [expr]
                            else:
                                if isinstance(expr, type(t)): return [expr]
            result = []
            #search for a suitable iterator
            if isinstance(expr, Basic):
                iter = expr.iter_basic_args()
            elif isinstance(expr, (tuple, list)):
                iter = expr.__iter__()
            else:
                iter = []

            for obj in iter:
                result.extend(_atoms(obj, typ))
            return result

        return set(_atoms(self, typ=types))

    def is_hypergeometric(self, k):
        from sympy.simplify import hypersimp
        return hypersimp(self, k) is not None

    @property
    def is_number(self):
        """Returns True if 'self' is a number.

           >>> from sympy import log
           >>> from sympy.abc import x, y

           >>> x.is_number
           False
           >>> (2*x).is_number
           False
           >>> (2 + log(2)).is_number
           True

        """
        result = False
        for obj in self.iter_basic_args():
            if obj.is_number:
                result = True
            else:
                return False
        return result

    @property
    def func(self):
        """
        The top-level function in an expression.

        The following should hold for all objects::

            >> x == x.func(*x.args)

        Example:

        >>> from sympy.abc import x
        >>> a = 2*x
        >>> a.func
        <class 'sympy.core.mul.Mul'>
        >>> a.args
        (2, x)
        >>> a.func(*a.args)
        2*x
        >>> a == a.func(*a.args)
        True

        """
        return self.__class__

    @property
    def args(self):
        """Returns a tuple of arguments of 'self'.

        Example:

        >>> from sympy import symbols, cot
        >>> from sympy.abc import x, y

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
        return self._args

    def iter_basic_args(self):
        """
        Iterates arguments of 'self'.

        Example:

        >>> from sympy.abc import x
        >>> a = 2*x
        >>> a.iter_basic_args()
        <tupleiterator object at 0x...>
        >>> list(a.iter_basic_args())
        [2, x]

        """
        return iter(self.args)

    def is_rational_function(self, *syms):
        """
        Test whether function is rational function - ratio of two polynomials.
        When no arguments are present, Basic.atoms(Symbol) is used instead.

        Example:

        >>> from sympy import symbols, sin
        >>> from sympy.abc import x, y

        >>> (x/y).is_rational_function()
        True

        >>> (x**2).is_rational_function()
        True

        >>> (x/sin(y)).is_rational_function(y)
        False

        """
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
            syms = list(self.atoms(C.Symbol))

        if not syms: # constant polynomial
            return True
        else:
            return self._eval_is_polynomial(syms)

    def as_poly(self, *gens, **args):
        """Converts `self` to a polynomial or returns `None`.

           >>> from sympy import Poly, sin
           >>> from sympy.abc import x, y

           >>> print (x**2 + x*y).as_poly()
           Poly(x**2 + x*y, x, y, domain='ZZ')

           >>> print (x**2 + x*y).as_poly(x, y)
           Poly(x**2 + x*y, x, y, domain='ZZ')

           >>> print (x**2 + sin(y)).as_poly(x, y)
           None

        """
        from sympy.polys import Poly, PolynomialError

        try:
            poly = Poly(self, *gens, **args)

            if not poly.is_Poly:
                return None
            else:
                return poly
        except PolynomialError:
            return None

    def as_basic(self):
        """Converts polynomial to a valid sympy expression.

           >>> from sympy import sin
           >>> from sympy.abc import x, y

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

        >>> from sympy import pi
        >>> from sympy.abc import x, y
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
            raise TypeError("subs accepts either 1 or 2 arguments")

    @cacheit
    def _subs_old_new(self, old, new):
        """Substitutes an expression old -> new."""
        old = sympify(old)
        new = sympify(new)
        return self._eval_subs(old, new)

    def _eval_subs(self, old, new):
        if self == old:
            return new
        else:
            return self.func(*[arg._eval_subs(old, new) for arg in self.args])

    def _subs_list(self, sequence):
        """
        Performs an order sensitive substitution from the
        input sequence list.

        Examples:

        >>> from sympy.abc import x, y
        >>> (x+y)._subs_list( [(x, 3),     (y, x**2)] )
        3 + x**2
        >>> (x+y)._subs_list( [(y, x**2),  (x, 3)   ] )
        12

        """
        if not isinstance(sequence, (list, tuple)):
            raise TypeError("Not an iterable container")
        result = self
        for old, new in sequence:
            if hasattr(result, 'subs'):
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

           >>> from sympy import sqrt, sin, cos, exp
           >>> from sympy.abc import x, y

           >>> from sympy.abc import a, b, c, d, e

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
        from function import FunctionClass
        if isinstance(self.func, FunctionClass):
            args = self.args
        else:
            args = (self.func,)+self
        return self.func(*[s.subs(old, new) for s in args])

    def __contains__(self, what):
        if self == what or self.is_Function and self.func is what: return True
        for x in self._args:
            # x is not necessarily of type Basic and so 'x in x == True'
            # may not hold.
            if x == what:
                return True

            # Not all arguments implement __contains__.
            try:
                if what in x:
                    return True
            except TypeError:
                continue
        return False

    @cacheit
    def has_any_symbols(self, *syms):
        """Return True if 'self' has any of the symbols.

           >>> from sympy import sin
           >>> from sympy.abc import x, y, z

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
                if type(expr) in (tuple, set, list):
                    for i in expr:
                        if search(i):
                            return True
                elif not isinstance(expr, Basic):
                    pass
                elif expr.is_Atom:
                    if expr.is_Symbol:
                        return expr in syms
                    else:
                        return False
                else:
                    for term in expr.iter_basic_args():
                        if search(term):
                            return True
                    else:
                        return False

            return search(self)

    @cacheit
    def has_all_symbols(self, *syms):
        """Return True if 'self' has all of the symbols.

           >>> from sympy import sin
           >>> from sympy.abc import x, y, z

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
                if type(expr) in (tuple, set, list):
                    for i in expr:
                        search(i)
                elif not isinstance(expr, Basic):
                    pass
                elif expr.is_Atom:
                    if expr.is_Symbol and expr in syms:
                        syms.remove(expr)
                else:
                    for term in expr.iter_basic_args():
                        if not syms:
                            break
                        else:
                            search(term)

            search(self)

            return not syms

    def has(self, *patterns):
        """
        Return True if self has any of the patterns.

        Example:
        >>> from sympy.abc import x
        >>> (2*x).has(x)
        True
        >>> (2*x/x).has(x)
        False

        """
        from sympy.utilities.iterables import flatten
        from sympy.core.symbol import Wild
        if len(patterns)>1:
            for p in patterns:
                if self.has(p):
                    return True
            return False
        elif not patterns:
            raise TypeError("has() requires at least 1 argument (got none)")
        p = sympify(patterns[0])
        if isinstance(p, BasicType):
            return bool(self.atoms(p))
        if p.is_Atom and not isinstance(p, Wild):
            return p in self.atoms(p.func)
        if p.matches(self) is not None:
            return True
        for e in flatten(self.args):
            if isinstance(e, Basic) and e.has(p):
                return True
        return False

    def matches(self, expr, repl_dict={}, evaluate=False):
        """
        Helper method for match() - switches the pattern and expr.

        Can be used to solve linear equations:
          >>> from sympy import Symbol, Wild, Integer
          >>> a,b = map(Symbol, 'ab')
          >>> x = Wild('x')
          >>> (a+b*x).matches(Integer(0))
          {x_: -a/b}

        """
        if evaluate:
            return self.subs(repl_dict).matches(expr, repl_dict)

        expr = sympify(expr)
        if not isinstance(expr, self.__class__):
            return None

        if self == expr:
            return repl_dict
        if len(self.args) != len(expr.args):
            return None

        d = repl_dict.copy()
        for arg, other_arg in zip(self.args, expr.args):
            if arg == other_arg:
                continue
            d = arg.subs(d).matches(other_arg, d)
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

        Example:

        >>> from sympy import symbols, Wild
        >>> from sympy.abc import x, y
        >>> p = Wild("p")
        >>> q = Wild("q")
        >>> r = Wild("r")
        >>> e = (x+y)**(x+y)
        >>> e.match(p**p)
        {p_: x + y}
        >>> e.match(p**q)
        {p_: x + y, q_: x + y}
        >>> e = (2*x)**2
        >>> e.match(p*q**r)
        {p_: 4, q_: x, r_: 2}
        >>> (p*q**r).subs(e.match(p*q**r))
        4*x**2

        """
        pattern = sympify(pattern)
        return pattern.matches(self, {})

    @cacheit
    def count_ops(self, symbolic=True):
        """ Return the number of operations in expressions.

        Examples:
        >>> from sympy.abc import a, b, x
        >>> from sympy import sin
        >>> (1+a+b**2).count_ops()
        POW + 2*ADD
        >>> (sin(x)*x+sin(x)**2).count_ops()
        2 + ADD + MUL + POW

        """
        return Integer(len(self)-1) + sum([t.count_ops(symbolic=symbolic) for t in self])

    def doit(self, **hints):
        """Evaluate objects that are not evaluated by default like limits,
           integrals, sums and products. All objects of this kind will be
           evaluated recursively, unless some species were excluded via 'hints'
           or unless the 'deep' hint was set to 'False'.

           >>> from sympy import Integral
           >>> from sympy.abc import x, y

           >>> 2*Integral(x, x)
           2*Integral(x, x)

           >>> (2*Integral(x, x)).doit()
           x**2

           >>> (2*Integral(x, x)).doit(deep = False)
           2*Integral(x, x)

        """
        if hints.get('deep', True):
            terms = [ term.doit(**hints) for term in self.args ]
            return self.new(*terms)
        else:
            return self

    def _eval_rewrite(self, pattern, rule, **hints):
        if self.is_Atom:
            return self
        sargs = self.args
        terms = [ t._eval_rewrite(pattern, rule, **hints) for t in sargs ]
        return self.new(*terms)

    def rewrite(self, *args, **hints):
        """Rewrites expression containing applications of functions
           of one kind in terms of functions of different kind. For
           example you can rewrite trigonometric functions as complex
           exponentials or combinatorial functions as gamma function.

           As a pattern this function accepts a list of functions to
           to rewrite (instances of DefinedFunction class). As rule
           you can use string or a destination function instance (in
           this case rewrite() will use the str() function).

           There is also possibility to pass hints on how to rewrite
           the given expressions. For now there is only one such hint
           defined called 'deep'. When 'deep' is set to False it will
           forbid functions to rewrite their contents.

           >>> from sympy import sin, exp, I
           >>> from sympy.abc import x, y

           >>> sin(x).rewrite(sin, exp)
           -I*(exp(I*x) - exp(-I*x))/2

        """
        if self.is_Atom or not args:
            return self
        else:
            pattern = args[:-1]
            rule = '_eval_rewrite_as_' + str(args[-1])

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

    def __call__(self, subsdict):
        """Use call as a shortcut for subs, but only support the dictionary version"""
        if not isinstance(subsdict, dict):
            raise TypeError("argument must be a dictionary")
        return self.subs(subsdict)

class Atom(Basic):
    """
    A parent class for atomic things. An atom is an expression with no subexpressions.

    Examples: Symbol, Number, Rational, Integer, ...
    But not: Add, Mul, Pow, ...
    """

    is_Atom = True

    __slots__ = []

    def _eval_derivative(self, s):
        if self==s: return S.One
        return S.Zero

    def matches(self, expr, repl_dict, evaluate=False):
        if self == expr:
            return repl_dict
        return None

    def _eval_subs(self, old, new):
        if self == old:
            return new
        else:
            return self

    def as_numer_denom(self):
        return self, S.One

    def count_ops(self, symbolic=True):
        return S.Zero

    def doit(self, **hints):
        return self

    def _eval_is_polynomial(self, syms):
        return True

    @property
    def is_number(self):
        return True

    def _eval_nseries(self, x, x0, n):
        return self



from sympify import _sympify, sympify, SympifyError
from singleton import S
from core import C
