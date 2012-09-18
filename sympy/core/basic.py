"""Base class for all the objects in SymPy"""
from copy import copy
from sympy.core.assumptions import ManagedProperties
from sympy.core.cache import cacheit
from sympy.core.core import BasicType, C
from sympy.core.sympify import _sympify, sympify, SympifyError
from sympy.core.compatibility import callable, reduce, cmp, iterable
from sympy.core.decorators import deprecated
from sympy.core.singleton import S

class Basic(object):
    """
    Base class for all objects in sympy.

    Conventions:

    1) Always use ``.args``, when accessing parameters of some instance:

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


    2) Never use internal methods or variables (the ones prefixed with ``_``):

        >>> cot(x)._args    # do not use this, use cot(x).args instead
        (x,)

    """
    __metaclass__ = ManagedProperties
    __slots__ = ['_mhash',              # hash value
                 '_args',               # arguments
                 '_assumptions'
                ]

    # To be overridden with True in the appropriate subclasses
    is_Atom = False
    is_Symbol = False
    is_Dummy = False
    is_Wild = False
    is_Function = False
    is_Add = False
    is_Mul = False
    is_Pow = False
    is_Number = False
    is_Float = False
    is_Rational = False
    is_Integer = False
    is_NumberSymbol = False
    is_Order = False
    is_Derivative = False
    is_Piecewise = False
    is_Poly = False
    is_AlgebraicNumber = False
    is_Relational = False
    is_Equality = False
    is_Boolean = False
    is_Not = False
    is_Matrix = False

    @property
    @deprecated(useinstead="is_Float", issue=1721, deprecated_since_version="0.7.0")
    def is_Real(self):  # pragma: no cover
        """Deprecated alias for ``is_Float``"""
        # When this is removed, remove the piece of code disabling the warning
        # from test_pickling.py
        return self.is_Float

    def __new__(cls, *args):
        obj = object.__new__(cls)
        obj._assumptions  = cls.default_assumptions
        obj._mhash = None # will be set by __hash__ method.

        obj._args = args  # all items in args must be Basic objects
        return obj

    def copy(self):
        return self.func(*self.args)

    def __reduce_ex__(self, proto):
        """ Pickling support."""
        return type(self), self.__getnewargs__(), self.__getstate__()

    def __getnewargs__(self):
        return self.args

    def __getstate__(self):
        return {}

    def __setstate__(self, state):
        for k, v in state.iteritems():
            setattr(self, k, v)

    def __hash__(self):
        # hash cannot be cached using cache_it because infinite recurrence
        # occurs as hash is needed for setting cache dictionary keys
        h = self._mhash
        if h is None:
            h = hash((type(self).__name__,) + self._hashable_content())
            self._mhash = h
        return h

    def _hashable_content(self):
        """Return a tuple of information about self that can be used to
        compute the hash. If a class defines additional attributes,
        like ``name`` in Symbol, then this method should be updated
        accordingly to return such relevent attributes.

        Defining more than _hashable_content is necessary if __eq__ has
        been defined by a class. See note about this in Basic.__eq__."""
        return self._args

    @property
    def assumptions0(self):
        """
        Return object `type` assumptions.

        For example:

          Symbol('x', real=True)
          Symbol('x', integer=True)

        are different objects. In other words, besides Python type (Symbol in
        this case), the initial assumptions are also forming their typeinfo.

        Examples
        ========

        >>> from sympy import Symbol
        >>> from sympy.abc import x
        >>> x.assumptions0
        {'commutative': True}
        >>> x = Symbol("x", positive=True)
        >>> x.assumptions0
        {'commutative': True, 'complex': True, 'hermitian': True,
        'imaginary': False, 'negative': False, 'nonnegative': True,
        'nonpositive': False, 'nonzero': True, 'positive': True, 'real': True,
        'zero': False}

        """
        return {}

    def compare(self, other):
        """
        Return -1, 0, 1 if the object is smaller, equal, or greater than other.

        Not in the mathematical sense. If the object is of a different type
        from the "other" then their classes are ordered according to
        the sorted_classes list.

        Examples
        ========

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
        if self is other:
            return 0
        c = cmp(self.__class__, other.__class__)
        if c:
            return c
        #
        st = self._hashable_content()
        ot = other._hashable_content()
        c = cmp(len(st), len(ot))
        if c:
            return c
        for l, r in zip(st, ot):
            if isinstance(l, Basic):
                c = l.compare(r)
            elif isinstance(l, frozenset):
                c = 0
            else:
                c = cmp(l, r)
            if c:
                return c
        return 0

    @staticmethod
    def _compare_pretty(a, b):
        from sympy.series.order import Order
        if isinstance(a, Order) and not isinstance(b, Order):
            return 1
        if not isinstance(a, Order) and isinstance(b, Order):
            return -1

        if a.is_Rational and b.is_Rational:
            return cmp(a.p*b.q, b.p*a.q)
        else:
            from sympy.core.symbol import Wild
            p1, p2, p3 = Wild("p1"), Wild("p2"), Wild("p3")
            r_a = a.match(p1 * p2**p3)
            if r_a and p3 in r_a:
                a3 = r_a[p3]
                r_b = b.match(p1 * p2**p3)
                if r_b and p3 in r_b:
                    b3 = r_b[p3]
                    c = Basic.compare(a3, b3)
                    if c != 0:
                        return c

        return Basic.compare(a,b)

    @staticmethod
    @deprecated(useinstead="default_sort_key", issue=1491, deprecated_since_version="0.7.2")
    def compare_pretty(a, b):
        """
        Is a > b in the sense of ordering in printing?

        THIS FUNCTION IS DEPRECATED.  Use ``default_sort_key`` instead.

        ::

          yes ..... return 1
          no ...... return -1
          equal ... return 0

        Strategy:

        It uses Basic.compare as a fallback, but improves it in many cases,
        like ``x**3``, ``x**4``, ``O(x**3)`` etc. In those simple cases, it just parses the
        expression and returns the "sane" ordering such as::

          1 < x < x**2 < x**3 < O(x**4) etc.

        Examples
        ========

        >>> from sympy.abc import x
        >>> from sympy import Basic, Number
        >>> Basic._compare_pretty(x, x**2)
        -1
        >>> Basic._compare_pretty(x**2, x**2)
        0
        >>> Basic._compare_pretty(x**3, x**2)
        1
        >>> Basic._compare_pretty(Number(1, 2), Number(1, 3))
        1
        >>> Basic._compare_pretty(Number(0), Number(-1))
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
        return cmp(a.sort_key(), b.sort_key())

    @classmethod
    def fromiter(cls, args, **assumptions):
        """
        Create a new object from an iterable.

        This is a convenience function that allows one to create objects from
        any iterable, without having to convert to a list or tuple first.

        Examples
        ========

        >>> from sympy import Tuple
        >>> Tuple.fromiter(i for i in xrange(5))
        (0, 1, 2, 3, 4)

        """
        return cls(*tuple(args), **assumptions)

    @classmethod
    def class_key(cls):
        """Nice order of classes. """
        return 5, 0, cls.__name__

    @cacheit
    def sort_key(self, order=None):
        """
        Return a sort key.

        Examples
        ========

        >>> from sympy.core import Basic, S, I
        >>> from sympy.abc import x

        >>> sorted([S(1)/2, I, -I], key=lambda x: x.sort_key())
        [1/2, -I, I]

        >>> S("[x, 1/x, 1/x**2, x**2, x**(1/2), x**(1/4), x**(3/2)]")
        [x, 1/x, x**(-2), x**2, sqrt(x), x**(1/4), x**(3/2)]
        >>> sorted(_, key=lambda x: x.sort_key())
        [x**(-2), 1/x, x**(1/4), sqrt(x), x, x**(3/2), x**2]

        """

        # XXX: remove this when issue #2070 is fixed
        def inner_key(arg):
            if isinstance(arg, Basic):
                return arg.sort_key(order)
            else:
                return arg

        args = self._sorted_args
        args = len(args), tuple([ inner_key(arg) for arg in args ])
        return self.class_key(), args, S.One.sort_key(), S.One

    def __eq__(self, other):
        """Return a boolean indicating whether a == b on the basis of
        their symbolic trees.

        This is the same as a.compare(b) == 0 but faster.

        Notes
        =====

        If a class that overrides __eq__() needs to retain the
        implementation of __hash__() from a parent class, the
        interpreter must be told this explicitly by setting __hash__ =
        <ParentClass>.__hash__. Otherwise the inheritance of __hash__()
        will be blocked, just as if __hash__ had been explicitly set to
        None.

        References
        ==========

        from http://docs.python.org/dev/reference/datamodel.html#object.__hash__
        """

        if type(self) is not type(other):
            # issue 3001 a**1.0 == a like a**2.0 == a**2
            while isinstance(self, C.Pow) and self.exp == 1:
                self = self.base
            while isinstance(other, C.Pow) and other.exp == 1:
                other = other.base
            try:
                other = _sympify(other)
            except SympifyError:
                return False    # sympy != other

            if type(self) is not type(other):
                return False

        return self._hashable_content() == other._hashable_content()


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

        return self._hashable_content() != other._hashable_content()

    def dummy_eq(self, other, symbol=None):
        """
        Compare two expressions and handle dummy symbols.

        Examples
        ========

        >>> from sympy import Dummy
        >>> from sympy.abc import x, y

        >>> u = Dummy('u')

        >>> (u**2 + 1).dummy_eq(x**2 + 1)
        True
        >>> (u**2 + 1) == (x**2 + 1)
        False

        >>> (u**2 + y).dummy_eq(x**2 + y, x)
        True
        >>> (u**2 + y).dummy_eq(x**2 + y, y)
        False

        """
        dummy_symbols = [ s for s in self.free_symbols if s.is_Dummy ]

        if not dummy_symbols:
            return self == other
        elif len(dummy_symbols) == 1:
            dummy = dummy_symbols.pop()
        else:
            raise ValueError("only one dummy symbol allowed on the left-hand side")

        if symbol is None:
            symbols = other.free_symbols

            if not symbols:
                return self == other
            elif len(symbols) == 1:
                symbol = symbols.pop()
            else:
                raise ValueError("specify a symbol in which expressions should be compared")

        tmp = dummy.__class__()

        return self.subs(dummy, tmp) == other.subs(symbol, tmp)

    # Note, we always use the default ordering (lex) in __str__ and __repr__,
    # regardless of the global setting.  See issue 2388.
    def __repr__(self):
        from sympy.printing import sstr
        return sstr(self, order=None)

    def __str__(self):
        from sympy.printing import sstr
        return sstr(self, order=None)

    def atoms(self, *types):
        """Returns the atoms that form the current object.

           By default, only objects that are truly atomic and can't
           be divided into smaller pieces are returned: symbols, numbers,
           and number symbols like I and pi. It is possible to request
           atoms of any type, however, as demonstrated below.

           Examples
           ========

           >>> from sympy import I, pi, sin
           >>> from sympy.abc import x, y
           >>> (1 + x + 2*sin(y + I*pi)).atoms()
           set([1, 2, I, pi, x, y])

           If one or more types are given, the results will contain only
           those types of atoms.

           Examples
           ========

           >>> from sympy import Number, NumberSymbol, Symbol
           >>> (1 + x + 2*sin(y + I*pi)).atoms(Symbol)
           set([x, y])

           >>> (1 + x + 2*sin(y + I*pi)).atoms(Number)
           set([1, 2])

           >>> (1 + x + 2*sin(y + I*pi)).atoms(Number, NumberSymbol)
           set([1, 2, pi])

           >>> (1 + x + 2*sin(y + I*pi)).atoms(Number, NumberSymbol, I)
           set([1, 2, I, pi])

           Note that I (imaginary unit) and zoo (complex infinity) are special
           types of number symbols and are not part of the NumberSymbol class.

           The type can be given implicitly, too:

           >>> (1 + x + 2*sin(y + I*pi)).atoms(x) # x is a Symbol
           set([x, y])

           Be careful to check your assumptions when using the implicit option
           since ``S(1).is_Integer = True`` but ``type(S(1))`` is ``One``, a special type
           of sympy atom, while ``type(S(2))`` is type ``Integer`` and will find all
           integers in an expression:

           >>> from sympy import S
           >>> (1 + x + 2*sin(y + I*pi)).atoms(S(1))
           set([1])

           >>> (1 + x + 2*sin(y + I*pi)).atoms(S(2))
           set([1, 2])

           Finally, arguments to atoms() can select more than atomic atoms: any
           sympy type (loaded in core/__init__.py) can be listed as an argument
           and those types of "atoms" as found in scanning the arguments of the
           expression recursively:

           >>> from sympy import Function, Mul
           >>> from sympy.core.function import AppliedUndef
           >>> f = Function('f')
           >>> (1 + f(x) + 2*sin(y + I*pi)).atoms(Function)
           set([f(x), sin(y + I*pi)])
           >>> (1 + f(x) + 2*sin(y + I*pi)).atoms(AppliedUndef)
           set([f(x)])

           >>> (1 + x + 2*sin(y + I*pi)).atoms(Mul)
           set([I*pi, 2*sin(y + I*pi)])

        """
        if types:
            types = tuple([t if isinstance(t, type) else type(t) for t in types])
        else:
            types = (Atom,)
        result = set()
        for expr in preorder_traversal(self):
            if isinstance(expr, types):
                result.add(expr)
        return result

    @property
    def free_symbols(self):
        """Return from the atoms of self those which are free symbols.

        For most expressions, all symbols are free symbols. For some classes
        this is not true. e.g. Integrals use Symbols for the dummy variables
        which are bound variables, so Integral has a method to return all symbols
        except those. Derivative keeps track of symbols with respect to which it
        will perform a derivative; those are bound variables, too, so it has
        its own symbols method.

        Any other method that uses bound variables should implement a symbols
        method."""
        union = set.union
        return reduce(union, [arg.free_symbols for arg in self.args], set())

    def is_hypergeometric(self, k):
        from sympy.simplify import hypersimp
        return hypersimp(self, k) is not None

    @property
    def is_number(self):
        """Returns ``True`` if 'self' is a number.

           >>> from sympy import log, Integral
           >>> from sympy.abc import x, y

           >>> x.is_number
           False
           >>> (2*x).is_number
           False
           >>> (2 + log(2)).is_number
           True
           >>> (2 + Integral(2, x)).is_number
           False
           >>> (2 + Integral(2, (x, 1, 2))).is_number
           True

        """
        # should be overriden by subclasses
        return False

    @property
    def is_comparable(self):
        is_real = self.is_real
        if is_real is False:
            return False
        is_number = self.is_number
        if is_number is False:
            return False
        if is_real and is_number:
            return True
        n, i = self.evalf(2).as_real_imag()
        if not i.is_Number or not n.is_Number:
            return False
        if i:
            if i._prec != 1:
                return False
        return n._prec != 1

    @property
    def func(self):
        """
        The top-level function in an expression.

        The following should hold for all objects::

            >> x == x.func(*x.args)

        Examples
        ========

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

        Examples
        ========

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

        Notes
        =====

        Never use self._args, always use self.args.
        Only use _args in __new__ when creating a new function.
        Don't override .args() from Basic (so that it's easy to
        change the interface in the future if needed).
        """
        return self._args

    @property
    def _sorted_args(self):
        """
        The same as ``args``.  Derived classes which don't fix an
        order on their arguments should override this method to
        produce the sorted representation.
        """
        return self.args

    def iter_basic_args(self):
        """
        Iterates arguments of ``self``.

        Examples
        ========

        >>> from sympy.abc import x
        >>> a = 2*x
        >>> a.iter_basic_args()
        <...iterator object at 0x...>
        >>> list(a.iter_basic_args())
        [2, x]

        """
        return iter(self.args)

    def as_poly(self, *gens, **args):
        """Converts ``self`` to a polynomial or returns ``None``.

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

    def as_content_primitive(self, radical=False):
        """A stub to allow Basic args (like Tuple) to be skipped when computing
        the content and primitive components of an expression.

        See docstring of Expr.as_content_primitive
        """
        return S.One, self

    def subs(self, *args, **kwargs):
        """
        Substitutes old for new in an expression after sympifying args.

        `args` is either:
          - two arguments, e.g. foo.subs(old, new)
          - one iterable argument, e.g. foo.subs(iterable). The iterable may be
             o an iterable container with (old, new) pairs. In this case the
               replacements are processed in the order given with successive
               patterns possibly affecting replacements already made.
             o a dict or set whose key/value items correspond to old/new pairs.
               In this case the old/new pairs will be sorted by op count and in
               case of a tie, by number of args and the default_sort_key. The
               resulting sorted list is then processed as an iterable container
               (see previous).

        If the keyword ``simultaneous`` is True, the subexpressions will not be
        evaluated until all the substitutions have been made.

        Examples
        ========

        >>> from sympy import pi, exp
        >>> from sympy.abc import x, y
        >>> (1 + x*y).subs(x, pi)
        pi*y + 1
        >>> (1 + x*y).subs({x:pi, y:2})
        1 + 2*pi
        >>> (1 + x*y).subs([(x, pi), (y, 2)])
        1 + 2*pi
        >>> reps = [(y, x**2), (x, 2)]
        >>> (x + y).subs(reps)
        6
        >>> (x + y).subs(reversed(reps))
        x**2 + 2

        >>> (x**2 + x**4).subs(x**2, y)
        y**2 + y

        To replace only the x**2 but not the x**4, use xreplace:

        >>> (x**2 + x**4).xreplace({x**2: y})
        x**4 + y

        To delay evaluation until all substitutions have been made,
        set the keyword ``simultaneous`` to True:

        >>> (x/y).subs([(x, 0), (y, 0)])
        0
        >>> (x/y).subs([(x, 0), (y, 0)], simultaneous=True)
        nan

        This has the added feature of not allowing subsequent substitutions
        to affect those already made:

        >>> ((x + y)/y).subs({x + y: y, y: x + y})
        1
        >>> ((x + y)/y).subs({x + y: y, y: x + y}, simultaneous=True)
        y/(x + y)

        In order to obtain a canonical result, unordered iterables are
        sorted by count_op length, number of arguments and by the
        default_sort_key to break any ties. All other iterables are left
        unsorted.

        >>> from sympy import sqrt, sin, cos, exp
        >>> from sympy.abc import a, b, c, d, e

        >>> A = (sqrt(sin(2*x)), a)
        >>> B = (sin(2*x), b)
        >>> C = (cos(2*x), c)
        >>> D = (x, d)
        >>> E = (exp(x), e)

        >>> expr = sqrt(sin(2*x))*sin(exp(x)*x)*cos(2*x) + sin(2*x)

        >>> expr.subs(dict([A,B,C,D,E]))
        a*c*sin(d*e) + b

        See Also
        ========
        replace: replacement capable of doing wildcard-like matching,
                 parsing of match, and conditional replacements
        xreplace: exact node replacement in expr tree; also capable of
                  using matching rules

        """
        from sympy.core.containers import Dict
        from sympy.utilities import default_sort_key

        unordered = False
        if len(args) == 1:
            sequence = args[0]
            if isinstance(sequence, set):
                unordered = True
            elif isinstance(sequence, (Dict, dict)):
                unordered = True
                sequence = sequence.items()
            elif not iterable(sequence):
                from sympy.utilities.misc import filldedent
                raise ValueError(filldedent("""
                   When a single argument is passed to subs
                   it should be a dictionary of old: new pairs or an iterable
                   of (old, new) tuples."""))
        elif len(args) == 2:
            sequence = [args]
        else:
            raise ValueError("subs accepts either 1 or 2 arguments")

        sequence = list(sequence)
        for i in range(len(sequence)):
            o, n = sequence[i]
            so, sn = sympify(o), sympify(n)
            if not isinstance(so, Basic):
                if type(o) is str:
                    so = C.Symbol(o)
            sequence[i] = (so, sn)
            if _aresame(so, sn):
                sequence[i] = None
                continue
        sequence = filter(None, sequence)

        if unordered:
            sequence = dict(sequence)
            if not all(k.is_Atom for k in sequence):
                d = {}
                for o, n in sequence.iteritems():
                    try:
                        ops = o.count_ops(), len(o.args)
                    except TypeError:
                        ops = (0, 0)
                    d.setdefault(ops, []).append((o, n))
                newseq = []
                for k in sorted(d.keys(), reverse=True):
                    newseq.extend(sorted([v[0] for v in d[k]], key=default_sort_key))
                sequence = [(k, sequence[k]) for k in newseq]
                del newseq, d
            else:
                sequence = sorted([(k, v) for (k, v) in sequence.iteritems()],
                                  key=default_sort_key)

        if kwargs.pop('simultaneous', False): # XXX should this be the default for dict subs?
            reps = {}
            rv = self
            for old, new in sequence:
                d = C.Dummy()
                rv = rv._subs(old, d)
                reps[d] = new
                if not isinstance(rv, Basic):
                    break
            return rv.xreplace(reps)
        else:
            rv = self
            for old, new in sequence:
                rv = rv._subs(old, new)
                if not isinstance(rv, Basic):
                    break
            return rv

    @cacheit
    def _subs(self, old, new, **hints):
        """Substitutes an expression old -> new.

        If self is not equal to old then _eval_subs is called.
        If _eval_subs doesn't want to make any special replacement
        then a None is received which indicates that the fallback
        should be applied wherein a search for replacements is made
        amongst the arguments of self.

        >>> from sympy import Basic, Add, Mul
        >>> from sympy.abc import x, y, z

        Examples
        ========

        Add's _eval_subs knows how to target x + y in the following
        so it makes the change:

            >>> (x + y + z).subs(x + y, 1)
            z + 1

        Add's _eval_subs doesn't need to know how to find x + y in
        the following:

            >>> Add._eval_subs(z*(x + y) + 3, x + y, 1) is None
            True

        The returned None will cause the fallback routine to traverse the args and
        pass the z*(x + y) arg to Mul where the change will take place and the
        substitution will succeed:

            >>> (z*(x + y) + 3).subs(x + y, 1)
            z + 3

        ** Developers Notes **

        An _eval_subs routine for a class should be written if:

            1) any arguments are not instances of Basic (e.g. bool, tuple);

            2) some arguments should not be targeted (as in integration
               variables);

            3) if there is something other than a literal replacement
               that should be attempted (as in Piecewise where the condition
               may be updated without doing a replacement).

        If it is overridden, here are some special cases that might arise:

            1) If it turns out that no special change was made and all
               the original sub-arguments should be checked for
               replacements then None should be returned.

            2) If it is necessary to do substitutions on a portion of
               the expression then _subs should be called. _subs will
               handle the case of any sub-expression being equal to old
               (which usually would not be the case) while its fallback
               will handle the recursion into the sub-arguments. For
               example, after Add's _eval_subs removes some matching terms
               it must process the remaining terms so it calls _subs
               on each of the un-matched terms and then adds them
               onto the terms previously obtained.

           3) If the initial expression should remain unchanged then
              the original expression should be returned. (Whenever an
              expression is returned, modified or not, no further
              substitution of old -> new is attempted.) Sum's _eval_subs
              routine uses this strategy when a substitution is attempted
              on any of its summation variables.
        """

        def fallback(self, old, new):
            """
            Try to replace old with new in any of self's arguments.
            """
            hit = False
            args = list(self.args)
            for i, arg in enumerate(args):
                if not hasattr(arg, '_eval_subs'):
                    continue
                arg = arg._subs(old, new, **hints)
                if arg is not args[i]:
                    hit = True
                    args[i] = arg
            if hit:
                return self.func(*args)
            return self

        if _aresame(self, old):
            return new

        rv = self._eval_subs(old, new)
        if rv is None:
            rv = fallback(self, old, new)
        return rv

    def _eval_subs(self, old, new):
        """Override this stub if you want to do anything more than
        attempt a replacement of old with new in the arguments of self.

        See also: _subs
        """
        return None

    def xreplace(self, rule):
        """
        Replace occurrences of objects within the expression.

        Parameters
        ==========
        rule : dict-like
            Expresses a replacement rule

        Returns
        =======
        xreplace : the result of the replacement

        Examples
        ========
        >>> from sympy import symbols, pi, exp
        >>> x, y, z = symbols('x y z')
        >>> (1 + x*y).xreplace({x: pi})
        pi*y + 1
        >>> (1 + x*y).xreplace({x:pi, y:2})
        1 + 2*pi

        Replacements occur only if an entire node in the expression tree is
        matched:

        >>> (x*y + z).xreplace({x*y: pi})
        z + pi
        >>> (x*y*z).xreplace({x*y: pi})
        x*y*z
        >>> (2*x).xreplace({2*x: y, x: z})
        y
        >>> (2*2*x).xreplace({2*x: y, x: z})
        4*z
        >>> (x + y + 2).xreplace({x + y: 2})
        x + y + 2
        >>> (x + 2 + exp(x + 2)).xreplace({x + 2: y})
        x + exp(y) + 2

        xreplace doesn't differentiate between free and bound symbols. In the
        following, subs(x, y) would not change x since it is a bound symbol,
        but xreplace does:

        >>> from sympy import Integral
        >>> Integral(x, (x, 1, 2*x)).xreplace({x: y})
        Integral(y, (y, 1, 2*y))

        Trying to replace x with an expression raises an error:

        >>> Integral(x, (x, 1, 2*x)).xreplace({x: 2*y}) #doctest: +SKIP
        ValueError: Invalid limits given: ((2*y, 1, 4*y),)

        See Also
        ========
        replace: replacement capable of doing wildcard-like matching,
                 parsing of match, and conditional replacements
        subs: substitution of subexpressions as defined by the objects
              themselves.

        """
        if self in rule:
            return rule[self]
        elif rule:
            args = tuple([arg.xreplace(rule) for arg in self.args])
            if args != self.args:
                return self.func(*args)
        return self

    @deprecated(useinstead="has", issue=2389, deprecated_since_version="0.7.2")
    def __contains__(self, obj):
        if self == obj:
            return True
        for arg in self.args:
            try:
                if obj in arg:
                    return True
            except TypeError:
                if obj == arg:
                    return True
        return False

    @cacheit
    def has(self, *patterns):
        """
        Test whether any subexpression matches any of the patterns.

        Examples
        ========

        >>> from sympy import sin, S
        >>> from sympy.abc import x, y, z
        >>> (x**2 + sin(x*y)).has(z)
        False
        >>> (x**2 + sin(x*y)).has(x, y, z)
        True
        >>> x.has(x)
        True

        Note that ``expr.has(*patterns)`` is exactly equivalent to
        ``any(expr.has(p) for p in patterns)``. In particular, ``False`` is
        returned when the list of patterns is empty.

        >>> x.has()
        False

        """
        return any(self._has(pattern) for pattern in patterns)

    def _has(self, pattern):
        """Helper for .has()"""
        from sympy.core.function import UndefinedFunction, Function
        if isinstance(pattern, UndefinedFunction):
            return any(f.func == pattern or f == pattern
            for f in self.atoms(Function, UndefinedFunction))

        pattern = sympify(pattern)
        if isinstance(pattern, BasicType):
            return any(isinstance(arg, pattern)
            for arg in preorder_traversal(self))

        try:
            match = pattern._has_matcher()
            return any(match(arg) for arg in preorder_traversal(self))
        except AttributeError:
            return any(arg == pattern for arg in preorder_traversal(self))


    def _has_matcher(self):
        """Helper for .has()"""
        return self.__eq__


    def replace(self, query, value, map=False):
        """
        Replace matching subexpressions of ``self`` with ``value``.

        If ``map = True`` then also return the mapping {old: new} where ``old``
        was a sub-expression found with query and ``new`` is the replacement
        value for it.

        Traverses an expression tree and performs replacement of matching
        subexpressions from the bottom to the top of the tree. The list of
        possible combinations of queries and replacement values is listed
        below:

        Examples
        ========

        Initial setup

            >>> from sympy import log, sin, cos, tan, Wild
            >>> from sympy.abc import x, y
            >>> f = log(sin(x)) + tan(sin(x**2))

        1.1. type -> type
            obj.replace(sin, tan)

            >>> f.replace(sin, cos)
            log(cos(x)) + tan(cos(x**2))
            >>> sin(x).replace(sin, cos, map=True)
            (cos(x), {sin(x): cos(x)})

        1.2. type -> func
            obj.replace(sin, lambda arg: ...)

            >>> f.replace(sin, lambda arg: sin(2*arg))
            log(sin(2*x)) + tan(sin(2*x**2))

        2.1. expr -> expr
            obj.replace(sin(a), tan(a))

            >>> a = Wild('a')
            >>> f.replace(sin(a), tan(a))
            log(tan(x)) + tan(tan(x**2))

        2.2. expr -> func
            obj.replace(sin(a), lambda a: ...)

            >>> f.replace(sin(a), cos(a))
            log(cos(x)) + tan(cos(x**2))
            >>> f.replace(sin(a), lambda a: sin(2*a))
            log(sin(2*x)) + tan(sin(2*x**2))

        3.1. func -> func
            obj.replace(lambda expr: ..., lambda expr: ...)

            >>> g = 2*sin(x**3)
            >>> g.replace(lambda expr: expr.is_Number, lambda expr: expr**2)
            4*sin(x**9)

        See Also
        ========
        subs: substitution of subexpressions as defined by the objects
              themselves.
        xreplace: exact node replacement in expr tree; also capable of
                  using matching rules

        """
        if isinstance(query, type):
            _query = lambda expr: isinstance(expr, query)

            if isinstance(value, type):
                _value = lambda expr, result: value(*expr.args)
            elif callable(value):
                _value = lambda expr, result: value(*expr.args)
            else:
                raise TypeError("given a type, replace() expects another type or a callable")
        elif isinstance(query, Basic):
            _query = lambda expr: expr.match(query)

            if isinstance(value, Basic):
                _value = lambda expr, result: value.subs(result)
            elif callable(value):
                _value = lambda expr, result: value(**dict([ (str(key)[:-1], val) for key, val in result.iteritems() ]))
            else:
                raise TypeError("given an expression, replace() expects another expression or a callable")
        elif callable(query):
            _query = query

            if callable(value):
                _value = lambda expr, result: value(expr)
            else:
                raise TypeError("given a callable, replace() expects another callable")
        else:
            raise TypeError("first argument to replace() must be a type, an expression or a callable")

        mapping = {}

        def rec_replace(expr):
            args, construct = [], False

            for arg in expr.args:
                result = rec_replace(arg)

                if result is not None:
                    construct = True
                else:
                    result = arg

                args.append(result)

            if construct:
                return expr.__class__(*args)
            else:
                result = _query(expr)

                if result:
                    value = _value(expr, result)

                    if map:
                        mapping[expr] = value

                    return value
                else:
                    return None

        result = rec_replace(self)

        if result is None:
            result = self

        if not map:
            return result
        else:
            return result, mapping

    def find(self, query, group=False):
        """Find all subexpressions matching a query. """
        query = _make_find_query(query)
        results = filter(query, preorder_traversal(self))

        if not group:
            return set(results)
        else:
            groups = {}

            for result in results:
                if result in groups:
                    groups[result] += 1
                else:
                    groups[result] = 1

            return groups

    def count(self, query):
        """Count the number of matching subexpressions. """
        query = _make_find_query(query)
        return sum(bool(query(sub)) for sub in preorder_traversal(self))

    def matches(self, expr, repl_dict={}):
        """
        Helper method for match() - switches the pattern and expr.

        Can be used to solve linear equations:

        >>> from sympy import Symbol, Wild, Integer
        >>> a,b = map(Symbol, 'ab')
        >>> x = Wild('x')
        >>> (a+b*x).matches(Integer(0))
        {x_: -a/b}

        """
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
            d = arg.xreplace(d).matches(other_arg, d)
            if d is None:
                return None
        return d

    def match(self, pattern):
        """
        Pattern matching.

        Wild symbols match all.

        Return ``None`` when expression (self) does not match
        with pattern. Otherwise return a dictionary such that::

          pattern.xreplace(self.match(pattern)) == self

        Examples
        ========

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
        >>> (p*q**r).xreplace(e.match(p*q**r))
        4*x**2

        """
        pattern = sympify(pattern)
        return pattern.matches(self)

    def count_ops(self, visual=None):
        """wrapper for count_ops that returns the operation count."""
        from sympy import count_ops
        return count_ops(self, visual)
        return sum(a.count_ops(visual) for a in self.args)

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
            return self.func(*terms)
        else:
            return self

    def _eval_rewrite(self, pattern, rule, **hints):
        if self.is_Atom:
            return self
        sargs = self.args
        terms = [ t._eval_rewrite(pattern, rule, **hints) for t in sargs ]
        return self.func(*terms)

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
            if isinstance(args[-1], basestring):
                rule = '_eval_rewrite_as_' + args[-1]
            else:
                rule = '_eval_rewrite_as_' + args[-1].__name__

            if not pattern:
                return self._eval_rewrite(None, rule, **hints)
            else:
                if iterable(pattern[0]):
                    pattern = pattern[0]

                pattern = [ p.__class__ for p in pattern if self.has(p) ]

                if pattern:
                    return self._eval_rewrite(tuple(pattern), rule, **hints)
                else:
                    return self

class Atom(Basic):
    """
    A parent class for atomic things. An atom is an expression with no subexpressions.

    Examples
    ========

    Symbol, Number, Rational, Integer, ...
    But not: Add, Mul, Pow, ...
    """

    is_Atom = True

    __slots__ = []

    def matches(self, expr, repl_dict={}):
        if self == expr:
            return repl_dict

    def xreplace(self, rule):
        return rule.get(self, self)

    def doit(self, **hints):
        return self

    @classmethod
    def class_key(cls):
        return 2, 0, cls.__name__

    @cacheit
    def sort_key(self, order=None):
        from sympy.core import S
        return self.class_key(), (1, (str(self),)), S.One.sort_key(), S.One

    @property
    def _sorted_args(self):
        # this is here as a safeguard against accidentally using _sorted_args
        # on Atoms -- they cannot be rebuilt as atom.func(*atom._sorted_args)
        # since there are no args. So the calling routine should be checking
        # to see that this property is not called for Atoms.
        raise AttributeError('Atoms have no args. It might be necessary'
        ' to make a check for Atoms in the calling code.')

def _aresame(a, b):
    """Return True if a and b are structurally the same, else False.

    Examples
    ========

    To SymPy, 2.0 == 2:

    >>> from sympy import S, Symbol, cos, sin
    >>> 2.0 == S(2)
    True

    Since a simple 'same or not' result is sometimes useful, this routine was
    written to provide that query:

    >>> from sympy.core.basic import _aresame
    >>> _aresame(S(2.0), S(2))
    False

    """
    from itertools import izip

    for i, j in izip(preorder_traversal(a), preorder_traversal(b)):
        if i != j or type(i) != type(j):
            return False
    else:
        return True

def _atomic(e):
    """Return atom-like quantities as far as substitution is
    concerned: Derivatives, Functions and Symbols. Don't
    return any 'atoms' that are inside such quantities unless
    they also appear outside, too.

    Examples
    ========
    >>> from sympy import Derivative, Function, cos
    >>> from sympy.abc import x, y
    >>> from sympy.core.basic import _atomic
    >>> f = Function('f')
    >>> _atomic(x + y)
    set([x, y])
    >>> _atomic(x + f(y))
    set([x, f(y)])
    >>> _atomic(Derivative(f(x), x) + cos(x) + y)
    set([y, cos(x), Derivative(f(x), x)])

    """
    from sympy import Derivative, Function, Symbol
    pot = preorder_traversal(e)
    seen = set()
    try:
        free = e.free_symbols
    except AttributeError:
        return set([e])
    atoms = set()
    for p in pot:
        if p in seen:
            pot.skip()
            continue
        seen.add(p)
        if isinstance(p, Symbol) and p in free:
            atoms.add(p)
        elif isinstance(p, (Derivative, Function)):
            pot.skip()
            atoms.add(p)
    return atoms

class preorder_traversal(object):
    """
    Do a pre-order traversal of a tree.

    This iterator recursively yields nodes that it has visited in a pre-order
    fashion. That is, it yields the current node then descends through the
    tree breadth-first to yield all of a node's children's pre-order
    traversal.


    For an expression, the order of the traversal depends on the order of
    .args, which in many cases can be arbitrary.

    Parameters
    ==========
    node : sympy expression
        The expression to traverse.
    key : (default None) sort key
        The key used to sort args of Basic objects. When None, args of Basic
        objects are processed in arbitrary order.

    Yields
    ======
    subtree : sympy expression
        All of the subtrees in the tree.

    Examples
    ========
    >>> from sympy import symbols
    >>> from sympy import symbols, default_sort_key
    >>> from sympy.core.basic import preorder_traversal
    >>> x, y, z = symbols('x y z')

    The nodes are returned in the order that they are encountered unless key
    is given.

    >>> list(preorder_traversal((x + y)*z, key=None)) # doctest: +SKIP
    [z*(x + y), z, x + y, y, x]
    >>> list(preorder_traversal((x + y)*z, key=default_sort_key))
    [z*(x + y), z, x + y, x, y]

    """
    def __init__(self, node, key=None):
        self._skip_flag = False
        self._pt = self._preorder_traversal(node, key)

    def _preorder_traversal(self, node, key):
        yield node
        if self._skip_flag:
            self._skip_flag = False
            return
        if isinstance(node, Basic):
            args = node.args
            if key:
                args = list(args)
                args.sort(key=key)
            for arg in args:
                for subtree in self._preorder_traversal(arg, key):
                    yield subtree
        elif iterable(node):
            for item in node:
                for subtree in self._preorder_traversal(item, key):
                    yield subtree


    def skip(self):
        """
        Skip yielding current node's (last yielded node's) subtrees.

        Examples
        --------
        >>> from sympy.core import symbols
        >>> from sympy.core.basic import preorder_traversal
        >>> x, y, z = symbols('x y z')
        >>> pt = preorder_traversal((x+y*z)*z)
        >>> for i in pt:
        ...     print i
        ...     if i == x+y*z:
        ...             pt.skip()
        z*(x + y*z)
        z
        x + y*z
        """
        self._skip_flag = True

    def next(self):
        return self._pt.next()

    def __iter__(self):
        return self

def _make_find_query(query):
    """Convert the argument of Basic.find() into a callable"""
    try:
        query = sympify(query)
    except SympifyError:
        pass
    if isinstance(query, type):
        return lambda expr: isinstance(expr, query)
    elif isinstance(query, Basic):
        return lambda expr: expr.match(query) is not None
    return query
