from __future__ import print_function, division

from sympy.core import S, C, sympify
from sympy.core.add import Add
from sympy.core.basic import Basic
from sympy.core.containers import Tuple
from sympy.core.numbers import Rational
from sympy.core.operations import LatticeOp, ShortCircuit
from sympy.core.function import Application, Lambda, ArgumentIndexError
from sympy.core.expr import Expr
from sympy.core.singleton import Singleton
from sympy.core.rules import Transform
from sympy.core.compatibility import as_int, with_metaclass, xrange
from sympy.core.logic import fuzzy_and


class IdentityFunction(with_metaclass(Singleton, Lambda)):
    """
    The identity function

    Examples
    ========

    >>> from sympy import Id, Symbol
    >>> x = Symbol('x')
    >>> Id(x)
    x

    """

    def __new__(cls):
        from sympy.sets.sets import FiniteSet
        x = C.Dummy('x')
        #construct "by hand" to avoid infinite loop
        obj = Expr.__new__(cls, Tuple(x), x)
        obj.nargs = FiniteSet(1)
        return obj

Id = S.IdentityFunction

###############################################################################
############################# ROOT and SQUARE ROOT FUNCTION ###################
###############################################################################


def sqrt(arg):
    """The square root function

    sqrt(x) -> Returns the principal square root of x.

    Examples
    ========

    >>> from sympy import sqrt, Symbol
    >>> x = Symbol('x')

    >>> sqrt(x)
    sqrt(x)

    >>> sqrt(x)**2
    x

    Note that sqrt(x**2) does not simplify to x.

    >>> sqrt(x**2)
    sqrt(x**2)

    This is because the two are not equal to each other in general.
    For example, consider x == -1:

    >>> from sympy import Eq
    >>> Eq(sqrt(x**2), x).subs(x, -1)
    False

    This is because sqrt computes the principal square root, so the square may
    put the argument in a different branch.  This identity does hold if x is
    positive:

    >>> y = Symbol('y', positive=True)
    >>> sqrt(y**2)
    y

    You can force this simplification by using the powdenest() function with
    the force option set to True:

    >>> from sympy import powdenest
    >>> sqrt(x**2)
    sqrt(x**2)
    >>> powdenest(sqrt(x**2), force=True)
    x

    To get both branches of the square root you can use the RootOf function:

    >>> from sympy import RootOf

    >>> [ RootOf(x**2-3,i) for i in (0,1) ]
    [-sqrt(3), sqrt(3)]

    See Also
    ========

    sympy.polys.rootoftools.RootOf, root, real_root

    References
    ==========

    * http://en.wikipedia.org/wiki/Square_root
    * http://en.wikipedia.org/wiki/Principal_value

    """
    # arg = sympify(arg) is handled by Pow
    return C.Pow(arg, S.Half)



def cbrt(arg):
    """This function computes the principial cube root of `arg`, so
    it's just a shortcut for `arg**Rational(1, 3)`.

    Examples
    ========

    >>> from sympy import cbrt, Symbol
    >>> x = Symbol('x')

    >>> cbrt(x)
    x**(1/3)

    >>> cbrt(x)**3
    x

    Note that cbrt(x**3) does not simplify to x.

    >>> cbrt(x**3)
    (x**3)**(1/3)

    This is because the two are not equal to each other in general.
    For example, consider `x == -1`:

    >>> from sympy import Eq
    >>> Eq(cbrt(x**3), x).subs(x, -1)
    False

    This is because cbrt computes the principal cube root, this
    identity does hold if `x` is positive:

    >>> y = Symbol('y', positive=True)
    >>> cbrt(y**3)
    y

    See Also
    ========

    sympy.polys.rootoftools.RootOf, root, real_root

    References
    ==========

    * http://en.wikipedia.org/wiki/Cube_root
    * http://en.wikipedia.org/wiki/Principal_value

    """
    return C.Pow(arg, C.Rational(1, 3))


def root(arg, n):
    """The n-th root function (a shortcut for ``arg**(1/n)``)

    root(x, n) -> Returns the principal n-th root of x.


    Examples
    ========

    >>> from sympy import root, Rational
    >>> from sympy.abc import x, n

    >>> root(x, 2)
    sqrt(x)

    >>> root(x, 3)
    x**(1/3)

    >>> root(x, n)
    x**(1/n)

    >>> root(x, -Rational(2, 3))
    x**(-3/2)


    To get all n n-th roots you can use the RootOf function.
    The following examples show the roots of unity for n
    equal 2, 3 and 4:

    >>> from sympy import RootOf, I

    >>> [ RootOf(x**2-1,i) for i in (0,1) ]
    [-1, 1]

    >>> [ RootOf(x**3-1,i) for i in (0,1,2) ]
    [1, -1/2 - sqrt(3)*I/2, -1/2 + sqrt(3)*I/2]

    >>> [ RootOf(x**4-1,i) for i in (0,1,2,3) ]
    [-1, 1, -I, I]

    SymPy, like other symbolic algebra systems, returns the
    complex root of negative numbers. This is the principal
    root and differs from the text-book result that one might
    be expecting. For example, the cube root of -8 does not
    come back as -2:

    >>> root(-8, 3)
    2*(-1)**(1/3)

    The real_root function can be used to either make such a result
    real or simply return the real root in the first place:

    >>> from sympy import real_root
    >>> real_root(_)
    -2
    >>> real_root(-32, 5)
    -2

    See Also
    ========

    sympy.polys.rootoftools.RootOf
    sympy.core.power.integer_nthroot
    sqrt, real_root

    References
    ==========

    * http://en.wikipedia.org/wiki/Square_root
    * http://en.wikipedia.org/wiki/Real_root
    * http://en.wikipedia.org/wiki/Root_of_unity
    * http://en.wikipedia.org/wiki/Principal_value
    * http://mathworld.wolfram.com/CubeRoot.html

    """
    n = sympify(n)
    return C.Pow(arg, 1/n)


def real_root(arg, n=None):
    """Return the real nth-root of arg if possible. If n is omitted then
    all instances of -1**(1/odd) will be changed to -1.

    Examples
    ========

    >>> from sympy import root, real_root, Rational
    >>> from sympy.abc import x, n

    >>> real_root(-8, 3)
    -2
    >>> root(-8, 3)
    2*(-1)**(1/3)
    >>> real_root(_)
    -2

    See Also
    ========

    sympy.polys.rootoftools.RootOf
    sympy.core.power.integer_nthroot
    root, sqrt
    """
    if n is not None:
        n = as_int(n)
        rv = C.Pow(arg, Rational(1, n))
        if n % 2 == 0:
            return rv
    else:
        rv = sympify(arg)
    n1pow = Transform(lambda x: S.NegativeOne,
                      lambda x:
                      x.is_Pow and
                      x.base is S.NegativeOne and
                      x.exp.is_Rational and
                      x.exp.p == 1 and x.exp.q % 2)
    return rv.xreplace(n1pow)

###############################################################################
############################# MINIMUM and MAXIMUM #############################
###############################################################################


class MinMaxBase(Expr, LatticeOp):
    def __new__(cls, *args, **assumptions):
        if not args:
            raise ValueError("The Max/Min functions must have arguments.")

        args = (sympify(arg) for arg in args)

        # first standard filter, for cls.zero and cls.identity
        # also reshape Max(a, Max(b, c)) to Max(a, b, c)
        try:
            _args = frozenset(cls._new_args_filter(args))
        except ShortCircuit:
            return cls.zero

        # second filter
        # variant I: remove ones which can be removed
        # args = cls._collapse_arguments(set(_args), **assumptions)

        # variant II: find local zeros
        args = cls._find_localzeros(set(_args), **assumptions)

        if not args:
            return cls.identity
        elif len(args) == 1:
            return args.pop()
        else:
            # base creation
            # XXX should _args be made canonical with sorting?
            _args = frozenset(args)
            obj = Expr.__new__(cls, _args, **assumptions)
            obj._argset = _args
            return obj

    @classmethod
    def _new_args_filter(cls, arg_sequence):
        """
        Generator filtering args.

        first standard filter, for cls.zero and cls.identity.
        Also reshape Max(a, Max(b, c)) to Max(a, b, c),
        and check arguments for comparability
        """
        for arg in arg_sequence:

            # pre-filter, checking comparability of arguments
            if (arg.is_real is False) or (arg is S.ComplexInfinity):
                raise ValueError("The argument '%s' is not comparable." % arg)

            if arg == cls.zero:
                raise ShortCircuit(arg)
            elif arg == cls.identity:
                continue
            elif arg.func == cls:
                for x in arg.iter_basic_args():
                    yield x
            else:
                yield arg

    @classmethod
    def _find_localzeros(cls, values, **options):
        """
        Sequentially allocate values to localzeros.

        When a value is identified as being more extreme than another member it
        replaces that member; if this is never true, then the value is simply
        appended to the localzeros.
        """
        localzeros = set()
        for v in values:
            is_newzero = True
            localzeros_ = list(localzeros)
            for z in localzeros_:
                if id(v) == id(z):
                    is_newzero = False
                elif cls._is_connected(v, z):
                    is_newzero = False
                    if cls._is_asneeded(v, z):
                        localzeros.remove(z)
                        localzeros.update([v])
            if is_newzero:
                localzeros.update([v])
        return localzeros

    @classmethod
    def _is_connected(cls, x, y):
        """
        Check if x and y are connected somehow.
        """
        xy = x > y
        yx = x < y
        if (x == y) or xy == True or xy == False or yx == True or yx == False:
            return True
        if x.is_Number and y.is_Number:
            return True
        return False

    @classmethod
    def _is_asneeded(cls, x, y):
        """
        Check if x and y satisfy relation condition.

        The relation condition for Max function is x > y,
        for Min function is x < y. They are defined in children Max and Min
        classes through the method _rel(cls, x, y)
        """
        if (x == y):
            return False
        if x.is_Number and y.is_Number:
            if cls._rel(x, y):
                return True
        xy = cls._rel(x, y)
        if xy == True or xy == False:
            return bool(xy)
        yx = cls._rel_inversed(x, y)
        if yx == True or yx == False:
            return not bool(yx)
        return False

    def _eval_derivative(self, s):
        # f(x).diff(s) -> x.diff(s) * f.fdiff(1)(s)
        i = 0
        l = []
        for a in self.args:
            i += 1
            da = a.diff(s)
            if da is S.Zero:
                continue
            try:
                df = self.fdiff(i)
            except ArgumentIndexError:
                df = Function.fdiff(self, i)
            l.append(df * da)
        return Add(*l)

    @property
    def is_real(self):
        return fuzzy_and(arg.is_real for arg in self.args)

class Max(MinMaxBase, Application):
    """
    Return, if possible, the maximum value of the list.

    When number of arguments is equal one, then
    return this argument.

    When number of arguments is equal two, then
    return, if possible, the value from (a, b) that is >= the other.

    In common case, when the length of list greater than 2, the task
    is more complicated. Return only the arguments, which are greater
    than others, if it is possible to determine directional relation.

    If is not possible to determine such a relation, return a partially
    evaluated result.

    Assumptions are used to make the decision too.

    Also, only comparable arguments are permitted.

    Examples
    ========

    >>> from sympy import Max, Symbol, oo
    >>> from sympy.abc import x, y
    >>> p = Symbol('p', positive=True)
    >>> n = Symbol('n', negative=True)

    >>> Max(x, -2)                  #doctest: +SKIP
    Max(x, -2)

    >>> Max(x, -2).subs(x, 3)
    3

    >>> Max(p, -2)
    p

    >>> Max(x, y)                   #doctest: +SKIP
    Max(x, y)

    >>> Max(x, y) == Max(y, x)
    True

    >>> Max(x, Max(y, z))           #doctest: +SKIP
    Max(x, y, z)

    >>> Max(n, 8, p, 7, -oo)        #doctest: +SKIP
    Max(8, p)

    >>> Max (1, x, oo)
    oo

    Algorithm

    The task can be considered as searching of supremums in the
    directed complete partial orders [1]_.

    The source values are sequentially allocated by the isolated subsets
    in which supremums are searched and result as Max arguments.

    If the resulted supremum is single, then it is returned.

    The isolated subsets are the sets of values which are only the comparable
    with each other in the current set. E.g. natural numbers are comparable with
    each other, but not comparable with the `x` symbol. Another example: the
    symbol `x` with negative assumption is comparable with a natural number.

    Also there are "least" elements, which are comparable with all others,
    and have a zero property (maximum or minimum for all elements). E.g. `oo`.
    In case of it the allocation operation is terminated and only this value is
    returned.

    Assumption:
       - if A > B > C then A > C
       - if A==B then B can be removed

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Directed_complete_partial_order
    .. [2] http://en.wikipedia.org/wiki/Lattice_%28order%29

    See Also
    ========

    Min : find minimum values
    """
    zero = S.Infinity
    identity = S.NegativeInfinity

    @classmethod
    def _rel(cls, x, y):
        """
        Check if x > y.
        """
        return (x > y)

    @classmethod
    def _rel_inversed(cls, x, y):
        """
        Check if x < y.
        """
        return (x < y)

    def fdiff( self, argindex ):
        from sympy.functions.special.delta_functions import Heaviside
        n = len(self.args)
        if 0 < argindex and argindex <= n:
            argindex -= 1
            if n == 2:
                return Heaviside( self.args[argindex] - self.args[1-argindex] )
            newargs = tuple([self.args[i] for i in xrange(n) if i != argindex])
            return Heaviside( self.args[argindex] - Max(*newargs) )
        else:
            raise ArgumentIndexError(self, argindex)


class Min(MinMaxBase, Application):
    """
    Return, if possible, the minimum value of the list.

    Examples
    ========

    >>> from sympy import Min, Symbol, oo
    >>> from sympy.abc import x, y
    >>> p = Symbol('p', positive=True)
    >>> n = Symbol('n', negative=True)

    >>> Min(x, -2)                  #doctest: +SKIP
    Min(x, -2)

    >>> Min(x, -2).subs(x, 3)
    -2

    >>> Min(p, -3)
    -3

    >>> Min(x, y)                   #doctest: +SKIP
    Min(x, y)

    >>> Min(n, 8, p, -7, p, oo)     #doctest: +SKIP
    Min(n, -7)

    See Also
    ========

    Max : find maximum values
    """
    zero = S.NegativeInfinity
    identity = S.Infinity

    @classmethod
    def _rel(cls, x, y):
        """
        Check if x < y.
        """
        return (x < y)

    @classmethod
    def _rel_inversed(cls, x, y):
        """
        Check if x > y.
        """
        return (x > y)

    def fdiff( self, argindex ):
        from sympy.functions.special.delta_functions import Heaviside
        n = len(self.args)
        if 0 < argindex and argindex <= n:
            argindex -= 1
            if n == 2:
                return Heaviside( self.args[1-argindex] - self.args[argindex] )
            newargs = tuple([ self.args[i] for i in xrange(n) if i != argindex])
            return Heaviside( Min(*newargs) - self.args[argindex] )
        else:
            raise ArgumentIndexError(self, argindex)
