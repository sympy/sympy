"""Hypergeometric and Meijer G-functions"""

from sympy import S
from sympy.core.function import Function, ArgumentIndexError
from sympy.core.containers import Tuple
from sympy.core.sympify import sympify
from sympy.core.mul import Mul

# TODO should __new__ accept **options?
# TODO should constructors should check if parameters are sensible?

# TODO when pull request #399 is in, this should be no longer necessary
def _make_tuple(v):
    """
    Turn an iterable argument V into a Tuple.

    Examples:
    >>> from sympy.functions.special.hyper import _make_tuple as mt
    >>> from sympy.core.containers import Tuple
    >>> mt([1, 2, 3])
    (1, 2, 3)
    >>> mt((4, 5))
    (4, 5)
    >>> mt((7, 8, 9))
    (7, 8, 9)
    """
    return Tuple(*[sympify(x) for x in v])

class TupleParametersBase(Function):
    """ Base class that takes care of differentiation, when some of
        the arguments are actually tuples. """
    def _eval_derivative(self, s):
        if self.args[0].has(s) or self.args[1].has(s):
            raise NotImplementedError('differentiation with respect to ' \
                                      'a parameter')
        return self.fdiff(3)*self.args[2].diff(s)

class hyper(TupleParametersBase):
    r"""
    The (generalized) hypergeometric function is defined by a series where
    the ratios of successive terms are a rational function of the summation
    index. When convergent, it is continued analytically to the largest
    possible domain.

    The hypergeometric function depends on two vectors of parameters, called
    the numerator parameters :math:`a_p`, and the denominator parameters
    :math:`b_q`. It also has an argument :math:`z`. The series definition is

    .. math ::
        {}_pF_q\left.\left(\begin{matrix} a_1, \dots, a_p \\ b_1, \dots, b_q \end{matrix}
                     \right| z \right)
        = \sum_{n=0}^\infty \frac{(a_1)_n \dots (a_p)_n}{(b_1)_n \dots (b_q)_n}
                            \frac{z^n}{n!},

    where :math:`(a)_n = (a)(a+1)\dots(a+n-1)` denotes the rising factorial.

    If one of the :math:`b_q` is a non-positive integer then the series is
    undefined unless one of the `a_p` is a larger (i.e. smaller in magnitude)
    non-positive integer. If none
    of the :math:`b_q` is a non-positive integer and one of the :math:`a_p` is
    a non-positive integer, then the series reduces to a polynomial. To simplify
    the following discussion, we assume that none of the :math:`a_p` or
    :math:`b_q` is a non-positive integer. For more details, see the references.

    The series converges for all :math:`z` if :math:`p \le q`, and thus
    defines an entire single-valued function in this case. If
    :math:`p = q+1` the series converges for :math:`|z| < 1`, and can be
    continued analytically into a half-plane. If :math:`p > q+1` the series is
    divergent for all :math:`z`.

    Note: The hypergeometric function constructor currently does *not* check if the
    parameters actually yield a well-defined function.


    **Examples**

    The parameters :math:`a_p` and :math:`b_q` can be passed as arbitrary
    iterables, for example:

    >>> from sympy.functions import hyper
    >>> from sympy.abc import x, n, a
    >>> hyper((1, 2, 3), [3, 4], x)
    hyper((1, 2, 3), (3, 4), x)

    There is also pretty printing (it looks better using unicode):

    >>> from sympy import pprint
    >>> pprint(hyper((1, 2, 3), [3, 4], x), use_unicode=False)
      _
     |_  /1, 2, 3 |  \
     |   |        | x|
    3  2 \  3, 4  |  /

    The parameters must always be iterables, even if they are vectors of
    length one or zero:

    >>> hyper((1, ), [], x)
    hyper((1,), (), x)

    But of course they may be variables (but if they depend on x then you
    should not expect much implemented functionality):

    >>> hyper((n, a), (n**2,), x)
    hyper((n, a), (n**2,), x)


    The hypergeometric function generalises many named special functions.
    The function hyperexpand() tries to express a hypergeometric function
    using named special functions.
    For example:

    >>> from sympy import hyperexpand
    >>> hyperexpand(hyper([], [], x))
    exp(x)

    You can also use expand_func:

    >>> from sympy import expand_func
    >>> expand_func(x*hyper([1, 1], [2], -x))
    log(x + 1)

    More examples:

    >>> from sympy import S
    >>> hyperexpand(hyper([], [S(1)/2], -x**2/4))
    cos(x)
    >>> hyperexpand(x*hyper([S(1)/2, S(1)/2], [S(3)/2], x**2))
    asin(x)

    We can also sometimes hyperexpand parametric functions:

    >>> from sympy.abc import a
    >>> hyperexpand(hyper([-a], [], x))
    (-x + 1)**a

    See Also:

    - :func:`sympy.simplify.hyperexpand`

    **References**

    - Luke, Y. L. (1969), The Special Functions and Their Approximations,
      Volume 1
    - http://en.wikipedia.org/wiki/Generalized_hypergeometric_function
    """

    nargs = 3

    def __new__(cls, ap, bq, z):
        # TODO should we check convergence conditions?
        return Function.__new__(cls, _make_tuple(ap), _make_tuple(bq), z)

    def fdiff(self, argindex=3):
        if argindex != 3:
            raise ArgumentIndexError(self, argindex)
        nap = Tuple(*[a + 1 for a in self.ap])
        nbq = Tuple(*[b + 1 for b in self.bq])
        fac = Mul(*self.ap)/Mul(*self.bq)
        return fac*hyper(nap, nbq, self.argument)

    def _eval_expand_func(self, deep=True, **hints):
        from sympy import gamma, hyperexpand
        if len(self.ap) == 2 and len(self.bq) == 1 and self.argument == 1:
            a, b = self.ap
            c    = self.bq[0]
            return gamma(c)*gamma(c - a - b)/gamma(c - a)/gamma(c - b)
        return hyperexpand(self)

    @property
    def argument(self):
        """ Argument of the hypergeometric function. """
        return self.args[2]

    @property
    def ap(self):
        """ Numerator parameters of the hypergeometric function. """
        return self.args[0]

    @property
    def bq(self):
        """ Denominator parameters of the hypergeometric function. """
        return self.args[1]

    @property
    def eta(self):
        """ A quantity related to the convergence of the series. """
        return sum(self.ap) - sum(self.bq)

    @property
    def radius_of_convergence(self):
        """
        Compute the radius of convergence of the defining series.

        Note that even if this is not oo, the function may still be evaluated
        outside of the radius of convergence by analytic continuation. But if
        this is zero, then the function is not actually defined anywhere else.

        >>> from sympy.functions import hyper
        >>> from sympy.abc import z
        >>> hyper((1, 2), [3], z).radius_of_convergence
        1
        >>> hyper((1, 2, 3), [4], z).radius_of_convergence
        0
        >>> hyper((1, 2), (3, 4), z).radius_of_convergence
        oo
        """
        from sympy import oo
        if any(a.is_integer and a <= 0 for a in self.ap + self.bq):
            aints = [a for a in self.ap if a.is_Integer and a <= 0]
            bints = [a for a in self.bq if a.is_Integer and a <= 0]
            if len(aints) < len(bints):
                return S(0)
            popped = False
            for b in bints:
                cancelled = False
                while aints:
                    a = aints.pop()
                    if a >= b:
                        cancelled = True
                        break
                    popped = True
                if not cancelled:
                    return S(0)
            if aints or popped:
                # There are still non-positive numerator parameters.
                # This is a polynomial.
                return oo
        if len(self.ap) == len(self.bq) + 1:
            return S(1)
        elif len(self.ap) <= len(self.bq):
            return oo
        else:
            return S(0)

    @property
    def convergence_statement(self):
        """ Return a condition on z under which the series converges. """
        from sympy import And, Or, re, Ne, oo
        R = self.radius_of_convergence
        if R == 0:
            return False
        if R == oo:
            return True
        # The special functions and their approximations, page 44
        e = self.eta
        z = self.argument
        c1 = And(re(e) < 0, abs(z) <= 1)
        c2 = And(0 <= re(e), re(e) < 1, abs(z) <= 1, Ne(z, 1))
        c3 = And(re(e) >= 1, abs(z) < 1)
        return Or(c1, c2, c3)

class meijerg(TupleParametersBase):
    r"""
    The Meijer G-function is defined by a Mellin-Barnes type integral that
    resembles an inverse Mellin transform. It generalises the hypergeometric
    functions.

    The Meijer G-function depends on four sets of parameters. There are
    "*numerator parameters*"
    :math:`a_1, \dots, a_n` and :math:`a_{n+1}, \dots, a_p`, and there are
    "*denominator parameters*"
    :math:`b_1, \dots, b_m` and :math:`b_{m+1}, \dots, b_q`.
    Confusingly, it is traditionally denoted as follows (note the position of
    `m`, `n`, `p`, `q`, and how they relate to the lengths of the four parameter
    vectors):

    .. math ::
        G_{p,q}^{m,n} \left.\left(\begin{matrix}a_1, \dots, a_n & a_{n+1}, \dots, a_p \\
                                        b_1, \dots, b_m & b_{m+1}, \dots, b_q
                          \end{matrix} \right| z \right).

    However, in sympy the four parameter vectors are always available
    separately (see examples), so that there is no need to keep track of the
    decorating sub- and super-scripts on the G symbol.

    The G function is defined as the following integral:

    .. math ::
         \frac{1}{2 \pi i} \int_L \frac{\prod_{j=1}^m \Gamma(b_j - s)
         \prod_{j=1}^n \Gamma(1 - a_j + s)}{\prod_{j=m+1}^q \Gamma(1- b_j +s)
         \prod_{j=n+1}^p \Gamma(a_j - s)} z^s \mathrm{d}s,

    where :math:`\Gamma(z)` is the gamma function. There are three possible
    contours which we will not describe in detail here (see the references).
    If the integral converges along more than one of them the definitions
    agree. The contours all separate the poles of :math:`\Gamma(1-a_j+s)`
    from the poles of :math:`\Gamma(b_k-s)`, so in particular the G function
    is undefined if :math:`a_j - b_k \in \mathbb{Z}_{>0}` for some
    :math:`j \le n` and :math:`k \le m`.

    The conditions under which one of the contours yields a convergent integral
    are complicated and we do not state them here, see the references.

    Note: Currently the Meijer G-function constructor does *not* check any
    convergence conditions.


    **Examples**

    You can pass the parameters either as four separate vectors:

    >>> from sympy.functions import meijerg
    >>> from sympy.abc import x, a
    >>> from sympy.core.containers import Tuple
    >>> from sympy import pprint
    >>> pprint(meijerg((1, 2), (a, 4), (5,), [], x), use_unicode=False)
     __1, 2 /1, 2  a, 4 |  \
    /__     |           | x|
    \_|4, 1 \ 5         |  /

    or as two nested vectors:

    >>> pprint(meijerg([(1, 2), (3, 4)], ([5], Tuple()), x), use_unicode=False)
     __1, 2 /1, 2  3, 4 |  \
    /__     |           | x|
    \_|4, 1 \ 5         |  /

    As with the hypergeometric function, the parameters may be passed as
    arbitrary iterables. Vectors of length zero and one also have to be
    passed as iterables. The parameters need not be constants, but if they
    depend on the argument then not much implemented functionality should be
    expected.

    All the subvectors of parameters are available:

    >>> from sympy import pprint
    >>> g = meijerg([1], [2], [3], [4], x)
    >>> pprint(g, use_unicode=False)
     __1, 1 /1  2 |  \
    /__     |     | x|
    \_|2, 2 \3  4 |  /
    >>> g.an
    (1,)
    >>> g.ap
    (1, 2)
    >>> g.aother
    (2,)
    >>> g.bm
    (3,)
    >>> g.bq
    (3, 4)
    >>> g.bother
    (4,)


    The Meijer G-function generalises the hypergeometric functions.
    In some cases it can be expressed in terms of hypergeometric functions,
    using Slater's theorem. For example:

    >>> from sympy import hyperexpand
    >>> from sympy.abc import a, b, c
    >>> hyperexpand(meijerg([a], [], [c], [b], x), allow_hyper=True)
    x**c*gamma(-a + c + 1)*hyper((-a + c + 1,), (-b + c + 1,), -x)/gamma(-b + c + 1)

    Thus the Meijer G-function also subsumes many named functions as special
    cases. You can use expand_func or hyperexpand to (try to) rewrite a
    Meijer G-function in terms of named special functions. For example:

    >>> from sympy import expand_func, S
    >>> expand_func(meijerg([[],[]], [[0],[]], -x))
    exp(x)
    >>> hyperexpand(meijerg([[],[]], [[S(1)/2],[0]], (x/2)**2))
    sin(x)/pi**(1/2)

    See Also:

    - :func:`sympy.simplify.hyperexpand`

    **References**

    - Luke, Y. L. (1969), The Special Functions and Their Approximations,
      Volume 1
    - http://en.wikipedia.org/wiki/Meijer_G-function
    """

    nargs = 3

    def __new__(cls, *args):
        if len(args) == 5:
            args = [(args[0], args[1]), (args[2], args[3]), args[4]]
        if len(args) != 3:
            raise TypeError("args must eiter be as, as', bs, bs', z or " \
                            "as, bs, z")
        def tr(p):
            if len(p) != 2:
                raise TypeError("wrong argument")
            return Tuple(_make_tuple(p[0]), _make_tuple(p[1]))

        # TODO should we check convergence conditions?
        return Function.__new__(cls, tr(args[0]), tr(args[1]), args[2])

    def fdiff(self, argindex=3):
        if argindex != 3:
            raise ArgumentIndexError(self, argindex)
        if len(self.an) >= 1:
            a = list(self.an)
            a[0] -= 1
            G = meijerg(a, self.aother, self.bm, self.bother, self.argument)
            return 1/self.argument * ((self.an[0]-1)*self + G)
        elif len(self.bm) >= 1:
            b = list(self.bm)
            b[0] += 1
            G = meijerg(self.an, self.aother, b, self.bother, self.argument)
            return 1/self.argument * (self.bm[0]*self - G)
        else:
            return S.Zero

    def _eval_expand_func(self, deep=True, **hints):
        from sympy import hyperexpand
        return hyperexpand(self)

    @property
    def argument(self):
        """ Argument of the Meijer G-function. """
        return self.args[2]

    @property
    def an(self):
        """ First set of numerator parameters. """
        return self.args[0][0]

    @property
    def ap(self):
        """ Combined numerator parameters. """
        return self.args[0][0] + self.args[0][1]

    @property
    def aother(self):
        """ Second set of numerator parameters. """
        return self.args[0][1]

    @property
    def bm(self):
        """ First set of denominator parameters. """
        return self.args[1][0]

    @property
    def bq(self):
        """ Combined denominator parameters. """
        return self.args[1][0] + self.args[1][1]

    @property
    def bother(self):
        """ Second set of denominator parameters. """
        return self.args[1][1]

    @property
    def nu(self):
        """ A quantity related to the convergence region of the integral,
            c.f. references. """
        return sum(self.bq) - sum(self.ap)

    @property
    def delta(self):
        """ A quantity related to the convergence region of the integral,
            c.f. references. """
        return len(self.bm) + len(self.an) - S(len(self.ap) + len(self.bq))/2
