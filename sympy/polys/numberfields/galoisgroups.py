"""Galois groups of polynomials. """

import random

from sympy.core.evalf import (
    evalf, fastlog, _evalf_with_bounded_error, quad_to_mpmath,
)
from sympy.core.numbers import Integer
from sympy.core.symbol import Dummy, symbols
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.polys.domains import ZZ
from sympy.polys.densebasic import dup_random
from sympy.polys.densetools import dup_eval
from sympy.polys.euclidtools import dup_discriminant
from sympy.polys.factortools import dup_factor_list
from sympy.polys.numberfields.utilities import coeff_search
from sympy.polys.polytools import Poly
from sympy.polys.sqfreetools import dup_sqf_p
from sympy.utilities.lambdify import lambdify

from mpmath import MPContext


class GaloisGroupException(Exception):
    ...


class ResolventException(GaloisGroupException):
    ...


class MaxTriesException(GaloisGroupException):
    ...


class Resolvent:
    r"""
    If $G$ is a subgroup of the symmetric group $S_n$,
    $F$ a multivariate polynomial in $\mathbb{Z}[X_1, \ldots, X_n]$,
    $H$ the stabilizer of $F$ in $G$ (i.e. the permutations $\sigma$ such that
    $F(X_{\sigma(1)}, \ldots, X_{\sigma(n)}) = F(X_1, \ldots, X_n)$), and $s$
    a set of left coset representatives of $H$ in $G$, then the resolvent
    polynomial $R(Y)$ is the product over $\sigma \in s$ of
    $Y - F(X_{\sigma(1)}, \ldots, X_{\sigma(n)})$.

    For example, consider the resolvent for the form
    $$F = X_0 X_2 + X_1 X_3$$
    and the group $G = S_4$. In this case, the stabilizer $H$ is the dihedral
    group $D4 = < (0123), (02) >$, and a set of representatives of $G/H$ is
    $\{I, (01), (03)\}$. The resolvent can be constructed as follows:

    >>> from sympy.combinatorics.permutations import Permutation
    >>> from sympy.core.symbol import symbols
    >>> from sympy.polys.numberfields.galoisgroups import Resolvent
    >>> X = symbols('X0 X1 X2 X3')
    >>> F = X[0]*X[2] + X[1]*X[3]
    >>> s = [Permutation([0, 1, 2, 3]), Permutation([1, 0, 2, 3]),
    ... Permutation([3, 1, 2, 0])]
    >>> R = Resolvent(F, X, s)

    This resolvent has three roots, which are the conjugates of ``F`` under the
    three permutations in ``s``:

    >>> R.root_lambdas[0](*X)
    X0*X2 + X1*X3
    >>> R.root_lambdas[1](*X)
    X0*X3 + X1*X2
    >>> R.root_lambdas[2](*X)
    X0*X1 + X2*X3

    Resolvents are useful for computing Galois groups. Given a polynomial $T$
    of degree $n$, we will use a resolvent $R$ where $Gal(T) \leq G \leq S_n$.
    We will then want to substitute the roots of $T$ for the variables $X_i$
    in $R$, and study things like the discriminant of $R$, and the way $R$
    factors over $\mathbb{Q}$.

    From the symmetry in $R$'s construction, and since $Gal(T) \leq G$, we know
    from Galois theory that the coefficients of $R$ must lie in $\mathbb{Z}$.
    This allows us to compute the coefficients of $R$ by approximating the
    roots of $T$ to sufficient precision, plugging these values in for the
    variables $X_i$ in the coefficient expressions of $R$, and then simply
    rounding to the nearest integer.

    In order to determine a sufficient precision for the roots of $T$, this
    ``Resolvent`` class imposes certain requirements on the form ``F``. It
    could be possible to design a different ``Resolvent`` class, that made
    different precision estimates, and different assumptions about ``F``.

    ``F`` must be homogeneous, and all terms must have unit coefficient.
    Furthermore, if $r$ is the number of terms in ``F``, and $t$ the total
    degree, and if $m$ is the number of conjugates of ``F``, i.e. the number
    of permutations in ``s``, then we require that $m < r 2^t$. Again, it is
    not impossible to work with forms ``F`` that violate these assumptions, but
    this ``Resolvent`` class requires them.

    Since determining the integer coefficients of the resolvent for a given
    polynomial $T$ is one of the main problems this class solves, we take some
    time to explain the precision bounds it uses.

    The general problem is:
    Given a multivariate polynomial $P \in \mathbb{Z}[X_1, \ldots, X_n]$, and a
    bound $M \in \mathbb{R}_+$, compute an $\varepsilon > 0$ such that for any
    complex numbers $a_1, \ldots, a_n$ with $|a_i| < M$, if the $a_i$ are
    approximated to within an accuracy of $\varepsilon$ by $b_i$, that is,
    $|a_i - b_i| < \varepsilon$ for $i = 1, \ldots, n$, then
    $|P(a_1, \ldots, a_n) - P(b_1, \ldots, b_n)| < 1/2$. In other words, if it
    is known that $P(a_1, \ldots, a_n) = c$ for some $c \in \mathbb{Z}$, then
    $P(b_1, \ldots, b_n)$ can be rounded to the nearest integer in order to
    determine $c$.

    To derive our error bound, consider the monomial $xyz$. Defining
    $d_i = b_i - a_i$, our error is
    $|(a_1 + d_1)(a_2 + d_2)(a_3 + d_3) - a_1 a_2 a_3|$, which is bounded
    above by $|(M + \varepsilon)^3 - M^3|$. Passing to a general monomial of
    total degree $t$, this expression is bounded by
    $M^{t-1}\varepsilon(t + 2^t\varepsilon/M)$ provided $\varepsilon < M$,
    and by $(t+1)M^{t-1}\varepsilon$ provided $\varepsilon < M/2^t$.
    But since our goal is to make the error less than $1/2$, we will choose
    $\varepsilon < 1/(2(t+1)M^{t-1})$, which implies the condition that
    $\varepsilon < M/2^t$, as long as $M \geq 2$.

    Passing from the general monomial to the general polynomial is easy, by
    scaling and summing error bounds.

    In our specific case, we are given a homogeneous polynomial $F$ of
    $r$ terms and total degree $t$, all of whose coefficients are $\pm 1$. We
    are given the $m$ permutations that make the conjugates of $F$, and
    we want to bound the error in the coefficients of the monic polynomial
    $R(Y)$ having $F$ and its conjugates as roots (i.e. the resolvent).

    For $j$ from $1$ to $m$, the coefficient of $Y^{m-j}$ in $R(Y)$ is the
    $j$th elementary symmetric polynomial in the conjugates of $F$. This sums
    the products of these conjugates, taken $j$ at a time, in all possible
    combinations. There are $\binom{m}{j}$ such combinations, and each product
    of $j$ conjugates of $F$ expands to a sum of $r^j$ terms, each of unit
    coefficient, and total degree $jt$. An error bound for the $j$th coeff of
    $R$ is therefore
    $$\binom{m}{j} r^j (jt + 1) M^{jt - 1} \varepsilon$$
    When our goal is to evaluate all the coefficients of $R$, we will want to
    use the maximum of these error bounds. It is clear that this bound is
    strictly increasing for $j$ up to the ceiling of $m/2$. After that point,
    the first factor $\binom{m}{j}$ begins to decrease, while the others
    continue to increase. However, the binomial coefficient never falls by more
    than a factor of $1/m$, so our assumptions that $M \geq 2$ and $m < r 2^t$
    are enough to tell us that the constant coefficient of $R$, i.e. that where
    $j = m$, has the largest error bound. Therefore we can use
    $$r^m (mt + 1) M^{mt - 1} \varepsilon$$
    as our error bound for all the coefficients.

    Note that this bound is also (more than) adequate to determine whether any
    of the roots of $R$ is an integer. Each of these roots is a single
    conjugate of $F$, which contains less error than the trace, i.e. the
    coefficient of $Y^{m - 1}$. By rounding the roots of $R$ to the nearest
    integers, we therefore get all the candidates for integer roots of $R$. By
    plugging these candidates into $R$, we can check whether any of them
    actually is a root.

    Note: We take the definition of resolvent from Cohen, but the error bound
    is ours.

    References
    ==========

    .. [1] Cohen, H. *A Course in Computational Algebraic Number Theory*.
       (Def 6.3.2)

    """

    def __init__(self, F, X, s):
        r"""
        Parameters
        ==========

        F : :py:class:`~.Expr`
            polynomial in the symbols in *X*
        X : list of :py:class:`~.Symbol`
        s : list of :py:class:`~.Permutation`
            representing the cosets of the stabilizer of *F* in
            some subgroup $G$ of $S_n$, where $n$ is the length of *X*.
        """
        self.F = F
        self.X = X
        self.s = s

        # Number of conjugates:
        self.m = len(s)
        # Total degree of F (computed below):
        self.t = None
        # Number of terms in F (computed below):
        self.r = 0

        for monom, coeff in Poly(F).terms():
            if abs(coeff) != 1:
                raise ResolventException('Resolvent class expects forms with unit coeffs')
            t = sum(monom)
            if t != self.t and self.t is not None:
                raise ResolventException('Resolvent class expects homogeneous forms')
            self.t = t
            self.r += 1

        m, t, r = self.m, self.t, self.r
        if not m < r * 2**t:
            raise ResolventException('Resolvent class expects m < r*2^t')
        M = symbols('M')
        self.coeff_prec_func = Poly(r**m*(m*t + 1)*M**(m*t - 1))

        # The conjugates of F are the roots of the resolvent.
        # For evaluating these to required numerical precisions, we need
        # lambdified versions.
        # Note: for a given permutation sigma, the conjugate (sigma F) is
        # equivalent to lambda [sigma^(-1) X]: F.
        self.root_lambdas = [
            lambdify((~s[j])(X), F)
            for j in range(self.m)
        ]

        # For evaluating the coeffs, we'll also need lambdified versions of
        # the elementary symmetric functions for degree m.
        Y = symbols('Y')
        R = symbols(' '.join(f'R{i}' for i in range(m)))
        f = 1
        for r in R:
            f *= (Y - r)
        C = Poly(f, Y).coeffs()
        self.esf_lambdas = [lambdify(R, c) for c in C]

    def get_prec(self, M):
        r"""
        For a given upper bound *M* on the magnitude of the complex numbers to
        be plugged in for this resolvent's symbols, compute a sufficient
        precision for evaluating those complex numbers, such that the
        coefficients of the resolvent can be determined.

        Parameters
        ==========

        M : real number
            Upper bound on magnitude of the complex numbers to be plugged in.

        Returns
        =======

        int $m$
            such that $2^{-m}$ is a sufficient upper bound on the
            error in approximating the complex numbers to be plugged in.

        """
        M = max(M, 2)
        f = self.coeff_prec_func
        r, _, _, _ = evalf(2*f(M), 1, {})
        return fastlog(r) + 1

    def approximate_roots_of_poly(self, T):
        """
        Approximate the roots of a given polynomial *T* to sufficient precision
        in order to evaluate this resolvent's coefficients, or determine
        whether the resolvent has an integer root.

        Parameters
        ==========

        T : :py:class:`~.Poly`

        Returns
        =======

        list of elements of :ref:`CC`

        """
        ctx = MPContext()
        # Since we're going to be approximating the roots of T anyway, we can
        # get a good upper bound on the magnitude of the roots by starting with
        # a very low precision approx.
        T_roots = T.all_roots()
        approx0 = [quad_to_mpmath(_evalf_with_bounded_error(r, m=0), ctx=ctx) for r in T_roots]
        # Here we add 1 to account for the possible error in our initial approximation.
        M = max(abs(b) for b in approx0) + 1
        m = self.get_prec(M)
        # Now we can actually approximate the roots to the precision needed.
        ctx.prec = m
        approx1 = [quad_to_mpmath(_evalf_with_bounded_error(r, m=m), ctx=ctx) for r in T_roots]
        return approx1

    def eval_for_poly(self, T, find_integer_root=False):
        r"""
        Compute the integer values of the coefficients of this resolvent, when
        plugging in the roots of a given polynomial.

        Parameters
        ==========

        T : :py:class:`~.Poly`

        find_integer_root : ``bool``, default ``False``
            If ``True``, then also determine whether the resolvent has an
            integer root, and return the first one found, along with its
            index, i.e. the index of the permutation ``self.s[i]`` it
            corresponds to.

        Returns
        =======

        Tuple ``(R, a, i)``

            ``R`` is this resolvent as a dense
            univariate polynomial over :ref:`ZZ`, i.e. a list of :ref:`ZZ`.

            If *find_integer_root* was ``True``, then ``a`` and ``i`` are the
            first integer root found, and its index, if one exists.
            Otherwise ``a`` and ``i`` are both ``None``.

        """
        approx_roots_of_T = self.approximate_roots_of_poly(T)
        approx_roots_of_self = [r(*approx_roots_of_T) for r in self.root_lambdas]
        approx_coeffs_of_self = [c(*approx_roots_of_self) for c in self.esf_lambdas]

        R = []
        for c in approx_coeffs_of_self:
            if round(c.imag) != 0:
                # If precision was enough, this should never happen.
                raise ResolventException(f"Got non-integer coeff for resolvent: {c}")
            R.append(ZZ(round(c.real)))

        a0, i0 = None, None

        if find_integer_root:
            for i, r in enumerate(approx_roots_of_self):
                if round(r.imag) != 0:
                    continue
                if not dup_eval(R, (a := ZZ(round(r.real))), ZZ):
                    a0, i0 = a, i
                    break

        return R, a0, i0


def tschirnhausen_transformation(T, max_coeff=10, max_tries=30, history=None,
                                 fixed_order=True):
    r"""
    Given a univariate, monic, irreducible polynomial over the integers, find
    another such polynomial defining the same number field.

    Parameters
    ==========

    T : Poly
        The given polynomial
    max_coeff : int
        When choosing a random polynomial as part of the process,
        keep the coeffs between plus and minus this.
    max_tries : int
        Consider at most this many random polynomials.
    history : set, None, optional (default=None)
        Pass a set of ``Poly.rep``'s in order to prevent any of these
        polynomials from being chosen. The given poly *T* will automatically be
        added to this set, before we try to find a new one.
    fixed_order : bool, default True
        If ``True``, work through candidate transformations A(x) in a fixed
        order, from small coeffs to large, resulting in deterministic behavior.
        If ``False``, the A(x) are chosen randomly, while still working our way
        up from small coefficients to larger ones.

    Returns
    =======

    Pair ``(A, U)``

        ``A`` and ``U`` are ``Poly``, ``A`` is the
        transformation, and ``U`` is the transformed polynomial that defines the
        same number field as *T*. The polynomial ``A`` maps the roots of *T* to the
        roots of ``U``.

    Raises
    ======

    MaxTriesException
        if could not find a polynomial before exceeding *max_tries*.

    References
    ==========

    .. [1] Cohen, H. *A Course in Computational Algebraic Number Theory*.
       (Alg 6.3.4)

    """
    X = Dummy('X')
    n = T.degree()
    if history is None:
        history = set()
    history.add(T.rep)

    if fixed_order:
        coeff_generators = {}
        deg_coeff_sum = 3
        current_degree = 2

    def get_coeff_generator(degree):
        gen = coeff_generators.get(degree, coeff_search(degree, 1))
        coeff_generators[degree] = gen
        return gen

    for i in range(max_tries):

        # We never use linear A(x), since applying a fixed linear transformation
        # to all roots will only multiply the discriminant of T by a square
        # integer. This will change nothing important. In particular, if disc(T)
        # was zero before, it will still be zero now, and typically we apply
        # the transformation in hopes of replacing T by a squarefree poly.

        if fixed_order:
            # If d is degree and c max coeff, we move through the dc-space
            # along lines of constant sum. First d + c = 3 with (d, c) = (2, 1).
            # Then d + c = 4 with (d, c) = (3, 1), (2, 2). Then d + c = 5 with
            # (d, c) = (4, 1), (3, 2), (2, 3), and so forth. For a given (d, c)
            # we go though all sets of coeffs where max = c, before moving on.
            gen = get_coeff_generator(current_degree)
            coeffs = next(gen)
            m = max(abs(c) for c in coeffs)
            if current_degree + m > deg_coeff_sum:
                if current_degree == 2:
                    deg_coeff_sum += 1
                    current_degree = deg_coeff_sum - 1
                else:
                    current_degree -= 1
                gen = get_coeff_generator(current_degree)
                coeffs = next(gen)
            a = [ZZ(1)] + [ZZ(c) for c in coeffs]

        else:
            # We use a progressive coeff bound, up to the max specified, since it
            # is preferable to succeed with smaller coeffs.
            # Give each coeff bound five tries, before incrementing.
            C = min(i//5 + 1, max_coeff)
            d = random.randint(2, n - 1)
            a = dup_random(d, -C, C, ZZ)

        A = Poly(a, T.gen)
        U = Poly(T.resultant(X - A), X)
        if U.rep not in history and dup_sqf_p(U.rep.rep, ZZ):
            return A, U
    raise MaxTriesException


def M20():
    """
    Return a representation of the group M20, a subgroup of S5 that is one of
    the possible Galois groups for polys of degree 5.
    """
    from sympy.combinatorics.perm_groups import PermutationGroup
    from sympy.combinatorics.permutations import Permutation
    return PermutationGroup(Permutation(0, 1, 2, 3, 4), Permutation(1, 2, 3, 4))


def G36minus(verbose=False):
    """
    Return a representation of the group G36-, a subgroup of S6 isomorphic
    to the semidirect product of C3^2 with C2^2.

    This is one of the possible Galois groups for polys of degree 6 (not yet
    supported!).
    """
    from sympy.combinatorics.perm_groups import PermutationGroup
    from sympy.combinatorics.permutations import Permutation
    # If verbose, conduct (and show) the search procedure whereby we found
    # this representation initially. Otherwise just form and return the group.
    if verbose:
        from sympy.combinatorics.named_groups import SymmetricGroup
        # For our normal subgroup N we take one obvious instance of C3^2 in S6:
        N = PermutationGroup(
            Permutation(5)(0, 1, 2),
            Permutation(5)(3, 4, 5)
        )
        # Now we need an instance H of the Klein Viergruppe, such that
        # H ^ N = {e}, and such that N will indeed be normal in N x H.
        # This means finding two elements a != b in S6 such that:
        #   order(a) == 2 == order(b)
        #   ab = ba
        #   a, b not in N
        #   aNa^1 subseteq N supseteq bNb^1
        cand = []
        for a in SymmetricGroup(6).elements:
            if a.order() == 2 and a not in N:
                if all(a * n * (~a) in N for n in N.elements):
                    cand.append(a)
        print(f'\n{len(cand)} candidates')
        print(cand)
        for a in cand:
            for b in cand:
                if b == a:
                    continue
                c = a*b
                if c != b*a:
                    continue
                if c not in cand:
                    continue
                print(f'{a} -- {b} -- {c}')
    # Form and return the group.
    G = PermutationGroup(
        Permutation(5)(0, 1, 2),
        Permutation(5)(3, 4, 5),
        Permutation(5)(0, 1),
        Permutation(5)(4, 5),
    )
    return G


def is_square(n):
    """Convenience to check if an Integer is a square. """
    return isinstance(sqrt(n), Integer)


def has_square_disc(T):
    """Convenience to check if a Poly has square discriminant. """
    return is_square(T.discriminant())


def _galois_group_degree_3(T, max_tries=30, randomize=False):
    r"""
    Compute the Galois group of a polynomial of degree 3.

    References
    ==========

    .. [1] Cohen, H. *A Course in Computational Algebraic Number Theory*.
    (Prop 6.3.5)
    """
    from sympy.combinatorics.named_groups import AlternatingGroup, SymmetricGroup
    return (AlternatingGroup(3), True) if has_square_disc(T) else (SymmetricGroup(3), False)


def _galois_group_degree_4_simple(T, max_tries=30, randomize=False):
    r"""
    Compute the Galois group of a polynomial of degree 4, using Alg 6.3.6
    of Cohen.

    References
    ==========

    .. [1] Cohen, H. *A Course in Computational Algebraic Number Theory*.

    """
    from sympy.combinatorics.permutations import Permutation
    from sympy.combinatorics.named_groups import (
        CyclicGroup, AbelianGroup, DihedralGroup, AlternatingGroup, SymmetricGroup
    )
    # Consider the resolvent for the form
    #   F = X0*X1^2 + X1*X2^2 + X2*X3^2 + X3*X0^2
    # and the group G = S4. In this case, the stabilizer H is C4 = < (0123) >,
    # and a set of representatives of G/H is {I, (01), (02), (03), (12), (23)}.
    X = symbols('X0 X1 X2 X3')
    F = X[0]*X[1]**2 + X[1]*X[2]**2 + X[2]*X[3]**2 + X[3]*X[0]**2
    s = [
        Permutation(3),
        Permutation(3)(0, 1),
        Permutation(3)(0, 2),
        Permutation(3)(0, 3),
        Permutation(3)(1, 2),
        Permutation(3)(2, 3),
    ]
    R = Resolvent(F, X, s)
    history = set()
    for i in range(max_tries):
        R_dup, _, _ = R.eval_for_poly(T)
        # If R is squarefree, we can proceed. Otherwise, apply a
        # Tschirnhausen transformation on T and try again.
        if dup_sqf_p(R_dup, ZZ):
            break
        _, T = tschirnhausen_transformation(T, max_tries=max_tries, history=history, fixed_order=not randomize)
    else:
        raise MaxTriesException

    # Compute list L of degrees of irreducible factors of R, in increasing order:
    fl = dup_factor_list(R_dup, ZZ)
    L = sorted(sum([
        [len(r) - 1] * e for r, e in fl[1]
    ], []))

    if L == [6]:
        return (AlternatingGroup(4), True) if has_square_disc(T) else (SymmetricGroup(4), False)

    if L == [1, 1, 4]:
        return (CyclicGroup(4), False)

    if L == [2, 2, 2]:
        return (AbelianGroup(2, 2), True)

    assert L == [2, 4]
    return (DihedralGroup(4), False)


def _galois_group_degree_4(T, max_tries=30, randomize=False):
    r"""
    Compute the Galois group of a polynomial of degree 4, following Alg 6.3.7
    of Cohen.

    References
    ==========

    .. [1] Cohen, H. *A Course in Computational Algebraic Number Theory*.

    """
    from sympy.combinatorics.permutations import Permutation
    from sympy.combinatorics.named_groups import (
        CyclicGroup, AbelianGroup, DihedralGroup, AlternatingGroup, SymmetricGroup
    )

    X = symbols('X0 X1 X2 X3')
    # We start by considering the resolvent for the form
    #   F = X0*X2 + X1*X3
    # and the group G = S4. In this case, the stabilizer H is D4 = < (0123), (02) >,
    # and a set of representatives of G/H is {I, (01), (03)}
    F1 = X[0] * X[2] + X[1] * X[3]
    s1 = [
        Permutation(3),
        Permutation(3)(0, 1),
        Permutation(3)(0, 3)
    ]
    R1 = Resolvent(F1, X, s1)

    # In the second half of the algorithm (if we reach it), we use another
    # form and set of coset representatives. However, we may need to permute
    # them first, so cannot form their resolvent now.
    F2_pre = X[0] * X[1] ** 2 + X[1] * X[2] ** 2 + X[2] * X[3] ** 2 + X[3] * X[0] ** 2
    s2_pre = [
        Permutation(3),
        Permutation(3)(0, 2)
    ]

    history = set()
    for i in range(max_tries):
        if i > 0:
            # If we're retrying, need a new polynomial T.
            _, T = tschirnhausen_transformation(T, max_tries=max_tries, history=history, fixed_order=not randomize)

        R_dup, _, i0 = R1.eval_for_poly(T, find_integer_root=True)
        # If R is not squarefree, must retry.
        if not dup_sqf_p(R_dup, ZZ):
            continue

        # By Prop 6.3.1 of [1], Gal(T) is contained in A4 iff disc(T) is square.
        sq_disc = has_square_disc(T)

        if i0 is None:
            # By Thm 6.3.3 of [1], Gal(T) is not conjugate to any subgroup of the
            # stabilizer H = D4 that we chose. This means Gal(T) is either A4 or S4.
            return (AlternatingGroup(4), True) if sq_disc else (SymmetricGroup(4), False)

        # Gal(T) is conjugate to a subgroup of H = D4, so it is either V, C4
        # or D4 itself.

        if sq_disc:
            # Neither C4 nor D4 is contained in A4, so Gal(T) must be V.
            return (AbelianGroup(2, 2), True)

        # Gal(T) can only be D4 or C4.
        # We will now use our second resolvent, with G being that conjugate of D4 that
        # Gal(T) is contained in. To determine the right conjugate, we will need
        # the permutation corresponding to the integer root we found.
        sigma = s1[i0]
        # Applying sigma means permuting the args of F, and
        # conjugating the set of coset representatives.
        F2 = F2_pre.subs(zip(X, sigma(X)), simultaneous=True)
        s2 = [sigma*tau*sigma for tau in s2_pre]
        R2 = Resolvent(F2, X, s2)
        R_dup, _, _ = R2.eval_for_poly(T)
        d = dup_discriminant(R_dup, ZZ)
        # If d is zero (R has a repeated root), must retry.
        if d == 0:
            continue
        if is_square(d):
            return (CyclicGroup(4), False)
        else:
            return (DihedralGroup(4), False)

    raise MaxTriesException


def _galois_group_degree_5(T, max_tries=30, randomize=False):
    r"""
    Compute the Galois group of a polynomial of degree 5, following Alg 6.3.9
    of Cohen.

    References
    ==========

    .. [1] Cohen, H. *A Course in Computational Algebraic Number Theory*.

    """
    from sympy.combinatorics.permutations import Permutation
    from sympy.combinatorics.named_groups import (
        CyclicGroup, DihedralGroup, AlternatingGroup, SymmetricGroup
    )

    # The ideas here are all the same as in the degree-4 method.
    # The specific resolvents we use, and how we interpret the results, are
    # adapted to the degree-5 case.

    X = symbols('X0 X1 X2 X3 X4')
    # For the first resolvent, we have G = S5,
    # and stabilizer H = M20 = < (01234), (1234) >.
    F1 = (X[0]**2*(X[1]*X[4] + X[2]*X[3])
        + X[1]**2*(X[2]*X[0] + X[3]*X[4])
        + X[2]**2*(X[3]*X[1] + X[4]*X[0])
        + X[3]**2*(X[4]*X[2] + X[0]*X[1])
        + X[4]**2*(X[0]*X[3] + X[1]*X[2]))
    s1 = [
        Permutation(4),
        Permutation(4)(0, 1),
        Permutation(4)(0, 2),
        Permutation(4)(0, 3),
        Permutation(4)(0, 4),
        Permutation(4)(1, 4)
    ]
    R1 = Resolvent(F1, X, s1)

    # For the second resolvent, we'll have G = D5, H = C5.
    F2_pre = X[0]*X[1]**2 + X[1]*X[2]**2 + X[2]*X[3]**2 + X[3]*X[4]**2 + X[4]*X[0]**2
    s2_pre = [
        Permutation(4),
        Permutation(4)(0, 1)(2, 4)
    ]

    history = set()
    for i in range(max_tries):
        if i > 0:
            _, T = tschirnhausen_transformation(T, max_tries=max_tries, history=history, fixed_order=not randomize)

        R_dup, _, i0 = R1.eval_for_poly(T, find_integer_root=True)
        if not dup_sqf_p(R_dup, ZZ):
            continue

        sq_disc = has_square_disc(T)

        if i0 is None:
            return (AlternatingGroup(5), True) if sq_disc else (SymmetricGroup(5), False)

        if not sq_disc:
            return (M20(), False)

        sigma = s1[i0]
        F2 = F2_pre.subs(zip(X, sigma(X)), simultaneous=True)
        s2 = [sigma*tau*sigma for tau in s2_pre]
        R2 = Resolvent(F2, X, s2)
        R_dup, _, _ = R2.eval_for_poly(T)
        d = dup_discriminant(R_dup, ZZ)
        if d == 0:
            continue
        if is_square(d):
            return (CyclicGroup(5), True)
        else:
            return (DihedralGroup(5), True)

    raise MaxTriesException



def galois_group(T, max_tries=30, randomize=False):
    r"""
    Compute the Galois group for polynomials *T* up to degree 5.

    Parameters
    ==========

    T : Poly
        Irreducible, monic polynomial over :ref:`ZZ`, whose Galois group
        is to be determined.
    max_tries : int, default 30
        Make at most this many attempts in those steps that involve
        generating Tschirnhausen transformations.
    randomize : bool, default False
        If ``True``, then use random coefficients when generating Tschirnhausen
        transformations. Otherwise try transformations in a fixed order,
        starting with small coefficients and degrees and working upward.

    Returns
    =======

    Pair ``(PermutationGroup, bool)``
        The first element is the Galois group, and the second says whether the
        group is contained in the alternating group $A_n$ ($n$ the degree of
        *T*).

    Raises
    ======

    ValueError
        if *T* is of an unsupported degree.

    MaxTriesException
        if could not complete before exceeding *max_tries* in those steps
        that involve generating Tschirnhausen transformations.

    """
    from sympy.combinatorics.named_groups import CyclicGroup
    gg = {
        3: _galois_group_degree_3,
        4: _galois_group_degree_4,
        5: _galois_group_degree_5,
    }
    max_supported = max(gg.keys())
    n = T.degree()
    if n > max_supported:
        raise ValueError(f"Only polynomials up to degree {max_supported} are supported.")
    if n < 1:
        raise ValueError("Constant polynomial has no Galois group.")
    if n < 3:
        return (CyclicGroup(n), n == 1)
    return gg[n](T, max_tries=max_tries, randomize=randomize)
