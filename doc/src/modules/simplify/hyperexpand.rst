Details on the Hypergeometric Function Expansion Module
#######################################################

This page describes how the function :func:`hyperexpand` and related code
work. For usage, see the documentation of the symplify module.

Hypergeometric Function Expansion Algorithm
*******************************************

This section describes the algorithm used to expand hypergeometric functions.
Most of it is based on the papers [Roach1996]_ and [Roach1997]_.

Recall that the hypergeometric function is (initially) defined as

.. math ::
    {}_pF_q\left(\begin{matrix} a_1, \dots, a_p \\ b_1, \dots, b_q \end{matrix}
                     \middle| z \right)
        = \sum_{n=0}^\infty \frac{(a_1)_n \dots (a_p)_n}{(b_1)_n \dots (b_q)_n}
                            \frac{z^n}{n!}.

It turns out that there are certain differential operators that can change the
`a_p` and `p_q` parameters by integers. If a sequence of such
operators is known that converts the set of indices `a_r^0` and
`b_s^0` into `a_p` and `b_q`, then we shall say the pair `a_p,
b_q` is reachable from `a_r^0, b_s^0`. Our general strategy is thus as
follows: given a set `a_p, b_q` of parameters, try to look up an origin
`a_r^0, b_s^0` for which we know an expression, and then apply the
sequence of differential operators to the known expression to find an
expression for the Hypergeometric function we are interested in.

Notation
========

In the following, the symbol `a` will always denote a numerator parameter
and the symbol `b` will always denote a denominator parameter. The
subscripts `p, q, r, s` denote vectors of that length, so e.g.
`a_p` denotes a vector of `p` numerator parameters. The subscripts
`i` and `j` denote "running indices", so they should usually be used in
conjuction with a "for all `i`". E.g. `a_i < 4` for all `i`.
Uppercase subscripts `I` and `J` denote a chosen, fixed index. So
for example `a_I > 0` is true if the inequality holds for the one index
`I` we are currently interested in.

Incrementing and decrementing indices
=====================================

Suppose `a_i \ne 0`. Set `A(a_i) =
\frac{z}{a_i}\frac{\mathrm{d}}{dz}+1`. It is then easy to show that
`A(a_i) {}_p F_q\left({a_p \atop b_q} \middle| z \right) = {}_p F_q\left({a_p +
e_i \atop b_q} \middle| z \right)`, where `e_i` is the i-th unit vector.
Similarly for `b_j \ne 1` we set `B(b_j) = \frac{z}{b_j-1}
\frac{\mathrm{d}}{dz}+1` and find `B(b_j) {}_p F_q\left({a_p \atop b_q}
\middle| z \right) = {}_p F_q\left({a_p \atop b_q - e_i} \middle| z \right)`.
Thus we can increment upper and decrement lower indices at will, as long as we
don't go through zero. The `A(a_i)` and `B(b_j)` are called shift
operators.

It is also easy to show that `\frac{\mathrm{d}}{dz} {}_p F_q\left({a_p
\atop b_q} \middle| z \right) = \frac{a_1 \dots a_p}{b_1 \dots b_q} {}_p
F_q\left({a_p + 1 \atop b_q + 1} \middle| z \right)`, where `a_p + 1` is
the vector `a_1 + 1, a_2 + 1, \dots` and similarly for `b_q + 1`.
Combining this with the shift operators,  we arrive at one form of the
Hypergeometric differential equation: `\left[ \frac{\mathrm{d}}{dz}
\prod_{j=1}^q B(b_j) - \frac{a_1 \dots a_p}{(b_1-1) \dots (b_q-1)}
\prod_{i=1}^p A(a_i) \right] {}_p F_q\left({a_p \atop b_q} \middle| z \right) =
0`. This holds if all shift operators are defined, i.e. if no `a_i = 0`
and no `b_j = 1`. Clearing denominators and multiplying through by z we
arrive at the following equation: `\left[ z\frac{\mathrm{d}}{dz}
\prod_{j=1}^q \left(z\frac{\mathrm{d}}{dz} + b_j-1 \right) - z \prod_{i=1}^p
\left( z\frac{\mathrm{d}}{\mathrm{d}z} + a_i \right) \right] {}_p F_q\left({a_p
\atop b_q} \middle| z\right) = 0`. Even though our derivation does not show it,
it can be checked that this equation holds whenever the `{}_p F_q` is
defined.

Notice that, under suitable conditions on `a_I, b_J`, each of the
operators `A(a_i)`, `B(b_j)` and
`z\frac{\mathrm{d}}{\mathrm{d}z}` can be expressed in terms of `A(a_I)` or
`B(b_J)`. Our next aim is to write the Hypergeometric differential
equation as follows: `[X A(a_I) - r] {}_p F_q\left({a_p \atop b_q}
\middle| z\right) = 0`, for some operator `X` and some constant `r`
to be determined. If
`r \ne 0`, then we can write this as `\frac{-1}{r} X {}_p F_q\left({a_p +
e_I \atop b_q} \middle| z\right) = {}_p F_q\left({a_p \atop b_q} \middle|
z\right)`, and so `\frac{-1}{r}X` undoes the shifting of `A(a_I)`,
whence it will be called an inverse-shift operator.

Now `A(a_I)` exists if `a_I \ne 0`, and then
`z\frac{\mathrm{d}}{\mathrm{d}z} = a_I A(a_I) - a_I`. Observe also that all the
operators `A(a_i)`, `B(b_j)` and
`z\frac{\mathrm{d}}{\mathrm{d}z}` commute. We have `\prod_{i=1}^p \left(
z\frac{\mathrm{d}}{\mathrm{d}z} + a_i \right) = \left(\prod_{i=1, i \ne I}^p
\left( z\frac{\mathrm{d}}{\mathrm{d}z} + a_i \right)\right) a_I A(a_I)`, so
this gives us the first half of `X`. The other half does not have such a
nice expression. We find `z\frac{\mathrm{d}}{dz} \prod_{j=1}^q
\left(z\frac{\mathrm{d}}{dz} + b_j-1 \right) = \left(a_I A(a_I) - a_I\right)
\prod_{j=1}^q \left(a_I A(a_I) - a_I + b_j - 1\right)`. Since the first
half had no constant term, we infer `r = -a_I\prod_{j=1}^q(b_j - 1
-a_I)`.

This tells us under which conditions we can "un-shift" `A(a_I)`, namely
when `a_I \ne 0` and `r \ne 0`. Substituting `a_I - 1` for
`a_I` then tells us under what conditions we can decrement the index
`a_I`. Doing a similar analysis for `B(a_J)`, we arrive at the
following rules:

* An index `a_I` can be decremented if `a_I \ne 1` and
  `a_I \ne b_j` for all `b_j`.
* An index `b_J` can be
  incremented if `b_J \ne -1` and `b_J \ne a_i` for all
  `a_i`.

Combined with the conditions (stated above) for the existence of shift
operators, we have thus established the rules of the game!


Reduction of Order
==================

Notice that, quite trivially, if `a_I = b_J`, we have `{}_p
F_q\left({a_p \atop b_q} \middle| z \right) = {}_{p-1} F_{q-1}\left({a_p^*
\atop b_q^*} \middle| z \right)`, where `a_p^*` means `a_p` with
`a_I` omitted, and similarly for `b_q^*`. We call this reduction of
order.

In fact, we can do even better. If `a_I - b_J \in \mathbb{Z}_{>0}`, then
it is easy to see that `\frac{(a_I)_n}{(b_J)_n}` is actually a polynomial
in `n`. It is also easy to see that
`(z\frac{\mathrm{d}}{\mathrm{d}z})^k z^n = n^k z^n`. Combining these two remarks
we find:

  If `a_I - b_J \in \mathbb{Z}_{>0}`, then there exists a polynomial
  `p(n) = p_0 + p_1 n + \dots` (of degree `a_I - b_J`) such
  that `\frac{(a_I)_n}{(b_J)_n} = p(n)` and `{}_p F_q\left({a_p
  \atop b_q} \middle| z \right) = \left(p_0 + p_1
  z\frac{\mathrm{d}}{\mathrm{d}z} + p_2
  \left(z\frac{\mathrm{d}}{\mathrm{d}z}\right)^2 + \dots \right) {}_{p-1}
  F_{q-1}\left({a_p^* \atop b_q^*} \middle| z \right)`.

Thus any set of parameters `a_p, b_q` is reachable from a set of
parameters `c_r, d_s` where `c_i - d_j \in \mathbb{Z}` implies
`c_i < d_j`. Such a set of parameters `c_r, d_s` is called
suitable. Our database of known formulae should only contain suitable origins.
The reasons are twofold: firstly, working from suitable origins is easier, and
secondly, a formula for a non-suitable origin can be deduced from a lower order
formula, and we should put this one into the database instead.

Moving Around in the Parameter Space
====================================

It remains to investigate the following question: suppose `a_p, b_q` and
`a_p^0, b_q^0` are both suitable, and also `a_i - a_i^0 \in
\mathbb{Z}`, `b_j - b_j^0 \in \mathbb{Z}`. When is `a_p, b_q`
reachable from `a_p^0, b_q^0`? It is clear that we can treat all
parameters independently that are incongruent mod 1. So assume that `a_i`
and `b_j` are congruent to `r` mod 1, for all `i` and
`j`. The same then follows for `a_i^0` and `b_j^0`.

If `r \ne 0`, then any such `a_p, b_q` is reachable from any
`a_p^0, b_q^0`. To see this notice that there exist constants `c, c^0`,
congruent mod 1, such that `a_i < c < b_j` for all `i` and
`j`, and similarly `a_i^0 < c^0 < b_j^0`. If `n = c - c^0 > 0` then
we first inverse-shift all the `b_j^0` `n` times up, and then
similarly shift shift up all the `a_i^0` `n` times. If `n <
0` then we first inverse-shift down the `a_i^0` and then shift down the
`b_j^0`. This reduces to the case `c = c^0`. But evidently we can
now shift or inverse-shift around the `a_i^0` arbitrarily so long as we
keep them less than `c`, and similarly for the `b_j^0` so long as
we keep them bigger than `c`. Thus `a_p, b_q` is reachable from
`a_p^0, b_q^0`.

If `r = 0` then the problem is slightly more involved. WLOG no parameter
is zero. We now have one additional complication: no parameter can ever move
through zero. Hence `a_p, b_q` is reachable from `a_p^0, b_q^0` if
and only if the number of `a_i < 0` equals the number of `a_i^0 <
0`, and similarly for the `b_i` and `b_i^0`. But in a suitable set
of parameters, all `b_j > 0`! This is because the Hypergeometric function
is undefined if one of the `b_j` is a non-positive integer and all
`a_i` are smaller than the `b_j`. Hence the number of `b_j \le 0` is
always zero.

We can thus associate to every suitable set of parameters `a_p, b_q`,
where no `a_i = 0`, the following invariants:

    * For every `r \in [0, 1)` the number `\alpha_r` of parameters
      `a_i \equiv r \pmod{1}`, and similarly the number `\beta_r`
      of parameters `b_i \equiv r \pmod{1}`.
    * The number `\gamma`
      of integers `a_i` with `a_i < 0`.

The above reasoning shows that `a_p, b_q` is reachable from `a_p^0,
b_q^0` if and only if the invariants `\alpha_r, \beta_r, \gamma` all
agree. Thus in particular "being reachable from" is a symmetric relation on
suitable parameters without zeros.


Applying the Operators
======================

If all goes well then for a given set of parameters we find an origin in our
database for which we have a nice formula. We now have to apply (potentially)
many differential operators to it. If we do this blindly then the result will
be very messy. This is because with Hypergeometric type functions, the
derivative is usually expressed as a sum of two contiguous functions. Hence if
we compute `N` derivatives, then the answer will involve `2N`
contiguous functions! This is clearly undesirable. In fact we know from the
Hypergeometric differential equation that we need at most `\max(p, q+1)`
contiguous functions to express all derivatives.

Hence instead of differentiating blindly, we will work with a
`\mathbb{C}(z)`-module basis: for an origin `a_r^0, b_s^0` we either store
(for particularly pretty answers) or compute a set of `N` functions
(typically `N = \max(r, s+1)`) with the property that the derivative of
any of them is a `\mathbb{C}(z)`-linear combination of them. In formulae,
we store a vector `B` of `N` functions, a matrix `M` and a
vector `C` (the latter two with entries in `\mathbb{C}(z)`), with
the following properties:

* `{}_r F_s\left({a_r^0 \atop b_s^0} \middle| z \right) = C B`
* `z\frac{\mathrm{d}}{\mathrm{d}z} B = M B`.

Then we can compute as many derivatives as we want and we will always end up
with `\mathbb{C}(z)`-linear combination of at most `N` special
functions.

As hinted above, `B`, `M` and `C` can either all be stored
(for particularly pretty answers) or computed from a single `{}_p F_q`
formula.


Loose Ends
==========

This describes the bulk of the hypergeometric function algorithm. There a few
further tricks, described in the hyperexpand.py source file. The extension to
Meijer G-functions is also described there.

Meijer G-Functions of Finite Confluence
***************************************

Slater's theorem essentially evaluates a `G`-function as a sum of residues.
If all poles are simple, the resulting series can be recognised as
hypergeometric series. Thus a `G`-function can be evaluated as a sum of
Hypergeometric functions.

If the poles are not simple, the resulting series are not hypergeometric. This
is known as the "confluent" or "logarithmic" case (the latter because the
resulting series tend to contain logarithms). The answer depends in a
complicated way on the multiplicities of various poles, and there is no
accepted notation for representing it (as far as I know).
However if there are only finitely many
multiple poles, we can evaluate the `G` function as a sum of hypergeometric
functions, plus finitely many extra terms. I could not find any good reference
for this, which is why I work it out here.

Recall the general setup. We define

.. math::
    G(z) = \frac{1}{2\pi i} \int_L \frac{\prod_{j=1}^m \Gamma(b_j - s)
      \prod_{j=1}^n \Gamma(1 - a_j + s)}{\prod_{j=m+1}^q \Gamma(1 - b_j + s)
      \prod_{j=n+1}^p \Gamma(a_j - s)} z^s \mathrm{d}s,

where `L` is a contour starting and ending at `+\infty`, enclosing all of the
poles of `\Gamma(b_j - s)` for `j = 1, \dots, n` once in the negative
direction, and no other poles. Also the integral is assumed absolutely
convergent.

In what follows, for any complex numbers `a, b`, we write `a \equiv b \pmod{1}` if
and only if there exists an integer `k` such that `a - b = k`. Thus there are
double poles iff `a_i \equiv a_j \pmod{1}` for some `i \ne j \le n`.

We now assume that whenever `b_j \equiv a_i \pmod{1}` for `i \le m`, `j > n`
then `b_j < a_i`. This means that no quotient of the relevant gamma functions
is a polynomial, and can always be achieved by "reduction of order". Fix a
complex number `c` such that `\{b_i | b_i \equiv c \pmod{1}, i \le  m\}` is
not empty. Enumerate this set as `b, b+k_1, \dots, b+k_u`, with `k_i`
non-negative integers. Enumerate similarly
`\{a_j | a_j \equiv c \pmod{1}, j > n\}` as `b + l_1, \dots, b + l_v`.
Then `l_i > k_j` for all `i, j`. For finite confluence, we need to assume
`v \ge u` for all such `c`.

Let `c_1, \dots, c_w` be distinct `\pmod{1}` and exhaust the congruence classes
of the `b_i`. I claim

.. math :: G(z) = -\sum_{j=1}^w (F_j(z) + R_j(z)),

where `F_j(z)` is a hypergeometric function and `R_j(z)` is a finite sum, both
to be specified later. Indeed corresponding to every `c_j` there is
a sequence of poles, at mostly finitely many of them multiple poles. This is where
the `j`-th term comes from.

Hence fix again `c`, enumerate the relevant `b_i` as
`b, b + k_1, \dots, b + k_u`. We will look at the `a_j` corresponding to
`a + l_1, \dots, a + l_u`. The other `a_i` are not treated specially. The
corresponding gamma functions have poles at (potentially) `s = b + r` for
`r = 0, 1, \dots`. For `r \ge l_u`, pole of the integrand is simple. We thus set

.. math :: R(z) = \sum_{r=0}^{l_u - 1} res_{s = r + b}.

We finally need to investigate the other poles. Set `r = l_u + t`, `t \ge 0`.
A computation shows

.. math ::
       \frac{\Gamma(k_i - l_u - t)}{\Gamma(l_i - l_u - t)}
            = \frac{1}{(k_i - l_u - t)_{l_i - k_i}}
            = \frac{(-1)^{\delta_i}}{(l_u - l_i + 1)_{\delta_i}}
              \frac{(l_u - l_i + 1)_t}{(l_u - k_i + 1)_t},

where `\delta_i = l_i - k_i`.

Also

.. math ::
    \Gamma(b_j - l_u - b - t) =
        \frac{\Gamma(b_j - l_u - b)}{(-1)^t(l_u + b + 1 - b_j)_t}, \\

    \Gamma(1 - a_j + l_u + b + t) =
        \Gamma(1 - a_j + l_u + b) (1 - a_j + l_u + b)_t

and

.. math ::
    res_{s = b + l_u + t} \Gamma(b - s) = -\frac{(-1)^{l_u + t}}{(l_u + t)!}
              = -\frac{(-1)^{l_u}}{l_u!} \frac{(-1)^t}{(l_u+1)_t}.

Hence

.. math ::
    res_{s = b + l_u + t} =& -z^{b + l_u}
       \frac{(-1)^{l_u}}{l_u!}
       \prod_{i=1}^{u} \frac{(-1)^{\delta_i}}{(l_u - k_i + 1)_{\delta_i}}
       \frac{\prod_{j=1}^n \Gamma(1 - a_j + l_u + b)
             \prod_{j=1}^m \Gamma(b_j - l_u - b)^*}
            {\prod_{j=n+1}^p \Gamma(a_j - l_u - b)^* \prod_{j=m+1}^q
             \Gamma(1 - b_j + l_u + b)}
       \\ &\times
       z^t
       \frac{(-1)^t}{(l_u+1)_t}
       \prod_{i=1}^{u} \frac{(l_u - l_i + 1)_t}{(l_u - k_i + 1)_t}
       \frac{\prod_{j=1}^n (1 - a_j + l_u + b)_t
             \prod_{j=n+1}^p (-1)^t (l_u + b + 1 - a_j)_t^*}
            {\prod_{j=1}^m (-1)^t (l_u + b + 1 - b_j)_t^*
             \prod_{j=m+1}^q (1 - b_j + l_u + b)_t},

where the `*` means to omit the terms we treated specially.

We thus arrive at

.. math ::
    F(z) = C \times {}_{p+1}F_{q}\left(
        \begin{matrix} 1, (1 + l_u - l_i), (1 + l_u + b - a_i)^* \\
                       1 + l_u, (1 + l_u - k_i), (1 + l_u + b - b_i)^*
        \end{matrix} \middle| (-1)^{p-m-n} z\right),

where `C` designates the factor in the residue independent of `t`.
(This result can also be written in slightly simpler form by converting
all the `l_u` etc back to `a_* - b_*`, but doing so is going to require more
notation still and is not helpful for computation.)

Extending The Hypergeometric Tables
***********************************

Adding new formulae to the tables is straightforward. At the top of the file
``sympy/simplify/hyperexpand.py``, there is a function called
:func:`add_formulae`. Nested in it are defined two helpers,
``add(ap, bq, res)`` and ``addb(ap, bq, B, C, M)``, as well as dummys
``a``, ``b``, ``c``, and ``z``.

The first step in adding a new formula is by using ``add(ap, bq, res)``. This
declares ``hyper(ap, bq, z) == res``. Here ``ap`` and ``bq`` may use the
dummys ``a``, ``b``, and ``c`` as free symbols. For example the well-known formula
`\sum_0^\infty \frac{(-a)_n z^n}{n!} = (1-z)^a` is declared by the following
line: ``add((-a, ), (), (1-z)**a)``.

From the information provided, the matrices `B`, `C` and `M` will be computed,
and the formula is now available when expanding hypergeometric functions.
Next the test file ``sympy/simplify/tests/test_hyperexpand.py`` should be run,
in particular the test :func:`test_formulae`. This will test the newly added
formula numerically. If it fails, there is (presumably) a typo in what was
entered.

Since all newly-added formulae are probably relatively complicated, chances
are that the automatically computed basis is rather suboptimal (there is no
good way of testing this, other than observing very messy output). In this
case the matrices `B`, `C` and `M` should be computed by hand. Then the helper
``addb`` can be used to declare a hypergeometric formula with hand-computed
basis.

An example
==========

Because this explanation so far might be very theoretical and difficult to
understand, we walk through an explicit example now. We take the Fresnel
function `C(z)` which obeys the following hypergeometric representation:

.. math ::
    C(z) = z \cdot {}_{1}F_{2}\left.\left(
        \begin{matrix} \frac{1}{4} \\
                       \frac{1}{2}, \frac{5}{4}
        \end{matrix} \right| -\frac{\pi^2 z^4}{16}\right) \,.

First we try to add this formula to the lookup table by using the
(simpler) function ``add(ap, bq, res)``. The first two arguments
are simply the lists containing the parameter sets of `{}_{1}F_{2}`.
The ``res`` argument is a little bit more complicated. We only know
`C(z)` in terms of `{}_{1}F_{2}(\ldots | f(z))` with `f`
a function of `z`, in our case

.. math ::
   f(z) = -\frac{\pi^2 z^4}{16} \,.

What we need is a formula where the hypergeometric function has
only `z` as argument `{}_{1}F_{2}(\ldots | z)`. We
introduce the new complex symbol `w` and search for a function
`g(w)` such that

.. math ::
   f(g(w)) = w

holds. Then we can replace every `z` in `C(z)` by `g(w)`.
In the case of our example the function `g` could look like

.. math ::
   g(w) = \frac{2}{\sqrt{\pi}} \exp\left(\frac{i \pi}{4}\right) w^{\frac{1}{4}} \,.

We get these functions mainly by guessing and testing the result. Hence
we proceed by computing `f(g(w))` (and simplifying naively)

.. math ::
   f(g(w)) &= -\frac{\pi^2 g(w)^4}{16} \\
           &= -\frac{\pi^2 g\left(\frac{2}{\sqrt{\pi}} \exp\left(\frac{i \pi}{4}\right) w^{\frac{1}{4}}\right)^4}{16} \\
           &= -\frac{\pi^2 \frac{2^4}{\sqrt{\pi}^4} \exp\left(\frac{i \pi}{4}\right)^4 {w^{\frac{1}{4}}}^4}{16} \\
           &= -\exp\left(i \pi\right) w \\
           &= w

and indeed get back `w`. (In case of branched functions we have to be
aware of branch cuts. In that case we take `w` to be a positive real
number and check the formula. If what we have found works for positive
`w`, then just replace :func:`exp` inside any branched function by
:func:`exp\_polar` and what we get is right for `all` `w`.) Hence
we can write the formula as

.. math ::
   C(g(w)) = g(w) \cdot {}_{1}F_{2}\left.\left(
        \begin{matrix} \frac{1}{4} \\
                       \frac{1}{2}, \frac{5}{4}
        \end{matrix} \right| w\right) \,.

and trivially

.. math ::
   {}_{1}F_{2}\left.\left(
   \begin{matrix} \frac{1}{4} \\
                  \frac{1}{2}, \frac{5}{4}
   \end{matrix} \right| w\right)
   = \frac{C(g(w))}{g(w)}
   = \frac{C\left(\frac{2}{\sqrt{\pi}} \exp\left(\frac{i \pi}{4}\right) w^{\frac{1}{4}}\right)}
          {\frac{2}{\sqrt{\pi}} \exp\left(\frac{i \pi}{4}\right) w^{\frac{1}{4}}}

which is exactly what is needed for the third paramenter,
``res``, in ``add``. Finally, the whole function call to add
this rule to the table looks like::

  add([S(1)/4],
      [S(1)/2, S(5)/4],
      fresnelc(exp(pi*I/4)*root(z,4)*2/sqrt(pi)) / (exp(pi*I/4)*root(z,4)*2/sqrt(pi))
     )

Using this rule we will find that it works but the results are not really nice
in terms of simplicity and number of special function instances included.
We can obtain much better results by adding the formula to the lookup table
in another way. For this we use the (more complicated) function ``addb(ap, bq, B, C, M)``.
The first two arguments are again the lists containing the parameter sets of
`{}_{1}F_{2}`. The remaining three are the matrices mentioned earlier
on this page.

We know that the `n = \max{\left(p, q+1\right)}`-th derivative can be
expressed as a linear combination of lower order derivatives. The matrix
`B` contains the basis `\{B_0, B_1, \ldots\}` and is of shape
`n \times 1`. The best way to get `B_i` is to take the first
`n = \max(p, q+1)` derivatives of the expression for `{}_p F_q`
and take out usefull pieces. In our case we find that
`n = \max{\left(1, 2+1\right)} = 3`. For computing the derivatives,
we have to use the operator `z\frac{\mathrm{d}}{\mathrm{d}z}`. The
first basis element `B_0` is set to the expression for `{}_1 F_2`
from above:

.. math ::
   B_0 = \frac{ \sqrt{\pi} \exp\left(-\frac{\mathbf{\imath}\pi}{4}\right)
   C\left( \frac{2}{\sqrt{\pi}} \exp\left(\frac{\mathbf{\imath}\pi}{4}\right) z^{\frac{1}{4}}\right)}
   {2 z^{\frac{1}{4}}}

Next we compute `z\frac{\mathrm{d}}{\mathrm{d}z} B_0`. For this we can
directly use SymPy!

   >>> from sympy import Symbol, sqrt, exp, I, pi, fresnelc, root, diff, expand
   >>> z = Symbol("z")
   >>> B0 = sqrt(pi)*exp(-I*pi/4)*fresnelc(2*root(z,4)*exp(I*pi/4)/sqrt(pi))/\
   ...          (2*root(z,4))
   >>> z * diff(B0, z)
   z*(cosh(2*sqrt(z))/(4*z) - sqrt(pi)*exp(-I*pi/4)*fresnelc(2*z**(1/4)*exp(I*pi/4)/sqrt(pi))/(8*z**(5/4)))
   >>> expand(_)
   cosh(2*sqrt(z))/4 - sqrt(pi)*exp(-I*pi/4)*fresnelc(2*z**(1/4)*exp(I*pi/4)/sqrt(pi))/(8*z**(1/4))

Formatting this result nicely we obtain

.. math ::
   B_1^\prime =
   - \frac{1}{4} \frac{
     \sqrt{\pi}
     \exp\left(-\frac{\mathbf{\imath}\pi}{4}\right)
     C\left( \frac{2}{\sqrt{\pi}} \exp\left(\frac{\mathbf{\imath}\pi}{4}\right) z^{\frac{1}{4}}\right)
   }
   {2 z^{\frac{1}{4}}}
   + \frac{1}{4} \cosh{\left( 2 \sqrt{z} \right )}

Computing the second derivative we find

   >>> from sympy import (Symbol, cosh, sqrt, pi, exp, I, fresnelc, root,
   ...                    diff, expand)
   >>> z = Symbol("z")
   >>> B1prime = cosh(2*sqrt(z))/4 - sqrt(pi)*exp(-I*pi/4)*\
   ...           fresnelc(2*root(z,4)*exp(I*pi/4)/sqrt(pi))/(8*root(z,4))
   >>> z * diff(B1prime, z)
   z*(-cosh(2*sqrt(z))/(16*z) + sinh(2*sqrt(z))/(4*sqrt(z)) + sqrt(pi)*exp(-I*pi/4)*fresnelc(2*z**(1/4)*exp(I*pi/4)/sqrt(pi))/(32*z**(5/4)))
   >>> expand(_)
   sqrt(z)*sinh(2*sqrt(z))/4 - cosh(2*sqrt(z))/16 + sqrt(pi)*exp(-I*pi/4)*fresnelc(2*z**(1/4)*exp(I*pi/4)/sqrt(pi))/(32*z**(1/4))

which can be printed as

.. math ::
   B_2^\prime =
   \frac{1}{16} \frac{
     \sqrt{\pi}
     \exp\left(-\frac{\mathbf{\imath}\pi}{4}\right)
     C\left( \frac{2}{\sqrt{\pi}} \exp\left(\frac{\mathbf{\imath}\pi}{4}\right) z^{\frac{1}{4}}\right)
   }
   {2 z^{\frac{1}{4}}}
   - \frac{1}{16} \cosh{\left(2\sqrt{z}\right)}
   + \frac{1}{4} \sinh{\left(2\sqrt{z}\right)} \sqrt{z}

We see the common pattern and can collect the pieces. Hence it makes sense to
choose `B_1` and `B_2` as follows

.. math ::
   B =
   \left( \begin{matrix}
     B_0 \\ B_1 \\ B_2
   \end{matrix} \right)
   =
   \left( \begin{matrix}
     \frac{
       \sqrt{\pi}
       \exp\left(-\frac{\mathbf{\imath}\pi}{4}\right)
       C\left( \frac{2}{\sqrt{\pi}} \exp\left(\frac{\mathbf{\imath}\pi}{4}\right) z^{\frac{1}{4}}\right)
     }{2 z^{\frac{1}{4}}} \\
     \cosh\left(2\sqrt{z}\right) \\
     \sinh\left(2\sqrt{z}\right) \sqrt{z}
   \end{matrix} \right)

(This is in contrast to the basis `B = \left(B_0, B_1^\prime, B_2^\prime\right)` that would
have been computed automatically if we used just ``add(ap, bq, res)``.)

Because it must hold that `{}_p F_q\left(\cdots \middle| z \right) = C B`
the entries of `C` are obviously

.. math ::
   C =
   \left( \begin{matrix}
     1 \\ 0 \\ 0
   \end{matrix} \right)

Finally we have to compute the entries of the `3 \times 3` matrix `M`
such that `z\frac{\mathrm{d}}{\mathrm{d}z} B = M B` holds. This is easy.
We already computed the first part `z\frac{\mathrm{d}}{\mathrm{d}z} B_0`
above. This gives us the first row of `M`. For the second row we have:

   >>> from sympy import Symbol, cosh, sqrt, diff
   >>> z = Symbol("z")
   >>> B1 = cosh(2*sqrt(z))
   >>> z * diff(B1, z)
   sqrt(z)*sinh(2*sqrt(z))

and for the third one

   >>> from sympy import Symbol, sinh, sqrt, expand, diff
   >>> z = Symbol("z")
   >>> B2 = sinh(2*sqrt(z))*sqrt(z)
   >>> expand(z * diff(B2, z))
   sqrt(z)*sinh(2*sqrt(z))/2 + z*cosh(2*sqrt(z))

Now we have computed the entries of this matrix to be

.. math ::
   M =
   \left( \begin{matrix}
     -\frac{1}{4} & \frac{1}{4} & 0 \\
     0            & 0           & 1 \\
     0            & z           & \frac{1}{2} \\
   \end{matrix} \right)

Note that the entries of `C` and `M` should typically be
rational functions in `z`, with rational coefficients. This is all
we need to do in order to add a new formula to the lookup table for
``hyperexpand``.

Implemented Hypergeometric Formulae
***********************************

A vital part of the algorithm is a relatively large table of hypergeometric
function representations. The following automatically generated list contains
all the representations implemented in SymPy (of course many more are
derived from them). These formulae are mostly taken from [Luke1969]_ and
[Prudnikov1990]_. They are all tested numerically.

.. automodule:: sympy.simplify.hyperexpand_doc

References
**********

.. [Roach1996] Kelly B. Roach.  Hypergeometric Function Representations.
      In: Proceedings of the 1996 International Symposium on Symbolic and
      Algebraic Computation, pages 301-308, New York, 1996. ACM.

.. [Roach1997] Kelly B. Roach.  Meijer G Function Representations.
      In: Proceedings of the 1997 International Symposium on Symbolic and
      Algebraic Computation, pages 205-211, New York, 1997. ACM.

.. [Luke1969] Luke, Y. L. (1969), The Special Functions and Their
              Approximations, Volume 1.

.. [Prudnikov1990] A. P. Prudnikov, Yu. A. Brychkov and O. I. Marichev (1990).
     Integrals and Series: More Special Functions, Vol. 3,
     Gordon and Breach Science Publisher.
