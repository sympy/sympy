Computing Integrals using Meijer G-Functions
********************************************

This text aims do describe in some detail the steps (and subtleties) involved
in using Meijer G-functions for computing definite and indefinite integrals.
We shall ignore proofs completely.

Overview
========

The algorithm to compute `\int f(x) \mathrm{d}x` or
`\int_0^\infty f(x) \mathrm{d}x` generally consists of three steps:

1. Rewrite the integrand using Meijer G-functions (one or sometimes two).
2. Apply an integration theorem, to get the answer (usually expressed as another
   G-function).
3. Expand the result in named special functions.

Step (3) is implemented in the function hyperexpand (q.v.). Steps (1) and (2)
are described below. Morever, G-functions are usually branched. Thus our treatment
of branched functions is described first.

Some other integrals (e.g. `\int_{-\infty}^\infty`) can also be computed by first
recasting them into one of the above forms. There is a lot of choice involved
here, and the algorithm is heuristic at best.

Polar Numbers and Branched Functions
====================================

Both Meijer G-Functions and Hypergeometric functions are typically branched
(possible branchpoints being `0`, `\pm 1`, `\infty`). This is not very important
when e.g. expanding a single hypergeometric function into named special functions,
since sorting out the branches can be left to the human user. However this
algorithm manipulates and transforms G-functions, and to do this correctly it needs
at least some crude understanding of the branchings involved.

To begin, we consider the set
`\mathcal{S} = \{(r, \theta) : r > 0, \theta \in \mathbb{R}\}`. We have a map
`p: \mathcal{S}: \rightarrow \mathbb{C}-\{0\}, (r, \theta) \mapsto r e^{i \theta}`.
Decreeing this to be a local biholomorphism gives `\mathcal{S}` both a topology
and a complex structure. This Riemann Surface is usually referred to as the
Riemann Surface of the logarithm, for the following reason:
We can define maps
`\operatorname{Exp}: \mathbb{C} \rightarrow \mathcal{S}, (x + i y) \mapsto (\exp(x), y)` and
`\operatorname{Log}: \mathcal{S} \rightarrow \mathbb{C}, (e^x, y) \mapsto x + iy`.
These can both be shown to be holomorphic, and are indeed mutual inverses.

We also sometimes formally attach a point "zero" (`0`) to `\mathcal{S}` and denote the
resulting object `\mathcal{S}_0`. Notably there is no complex structure
defined near `0`. A fundamental system of neighbourhoods is given by
`\{\operatorname{Exp}(z) : \Re(z) < k\}`, which at least defines a topology. Elements of
`\mathcal{S}_0` shall be called polar numbers.
We further define functions
`\operatorname{Arg}: \mathcal{S} \rightarrow \mathbb{R}, (r, \theta) \mapsto \theta` and
`|.|: \mathcal{S}_0 \rightarrow \mathbb{R}_{>0}, (r, \theta) \mapsto r`.
These have evident meaning and are both continuous everywhere.

Using these maps many operations can be extended from `\mathbb{C}` to
`\mathcal{S}`. We define `\operatorname{Exp}(a) \operatorname{Exp}(b) = \operatorname{Exp}(a + b)` for `a, b \in \mathbb{C}`,
also for `a \in \mathcal{S}` and `b \in \mathbb{C}` we define
`a^b = \operatorname{Exp}(b \operatorname{Log}(a))`.
It can be checked easily that using these definitions, many algebraic properties
holding for positive reals (e.g. `(ab)^c = a^c b^c`) which hold in `\mathbb{C}`
only for some numbers (because of branch cuts) hold indeed for all polar numbers.

As one peculiarity it should be mentioned that addition of polar numbers is not
usually defined. However, formal sums of polar numbers can be used to express
branching behaviour. For example, consider the functions `F(z) = \sqrt{1 + z}`
and `G(a, b) = \sqrt{a + b}`, where `a, b, z` are polar numbers.
The general rule is that functions of a single polar variable are defined in
such a way that they are continuous on circles, and agree with the usual
definition for positive reals. Thus if `S(z)` denotes the standard branch of
the square root function on `\mathbb{C}`, we are forced to define

.. math:: F(z) = \begin{cases}
    S(p(z)) &: |z| < 1 \\
    S(p(z)) &: -\pi < \operatorname{Arg}(z) + 4\pi n \le \pi \text{ for some } n \in \mathbb{Z} \\
    -S(p(z)) &: \text{else}
   \end{cases}.

(We are omitting `|z| = 1` here, this does not matter for integration.)
Finally we define `G(a, b) = \sqrt{a}F(b/a)`.

Representing Branched Functions on the Argand Plane
===================================================

Suppose `f: \mathcal{S} \to \mathbb{C}` is a holomorphic function. We wish to
define a function `F` on (part of) the complex numbers `\mathbb{C}` that
represents `f` as closely as possible. This process is knows as "introducing
branch cuts". In our situation, there is actually a canonical way of doing this
(which is adhered to in all of SymPy), as follows: Introduce the "cut complex
plane"
`C = \mathbb{C} \setminus \mathbb{R}_{\le 0}`. Define a function
`l: C \to \mathcal{S}` via `re^{i\theta} \mapsto r \operatorname{Exp}(i\theta)`. Here `r > 0`
and `-\pi < \theta \le \pi`. Then `l` is holomorphic, and we define
`G = f \circ l`. This called "lifting to the principal branch" throughout the
SymPy documentation.

Table Lookups and Inverse Mellin Transforms
===========================================

Suppose we are given an integrand `f(x)` and are trying to rewrite it as a
single G-function. To do this, we first split `f(x)` into the form `x^s g(x)`
(where `g(x)` is supposed to be simpler than `f(x)`). This is because multiplicative
powers can be absorbed into the G-function later. This splitting is done by
``_split_mul(f, x)``. Then we assemble a tuple of functions that occur in
`f` (e.g. if `f(x) = e^x \cos{x}`, we would assemble the tuple `(\cos, \exp)`).
This is done by the function ``_mytype(f, x)``. Next we index a lookup table
(created using ``_create_lookup_table()``) with this tuple. This (hopefully)
yields a list of Meijer G-function formulae involving these functions, we then
pattern-match all of them. If one fits, we were successful, otherwise not and we
have to try something else.

Suppose now we want to rewrite as a product of two G-functions. To do this,
we (try to) find all inequivalent ways of splitting `f(x)` into a product
`f_1(x) f_2(x)`.
We could try these splittings in any order, but it is often a good idea to
minimise (a) the number of powers occuring in `f_i(x)` and (b) the number of
different functions occuring in `f_i(x)`. Thus given e.g.
`f(x) = \sin{x}\, e^{x} \sin{2x}` we should try `f_1(x) = \sin{x}\, \sin{2x}`,
`f_2(x) = e^{x}` first.
All of this is done by the function ``_mul_as_two_parts(f)``.

Finally, we can try a recursive Mellin transform technique. Since the Meijer
G-function is defined essentially as a certain inverse mellin transform,
if we want to write a function `f(x)` as a G-function, we can compute its mellin
transform `F(s)`. If `F(s)` is in the right form, the G-function expression
can be read off. This technique generalises many standard rewritings, e.g.
`e^{ax} e^{bx} = e^{(a + b) x}`.

One twist is that some functions don't have mellin transforms, even though they
can be written as G-functions. This is true for example for `f(x) = e^x \sin{x}`
(the function grows too rapidly to have a mellin transform). However if the function
is recognised to be analytic, then we can try to compute the mellin-transform of
`f(ax)` for a parameter `a`, and deduce the G-function expression by analytic
continuation. (Checking for analyticity is easy. Since we can only deal with a
certain subset of functions anyway, we only have to filter out those which are
not analyitc.)

The function ``_rewrite_single`` does the table lookup and recursive mellin
transform. The functions ``_rewrite1`` and ``_rewrite2`` respectively use
above-mentioned helpers and ``_rewrite_single`` to rewrite their argument as
respectively one or two G-functions.

Applying the Integral Theorems
==============================

If the integrand has been recast into G-functions, evaluating the integral is
relatively easy. We first do some substitutions to reduce e.g. the exponent
of the argument of the G-function to unity (see ``_rewrite_saxena_1`` and
``_rewrite_saxena``, respectively, for one or two G-functions). Next we go through
a list of conditions under which the integral theorem applies. It can fail for
basically two reasons: either the integral does not exist, or the manipulations
in deriving the theorem may not be allowed (for more details, see this [BlogPost]_).

Sometimes this can be remedied by reducing the argument of the G-functions
involved. For example it is clear that the G-function representing `e^z`
is satisfies `G(\operatorname{Exp}(2 \pi i)z) = G(z)` for all `z \in \mathcal{S}`. The function
``meijerg.get_period()`` can be used to discover this, and the function
``principal_branch(z, period)`` in ``functions/elementary/complexes.py`` can
be used to exploit the information. This is done transparently by the
integration code.

.. [BlogPost] http://nessgrh.wordpress.com/2011/07/07/tricky-branch-cuts/

The G-Function Integration Theorems
***********************************

This section intends to display in detail the definite integration theorems
used in the code. The following two formulae go back to Meijer (In fact he
proved more general formulae; indeed in the literature formulae are usually
staded in more general form. However it is very easy to deduce the general
formulae from the ones we give here. It seemed best to keep the theorems as
simple as possible, since they are very complicated anyway.):

1. .. math:: \int_0^\infty
    G_{p, q}^{m, n} \left.\left(\begin{matrix} a_1, \dots, a_p \\
                                               b_1, \dots, b_q \end{matrix}
            \right| \eta x \right) \mathrm{d}x =
     \frac{\prod_{j=1}^m \Gamma(b_j + 1) \prod_{j=1}^n \Gamma(-a_j)}{\eta
           \prod_{j=m+1}^q \Gamma(-b_j) \prod_{j=n+1}^p \Gamma(a_j + 1)}

2. .. math:: \int_0^\infty
    G_{u, v}^{s, t} \left.\left(\begin{matrix} c_1, \dots, c_u \\
                                               d_1, \dots, d_v \end{matrix}
            \right| \sigma x \right)
    G_{p, q}^{m, n} \left.\left(\begin{matrix} a_1, \dots, a_p \\
                                               b_1, \dots, b_q \end{matrix}
            \right| \omega x \right)
    \mathrm{d}x =
    G_{v+p, u+q}^{m+t, n+s} \left.\left(
          \begin{matrix} a_1, \dots, a_n, -d_1, \dots, -d_v, a_{n+1}, \dots, a_p \\
                         b_1, \dots, b_m, -c_1, \dots, -c_u, b_{m+1}, \dots, b_q
          \end{matrix}
            \right| \frac{\omega}{\sigma} \right)

The more interesting question is under what conditions these formulae are
valid. Below we detail the conditions implemented in SymPy. They are an
amalgamation of conditions found in [Prudnikov1990]_ and [Luke1969]_; please
let us know if you find any errors.

Conditions of Convergence for Integral (1)
==========================================
.. TODO: Formatting could be improved.

We can without loss of generality assume `p \le q`, since the G-functions
of indices `m, n, p, q` and of indices `n, m, q, p` can be related easily
(see e.g. [Luke1969]_, section 5.3). We introduce the following notation:

.. math:: \xi = m + n - p \\
          \delta = m + n - \frac{p + q}{2}

.. math:: C_3: -\Re(b_j) < 1 \text{ for } j=1, \dots, m \\
               0 < -\Re(a_j) \text{ for } j=1, \dots, n

.. math:: C_3^*: -\Re(b_j) < 1 \text{ for } j=1, \dots, q \\
               0 < -\Re(a_j) \text{ for } j=1, \dots, p

.. math:: C_4: -\Re(\delta) + \frac{q + 1 - p}{2} > q - p

The convergence conditions will be detailed in several "cases", numbered one
to five. For later use it will be helpful to separate conditions "at infinity"
from conditions "at zero". By conditions "at infinity" we mean conditions that
only depend on the behaviour of the integrand for large, positive values
of `x`, whereas by conditions "at zero" we mean conditions that only depend on
the behaviour of the integrand on `(0, \epsilon)` for any `\epsilon > 0`.
Since all our conditions are specified in terms of parameters of the
G-functions, this distinction is not immediately visible. They are, however, of
very distinct character mathematically; the conditions at infinity being in
particular much harder to control.

In order for the integral theorem to be valid, conditions
`n` "at zero" and "at infinity" both have to be fulfilled, for some `n`.

These are the conditions "at infinity":

1. .. math:: \delta > 0 \wedge |\arg(\eta)| < \delta \pi \wedge (A \vee B \vee C),

   where

   .. math::
      A = 1 \le n \wedge p < q \wedge 1 \le m

   .. math::
      B = 1 \le p \wedge 1 \le m \wedge q = p+1 \wedge
                  \neg (n = 0 \wedge m = p + 1 )

   .. math::
      C = 1 \le n \wedge q = p \wedge |\arg(\eta)| \ne (\delta - 2k)\pi
             \text{ for } k = 0, 1, \dots
               \left\lceil \frac{\delta}{2} \right\rceil.
2. .. math:: n = 0 \wedge p + 1 \le m \wedge |\arg(\eta)| < \delta \pi
3. .. math:: (p < q \wedge 1 \le m \wedge \delta > 0 \wedge |\arg(\eta)| = \delta \pi)
              \vee (p \le q - 2 \wedge \delta = 0 \wedge \arg(\eta) = 0)
4. .. math:: p = q \wedge \delta = 0 \wedge \arg(\eta) = 0 \wedge \eta \ne 0
        \wedge \Re\left(\sum_{j=1}^p b_j - a_j \right) < 0
5. .. math:: \delta > 0 \wedge |\arg(\eta)| < \delta \pi

And these are the conditions "at zero":

1. .. math:: \eta \ne 0 \wedge C_3
2. .. math:: C_3
3. .. math:: C_3 \wedge C_4
4. .. math:: C_3
5. .. math:: C_3

Conditions of Convergence for Integral (2)
==========================================

We introduce the following notation:

.. many of the latex expressions below were generated semi-automatically

.. math:: b^* = s + t - \frac{u + v}{2}
.. math:: c^* = m + n - \frac{p + q}{2}
.. math:: \rho = \sum_{j=1}^v d_j - \sum_{j=1}^u c_j + \frac{u - v}{2} + 1
.. math:: \mu = \sum_{j=1}^q b_j - \sum_{j=1}^p a_j + \frac{p - q}{2} + 1
.. math:: \phi = q - p - \frac{u - v}{2} + 1
.. math:: \eta = 1 - (v - u) - \mu - \rho
.. math:: \psi = \frac{\pi(q - m - n) + |\arg(\omega)|}{q - p}
.. math:: \theta = \frac{\pi(v - s - t) + |\arg(\sigma)|)}{v - u}
.. math:: \lambda_c = (q - p)|\omega|^{1/(q - p)} \cos{\psi}
                    + (v - u)|\sigma|^{1/(v - u)} \cos{\theta}
.. math:: \lambda_{s0}(c_1, c_2) = c_1 (q - p)|\omega|^{1/(q - p)} \sin{\psi}
                    + c_2 (v - u)|\sigma|^{1/(v - u)} \sin{\theta}
.. math:: \lambda_s =
  \begin{cases} \operatorname{\lambda_{s0}}\left(-1,-1\right) \operatorname{\lambda_{s0}}\left(1,1\right) & \text{for}\: \arg(\omega) = 0 \wedge \arg(\sigma) = 0 \\\operatorname{\lambda_{s0}}\left(\operatorname{sign}\left(\operatorname{\arg}\left(\omega\right)\right),-1\right) \operatorname{\lambda_{s0}}\left(\operatorname{sign}\left(\operatorname{\arg}\left(\omega\right)\right),1\right) & \text{for}\: \arg(\omega) \ne 0 \wedge \arg(\sigma) = 0 \\\operatorname{\lambda_{s0}}\left(-1,\operatorname{sign}\left(\operatorname{\arg}\left(\sigma\right)\right)\right) \operatorname{\lambda_{s0}}\left(1,\operatorname{sign}\left(\operatorname{\arg}\left(\sigma\right)\right)\right) & \text{for}\: \arg(\omega) = 0 \wedge \arg(\sigma) \ne 0) \\\operatorname{\lambda_{s0}}\left(\operatorname{sign}\left(\operatorname{\arg}\left(\omega\right)\right),\operatorname{sign}\left(\operatorname{\arg}\left(\sigma\right)\right)\right) & \text{otherwise} \end{cases}
.. math:: z_0 = \frac{\omega}{\sigma} e^{-i\pi (b^* + c^*)}
.. math:: z_1 = \frac{\sigma}{\omega} e^{-i\pi (b^* + c^*)}

The following conditions will be helpful:

.. math:: C_1: (a_i - b_j \notin \mathbb{Z}_{>0} \text{ for } i = 1, \dots, n, j = 1, \dots, m) \\
               \wedge
               (c_i - d_j \notin \mathbb{Z}_{>0} \text{ for } i = 1, \dots, t, j = 1, \dots, s)
.. math:: C_2:
    \Re(1 + b_i + d_j) > 0 \text{ for } i = 1, \dots, m, j = 1, \dots, s
.. math:: C_3:
    \Re(a_i + c_j) < 1 \text{ for } i = 1, \dots, n, j = 1, \dots, t
.. math:: C_4:
    (p - q)\Re(c_i) - \Re(\mu) > -\frac{3}{2} \text{ for } i=1, \dots, t
.. math:: C_5:
    (p - q)\Re(1 + d_i) - \Re(\mu) > -\frac{3}{2} \text{ for } i=1, \dots, s
.. math:: C_6:
    (u - v)\Re(a_i) - \Re(\rho) > -\frac{3}{2} \text{ for } i=1, \dots, n
.. math:: C_7:
    (u - v)\Re(1 + b_i) - \Re(\rho) > -\frac{3}{2} \text{ for } i=1, \dots, m
.. math:: C_8:
    0 < \lvert{\phi}\rvert + 2 \Re\left(\left(\mu -1\right) \left(- u + v\right) + \left(- p + q\right) \left(\rho -1\right) + \left(- p + q\right) \left(- u + v\right)\right)
.. math:: C_9:
    0 < \lvert{\phi}\rvert - 2 \Re\left(\left(\mu -1\right) \left(- u + v\right) + \left(- p + q\right) \left(\rho -1\right) + \left(- p + q\right) \left(- u + v\right)\right)
.. math:: C_{10}:
    \lvert{\operatorname{arg}\left(\sigma\right)}\rvert < \pi b^{*}
.. math:: C_{11}:
    \lvert{\operatorname{arg}\left(\sigma\right)}\rvert = \pi b^{*}
.. math:: C_{12}:
    |\arg(\omega)| < c^*\pi
.. math:: C_{13}:
    |\arg(\omega)| = c^*\pi
.. math:: C_{14}^1:
    \left(z_0 \ne 1 \wedge |\arg(1 - z_0)| < \pi \right) \vee
    \left(z_0 = 1 \wedge \Re(\mu + \rho - u + v) < 1 \right)
.. math:: C_{14}^2:
    \left(z_1 \ne 1 \wedge |\arg(1 - z_1)| < \pi \right) \vee
    \left(z_1 = 1 \wedge \Re(\mu + \rho - p + q) < 1 \right)
.. math:: C_{14}:
    \phi = 0 \wedge b^* + c^* \le 1 \wedge (C_{14}^1 \vee C_{14}^2)
.. math:: C_{15}:
    \lambda_c > 0 \vee (\lambda_c = 0 \wedge \lambda_s \ne 0 \wedge \Re(\eta) > -1)
                  \vee (\lambda_c = 0 \wedge \lambda_s = 0 \wedge \Re(\eta) > 0)
.. math:: C_{16}: \int_0^\infty G_{u, v}^{s, t}(\sigma x) \mathrm{d} x
                     \text{ converges at infinity }
.. math:: C_{17}: \int_0^\infty G_{p, q}^{m, n}(\omega x) \mathrm{d} x
                     \text{ converges at infinity }

Note that `C_{16}` and `C_{17}` are the reason we split the convergence conditions for
integral (1).

With this notation established, the implemented convergence conditions can be enumerated
as follows:

1. .. math:: m n s t \neq 0 \wedge 0 < b^{*} \wedge 0 < c^{*} \wedge C_{1} \wedge C_{2} \wedge C_{3} \wedge C_{10} \wedge C_{12}
2. .. math:: u = v \wedge b^{*} = 0 \wedge 0 < c^{*} \wedge 0 < \sigma \wedge \Re{\rho} < 1 \wedge C_{1} \wedge C_{2} \wedge C_{3} \wedge C_{12}
3. .. math:: p = q \wedge u = v \wedge b^{*} = 0 \wedge c^{*} = 0 \wedge 0 < \sigma \wedge 0 < \omega \wedge \Re{\mu} < 1 \wedge \Re{\rho} < 1 \wedge \sigma \neq \omega \wedge C_{1} \wedge C_{2} \wedge C_{3}
4. .. math:: p = q \wedge u = v \wedge b^{*} = 0 \wedge c^{*} = 0 \wedge 0 < \sigma \wedge 0 < \omega \wedge \Re\left(\mu + \rho\right) < 1 \wedge \omega \neq \sigma \wedge C_{1} \wedge C_{2} \wedge C_{3}
5. .. math:: p = q \wedge u = v \wedge b^{*} = 0 \wedge c^{*} = 0 \wedge 0 < \sigma \wedge 0 < \omega \wedge \Re\left(\mu + \rho\right) < 1 \wedge \omega \neq \sigma \wedge C_{1} \wedge C_{2} \wedge C_{3}
6. .. math:: q < p \wedge 0 < s \wedge 0 < b^{*} \wedge 0 \leq c^{*} \wedge C_{1} \wedge C_{2} \wedge C_{3} \wedge C_{5} \wedge C_{10} \wedge C_{13}
7. .. math:: p < q \wedge 0 < t \wedge 0 < b^{*} \wedge 0 \leq c^{*} \wedge C_{1} \wedge C_{2} \wedge C_{3} \wedge C_{4} \wedge C_{10} \wedge C_{13}
8. .. math:: v < u \wedge 0 < m \wedge 0 < c^{*} \wedge 0 \leq b^{*} \wedge C_{1} \wedge C_{2} \wedge C_{3} \wedge C_{7} \wedge C_{11} \wedge C_{12}
9. .. math:: u < v \wedge 0 < n \wedge 0 < c^{*} \wedge 0 \leq b^{*} \wedge C_{1} \wedge C_{2} \wedge C_{3} \wedge C_{6} \wedge C_{11} \wedge C_{12}
10. .. math:: q < p \wedge u = v \wedge b^{*} = 0 \wedge 0 \leq c^{*} \wedge 0 < \sigma \wedge \Re{\rho} < 1 \wedge C_{1} \wedge C_{2} \wedge C_{3} \wedge C_{5} \wedge C_{13}
11. .. math:: p < q \wedge u = v \wedge b^{*} = 0 \wedge 0 \leq c^{*} \wedge 0 < \sigma \wedge \Re{\rho} < 1 \wedge C_{1} \wedge C_{2} \wedge C_{3} \wedge C_{4} \wedge C_{13}
12. .. math:: p = q \wedge v < u \wedge 0 \leq b^{*} \wedge c^{*} = 0 \wedge 0 < \omega \wedge \Re{\mu} < 1 \wedge C_{1} \wedge C_{2} \wedge C_{3} \wedge C_{7} \wedge C_{11}
13. .. math:: p = q \wedge u < v \wedge 0 \leq b^{*} \wedge c^{*} = 0 \wedge 0 < \omega \wedge \Re{\mu} < 1 \wedge C_{1} \wedge C_{2} \wedge C_{3} \wedge C_{6} \wedge C_{11}
14. .. math:: p < q \wedge v < u \wedge 0 \leq b^{*} \wedge 0 \leq c^{*} \wedge C_{1} \wedge C_{2} \wedge C_{3} \wedge C_{4} \wedge C_{7} \wedge C_{11} \wedge C_{13}
15. .. math:: q < p \wedge u < v \wedge 0 \leq b^{*} \wedge 0 \leq c^{*} \wedge C_{1} \wedge C_{2} \wedge C_{3} \wedge C_{5} \wedge C_{6} \wedge C_{11} \wedge C_{13}
16. .. math:: q < p \wedge v < u \wedge 0 \leq b^{*} \wedge 0 \leq c^{*} \wedge C_{1} \wedge C_{2} \wedge C_{3} \wedge C_{5} \wedge C_{7} \wedge C_{8} \wedge C_{11} \wedge C_{13} \wedge C_{14}
17. .. math:: p < q \wedge u < v \wedge 0 \leq b^{*} \wedge 0 \leq c^{*} \wedge C_{1} \wedge C_{2} \wedge C_{3} \wedge C_{4} \wedge C_{6} \wedge C_{9} \wedge C_{11} \wedge C_{13} \wedge C_{14}
18. .. math:: t = 0 \wedge 0 < s \wedge 0 < b^{*} \wedge 0 < \phi \wedge C_{1} \wedge C_{2} \wedge C_{10}
19. .. math:: s = 0 \wedge 0 < t \wedge 0 < b^{*} \wedge \phi < 0 \wedge C_{1} \wedge C_{3} \wedge C_{10}
20. .. math::  n = 0 \wedge 0 < m \wedge 0 < c^{*} \wedge \phi < 0 \wedge C_{1} \wedge C_{2} \wedge C_{12}
21. .. math:: m = 0 \wedge 0 < n \wedge 0 < c^{*} \wedge 0 < \phi \wedge C_{1} \wedge C_{3} \wedge C_{12}
22. .. math:: s t = 0 \wedge 0 < b^{*} \wedge 0 < c^{*} \wedge C_{1} \wedge C_{2} \wedge C_{3} \wedge C_{10} \wedge C_{12}
23. .. math:: m n = 0 \wedge 0 < b^{*} \wedge 0 < c^{*} \wedge C_{1} \wedge C_{2} \wedge C_{3} \wedge C_{10} \wedge C_{12}
24. .. math:: p < m + n \wedge t = 0 \wedge \phi = 0 \wedge 0 < s \wedge 0 < b^{*} \wedge c^{*} < 0 \wedge \lvert{\operatorname{arg}\left(\omega\right)}\rvert < \pi \left(m + n - p + 1\right) \wedge C_{1} \wedge C_{2} \wedge C_{10} \wedge C_{14} \wedge C_{15}
25. .. math:: q < m + n \wedge s = 0 \wedge \phi = 0 \wedge 0 < t \wedge 0 < b^{*} \wedge c^{*} < 0 \wedge \lvert{\operatorname{arg}\left(\omega\right)}\rvert < \pi \left(m + n - q + 1\right) \wedge C_{1} \wedge C_{3} \wedge C_{10} \wedge C_{14} \wedge C_{15}
26. .. math:: p = q -1 \wedge t = 0 \wedge \phi = 0 \wedge 0 < s \wedge 0 < b^{*} \wedge 0 \leq c^{*} \wedge \pi c^{*} < \lvert{\operatorname{arg}\left(\omega\right)}\rvert \wedge C_{1} \wedge C_{2} \wedge C_{10} \wedge C_{14} \wedge C_{15}
27. .. math:: p = q + 1 \wedge s = 0 \wedge \phi = 0 \wedge 0 < t \wedge 0 < b^{*} \wedge 0 \leq c^{*} \wedge \pi c^{*} < \lvert{\operatorname{arg}\left(\omega\right)}\rvert \wedge C_{1} \wedge C_{3} \wedge C_{10} \wedge C_{14} \wedge C_{15}
28. .. math:: p < q -1 \wedge t = 0 \wedge \phi = 0 \wedge 0 < s \wedge 0 < b^{*} \wedge 0 \leq c^{*} \wedge \pi c^{*} < \lvert{\operatorname{arg}\left(\omega\right)}\rvert \wedge \lvert{\operatorname{arg}\left(\omega\right)}\rvert < \pi \left(m + n - p + 1\right) \wedge C_{1} \wedge C_{2} \wedge C_{10} \wedge C_{14} \wedge C_{15}
29. .. math:: q + 1 < p \wedge s = 0 \wedge \phi = 0 \wedge 0 < t \wedge 0 < b^{*} \wedge 0 \leq c^{*} \wedge \pi c^{*} < \lvert{\operatorname{arg}\left(\omega\right)}\rvert \wedge \lvert{\operatorname{arg}\left(\omega\right)}\rvert < \pi \left(m + n - q + 1 \right) \wedge C_{1} \wedge C_{3} \wedge C_{10} \wedge C_{14} \wedge C_{15}
30. .. math:: n = 0 \wedge \phi = 0 \wedge 0 < s + t \wedge 0 < m \wedge 0 < c^{*} \wedge b^{*} < 0 \wedge \lvert{\operatorname{arg}\left(\sigma\right)}\rvert < \pi \left(s + t - u + 1\right) \wedge C_{1} \wedge C_{2} \wedge C_{12} \wedge C_{14} \wedge C_{15}
31. .. math:: m = 0 \wedge \phi = 0 \wedge v < s + t \wedge 0 < n \wedge 0 < c^{*} \wedge b^{*} < 0 \wedge \lvert{\operatorname{arg}\left(\sigma\right)}\rvert < \pi \left(s + t - v + 1\right) \wedge C_{1} \wedge C_{3} \wedge C_{12} \wedge C_{14} \wedge C_{15}
32. .. math:: n = 0 \wedge \phi = 0 \wedge u = v -1 \wedge 0 < m \wedge 0 < c^{*} \wedge 0 \leq b^{*} \wedge \pi b^{*} < \lvert{\operatorname{arg}\left(\sigma\right)}\rvert \wedge \lvert{\operatorname{arg}\left(\sigma\right)}\rvert < \pi \left(b^{*} + 1\right) \wedge C_{1} \wedge C_{2} \wedge C_{12} \wedge C_{14} \wedge C_{15}
33. .. math:: m = 0 \wedge \phi = 0 \wedge u = v + 1 \wedge 0 < n \wedge 0 < c^{*} \wedge 0 \leq b^{*} \wedge \pi b^{*} < \lvert{\operatorname{arg}\left(\sigma\right)}\rvert \wedge \lvert{\operatorname{arg}\left(\sigma\right)}\rvert < \pi \left(b^{*} + 1\right) \wedge C_{1} \wedge C_{3} \wedge C_{12} \wedge C_{14} \wedge C_{15}
34. .. math:: n = 0 \wedge \phi = 0 \wedge u < v -1 \wedge 0 < m \wedge 0 < c^{*} \wedge 0 \leq b^{*} \wedge \pi b^{*} < \lvert{\operatorname{arg}\left(\sigma\right)}\rvert \wedge \lvert{\operatorname{arg}\left(\sigma\right)}\rvert < \pi \left(s + t - u + 1\right) \wedge C_{1} \wedge C_{2} \wedge C_{12} \wedge C_{14} \wedge C_{15}
35. .. math:: m = 0 \wedge \phi = 0 \wedge v + 1 < u \wedge 0 < n \wedge 0 < c^{*} \wedge 0 \leq b^{*} \wedge \pi b^{*} < \lvert{\operatorname{arg}\left(\sigma\right)}\rvert \wedge \lvert{\operatorname{arg}\left(\sigma\right)}\rvert < \pi \left(s + t - v + 1 \right) \wedge C_{1} \wedge C_{3} \wedge C_{12} \wedge C_{14} \wedge C_{15}
36. .. math:: C_{17} \wedge t = 0 \wedge u < s \wedge 0 < b^{*} \wedge C_{10} \wedge C_{1} \wedge C_{2} \wedge C_{3}
37. .. math:: C_{17} \wedge s = 0 \wedge v < t \wedge 0 < b^{*} \wedge C_{10} \wedge C_{1} \wedge C_{2} \wedge C_{3}
38. .. math:: C_{16} \wedge n = 0 \wedge p < m \wedge 0 < c^{*} \wedge C_{12} \wedge C_{1} \wedge C_{2} \wedge C_{3}
39. .. math:: C_{16} \wedge m = 0 \wedge q < n \wedge 0 < c^{*} \wedge C_{12} \wedge C_{1} \wedge C_{2} \wedge C_{3}


The Inverse Laplace Transform of a G-function
*********************************************

The inverse laplace transform of a Meijer G-function can be expressed as
another G-function. This is a fairly versatile method for computing this
transform. However, I could not find the details in the literature, so I work
them out here. In [Luke1969]_, section 5.6.3, there is a formula for the inverse
Laplace transform of a G-function of argument `bz`, and convergence conditions
are also given. However, we need a formula for argument `bz^a` for rational `a`.

We are asked to compute

.. math ::
    f(t) = \frac{1}{2\pi i} \int_{c-i\infty}^{c+i\infty} e^{zt} G(bz^a) \mathrm{d}z,

for positive real `t`. Three questions arise:

1. When does this integral converge?
2. How can we compute the integral?
3. When is our computation valid?


How to compute the integral
===========================

We shall work formally for now. Denote by `\Delta(s)` the product of gamma
functions appearing in the definition of `G`, so that

.. math :: G(z) = \frac{1}{2\pi i} \int_L \Delta(s) z^s \mathrm{d}s.

Thus

.. math ::
    f(t) = \frac{1}{(2\pi i)^2} \int_{c - i\infty}^{c + i\infty} \int_L
                  e^{zt} \Delta(s) b^s z^{as} \mathrm{d}s \mathrm{d}z.

We interchange the order of integration to get

.. math ::
    f(t) = \frac{1}{2\pi i} \int_L b^s \Delta(s)
          \int_{c-i\infty}^{c+i\infty} e^{zt} z^{as} \frac{\mathrm{d}z}{2\pi i}
                \mathrm{d}s.

The inner integral is easily seen to be
`\frac{1}{\Gamma(-as)} \frac{1}{t^{1+as}}`. (Using Cauchy's theorem and Jordan's
lemma deform the contour to run from `-\infty` to `-\infty`, encircling `0` once
in the negative sense. For `as` real and greater than one,
this contour can be pushed onto
the negative real axis and the integral is recognised as a product of a sine and
a gamma function. The formula is then proved using the functional equation of the
gamma function, and extended to the entire domain of convergence of the original
integral by appealing to analytic continuation.)
Hence we find

.. math ::
  f(t) = \frac{1}{t} \frac{1}{2\pi i} \int_L \Delta(s) \frac{1}{\Gamma(-as)}
                \left(\frac{b}{t^a}\right)^s \mathrm{d}s,

which is a so-called Fox H function (of argument `\frac{b}{t^a}`). For rational
`a`, this can be expressed as a Meijer G-function using the gamma function
multiplication theorem.

When this computation is valid
==============================

There are a number of obstacles in this computation. Interchange of integrals
is only valid if all integrals involved are absolutely convergent. In
particular the inner integral has to converge. Also, for our identification of
the final integral as a Fox H / Meijer G-function to be correct, the poles of
the newly obtained gamma function must be separated properly.

It is easy to check that the inner integal converges absolutely for
`\Re(as) < -1`. Thus the contour `L` has to run left of the line `\Re(as) = -1`.
Under this condition, the poles of the newly-introduced gamma function are
separated properly.

It remains to observe that the Meijer G-function is an analytic, unbranched
function of its parameters, and of the coefficient `b`. Hence so is `f(t)`.
Thus the final computation remains valid as long as the initial integral
converges, and if there exists a changed set of parameters where the computation
is valid. If we assume w.l.o.g. that `a > 0`, then the latter condition is
fulfilled if `G` converges along contours (2) or (3) of [Luke1969]_,
section 5.2, i.e. either `\delta >= \frac{a}{2}` or `p \ge 1, p \ge q`.

When the integral exists
========================

Using [Luke1969]_, section 5.10, for any given meijer G-function we can find a
dominant term of the form `z^a e^{bz^c}` (although this expression might not be
the best possible, because of cancellation).

We must thus investigate

.. math :: \lim_{T \to \infty} \int_{c-iT}^{c+iT}
                     e^{zt} z^a e^{bz^c} \mathrm{d}z.

(This principal value integral is the exact statement used in the Laplace
inversion theorem.) We write `z = c + i \tau`. Then
`arg(z) \to \pm \frac{\pi}{2}`, and so `e^{zt} \sim e^{it \tau}` (where `\sim`
shall always mean "asymptotically equivalent up to a positive real
multiplicative constant"). Also
`z^{x + iy} \sim |\tau|^x e^{i y \log{|\tau|}} e^{\pm x i \frac{\pi}{2}}.`

Set `\omega_{\pm} = b e^{\pm i \Re(c) \frac{\pi}{2}}`. We have three cases:

1. `b=0` or `\Re(c) \le 0`.
   In this case the integral converges if `\Re(a) \le -1`.
2. `b \ne 0`, `\Im(c) = 0`, `\Re(c) > 0`.
   In this case the integral converges if `\Re(\omega_{\pm}) < 0`.
3. `b \ne 0`, `\Im(c) = 0`, `\Re(c) > 0`, `\Re(\omega_{\pm}) \le 0`, and at least
   one of `\Re(\omega_{\pm}) = 0`.
   Here the same condition as in (1) applies.

Implemented G-Function Formulae
*******************************

An important part of the algorithm is a table expressing various functions
as Meijer G-functions. This is essentially a table of Mellin Transforms in
disguise. The following automatically generated table shows the formulae
currently implemented in SymPy. An entry "generated" means that the
corresponding G-function has a variable number of parameters.
This table is intended to shrink in future, when the algorithm's capabilities
of deriving new formulae improve. Of course it has to grow whenever a new class
of special functions is to be dealt with.

.. automodule:: sympy.integrals.meijerint_doc
