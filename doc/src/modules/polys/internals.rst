.. _polys-internals:

===============================================
Internals of the Polynomial Manipulation Module
===============================================

The implementation of the polynomials module is structured internally in
"levels". There are four levels, called L0, L1, L2 and L3. The levels three
and four contain the user-facing functionality and were described in the
previous section. This section focuses on levels zero and one.

Level zero provides core polynomial manipulation functionality with C-like,
low-level interfaces. Level one wraps this low-level functionality into object
oriented structures. These are *not* the classes seen by the user, but rather
classes used internally throughout the polys module.

There is one additional complication in the implementation. This comes from the
fact that all polynomial manipulations are relative to a *ground domain*. For
example, when factoring a polynomial like `x^{10} - 1`, one has to decide what
ring the coefficients are supposed to belong to, or less trivially, what
coefficients are allowed to appear in the factorization. This choice of
coefficients is called a ground domain. Typical choices include the integers
`\mathbb{Z}`, the rational numbers `\mathbb{Q}` or various related rings and
fields. But it is perfectly legitimate (although in this case uninteresting)
to factorize over polynomial rings such as `k[Y]`, where `k` is some fixed
field.

Thus the polynomial manipulation algorithms (both
complicated ones like factoring, and simpler ones like addition or
multiplication) have to rely on other code to manipulate the coefficients.
In the polynomial manipulation module, such code is encapsulated in so-called
"domains". A domain is basically a factory object: it takes various
representations of data, and converts them into objects with unified interface.
Every object created by a domain has to implement the arithmetic operations
`+`, `-` and `\times`. Other operations are accessed through the domain, e.g.
as in ``ZZ.quo(ZZ(4), ZZ(2))``.

Note that there is some amount of *circularity*: the polynomial ring domains
use the level one classes, the level one classes use the level zero functions,
and level zero functions use domains. It is possible, in principle, but not in
the current implementation, to work in rings like `k[X][Y]`. This would create
even more layers. For this reason, working in the isomorphic ring `k[X, Y]`
is preferred.

Domains
=======

.. currentmodule:: sympy.polys.domains

Here we document the various implemented ground domains. There are three
types: abstract domains, concrete domains, and "implementation domains".
Abstract domains cannot be (usefully) instantiated at all, and just collect
together functionality shared by many other domains. Concrete domains are
those meant to be instantiated and used in the polynomial manipulation
algorithms. In some cases, there are various possible ways to implement the
data type the domain provides. For example, depending on what libraries are
available on the system, the integers are implemented either using the python
built-in integers, or using gmpy. Note that various aliases are created
automatically depending on the libraries available. As such e.g. ``ZZ`` always
refers to the most efficient implementation of the integer ring available.

Abstract Domains
****************

.. autoclass:: sympy.polys.domains.domain.Domain
   :members:

.. autoclass:: sympy.polys.domains.field.Field
   :members:

.. autoclass:: sympy.polys.domains.ring.Ring
   :members:

.. autoclass:: sympy.polys.domains.simpledomain.SimpleDomain
   :members:

.. autoclass:: sympy.polys.domains.compositedomain.CompositeDomain
   :members:

Concrete Domains
****************

.. autoclass:: FiniteField
   :members:

.. autoclass:: IntegerRing
   :members:

.. autoclass:: PolynomialRing
   :members:

.. autoclass:: RationalField
   :members:

.. autoclass:: AlgebraicField
   :members:

.. autoclass:: FractionField
   :members:

.. autoclass:: RealField
   :members:

.. autoclass:: ExpressionDomain
   :members:

Implementation Domains
**********************

.. autoclass:: PythonFiniteField
.. autoclass:: GMPYFiniteField

.. autoclass:: PythonIntegerRing
.. autoclass:: GMPYIntegerRing

.. autoclass:: PythonRationalField
.. autoclass:: GMPYRationalField

Level One
=========

.. currentmodule:: sympy.polys.polyclasses

.. autoclass:: DMP
   :members:

.. autoclass:: DMF
   :members:

.. autoclass:: ANP
   :members:

Level Zero
==========

Level zero contains the bulk code of the polynomial manipulation module.

Manipulation of dense, multivariate polynomials
***********************************************

These functions can be used to manipulate polynomials in `K[X_0, \dots, X_u]`.
Functions for manipulating multivariate polynomials in the dense representation
have the prefix ``dmp_``. Functions which only apply to univariate polynomials
(i.e. `u = 0`)
have the prefix ``dup__``. The ground domain `K` has to be passed explicitly.
For many multivariate polynomial manipulation functions also the level `u`,
i.e. the number of generators minus one, has to be passed.
(Note that, in many cases, ``dup_`` versions of functions are available, which
may be slightly more efficient.)

**Basic manipulation:**

.. currentmodule:: sympy.polys.densebasic

.. autofunction:: dmp_LC
.. autofunction:: dmp_TC
.. autofunction:: dmp_ground_LC
.. autofunction:: dmp_ground_TC
.. autofunction:: dmp_true_LT
.. autofunction:: dmp_degree
.. autofunction:: dmp_degree_in
.. autofunction:: dmp_degree_list
.. autofunction:: dmp_strip
.. autofunction:: dmp_validate
.. autofunction:: dup_reverse
.. autofunction:: dmp_copy
.. autofunction:: dmp_to_tuple
.. autofunction:: dmp_normal
.. autofunction:: dmp_convert
.. autofunction:: dmp_from_sympy
.. autofunction:: dmp_nth
.. autofunction:: dmp_ground_nth
.. autofunction:: dmp_zero_p
.. autofunction:: dmp_zero
.. autofunction:: dmp_one_p
.. autofunction:: dmp_one
.. autofunction:: dmp_ground_p
.. autofunction:: dmp_ground
.. autofunction:: dmp_zeros
.. autofunction:: dmp_grounds
.. autofunction:: dmp_negative_p
.. autofunction:: dmp_positive_p
.. autofunction:: dmp_from_dict
.. autofunction:: dmp_to_dict
.. autofunction:: dmp_swap
.. autofunction:: dmp_permute
.. autofunction:: dmp_nest
.. autofunction:: dmp_raise
.. autofunction:: dmp_deflate
.. autofunction:: dmp_multi_deflate
.. autofunction:: dmp_inflate
.. autofunction:: dmp_exclude
.. autofunction:: dmp_include
.. autofunction:: dmp_inject
.. autofunction:: dmp_eject
.. autofunction:: dmp_terms_gcd
.. autofunction:: dmp_list_terms
.. autofunction:: dmp_apply_pairs
.. autofunction:: dmp_slice
.. autofunction:: dup_random

**Arithmetic operations:**

.. currentmodule:: sympy.polys.densearith

.. autofunction:: dmp_add_term
.. autofunction:: dmp_sub_term
.. autofunction:: dmp_mul_term
.. autofunction:: dmp_add_ground
.. autofunction:: dmp_sub_ground
.. autofunction:: dmp_mul_ground
.. autofunction:: dmp_quo_ground
.. autofunction:: dmp_exquo_ground
.. autofunction:: dup_lshift
.. autofunction:: dup_rshift
.. autofunction:: dmp_abs
.. autofunction:: dmp_neg
.. autofunction:: dmp_add
.. autofunction:: dmp_sub
.. autofunction:: dmp_add_mul
.. autofunction:: dmp_sub_mul
.. autofunction:: dmp_mul
.. autofunction:: dmp_sqr
.. autofunction:: dmp_pow
.. autofunction:: dmp_pdiv
.. autofunction:: dmp_prem
.. autofunction:: dmp_pquo
.. autofunction:: dmp_pexquo
.. autofunction:: dmp_rr_div
.. autofunction:: dmp_ff_div
.. autofunction:: dmp_div
.. autofunction:: dmp_rem
.. autofunction:: dmp_quo
.. autofunction:: dmp_exquo
.. autofunction:: dmp_max_norm
.. autofunction:: dmp_l1_norm
.. autofunction:: dmp_expand

**Further tools:**

.. currentmodule:: sympy.polys.densetools

.. autofunction:: dmp_integrate
.. autofunction:: dmp_integrate_in
.. autofunction:: dmp_diff
.. autofunction:: dmp_diff_in
.. autofunction:: dmp_eval
.. autofunction:: dmp_eval_in
.. autofunction:: dmp_eval_tail
.. autofunction:: dmp_diff_eval_in
.. autofunction:: dmp_trunc
.. autofunction:: dmp_ground_trunc
.. autofunction:: dup_monic
.. autofunction:: dmp_ground_monic
.. autofunction:: dup_content
.. autofunction:: dmp_ground_content
.. autofunction:: dup_primitive
.. autofunction:: dmp_ground_primitive
.. autofunction:: dup_extract
.. autofunction:: dmp_ground_extract
.. autofunction:: dup_real_imag
.. autofunction:: dup_mirror
.. autofunction:: dup_scale
.. autofunction:: dup_shift
.. autofunction:: dup_transform
.. autofunction:: dmp_compose
.. autofunction:: dup_decompose
.. autofunction:: dmp_lift
.. autofunction:: dup_sign_variations
.. autofunction:: dmp_clear_denoms
.. autofunction:: dmp_revert

Manipulation of dense, univariate polynomials with finite field coefficients
****************************************************************************
.. currentmodule:: sympy.polys.galoistools

Functions in this module carry the prefix ``gf_``, referring to the classical
name "Galois Fields" for finite fields. Note that many polynomial
factorization algorithms work by reduction to the finite field case, so having
special implementations for this case is justified both by performance, and by
the necessity of certain methods which do not even make sense over general
fields.

.. autofunction:: gf_crt
.. autofunction:: gf_crt1
.. autofunction:: gf_crt2
.. autofunction:: gf_int
.. autofunction:: gf_degree
.. autofunction:: gf_LC
.. autofunction:: gf_TC
.. autofunction:: gf_strip
.. autofunction:: gf_trunc
.. autofunction:: gf_normal
.. autofunction:: gf_from_dict
.. autofunction:: gf_to_dict
.. autofunction:: gf_from_int_poly
.. autofunction:: gf_to_int_poly
.. autofunction:: gf_neg
.. autofunction:: gf_add_ground
.. autofunction:: gf_sub_ground
.. autofunction:: gf_mul_ground
.. autofunction:: gf_quo_ground
.. autofunction:: gf_add
.. autofunction:: gf_sub
.. autofunction:: gf_mul
.. autofunction:: gf_sqr
.. autofunction:: gf_add_mul
.. autofunction:: gf_sub_mul
.. autofunction:: gf_expand
.. autofunction:: gf_div
.. autofunction:: gf_rem
.. autofunction:: gf_quo
.. autofunction:: gf_exquo
.. autofunction:: gf_lshift
.. autofunction:: gf_rshift
.. autofunction:: gf_pow
.. autofunction:: gf_pow_mod
.. autofunction:: gf_gcd
.. autofunction:: gf_lcm
.. autofunction:: gf_cofactors
.. autofunction:: gf_gcdex
.. autofunction:: gf_monic
.. autofunction:: gf_diff
.. autofunction:: gf_eval
.. autofunction:: gf_multi_eval
.. autofunction:: gf_compose
.. autofunction:: gf_compose_mod
.. autofunction:: gf_trace_map
.. autofunction:: gf_random
.. autofunction:: gf_irreducible
.. autofunction:: gf_irreducible_p
.. autofunction:: gf_sqf_p
.. autofunction:: gf_sqf_part
.. autofunction:: gf_sqf_list
.. autofunction:: gf_Qmatrix
.. autofunction:: gf_Qbasis
.. autofunction:: gf_berlekamp
.. autofunction:: gf_zassenhaus
.. autofunction:: gf_shoup
.. autofunction:: gf_factor_sqf
.. autofunction:: gf_factor
.. autofunction:: gf_value
.. autofunction:: gf_csolve

Manipulation of sparse, distributed polynomials and vectors
***********************************************************

Dense representations quickly require infeasible amounts of storage and
computation time if the number of variables increases. For this reason,
there is code to manipulate polynomials in a *sparse* representation.



.. currentmodule:: sympy.polys.rings

Sparse polynomials are represented as dictionaries.

.. autofunction:: ring
.. autofunction:: xring
.. autofunction:: vring
.. autofunction:: sring

.. autoclass:: PolyRing
   :members:

.. autoclass:: PolyElement
   :members:

In commutative algebra, one often studies not only polynomials, but also
*modules* over polynomial rings. The polynomial manipulation module provides
rudimentary low-level support for finitely generated free modules. This is
mainly used for Groebner basis computations (see there), so manipulation
functions are only provided to the extend needed. They carry the prefix
``sdm_``. Note that in examples, the generators of the free module are called
`f_1, f_2, \dots`.

.. currentmodule:: sympy.polys.distributedmodules

.. autofunction:: sdm_monomial_mul
.. autofunction:: sdm_monomial_deg
.. autofunction:: sdm_monomial_divides
.. autofunction:: sdm_LC
.. autofunction:: sdm_to_dict
.. autofunction:: sdm_from_dict
.. autofunction:: sdm_add
.. autofunction:: sdm_LM
.. autofunction:: sdm_LT
.. autofunction:: sdm_mul_term
.. autofunction:: sdm_zero
.. autofunction:: sdm_deg
.. autofunction:: sdm_from_vector
.. autofunction:: sdm_to_vector

Polynomial factorization algorithms
***********************************

Many variants of Euclid's algorithm:

.. currentmodule:: sympy.polys.euclidtools

Classical remainder sequence
----------------------------

Let `K` be a field, and consider the ring `K[X]` of polynomials in a single
indeterminate `X` with coefficients in `K`. Given two elements `f` and `g`
of `K[X]` with `g\neq 0` there are unique polynomials `q` and `r` such that
`f = qg + r` and `\deg(r) < \deg(g)` or `r = 0`. 
They are denoted by `\mathrm{quo}(f,g)` 
(*quotient*) and `\mathrm{rem}(f,g)` (*remainder*), so we have
the *division identity*

.. math::

  f = \mathrm{quo}(f,g)g + \mathrm{rem}(f,g).

It follows that every ideal `I` of `K[X]` is a principal ideal, generated by
any element `\neq 0` of minimum degree (assuming `I` non-zero). In fact,
if `g` is such a polynomial and `f` is any element of `I`, 
`\mathrm{rem}(f,g)` belongs to `I` as a linear combination of `f` and `g`, 
hence must be zero; therefore `f` is a multiple of `g`. 

Using this result it is possible to find a `greatest common 
divisor <http://en.wikipedia.org/wiki/Greatest_common_divisor>`_
(gcd) of any polynomials `f,g,\ldots` in `K[X]`.
If `I` is the ideal formed by all linear combinations of the given polynomials
with coefficients in `K[X]`, and `d` is its generator,
then every common divisor of the polynomials also divides `d`.
On the other hand, the given polynomials are multiples of the generator `d`;
hence `d` is a gcd of the polynomials, denoted `\mathrm{gcd}(f,g,\ldots)`. 

An algorithm for the gcd of two polynomials `f` and `g` in `K[X]` can
now be obtained as follows.
By the division identity, `r = \mathrm{rem}(f,g)` is in the ideal generated
by `f` and `g`, as well as `f` is in the ideal generated by `g` and `r`.
Hence the ideals generated by the pairs `(f,g)` and `(g,r)` are the same.
Set `f_0 = f`, `f_1 = g`, and define recursively
`f_i = \mathrm{rem}(f_{i-2},f_{i-1})` for `i\ge 2`.
The recursion ends after a finite number of steps with `f_{k+1}=0`,
since the degrees of the polynomials are strictly decreasing. 
By the above remark, all the pairs `(f_{i-1},f_i)` generate the same ideal.
In particular, the ideal generated by `f` and `g` is generated by `f_k`
alone as `f_{k+1} = 0`. Hence `d = f_k` is a gcd of `f` and `g`.

The sequence of polynomials `f_0`, `f_1,\ldots, f_k` is called the
*Euclidean polynomial remainder sequence* determined by `(f,g)` because
of the analogy with the classical `Euclidean algorithm 
<http://en.wikipedia.org/wiki/Euclidean_algorithm>`_ for the gcd of
natural numbers.

The algorithm may be extended to obtain an expression for `d` in terms of
`f` and `g` by using the full division identities
to write recursively each `f_i` as a linear combination of `f` and `g`.
This leads to an equation

.. math::

   d = uf + vg\qquad (u,v \in K[X])

analogous to `Bézout's identity 
<http://en.wikipedia.org/wiki/B%C3%A9zout%27s_identity>`_
in the case of integers.
 
.. autofunction:: dmp_half_gcdex
.. autofunction:: dmp_gcdex
.. autofunction:: dmp_invert
.. autofunction:: dmp_euclidean_prs

Simplified remainder sequences
------------------------------

Assume, as is usual, that the coefficient field `K` is 
the field of fractions of an integral domain `A`.
In this case the coefficients (numerators and denominators)
of the polynomials in the Euclidean remainder sequence
tend to grow very fast.

If `A` is a unique factorization domain, the coefficients may be
reduced by cancelling common factors of numerators and denominators.
Further reduction is possible noting that a gcd of polynomials in 
`K[X]` is not unique:
it may be multiplied by any (non-zero) constant factor.

Any polynomial `f` in `K[X]` can be simplified by extracting
the denominators and common factors of the numerators of its coefficients.
This yields the representation `f = cF` where `c\in K` is
the *content* of `f` and `F` is a *primitive* polynomial, i.e.,
a polynomial in `A[X]` with coprime coefficients.

It is possible to start the algorithm by replacing the given polynomials
`f` and `g` with their primitive parts. This will only modify
`\mathrm{rem}(f,g)` by a constant factor.
Replacing it with its primitive part and continuing recursively
we obtain all the primitive parts of the polynomials in
the Euclidean remainder sequence, including the primitive
`\mathrm{gcd}(f,g)`.

This sequence is the *primitive polynomial remainder sequence*.
It is an example of *general polynomial remainder sequences* where
the computed remainders are modified by constant multipliers (or divisors)
in order to simplify the results.

.. autofunction:: dmp_primitive_prs

Subresultant sequence
---------------------

The coefficients of the primitive polynomial sequence do not grow
exceedingly, but the computation of the primitive parts requires
extra processing effort. Besides, the method only works with fraction fields
of unique factorization domains, excluding, for example, the general number
fields.

Collins [Collins67] realized that the so-called *subresultant polynomials*
of a pair of polynomials also form a generalized remainder sequence.
The coefficients of these polynomials
are expressible as determinants in the coefficients of the given
polynomials. Hence (the logarithm of) their size only grows linearly.
In addition, if the coefficients of the given polynomials
are in the subdomain `A`, so are those
of the subresultant polynomials. This means that the subresultant
sequence is comparable to the primitive remainder sequence without
relying on unique factorization in `A`.

To see how subresultants are associated with remainder sequences
recall that all polynomials `h` in the sequence are linear combinations of
the given polynomials `f` and `g`

.. math::

   h = uf+vg

with polynomials `u` and `v` in `K[X]`. Moreover, as is seen from the
extended Euclidean algorithm, the degrees of `u` and `v` are relatively
low, with limited growth from step to step.

Let `n = \deg(f)`, and `m = \deg(g)`, and assume `n\ge m`.
If `\deg(h) = j < m`, the coefficients of the powers `X^k` (`k > j`)
in the products `uf` and `vg` cancel each other. In particular, the
products must have the same degree, say, `l`.
Then `\deg(u) = l - n` and `\deg(v) = l - m` with a total of `2l -n - m + 2`
coefficients to be determined.

On the other hand, the equality `h = uf + vg` implies that `l - j`
linear combinations of the coefficients are zero, those associated with
the powers `X^i` (`j < i \leq l`), and one has a given non-zero value,
namely the leading coefficient of `h`.

To satisfy these `l - j + 1` linear equations the total number of
coefficients to be determined cannot be lower than `l - j + 1`, in general.
This leads to the inequality `l \ge n + m - j - 1`.
Taking `l = n + m - j - 1`, we obtain `\deg(u) = m - j - 1` and
`\deg(v) = n - j - 1`.

In the case `j = 0` the matrix of the resulting system of linear equations
is the `Sylvester matrix <http://en.wikipedia.org/wiki/Sylvester_matrix>`_ 
`S(f,g)` associated to `f` and `g`,
an `(n+m)\times (n+m)` matrix with coefficients of `f` and `g` as entries. 
Its determinant is the `resultant <http://en.wikipedia.org/wiki/Resultant>`_
`\mathrm{res}(f,g)` of the pair `(f,g)`. 
It is non-zero if and only if `f` and `g` are relatively prime.

For any `j` in the interval from `0` to `m` the matrix of the linear system is
an `(n+m-2j)\times (n+m-2j)` submatrix of the Sylvester matrix.
Its determinant `s_j(f,g)`
is called the `j` th *scalar subresultant* of `f` and `g`.

If `s_j(f,g)` is not zero, the associated equation `h = uf + vg` has
a unique solution where `\deg(h) = j` and the leading coefficient
of `h` has any given value; the one with leading coefficient
`s_j(f,g)` is the `j` th *subresultant polynomial* or, briefly,
*subresultant* of the pair `(f,g)`, and denoted `S_j(f,g)`.
This choice guarantees that the remainining coefficients
are also certain subdeterminants of the Sylvester matrix.
In particular, if `f` and `g` are in `A[X]`, so is `S_j(f,g)` as well.
This construction of subresultants applies to any `j` between
`0` and `m` regardless of the value of `s_j(f,g)`; if it is zero, then
`\deg(S_j(f,g)) < j`.

The properties of subresultants are as follows. Let `n_0 = \deg(f)`,
`n_1 = \deg(g)`, `n_2, \dots, n_k` be the decreasing sequence of
degrees of polynomials in a remainder sequence.
Let `0 \le j \le n_1`; then

- `s_j(f,g)\ne 0` if and only if `j = n_i` for some `i`.

- `S_j(f,g)\ne 0` if and only if `j = n_i` or `j = n_i - 1` for some `i`.

Normally, `n_{i-1} - n_i = 1` for `1 < i \le k`. If `n_{i-1} - n_i > 1`
for some `i` (the *abnormal* case), then `S_{n_{i-1}-1}(f,g)` and
`S_{n_i}(f,g)` are constant multiples of each other.
Hence either one could be included in the polynomial remainder sequence.
The former is given by smaller determinants,
so it is expected to have smaller coefficients.

Collins defined the *subresultant remainder sequence* by setting

.. math::

   f_i = S_{n_{i-1}-1}(f,g) \qquad (2\le i \le k).

In the normal case, these are the same as the `S_{n_i}(f,g)`. He also
derived expressions for the constants `\gamma_i` in the remainder
formulas

.. math::

   \gamma_i f_i = \mathrm{rem}(f_{i-2},f_{i-1})

in terms of the leading coefficients of `f_1,\ldots,f_{i-1}`, working
in the field `K`.

Brown and Traub [BrownTraub71] later developed a recursive procedure
for computing the coefficients `\gamma_i`. Their algorithm deals with elements
of the domain `A` exclusively (assuming `f,g\in A[X]`). However, in the
abnormal case there was a problem, a division in `A`
which could only be conjectured to be exact.

This was subsequently justified by Brown [Brown78] who showed that
the result of the division is, in fact, a scalar subresultant.
More specifically, the constant appearing in the computation of `f_i` is
`s_{n_{i-2}}(f,g)` (Theorem 3).
The implication of this discovery is that the scalar subresultants
are computed as by-products of the algorithm, all but `s_{n_k}(f,g)`
which is not needed after finding `f_{k+1} = 0`.
Completing the last step we obtain all non-zero scalar subresultants,
including the last one which is the resultant if this does not vanish.

.. autofunction:: dmp_inner_subresultants
.. autofunction:: dmp_subresultants
.. autofunction:: dmp_prs_resultant
.. autofunction:: dmp_zz_modular_resultant
.. autofunction:: dmp_zz_collins_resultant
.. autofunction:: dmp_qq_collins_resultant
.. autofunction:: dmp_resultant
.. autofunction:: dmp_discriminant
.. autofunction:: dmp_rr_prs_gcd
.. autofunction:: dmp_ff_prs_gcd
.. autofunction:: dmp_zz_heu_gcd
.. autofunction:: dmp_qq_heu_gcd
.. autofunction:: dmp_inner_gcd
.. autofunction:: dmp_gcd
.. autofunction:: dmp_lcm
.. autofunction:: dmp_content
.. autofunction:: dmp_primitive
.. autofunction:: dmp_cancel

Polynomial factorization in characteristic zero:

.. currentmodule:: sympy.polys.factortools

.. autofunction:: dmp_trial_division
.. autofunction:: dmp_zz_mignotte_bound
.. autofunction:: dup_zz_hensel_step
.. autofunction:: dup_zz_hensel_lift
.. autofunction:: dup_zz_zassenhaus
.. autofunction:: dup_zz_irreducible_p
.. autofunction:: dup_cyclotomic_p
.. autofunction:: dup_zz_cyclotomic_poly
.. autofunction:: dup_zz_cyclotomic_factor
.. autofunction:: dup_zz_factor_sqf
.. autofunction:: dup_zz_factor
.. autofunction:: dmp_zz_wang_non_divisors
.. autofunction:: dmp_zz_wang_test_points
.. autofunction:: dmp_zz_wang_lead_coeffs
.. autofunction:: dmp_zz_diophantine
.. autofunction:: dmp_zz_wang_hensel_lifting
.. autofunction:: dmp_zz_wang
.. autofunction:: dmp_zz_factor
.. autofunction:: dmp_ext_factor
.. autofunction:: dup_gf_factor
.. autofunction:: dmp_factor_list
.. autofunction:: dmp_factor_list_include
.. autofunction:: dmp_irreducible_p

Groebner basis algorithms
*************************

Groebner bases can be used to answer many problems in computational
commutative algebra. Their computation in rather complicated, and very
performance-sensitive. We present here various low-level implementations of
Groebner basis computation algorithms; please see the previous section of the
manual for usage.

.. currentmodule:: sympy.polys.groebnertools

.. autofunction:: groebner
.. autofunction:: spoly
.. autofunction:: red_groebner
.. autofunction:: is_groebner
.. autofunction:: is_minimal
.. autofunction:: is_reduced

.. currentmodule:: sympy.polys.fglmtools

.. autofunction:: matrix_fglm

Groebner basis algorithms for modules are also provided:

.. currentmodule:: sympy.polys.distributedmodules

.. autofunction:: sdm_spoly
.. autofunction:: sdm_ecart
.. autofunction:: sdm_nf_mora
.. autofunction:: sdm_groebner

Exceptions
==========

These are exceptions defined by the polynomials module.

TODO sort and explain

.. currentmodule:: sympy.polys.polyerrors

.. autoclass:: BasePolynomialError

.. autoclass:: ExactQuotientFailed
.. autoclass:: OperationNotSupported
.. autoclass:: HeuristicGCDFailed
.. autoclass:: HomomorphismFailed
.. autoclass:: IsomorphismFailed
.. autoclass:: ExtraneousFactors
.. autoclass:: EvaluationFailed
.. autoclass:: RefinementFailed
.. autoclass:: CoercionFailed
.. autoclass:: NotInvertible
.. autoclass:: NotReversible
.. autoclass:: NotAlgebraic
.. autoclass:: DomainError
.. autoclass:: PolynomialError
.. autoclass:: UnificationFailed
.. autoclass:: GeneratorsNeeded
.. autoclass:: ComputationFailed
.. autoclass:: GeneratorsError
.. autoclass:: UnivariatePolynomialError
.. autoclass:: MultivariatePolynomialError
.. autoclass:: PolificationFailed
.. autoclass:: OptionError
.. autoclass:: FlagError

Reference
=========

Modular GCD
***********

.. currentmodule:: sympy.polys.modulargcd

.. autofunction:: modgcd_univariate
.. autofunction:: modgcd_bivariate
.. autofunction:: modgcd_multivariate
.. autofunction:: func_field_modgcd

Manipulation of power series
****************************************************************************
.. currentmodule:: sympy.polys.ring_series

Functions in this module carry the prefix ``rs_``, standing for "ring series".
They manipulate finite power series in the sparse representation provided
by ``polys.ring.ring``.


.. autofunction:: rs_trunc
.. autofunction:: rs_mul
.. autofunction:: rs_square
.. autofunction:: rs_pow
.. autofunction:: rs_series_inversion
.. autofunction:: rs_series_from_list
.. autofunction:: rs_integrate
.. autofunction:: rs_log
.. autofunction:: rs_exp
.. autofunction:: rs_newton
.. autofunction:: rs_hadamard_exp
.. autofunction:: rs_compose_add


Undocumented
============

Many parts of the polys module are still undocumented, and even where there is
documentation it is scarce. Please contribute!
