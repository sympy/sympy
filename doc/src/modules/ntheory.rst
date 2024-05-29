.. _ntheory-module:

=============
Number Theory
=============

.. module:: sympy.ntheory.generate

Ntheory Class Reference
=======================

.. autoclass:: Sieve
   :members:

Ntheory Functions Reference
===========================

.. autofunction:: prime

.. autofunction:: primepi

.. autofunction:: nextprime

.. autofunction:: prevprime

.. autofunction:: primerange

.. autofunction:: randprime

.. autofunction:: primorial

.. autofunction:: cycle_length

.. autofunction:: composite

.. autofunction:: compositepi

.. module:: sympy.ntheory.factor_

.. autofunction:: smoothness

.. autofunction:: smoothness_p

.. autofunction:: multiplicity

.. autofunction:: perfect_power

.. autofunction:: pollard_rho

.. autofunction:: pollard_pm1

.. autofunction:: factorint

.. autofunction:: factorrat

.. autofunction:: primefactors

.. autofunction:: divisors

.. autofunction:: proper_divisors

.. autofunction:: divisor_count

.. autofunction:: proper_divisor_count

.. autofunction:: udivisors

.. autofunction:: udivisor_count

.. autofunction:: antidivisors

.. autofunction:: antidivisor_count

.. autofunction:: totient

.. autofunction:: reduced_totient

.. autofunction:: divisor_sigma

.. autofunction:: udivisor_sigma

.. autofunction:: core

.. autofunction:: digits

.. autofunction:: primenu

.. autofunction:: primeomega

.. autofunction:: mersenne_prime_exponent

.. autofunction:: is_perfect

.. autofunction:: abundance

.. autofunction:: is_abundant

.. autofunction:: is_deficient

.. autofunction:: is_amicable

.. autofunction:: is_carmichael

.. autofunction:: find_carmichael_numbers_in_range

.. autofunction:: find_first_n_carmichaels

.. module:: sympy.ntheory.modular

.. autofunction:: symmetric_residue

.. autofunction:: crt

.. autofunction:: crt1

.. autofunction:: crt2

.. autofunction:: solve_congruence

.. module:: sympy.ntheory.multinomial

.. autofunction:: binomial_coefficients

.. autofunction:: binomial_coefficients_list

.. autofunction:: multinomial_coefficients

.. autofunction:: multinomial_coefficients_iterator

.. module:: sympy.ntheory.partitions_

.. autofunction:: npartitions

.. module:: sympy.ntheory.primetest

.. autofunction:: is_fermat_pseudoprime

.. autofunction:: is_euler_pseudoprime

.. autofunction:: is_euler_jacobi_pseudoprime

.. autofunction:: is_square

.. autofunction:: mr

.. autofunction:: is_lucas_prp

.. autofunction:: is_strong_lucas_prp

.. autofunction:: is_extra_strong_lucas_prp

.. autofunction:: proth_test

.. autofunction:: is_mersenne_prime

.. autofunction:: isprime

.. autofunction:: is_gaussian_prime

.. module:: sympy.ntheory.residue_ntheory

.. autofunction:: n_order

.. autofunction:: is_primitive_root

.. autofunction:: primitive_root

.. autofunction:: sqrt_mod

.. autofunction:: sqrt_mod_iter

.. autofunction:: quadratic_residues

.. autofunction:: nthroot_mod

.. autofunction:: is_nthpow_residue

.. autofunction:: is_quad_residue

.. autofunction:: legendre_symbol

.. autofunction:: jacobi_symbol

.. autofunction:: mobius

.. autofunction:: discrete_log

.. autofunction:: quadratic_congruence

.. autofunction:: polynomial_congruence

.. autofunction:: binomial_mod

.. automodule:: sympy.ntheory.continued_fraction
   :members:

.. automodule:: sympy.ntheory.digits
   :members:

.. module:: sympy.ntheory.egyptian_fraction

.. autofunction:: egyptian_fraction

.. module:: sympy.ntheory.bbp_pi

.. autofunction:: pi_hex_digits

.. module:: sympy.ntheory.ecm

ECM function
============

The `ecm` function is a subexponential factoring algorithm capable of factoring
numbers of around ~35 digits comfortably within few seconds. The time complexity
of `ecm` is dependent on the smallest proper factor of the number. So even if the
number is really large but its factors are comparatively smaller then `ecm`
can easily factor them. For example we take `N` with 15 digit factors
`15154262241479`, `15423094826093`, `799333555511111`, `809709509409109`,
`888888877777777`, `914148152112161`. Now N is a 87 digit number. `ECM` takes
under around 47s to factorise this.

.. autofunction:: ecm

Examples
--------

 >>> from sympy.ntheory import ecm
 >>> ecm(7060005655815754299976961394452809, B1=100000, B2=1000000)
 {6988699669998001, 1010203040506070809}
 >>> ecm(122921448543883967430908091422761898618349713604256384403202282756086473494959648313841, B1=100000, B2=1000000)
 {15154262241479,
 15423094826093,
 799333555511111,
 809709509409109,
 888888877777777,
 914148152112161}

.. module:: sympy.ntheory.qs

QS function
===========

The `qs` function is a subexponential factoring algorithm, the fastest
factoring algorithm for numbers within 100 digits. The time complexity of
`qs` is dependent on the size of the number so it is used if the number contains
large factors. Due to this while factoring numbers first `ecm` is used to get
smaller factors of around ~15 digits then `qs` is used to get larger factors.

For factoring `2709077133180915240135586837960864768806330782747` which is a semi-prime number
with two 25 digit factors. `qs` is able to factorize this in around 248s.

.. autofunction:: qs

Examples
--------

 >>> from sympy.ntheory import qs
 >>> qs(5915587277*3267000013, 1000, 10000)
 {3267000013, 5915587277}
