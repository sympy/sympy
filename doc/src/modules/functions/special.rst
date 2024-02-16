.. _special-functions:

Special
=======

Dirac Delta and Related Discontinuous Functions
-----------------------------------------------

.. autoclass:: sympy.functions.special.delta_functions.DiracDelta
   :members:

.. autoclass:: sympy.functions.special.delta_functions.Heaviside
   :members:

.. module:: sympy.functions.special.singularity_functions

.. autoclass:: sympy.functions.special.singularity_functions.SingularityFunction
   :members:

Gamma, Beta and Related Functions
---------------------------------

.. module:: sympy.functions.special.gamma_functions

.. autoclass:: sympy.functions.special.gamma_functions.gamma
   :members:
.. autoclass:: sympy.functions.special.gamma_functions.loggamma
   :members:
.. autoclass:: sympy.functions.special.gamma_functions.polygamma
   :members:
.. autoclass:: sympy.functions.special.gamma_functions.digamma
   :members:
.. autoclass:: sympy.functions.special.gamma_functions.trigamma
   :members:
.. autoclass:: sympy.functions.special.gamma_functions.uppergamma
   :members:
.. autoclass:: sympy.functions.special.gamma_functions.lowergamma
   :members:
.. autoclass:: sympy.functions.special.gamma_functions.multigamma
   :members:
.. module:: sympy.functions.special.beta_functions
.. autoclass:: sympy.functions.special.beta_functions.beta
   :members:

Error Functions and Fresnel Integrals
-------------------------------------

.. module:: sympy.functions.special.error_functions

.. autoclass:: sympy.functions.special.error_functions.erf
   :members:
.. autoclass:: sympy.functions.special.error_functions.erfc
   :members:
.. autoclass:: sympy.functions.special.error_functions.erfi
   :members:
.. autoclass:: sympy.functions.special.error_functions.erf2
   :members:
.. autoclass:: sympy.functions.special.error_functions.erfinv
   :members:
.. autoclass:: sympy.functions.special.error_functions.erfcinv
   :members:
.. autoclass:: sympy.functions.special.error_functions.erf2inv
   :members:

.. autoclass:: sympy.functions.special.error_functions.FresnelIntegral
   :members:

.. autoclass:: fresnels
   :members:
.. autoclass:: fresnelc
   :members:

Exponential, Logarithmic and Trigonometric Integrals
----------------------------------------------------

.. autoclass:: Ei
   :members:
.. autoclass:: expint
   :members:
.. autofunction:: E1
.. autoclass:: li
   :members:
.. autoclass:: Li
   :members:
.. autoclass:: Si
   :members:
.. autoclass:: Ci
   :members:
.. autoclass:: Shi
   :members:
.. autoclass:: Chi
   :members:

Bessel Type Functions
---------------------

.. module:: sympy.functions.special.bessel

.. autoclass:: sympy.functions.special.bessel.BesselBase
   :members:

.. autoclass:: sympy.functions.special.bessel.besselj
   :members:
.. autoclass:: sympy.functions.special.bessel.bessely
   :members:
.. _besseli:
.. autoclass:: sympy.functions.special.bessel.besseli
   :members:
.. autoclass:: sympy.functions.special.bessel.besselk
   :members:
.. autoclass:: sympy.functions.special.bessel.hankel1
   :members:
.. autoclass:: sympy.functions.special.bessel.hankel2
   :members:
.. autoclass:: sympy.functions.special.bessel.jn
   :members:
.. autoclass:: sympy.functions.special.bessel.yn
   :members:

.. autofunction:: sympy.functions.special.bessel.jn_zeros

.. autoclass:: sympy.functions.special.bessel.marcumq
   :members:

Airy Functions
--------------

.. autoclass:: sympy.functions.special.bessel.AiryBase
   :members:

.. autoclass:: sympy.functions.special.bessel.airyai
   :members:
.. autoclass:: sympy.functions.special.bessel.airybi
   :members:
.. autoclass:: sympy.functions.special.bessel.airyaiprime
   :members:
.. autoclass:: sympy.functions.special.bessel.airybiprime
   :members:

B-Splines
---------

.. autofunction:: sympy.functions.special.bsplines.bspline_basis
.. autofunction:: sympy.functions.special.bsplines.bspline_basis_set
.. autofunction:: sympy.functions.special.bsplines.interpolating_spline

Riemann Zeta and Related Functions
----------------------------------
.. module:: sympy.functions.special.zeta_functions

.. autoclass:: zeta
   :members:
.. autoclass:: dirichlet_eta
   :members:
.. autoclass:: polylog
   :members:
.. autoclass:: lerchphi
   :members:
.. autoclass:: stieltjes
   :members:

Hypergeometric Functions
------------------------
.. autoclass:: sympy.functions.special.hyper.hyper
   :members:

.. autoclass:: sympy.functions.special.hyper.meijerg
   :members:

.. autoclass:: sympy.functions.special.hyper.appellf1
   :members:

Elliptic integrals
------------------
.. module:: sympy.functions.special.elliptic_integrals

.. autoclass:: elliptic_k
   :members:
.. autoclass:: elliptic_f
   :members:
.. autoclass:: elliptic_e
   :members:
.. autoclass:: elliptic_pi
   :members:

Mathieu Functions
-----------------
.. module:: sympy.functions.special.mathieu_functions

.. autoclass:: sympy.functions.special.mathieu_functions.MathieuBase
   :members:

.. autoclass:: sympy.functions.special.mathieu_functions.mathieus
   :members:
.. autoclass:: sympy.functions.special.mathieu_functions.mathieuc
   :members:
.. autoclass:: sympy.functions.special.mathieu_functions.mathieusprime
   :members:
.. autoclass:: sympy.functions.special.mathieu_functions.mathieucprime
   :members:

Orthogonal Polynomials
----------------------

.. automodule:: sympy.functions.special.polynomials

Jacobi Polynomials
++++++++++++++++++

.. autoclass:: sympy.functions.special.polynomials.jacobi
   :members:

.. autofunction:: sympy.functions.special.polynomials.jacobi_normalized

Gegenbauer Polynomials
++++++++++++++++++++++

.. autoclass:: sympy.functions.special.polynomials.gegenbauer
   :members:

Chebyshev Polynomials
+++++++++++++++++++++

.. autoclass:: sympy.functions.special.polynomials.chebyshevt
   :members:

.. autoclass:: sympy.functions.special.polynomials.chebyshevu
   :members:

.. autoclass:: sympy.functions.special.polynomials.chebyshevt_root
   :members:

.. autoclass:: sympy.functions.special.polynomials.chebyshevu_root
   :members:

Legendre Polynomials
++++++++++++++++++++

.. autoclass:: sympy.functions.special.polynomials.legendre
   :members:

.. autoclass:: sympy.functions.special.polynomials.assoc_legendre
   :members:

Hermite Polynomials
+++++++++++++++++++

.. autoclass:: sympy.functions.special.polynomials.hermite
   :members:

.. autoclass:: sympy.functions.special.polynomials.hermite_prob
   :members:

Laguerre Polynomials
++++++++++++++++++++

.. autoclass:: sympy.functions.special.polynomials.laguerre
   :members:
.. autoclass:: sympy.functions.special.polynomials.assoc_laguerre
   :members:

Spherical Harmonics
-------------------

.. autoclass:: sympy.functions.special.spherical_harmonics.Ynm
   :members:

.. autofunction:: sympy.functions.special.spherical_harmonics.Ynm_c

.. autoclass:: sympy.functions.special.spherical_harmonics.Znm
   :members:

Tensor Functions
----------------

.. autofunction:: sympy.functions.special.tensor_functions.Eijk

.. autofunction:: sympy.functions.special.tensor_functions.eval_levicivita

.. autoclass:: sympy.functions.special.tensor_functions.LeviCivita
   :members:

.. autoclass:: sympy.functions.special.tensor_functions.KroneckerDelta
   :members:
