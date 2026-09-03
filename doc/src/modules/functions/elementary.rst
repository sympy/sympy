.. _elementary-functions:

Elementary
==========

This module implements elementary functions such as trigonometric, hyperbolic, and
sqrt, as well as functions like ``Abs``, ``Max``, ``Min`` etc.


Complex Functions
-----------------

.. autoclass:: sympy.functions.elementary.complexes.re
   :members:

.. autoclass:: sympy.functions.elementary.complexes.im
   :members:

.. autoclass:: sympy.functions.elementary.complexes.sign
   :members:

.. autoclass:: sympy.functions.elementary.complexes.Abs
   :members:

.. autoclass:: sympy.functions.elementary.complexes.arg
   :members:

.. autoclass:: sympy.functions.elementary.complexes.conjugate
   :members:

.. autoclass:: sympy.functions.elementary.complexes.polar_lift
   :members:

.. autoclass:: sympy.functions.elementary.complexes.periodic_argument
   :members:

.. autoclass:: sympy.functions.elementary.complexes.principal_branch
   :members:


Trigonometric
-------------

.. _trionometric functions:

Trigonometric Functions
~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: sympy.functions.elementary.trigonometric.sin
   :members:

.. autoclass:: sympy.functions.elementary.trigonometric.cos
   :members:

.. autoclass:: sympy.functions.elementary.trigonometric.tan
   :members:

.. autoclass:: sympy.functions.elementary.trigonometric.cot
   :members:

.. autoclass:: sympy.functions.elementary.trigonometric.sec
   :members:

.. autoclass:: sympy.functions.elementary.trigonometric.csc
   :members:

.. autoclass:: sympy.functions.elementary.trigonometric.sinc
   :members:


Trigonometric Inverses
~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: sympy.functions.elementary.trigonometric.asin
   :members:

.. autoclass:: sympy.functions.elementary.trigonometric.acos
   :members:

.. autoclass:: sympy.functions.elementary.trigonometric.atan
   :members:

.. autoclass:: sympy.functions.elementary.trigonometric.acot
   :members:

.. autoclass:: sympy.functions.elementary.trigonometric.asec
   :members:

.. autoclass:: sympy.functions.elementary.trigonometric.acsc
   :members:

.. autoclass:: sympy.functions.elementary.trigonometric.atan2
   :members:


Hyperbolic
----------

Hyperbolic Functions
~~~~~~~~~~~~~~~~~~~~


.. autoclass:: sympy.functions.elementary.hyperbolic.HyperbolicFunction
   :members:

.. autoclass:: sympy.functions.elementary.hyperbolic.sinh
   :members:

.. autoclass:: sympy.functions.elementary.hyperbolic.cosh
   :members:

.. autoclass:: sympy.functions.elementary.hyperbolic.tanh
   :members:

.. autoclass:: sympy.functions.elementary.hyperbolic.coth
   :members:

.. autoclass:: sympy.functions.elementary.hyperbolic.sech
   :members:

.. autoclass:: sympy.functions.elementary.hyperbolic.csch
   :members:


Hyperbolic Inverses
~~~~~~~~~~~~~~~~~~~

.. autoclass:: sympy.functions.elementary.hyperbolic.asinh
   :members:

.. autoclass:: sympy.functions.elementary.hyperbolic.acosh
   :members:

.. autoclass:: sympy.functions.elementary.hyperbolic.atanh
   :members:

.. autoclass:: sympy.functions.elementary.hyperbolic.acoth
   :members:

.. autoclass:: sympy.functions.elementary.hyperbolic.asech
   :members:

.. autoclass:: sympy.functions.elementary.hyperbolic.acsch
   :members:

Integer Functions
-----------------

.. autoclass:: sympy.functions.elementary.integers.ceiling
   :members:

.. autoclass:: sympy.functions.elementary.integers.floor
   :members:

.. autoclass:: sympy.functions.elementary.integers.RoundFunction
   :members:

.. autoclass:: sympy.functions.elementary.integers.frac
   :members:

Exponential
-----------

.. autoclass:: sympy.functions.elementary.exponential.exp
   :members:

.. autoclass:: sympy.functions.elementary.exponential.LambertW
   :members:

.. autoclass:: sympy.functions.elementary.exponential.log
   :members:

.. autoclass:: sympy.functions.elementary.exponential.exp_polar
   :members:

.. autoclass:: sympy.functions.elementary.exponential.EML
   :members:

The binary primitive :class:`~sympy.functions.elementary.exponential.EML`,
defined by :math:`\operatorname{EML}(x, y) = e^{x} - \log(y)`, follows
Odrzywolek (arXiv:2603.21852): together with the constant ``1`` and the free
symbols, it can express every elementary (calculator) function through the
grammar :math:`S \to 1 \mid x \mid \operatorname{EML}(S, S)`.  Use
:func:`~sympy.functions.elementary.exponential.to_eml` to rewrite the
elementary functions of an expression into ``EML`` form and
:func:`~sympy.functions.elementary.exponential.from_eml` to expand them back
to ``exp``/``log``.  Note that ``EML`` generates ``exp``, ``log`` and
subtraction but not the addition or multiplication of independent terms, so
the surrounding arithmetic of an expression is preserved rather than being
folded into a single ``EML`` tree.

.. autofunction:: sympy.functions.elementary.exponential.to_eml

.. autofunction:: sympy.functions.elementary.exponential.from_eml


Piecewise
---------

.. autoclass:: sympy.functions.elementary.piecewise.ExprCondPair
   :members:

.. autoclass:: sympy.functions.elementary.piecewise.Piecewise
   :members:

   .. automethod:: sympy.functions.elementary.piecewise.Piecewise._eval_integral

.. autofunction:: sympy.functions.elementary.piecewise.piecewise_exclusive

.. autofunction:: sympy.functions.elementary.piecewise.piecewise_fold


Miscellaneous
-------------

.. autoclass:: sympy.functions.elementary.miscellaneous.IdentityFunction
   :members:

.. autoclass:: sympy.functions.elementary.miscellaneous.Min
   :members:

.. autoclass:: sympy.functions.elementary.miscellaneous.Max
   :members:

.. autofunction:: sympy.functions.elementary.miscellaneous.root

.. autofunction:: sympy.functions.elementary.miscellaneous.sqrt

.. autofunction:: sympy.functions.elementary.miscellaneous.cbrt

.. autofunction:: sympy.functions.elementary.miscellaneous.real_root
