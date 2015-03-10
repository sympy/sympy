Series Expansions
=================

.. module:: sympy.series

The series module implements series expansions as a function and many related
functions.

Limits
------

The main purpose of this module is the computation of limits.

.. autofunction:: sympy.series.limits.limit

.. autoclass:: sympy.series.limits.Limit
   :members:

As is explained above, the workhorse for limit computations is the
function gruntz() which implements Gruntz' algorithm for computing limits.

The Gruntz Algorithm
^^^^^^^^^^^^^^^^^^^^

This section explains the basics of the algorithm used for computing limits.
Most of the time the limit() function should just work. However it is still
useful to keep in mind how it is implemented in case something does not work
as expected.

First we define an ordering on functions. Suppose `f(x)` and `g(x)` are two
real-valued functions such that `\lim_{x \to \infty} f(x) = \infty` and
similarly `\lim_{x \to \infty} g(x) = \infty`. We shall say that `f(x)`
*dominates*
`g(x)`, written `f(x) \succ g(x)`, if for all `a, b \in \mathbb{R}_{>0}` we have
`\lim_{x \to \infty} \frac{f(x)^a}{g(x)^b} = \infty`.
We also say that `f(x)` and
`g(x)` are *of the same comparability class* if neither `f(x) \succ g(x)` nor
`g(x) \succ f(x)` and shall denote it as `f(x) \asymp g(x)`.

Note that whenever `a, b \in \mathbb{R}_{>0}` then
`a f(x)^b \asymp f(x)`, and we shall use this to extend the definition of
`\succ` to all functions which tend to `0` or `\pm \infty` as `x \to \infty`.
Thus we declare that `f(x) \asymp 1/f(x)` and `f(x) \asymp -f(x)`.

It is easy to show the following examples:

* `e^x \succ x^m`
* `e^{x^2} \succ e^{mx}`
* `e^{e^x} \succ e^{x^m}`
* `x^m \asymp x^n`
* `e^{x + \frac{1}{x}} \asymp e^{x + \log{x}} \asymp e^x`.

From the above definition, it is possible to prove the following property:

    Suppose `\omega`, `g_1, g_2, \dots` are functions of `x`,
    `\lim_{x \to \infty} \omega = 0` and `\omega \succ g_i` for
    all `i`. Let `c_1, c_2, \dots \in \mathbb{R}` with `c_1 < c_2 < \dots`.

    Then `\lim_{x \to \infty} \sum_i g_i \omega^{c_i} = \lim_{x \to \infty} g_1 \omega^{c_1}`.

For `g_1 = g` and `\omega` as above we also have the following easy result:

    * `\lim_{x \to \infty} g \omega^c = 0` for `c > 0`
    * `\lim_{x \to \infty} g \omega^c = \pm \infty` for `c < 0`,
      where the sign is determined by the (eventual) sign of `g`
    * `\lim_{x \to \infty} g \omega^0 = \lim_{x \to \infty} g`.


Using these results yields the following strategy for computing
`\lim_{x \to \infty} f(x)`:

1. Find the set of *most rapidly varying subexpressions* (MRV set) of `f(x)`.
   That is, from the set of all subexpressions of `f(x)`, find the elements that
   are maximal under the relation `\succ`.
2. Choose a function `\omega` that is in the same comparability class as
   the elements in the MRV set, such that `\lim_{x \to \infty} \omega = 0`.
3. Expand `f(x)` as a series in `\omega` in such a way that the antecedents of
   the above theorem are satisfied.
4. Apply the theorem and conclude the computation of
   `\lim_{x \to \infty} f(x)`, possibly by recursively working on `g_1(x)`.


Notes
"""""

This exposition glossed over several details. Many are described in the file
gruntz.py, and all can be found in Gruntz' very readable thesis. The most
important points that have not been explained are:

1. Given f(x) and g(x), how do we determine if `f(x) \succ g(x)`,
   `g(x) \succ f(x)` or `g(x) \asymp f(x)`?
2. How do we find the MRV set of an expression?
3. How do we compute series expansions?
4. Why does the algorithm terminate?

If you are interested, be sure to take a look at
`Gruntz Thesis <http://www.cybertester.com/data/gruntz.pdf>`_.

Reference
"""""""""

.. autofunction:: sympy.series.gruntz.gruntz

.. autofunction:: sympy.series.gruntz.compare

.. autofunction:: sympy.series.gruntz.rewrite

.. autofunction:: sympy.series.gruntz.build_expression_tree

.. autofunction:: sympy.series.gruntz.mrv_leadterm

.. autofunction:: sympy.series.gruntz.calculate_series

.. autofunction:: sympy.series.gruntz.limitinf

.. autofunction:: sympy.series.gruntz.sign

.. autofunction:: sympy.series.gruntz.mrv

.. autofunction:: sympy.series.gruntz.mrv_max1

.. autofunction:: sympy.series.gruntz.mrv_max3

.. autoclass:: sympy.series.gruntz.SubsSet
   :members:


More Intuitive Series Expansion
-------------------------------

This is achieved
by creating a wrapper around Basic.series(). This allows for the use of
series(x*cos(x),x), which is possibly more intuative than (x*cos(x)).series(x).

Examples
^^^^^^^^
    >>> from sympy import Symbol, cos, series
    >>> x = Symbol('x')
    >>> series(cos(x),x)
    1 - x**2/2 + x**4/24 + O(x**6)

Reference
^^^^^^^^^

.. autofunction:: sympy.series.series.series

Order Terms
-----------

This module also implements automatic keeping track of the order of your
expansion.

Examples
^^^^^^^^
     >>> from sympy import Symbol, Order
     >>> x = Symbol('x')
     >>> Order(x) + x**2
     O(x)
     >>> Order(x) + 1
     1 + O(x)

Reference
^^^^^^^^^

.. autoclass:: sympy.series.order.Order
   :members:

Series Acceleration
-------------------

TODO

Reference
^^^^^^^^^

.. autofunction:: sympy.series.acceleration.richardson

.. autofunction:: sympy.series.acceleration.shanks

Residues
--------

TODO

Reference
^^^^^^^^^

.. autofunction:: sympy.series.residues.residue
