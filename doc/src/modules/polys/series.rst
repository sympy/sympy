.. _polys-series:

=================
Power Series Ring
=================

This module provides tools for creating and performing arithmetic on univariate
power series. It supports rings over the integer :ref:`ZZ` and rational :ref:`QQ`
domains.

A power series is represented as a finite sequence of coefficients up to a
specified precision. For a ring with precision `prec`, a series is represented by
its terms up to the degree ``prec - 1``. Arithmetic operations that would result in
terms of degree `prec` or higher are truncated, and this is denoted using the
`O(x**prec)` notation.

A key feature of this implementation is its handling of precision. If an
operation performed using two polynomials results in a polynomial whose degree is
strictly less than the ring's precision, the result is returned as an exact polynomial
without any truncation. This ensures that calculations with finite polynomials
are mathematically exact, preserving accuracy while leveraging the series
framework.

The architecture includes two backends: a pure Python implementation and a
high-performance backend using the ``FLINT`` library. The system automatically
selects the ``FLINT`` backend if it is available, falling back to the Python
implementation otherwise. The primary factory function,
:func:`~sympy.polys.series.power_series_ring`, handles this selection
transparently.

To create a power series ring, use the :func:`~sympy.polys.series.power_series_ring`
function, specifying a domain and the desired precision. Ring elements can then be
created by passing a list of coefficients.

For example, let's create a ring with precision 8 over the integers:

    >>> from sympy.polys.series import power_series_ring
    >>> from sympy import ZZ
    >>> R = power_series_ring(ZZ, prec=8)
    >>> f = R([1, 2, 3])
    >>> g = R([4, 1])

Arithmetic operations are performed using the ring's methods. The print method
provides a readable string representation.

    >>> R.print(R.add(f, g))
    5 + 3*x + 3*x**2
    >>> R.print(R.multiply(f, g))
    4 + 9*x + 14*x**2 + 3*x**3

As shown below, when the result of an operation is a polynomial with a degree
less than the ring's precision, the exact result is returned.

    >>> p = R([2, -3])
    >>> q = R([0, 7, 6, 1])
    >>> R.print(R.multiply(p, q))
    14*x - 9*x**2 - 16*x**3 - 3*x**4

However, if an operation produces a result that exceeds the precision threshold,
it is automatically truncated.

    >>> r = R([1, 2, 3, 4, 5, 6, 7, 8])
    >>> s = R([0, 1, 1])
    >>> R.print(R.multiply(r, s))
    x + 3*x**2 + 5*x**3 + 7*x**4 + 9*x**5 + 11*x**6 + 13*x**7 + O(x**8)

.. autofunction:: sympy.polys.series.power_series_ring

.. py:class:: sympy.polys.series.ring.TSeries


Protocols for Power Series Rings
================================

.. autoclass:: sympy.polys.series.base.PowerSeriesRingProto

.. py:class:: TSeries


Python Implementation
=====================

.. currentmodule:: sympy.polys.series.ringpython

.. autoclass:: PythonPowerSeriesRingZZ
    :members:

.. autoclass:: PythonPowerSeriesRingQQ
    :members:

.. py:class:: USeries

Flint Implementation
====================

.. currentmodule:: sympy.polys.series.ringflint

.. autoclass:: FlintPowerSeriesRingZZ
    :members:

.. autoclass:: FlintPowerSeriesRingQQ
    :members:

.. py:class:: ZZSeries
.. py:class:: QQSeries

.. py:class:: flint.types.fmpz_poly.fmpz_poly
.. py:class:: flint.types.fmpq_poly.fmpq_poly
.. py:class:: flint.types.fmpz_series.fmpz_series
.. py:class:: flint.types.fmpq_series.fmpq_series
