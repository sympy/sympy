Discrete Module
===============

The ``discrete`` module in SymPy implements methods to compute discrete
transforms and convolutions of finite sequences.

.. automodule:: sympy.discrete

Since the discrete transforms can be used to reduce the computational complexity
of the discrete convolutions, the ``convolutions`` module makes use of the
``transforms`` module for efficient computation (notable for long input sequences).

Transforms
----------

.. module:: sympy.discrete.transforms

This section lists the methods which implement the basic transforms
for discrete sequences.

Fast Fourier Transform
^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: fft

.. autofunction:: ifft

Number Theoretic Transform
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: ntt

.. autofunction:: intt

Fast Walsh Hadamard Transform
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: fwht

.. autofunction:: ifwht

Möbius Transform
^^^^^^^^^^^^^^^^

.. autofunction:: mobius_transform

.. autofunction:: inverse_mobius_transform


Convolutions
------------

.. module:: sympy.discrete.convolutions

This section lists the methods which implement the basic convolutions
for discrete sequences.

Convolution
^^^^^^^^^^^

This is a general method for calculating the convolution of discrete
sequences, which internally calls one of the methods ``convolution_fft``,
``convolution_ntt``, ``convolution_fwht``, or ``convolution_subset``.

.. autofunction:: convolution

Convolution using Fast Fourier Transform
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: convolution_fft

Convolution using Number Theoretic Transform
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: convolution_ntt

Convolution using Fast Walsh Hadamard Transform
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: convolution_fwht

Subset Convolution
^^^^^^^^^^^^^^^^^^

.. autofunction:: convolution_subset

Covering Product
^^^^^^^^^^^^^^^^

.. autofunction:: covering_product

Intersecting Product
^^^^^^^^^^^^^^^^^^^^

.. autofunction:: intersecting_product
