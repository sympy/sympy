===========
SymPy Logo
===========

We would like to make it easy for you to include the SymPy project identity in
your next academic paper, course materials, or presentation.

.. image::  ../logo/sympy.svg
  :width: 600
  :align: center
  :alt: SymPy Logo

The above image logo is an SVG version of the logo. We also have a PNG version of the logo:

.. image::  ../../_build/logo/sympy-500px.png
  :width: 500
  :align: center
  :alt: SymPy Logo

If you would like one without the "SymPy" text, we have that too:

.. image::  ../../_build/logo/sympy-notext-500px.png
  :width: 500
  :align: center
  :alt: SymPy Logo

Note: The text version should be preferred unless the "SymPy" name is already present separately.

If you would like to generate SymPy's collection of official logos yourself,
you can do so by first :ref:`installing the required dependencies <build-the-documentation>`, and then running:

.. code-block:: none

    $ cd doc

    $ make logo # will be stored in the _build/logo subdirectory

which will generate the logos by using the ``sympy.svg`` file in your local
copy of SymPy.

There is also a ``sympy/doc/generate_logos.py`` script that allows for a wider
variety of options while generating the logo.

The license of all the logos is the same as SymPy: BSD. See the
`LICENSE file <https://github.com/sympy/sympy/blob/master/LICENSE>`_ for more information.
