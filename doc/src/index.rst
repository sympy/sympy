.. _documentation:

.. SymPy documentation master file, created by sphinx-quickstart.py on Sat Mar 22 19:34:32 2008.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SymPy's documentation!
----------------------------------

A PDF version of these docs can be found `here <https://github.com/sympy/sympy/releases>`_.

`SymPy <https://sympy.org>`_ is a Python library for symbolic mathematics.
If you are new to SymPy, start with the :ref:`Tutorial <tutorial>`.

This is the central page for all of SymPy's documentation.

==================================================================

..  rst-class:: column

:ref:`SymPy Tutorials <tutorial>`
-------------------------------------------------------

*Master SymPy concepts and features*


..  rst-class:: column

:ref:`How-to Guides <guides>`
-------------------------------------------------------

*Step-by-step intructions on how to do different key developer tasks*


..  rst-class:: column clearfix

:ref:`Explanation <explanation>`
-------------------------------------------------------

*Common pitfalls and advanced topics*


..  rst-class:: column

:ref:`SymPy API Reference <reference>`
-------------------------------------------------------

*SymPy internal and public modules*


..  rst-class:: clearfix row custom-headings

.. autosummary::
   :toctree: reference/private
   :template: custom-module-template.rst
   :recursive:

   sympy.abc

.. toctree::
   :hidden:

   tutorial/index.rst
   explanation/index.rst
   reference/index.rst
   guides/index.rst
   miscellaneous/index.rst
