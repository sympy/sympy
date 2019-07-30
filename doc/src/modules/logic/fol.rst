First Order Logic Module
========================

.. module:: sympy.logic.FOL

Introduction
------------

The First Order Logic module for SymPy facilitates creating, manipulating and inferring using First Order Predicate Calculus.

First Order Logic Classes
-------------------------

.. autoclass:: sympy.logic.FOL.Predicate

.. autoclass:: sympy.logic.FOL.Function

.. autoclass:: sympy.logic.FOL.Constant

.. autoclass:: sympy.logic.FOL.ForAll

.. autoclass:: sympy.logic.FOL.Exists


First Order Logic Functions
---------------------------

The following function can be used for obtaining the truth value under a given interpretation.

.. autofunction:: sympy.logic.FOL.fol_true

The following functions can be used to convert expressions to normal forms.

.. autofunction:: sympy.logic.FOL.to_pnf

.. autofunction:: sympy.logic.FOL.to_snf

The following functions can be used for Inference

.. autofunction:: sympy.logic.FOL.mgu

.. autofunction:: sympy.logic.FOL.resolve

.. autofunction:: sympy.logic.FOL.entails

.. autoclass:: sympy.logic.FOL.FOL_KB
