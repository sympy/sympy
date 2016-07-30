.. _combinatorics-fp_groups

Finitely Presented Group
========================

.. module:: sympy.combinatorics.fp_groups

This module contains the implementation of Finitely Presented Groups and the algorithms for computation
of its attributes. Algorithms include Coset Enumeration using the Todd Coxeter algorithm, Low Index Subgroups
algorithm, Subgroup Presentations using Reidemeister algorithm.

.. automethod:: fp_group

.. automethod:: xfp_group

.. automethod:: vfp_group

.. autoclass:: FpGroup

.. autoclass:: CosetTable

.. automethod:: coset_enumeration_r

.. automethod:: coset_enumeration_c

.. automethod:: low_index_subgroups

.. automethod:: descendant_subgroups

.. automethod:: try_descendant

.. automethod:: first_in_class

.. automethod:: define_schreier_generators

.. automethod:: reidemeister_relators

.. automethod:: rewrite

.. automethod:: elimination_technique_1

.. automethod:: elimination_technique_2

.. automethod:: simplify_presentation

.. automethod:: reidemeister_presentation
