.. _manualintegrate-rules:

=========================
Rules of manualintegrate
=========================

:func:`~sympy.integrals.manualintegrate.integral_steps` returns a tree of
rules. Each rule rewrites an integral into a result, or into simpler
integrals that are solved by its substeps; ``eval()`` evaluates the tree into
the antiderivative. The table is generated from the docstrings of the rule
classes, which are listed below it.

.. manualintegrate-rules::
