Category Theory Module
======================

Introduction
------------

The category theory module for SymPy will allow manipulating diagrams
within a single category, including drawing them in TikZ and deciding
whether they are commutative or not.

The general reference work this module tries to follow is

  [JoyOfCats] J. Adamek, H. Herrlich. G. E. Strecker: Abstract and
              Concrete Categories. The Joy of Cats.

The latest version of this book should be available for free download
from

   katmat.math.uni-bremen.de/acc/acc.pdf

The module is still in its pre-embryonic stage.

Base Class Reference
--------------------

.. module:: sympy.categories

This section lists the classes which implement some of the basic
notions in category theory: objects, morphisms, categories, and
diagrams.

.. autoclass:: Object
   :members:

.. autoclass:: Morphism
   :members:

.. autoclass:: NamedMorphism
   :members:

.. autoclass:: CompositeMorphism
   :members:

.. autoclass:: IdentityMorphism
   :members:

.. autoclass:: Category
   :members:

.. autoclass:: Diagram
   :members:

Diagram Drawing
---------------

.. module:: sympy.categories.diagram_drawing

This section lists the classes which allow automatic drawing of
diagrams.

.. autoclass:: DiagramGrid
   :members:

.. autoclass:: ArrowStringDescription

.. autoclass:: XypicDiagramDrawer
   :members:

.. autofunction:: xypic_draw_diagram

.. autofunction:: preview_diagram
