=======
Control
=======

Currently, :mod:`sympy.physics.control` is able to deal with LTI
(Linear, time-invariant) systems. The ``TransferFunction`` class is used to
represent Continuous-time Transfer functions in the Laplace domain; where Transfer
functions are input to output representations of dynamic systems. The additive
property is used for transfer functions in the ``Parallel`` class, and the
multiplicative property is used for transfer functions in the ``Series`` class.
Also, there is a ``Feedback`` class which is used to represent negative feedback
interconnection between two input/output systems.

The advantage of this symbolic Control system package is that the solutions obtained
from it are highly accurate and do not rely on numerical methods to approximate the
solutions. Symbolic solutions obtained are also in a compact form that can be used for
further analysis.

Transfer Function
=================

In :mod:`sympy.physics.control`, a transfer function of a linear, time-invariant (LTI) system is
defined as the ratio of the Laplace transform of the output (response function) to the
Laplace transform of the input (driving function) under the assumption that all initial
conditions are zero. The transfer function is a property of a system itself, independent of
the magnitude and nature of the input or driving function, and if it is known, the output
or response can be studied for various forms of inputs with a view toward understanding the
nature of the system. Here, the LTI systems can be strictly described by ratio of polynomials
in the Laplace Transform complex variable. To form a transfer function from given numerator and
denominator polynomials, we can do: ::

  >>> from sympy.physics.control import *
  >>> from sympy import Symbol
  >>> s = Symbol('s')  # Laplace-transform complex variable
  >>> num_poly = s - 1
  >>> den_poly = s**2 + s + 1
  >>> G = TransferFunction(num_poly, den_poly, s)
  >>> pprint(G)
    s - 1
  ──────────
   2
  s  + s + 1

Block Diagram Algebra
=====================

A block diagram of a system is a pictorial representation of the functions performed by
each component and of the flow of signals. It has the advantage of indicating more
realistically the signal flows of the actual system. In a block diagram all system variables
are linked to each other through blocks; where a block is basically a symbol for the mathematical
operation on the input signal to the block which produces the output. Here, the transfer functions
of the components are usually entered in the corresponding blocks, which are then connected by arrows
to indicate the direction of the flow of signals.

And block diagram algebra is the algebra involved with the basic elements of the block diagram, and
deals with the pictorial representation of algebraic equations.

There are three basic types of interconnection between two blocks: Series (Cascade), Parallel, and
Feedback.
