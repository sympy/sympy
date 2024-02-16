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
interconnection between two input/output systems. MIMO systems are also supported
with ``TransferFunctionMatrix`` as the base class for representing one. ``MIMOSeries``,
``MIMOParallel``  and ``MIMOFeedback`` are MIMO equivalent of ``Series``, ``Parallel``
and ``Feedback`` classes.

The advantage of this symbolic Control system package is that the solutions obtained
from it are highly accurate and do not rely on numerical methods to approximate the
solutions. Symbolic solutions obtained are also in a compact form that can be used for
further analysis.
