.. _control_tutorial:

=============================
Control Tutorials
=============================

The control module in SymPy provides essential tools for symbolic control system
analysis. The TransferFunction class allows for creating transfer functions and
analyzing their properties, such as stability is_stable(), poles, and zeros.
Series, Parallel, and Feedback classes enable constructing and simplifying
system interconnections. The TransferFunctionMatrix handles multi-input, multi-output
(MIMO) systems, while MIMOSeries, MIMOParallel, and MIMOFeedback extend these
functionalities for complex systems.

Additionally, the module includes the StateSpace class, which allows for modeling
control systems using state variables, inputs, and outputs in matrix form. This
representation is particularly useful for time-domain analysis and handling complex
MIMO systems.

This tutorial contains a breif guide on how to solve Control Problems using
`TransferFunction` and `StateSpace`.

.. toctree::
   :maxdepth: 1

   control_problems.rst
   electrical_problems.rst
   mechanics_problems.rst
