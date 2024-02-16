====================
Control System Plots
====================

This module contains plotting functions for some of the common plots
used in control system. Matplotlib will be required as an external dependency
if the user wants the plots. To get only the numerical data of the plots,
NumPy will be required as external dependency.

Pole-Zero Plot
--------------

.. automethod:: sympy.physics.control.control_plots.pole_zero_plot

.. automethod:: sympy.physics.control.control_plots.pole_zero_numerical_data

Bode Plot
---------

.. automethod:: sympy.physics.control.control_plots.bode_plot

.. automethod:: sympy.physics.control.control_plots.bode_magnitude_plot

.. automethod:: sympy.physics.control.control_plots.bode_phase_plot

.. automethod:: sympy.physics.control.control_plots.bode_magnitude_numerical_data

.. automethod:: sympy.physics.control.control_plots.bode_phase_numerical_data

Impulse-Response Plot
---------------------

.. automethod:: sympy.physics.control.control_plots.impulse_response_plot

.. automethod:: sympy.physics.control.control_plots.impulse_response_numerical_data

Step-Response Plot
------------------

.. automethod:: sympy.physics.control.control_plots.step_response_plot

.. automethod:: sympy.physics.control.control_plots.step_response_numerical_data

Ramp-Response Plot
------------------

.. automethod:: sympy.physics.control.control_plots.ramp_response_plot

.. automethod:: sympy.physics.control.control_plots.ramp_response_numerical_data
