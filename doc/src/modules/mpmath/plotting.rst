Plotting
========

If `matplotlib <http://matplotlib.sourceforge.net/>`_ is available, the functions ``plot`` and ``cplot`` in mpmath can be used to plot functions respectively as x-y graphs and in the complex plane. Also, ``splot`` can be used to produce 3D surface plots.

Function curve plots
-----------------------

.. image:: plot.png

Output of ``plot([cos, sin], [-4, 4])``

.. autofunction:: mpmath.plot

Complex function plots
-------------------------

.. image:: cplot.png

Output of ``fp.cplot(fp.gamma, points=100000)``

.. autofunction:: mpmath.cplot

3D surface plots
----------------

.. image:: splot.png

Output of ``splot`` for the donut example.

.. autofunction:: mpmath.splot

