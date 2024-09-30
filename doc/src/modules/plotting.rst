========
Plotting
========

.. module:: sympy.plotting.plot

Introduction
------------

The plotting module allows you to make 2-dimensional and 3-dimensional plots.
Presently the plots are rendered using :external:mod:`matplotlib` as a
backend. It is also possible to plot 2-dimensional plots using a
:class:`~.TextBackend` if you do not have :external:mod:`matplotlib`.

The plotting module has the following functions:

* :func:`~.plot`: Plots 2D line plots.
* :func:`~.plot_parametric`: Plots 2D parametric plots.
* :func:`~.plot_implicit`: Plots 2D implicit and region plots.
* :func:`~.plot3d`: Plots 3D plots of functions in two variables.
* :func:`~.plot3d_parametric_line`: Plots 3D line plots, defined by a parameter.
* :func:`~.plot3d_parametric_surface`: Plots 3D parametric surface plots.

The above functions are only for convenience and ease of use. It is possible to
plot any plot by passing the corresponding ``Series`` class to :class:`~.Plot` as
argument.

Plot Class
----------

.. autoclass:: sympy.plotting.plot::Plot
   :members:

Plotting Function Reference
---------------------------

.. autofunction:: plot

.. autofunction:: plot_parametric

.. autofunction:: plot3d

.. autofunction:: plot3d_parametric_line

.. autofunction:: plot3d_parametric_surface

.. autofunction:: sympy.plotting.plot_implicit::plot_implicit

PlotGrid Class
--------------

.. autoclass:: sympy.plotting.plot::PlotGrid
   :members:

Series Classes
--------------

.. autoclass:: sympy.plotting.series::BaseSeries
   :members:

.. autoclass:: sympy.plotting.series::Line2DBaseSeries
   :members:

.. autoclass:: sympy.plotting.series::LineOver1DRangeSeries
   :members:

.. autoclass:: sympy.plotting.series::Parametric2DLineSeries
   :members:

.. autoclass:: sympy.plotting.series::Line3DBaseSeries
   :members:

.. autoclass:: sympy.plotting.series::Parametric3DLineSeries
   :members:

.. autoclass:: sympy.plotting.series::SurfaceBaseSeries
   :members:

.. autoclass:: sympy.plotting.series::SurfaceOver2DRangeSeries
   :members:

.. autoclass:: sympy.plotting.series::ParametricSurfaceSeries
   :members:

.. autoclass:: sympy.plotting.series::ImplicitSeries
   :members:

Backends
--------

.. autoclass:: sympy.plotting.plot::MatplotlibBackend
   :members:

.. autoclass:: sympy.plotting.plot::TextBackend
   :members:

Pyglet Plotting
---------------

.. module:: sympy.plotting.pygletplot

This is the documentation for the old plotting module that uses pyglet.
This module has some limitations and is not actively developed anymore.
For an alternative you can look at the new plotting module.

The pyglet plotting module can do nice 2D and 3D plots that can be
controlled by console commands as well as keyboard and mouse, with
the only dependency being `pyglet <https://pyglet.org/>`_.

Here is the simplest usage:

    >>> from sympy import var
    >>> from sympy.plotting.pygletplot import PygletPlot as Plot
    >>> var('x y z')
    >>> Plot(x*y**3-y*x**3)

To see lots of plotting examples, see ``examples/pyglet_plotting.py`` and try running
it in interactive mode (``python -i plotting.py``)::

    $ python -i examples/pyglet_plotting.py

And type for instance ``example(7)`` or ``example(11)``.

See also the `Plotting Module <https://github.com/sympy/sympy/wiki/Plotting-capabilities>`_
wiki page for screenshots.


Plot Window Controls
--------------------

======================   ========
Camera                   Keys
======================   ========
Sensitivity Modifier     SHIFT
Zoom                     R and F, Page Up and Down, Numpad + and -
Rotate View X,Y axis     Arrow Keys, A,S,D,W, Numpad 4,6,8,2
Rotate View Z axis       Q and E, Numpad 7 and 9
Rotate Ordinate Z axis   Z and C, Numpad 1 and 3
View XY                  F1
View XZ                  F2
View YZ                  F3
View Perspective         F4
Reset                    X, Numpad 5
======================   ========

======================   ========
Axes                     Keys
======================   ========
Toggle Visible           F5
Toggle Colors            F6
======================   ========

======================   ========
Window                   Keys
======================   ========
Close                    ESCAPE
Screenshot               F8
======================   ========

The mouse can be used to rotate, zoom, and translate by dragging the left, middle,
and right mouse buttons respectively.

Coordinate Modes
----------------

``Plot`` supports several curvilinear coordinate modes, and they are independent
for each plotted function. You can specify a coordinate mode explicitly with
the 'mode' named argument, but it can be automatically determined for cartesian
or parametric plots, and therefore must only be specified for polar,
cylindrical, and spherical modes.

Specifically, ``Plot(function arguments)`` and ``Plot.__setitem__(i, function
arguments)`` (accessed using array-index syntax on the ``Plot`` instance) will
interpret your arguments as a cartesian plot if you provide one function and a
parametric plot if you provide two or three functions. Similarly, the arguments
will be interpreted as a curve is one variable is used, and a surface if two
are used.

Supported mode names by number of variables:

* 1 (curves): parametric, cartesian, polar
* 2 (surfaces): parametric, cartesian, cylindrical, spherical

::

    >>> Plot(1, 'mode=spherical; color=zfade4')

Note that function parameters are given as option strings of the form
``"key1=value1; key2 = value2"`` (spaces are truncated). Keyword arguments given
directly to plot apply to the plot itself.

Specifying Intervals for Variables
----------------------------------

The basic format for variable intervals is [var, min, max, steps]. However, the
syntax is quite flexible, and arguments not specified are taken from the
defaults for the current coordinate mode:

    >>> Plot(x**2) # implies [x,-5,5,100]
    >>> Plot(x**2, [], []) # [x,-1,1,40], [y,-1,1,40]
    >>> Plot(x**2-y**2, [100], [100]) # [x,-1,1,100], [y,-1,1,100]
    >>> Plot(x**2, [x,-13,13,100])
    >>> Plot(x**2, [-13,13]) # [x,-13,13,100]
    >>> Plot(x**2, [x,-13,13]) # [x,-13,13,100]
    >>> Plot(1*x, [], [x], 'mode=cylindrical') # [unbound_theta,0,2*Pi,40], [x,-1,1,20]

Using the Interactive Interface
-------------------------------
::

    >>> p = Plot(visible=False)
    >>> f = x**2
    >>> p[1] = f
    >>> p[2] = f.diff(x)
    >>> p[3] = f.diff(x).diff(x)
    >>> p
    [1]: x**2, 'mode=cartesian'
    [2]: 2*x, 'mode=cartesian'
    [3]: 2, 'mode=cartesian'
    >>> p.show()
    >>> p.clear()
    >>> p
    <blank plot>
    >>> p[1] =  x**2+y**2
    >>> p[1].style = 'solid'
    >>> p[2] = -x**2-y**2
    >>> p[2].style = 'wireframe'
    >>> p[1].color = z, (0.4,0.4,0.9), (0.9,0.4,0.4)
    >>> p[1].style = 'both'
    >>> p[2].style = 'both'
    >>> p.close()

Using Custom Color Functions
----------------------------

The following code plots a saddle and color it by the magnitude of its gradient:

    >>> fz = x**2-y**2
    >>> Fx, Fy, Fz = fz.diff(x), fz.diff(y), 0
    >>> p[1] = fz, 'style=solid'
    >>> p[1].color = (Fx**2 + Fy**2 + Fz**2)**(0.5)

The coloring algorithm works like this:

#. Evaluate the color function(s) across the curve or surface.
#. Find the minimum and maximum value of each component.
#. Scale each component to the color gradient.

When not specified explicitly, the default color gradient is
$f(0.0)=(0.4,0.4,0.4) \rightarrow f(1.0)=(0.9,0.9,0.9)$. In our case, everything is
gray-scale because we have applied the default color gradient uniformly for
each color component. When defining a color scheme in this way, you might want
to supply a color gradient as well:

    >>> p[1].color = (Fx**2 + Fy**2 + Fz**2)**(0.5), (0.1,0.1,0.9), (0.9,0.1,0.1)

Here's a color gradient with four steps:

    >>> gradient = [ 0.0, (0.1,0.1,0.9), 0.3, (0.1,0.9,0.1),
    ...              0.7, (0.9,0.9,0.1), 1.0, (1.0,0.0,0.0) ]
    >>> p[1].color = (Fx**2 + Fy**2 + Fz**2)**(0.5), gradient

The other way to specify a color scheme is to give a separate function for each
component r, g, b. With this syntax, the default color scheme is defined:

    >>> p[1].color = z,y,x, (0.4,0.4,0.4), (0.9,0.9,0.9)

This maps z->red, y->green, and x->blue. In some cases, you might prefer to use
the following alternative syntax:

    >>> p[1].color = z,(0.4,0.9), y,(0.4,0.9), x,(0.4,0.9)

You can still use multi-step gradients with three-function color schemes.

.. _plot_geom:

Plotting Geometric Entities
---------------------------

The plotting module is capable of plotting some 2D geometric entities like
line, circle and ellipse. The following example plots a circle centred at
origin and of radius 2 units.
::

    >>> from sympy import *
    >>> x,y = symbols('x y')
    >>> plot_implicit(Eq(x**2+y**2, 4))

Similarly, :func:`~.plot_implicit` may be used to plot any 2-D geometric structure from
its implicit equation.

Plotting polygons (Polygon, RegularPolygon, Triangle) are not supported
directly.

Plotting with ASCII art
-----------------------

.. autofunction:: sympy.plotting.textplot::textplot
