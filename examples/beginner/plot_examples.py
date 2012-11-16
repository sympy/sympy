#! /usr/bin/env python
# Check the plot docstring

from sympy import Symbol, exp, sin, cos
from sympy.plotting import (plot, plot_parametric,
                            plot3d_parametric_surface, plot3d_parametric_line,
                            plot3d)

lx = range(5)
ly = [i**2 for i in lx]

x = Symbol('x')
y = Symbol('y')
u = Symbol('u')
v = Symbol('v')
expr = x**2 - 1

b = plot(expr, (x, 2, 4), show=False)  # cartesian plot
e = plot(exp(-x), (x, 0, 4), show=False)  # cartesian plot (and coloring, see below)
f = plot3d_parametric_line(sin(x), cos(x), x, (x, 0, 10), show=False)  # 3d parametric line plot
g = plot3d(sin(x)*cos(y), (x, -5, 5), (y, -10, 10), show=False)  # 3d surface cartesian plot
h = plot3d_parametric_surface(cos(u)*v, sin(u)*v, u, (u, 0, 10), (v, -2, 2), show=False)  # 3d parametric surface plot

# Some aesthetics
e[0].line_color = lambda x: x / 4
f[0].line_color = lambda x, y, z: z / 10
g[0].surface_color = lambda x, y: sin(x)

# Some more stuff on aesthetics - coloring wrt coordinates or parameters
param_line_2d = plot_parametric((x*cos(x), x*sin(x), (x, 0, 15)), (1.1*x*cos(x), 1.1*x*sin(x), (x, 0, 15)), show=False)
param_line_2d[0].line_color = lambda u: sin(u)  # parametric
param_line_2d[1].line_color = lambda u, v: u**2 + v**2  # coordinates
param_line_2d.title = 'The inner one is colored by parameter and the outher one by coordinates'

param_line_3d = plot3d_parametric_line((x*cos(x), x*sin(x), x, (x, 0, 15)),
                     (1.5*x*cos(x), 1.5*x*sin(x), x, (x, 0, 15)),
                     (2*x*cos(x), 2*x*sin(x), x, (x, 0, 15)), show=False)
param_line_3d[0].line_color = lambda u: u  # parametric
param_line_3d[1].line_color = lambda u, v: u*v  # first and second coordinates
param_line_3d[2].line_color = lambda u, v, w: u*v*w  # all coordinates


if __name__ == '__main__':
    for p in [b, e, f, g, h, param_line_2d, param_line_3d]:
        p.show()
