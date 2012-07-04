"""Predefined R^n manifolds together with common coord. systems.

Coordinate systems are predefined as well as the transformation laws between
them.

Coordinate functions can be accessed as attributes of the manifold (eg `R2.x`),
as attributes of the coordinate systems (eg `R2_r.x` and `R2_p.theta`), or by
using the usual `coord_sys.coord_function(index, name)` interface.
"""

from differential_geometry import Manifold, Patch, CoordSystem
from sympy import sqrt, atan2, sin, cos, Dummy

###############################################################################
# R2
###############################################################################
R2 = Manifold('R^2', 2)
# Patch and coordinate systems.
R2_origin = Patch('R^2_o', R2)
R2_r = CoordSystem('R^2_r', R2_origin)
R2_p = CoordSystem('R^2_p', R2_origin)

# Connecting the coordinate charts.
x, y, r, theta = [Dummy(s) for s in ['x', 'y', 'r', 'theta']]
R2_r.connect_to(R2_p, [x, y],
                      [sqrt(x**2 + y**2), atan2(y, x)],
                inverse=False, fill_in_gaps=False)
R2_p.connect_to(R2_r, [r, theta],
                      [r*cos(theta), r*sin(theta)],
                inverse=False, fill_in_gaps=False)
del x, y, r, theta

# Defining the basis coordinate functions and adding shortcuts for them to the
# manifold and the patch.
## for rectangular chart
R2_r.x = R2_r.coord_function(0)
R2_origin.x = R2_r.x
R2.x = R2_r.x
R2_r.y = R2_r.coord_function(1)
R2_origin.y = R2_r.y
R2.y = R2_r.y
## for polar chart
R2_p.r = R2_p.coord_function(0)
R2_origin.r = R2_p.r
R2.r = R2_p.r
R2_p.theta = R2_p.coord_function(1)
R2_origin.theta = R2_p.theta
R2.theta = R2_p.theta

# Defining the basis vector fields and adding shortcuts for them to the
# manifold and the patch.
## for rectangular chart
R2_r.d_dx = R2_r.base_vector(0)
R2_origin.d_dx = R2_r.d_dx
R2.d_dx = R2_r.d_dx
R2_r.d_dy = R2_r.base_vector(1)
R2_origin.d_dy = R2_r.d_dy
R2.d_dy = R2_r.d_dy
## for polar chart
R2_p.d_dr = R2_p.base_vector(0)
R2_origin.d_dr = R2_p.d_dr
R2.d_dr = R2_p.d_dr
R2_p.d_dtheta = R2_p.base_vector(1)
R2_origin.d_dtheta = R2_p.d_dtheta
R2.d_dtheta = R2_p.d_dtheta
