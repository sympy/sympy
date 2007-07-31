"""
Plotting Examples

Note: In Python < 2.5, you will need the ctypes library
to use plotting. It is included with Python 2.5 and later.
"""
import sys
sys.path.append("..")

from sympy import symbols
from sympy import Plot
from sympy import sin, cos, pi, sqrt

import cProfile
from pstats import Stats

if __name__ == "__main__":
    x, y = symbols('xy')
    p = Plot(width=600, height=500, ortho=False)
    a,b = 1, 0.5 # radius, thickness for p[10] (a torus
                 # which takes forever to render)

    #p[1] = x**2
    #p[2] = x**2+y**2, [x,-1,1,20], [y,-1,1,20]
    #p[3] =-x**2-y**2, [x,-1,1,20], [y,-1,1,20]
    #p[4] = x**2-y**2
    #p[5] = y**2-x**2
    #p[5] = 1, 'mode=polar'
    p[6] = sin(4*x), 'mode=polar;color=0.5+t/2.0, 0.5+x/2.0, 0.5+y/2.0'
    #p[7] = 1, 'mode=cylindrical'
    #p[8] = 1/y, 'mode=cylindrical;style=wireframe', [x,0,2*pi,20], [y,-1,1,10]
    #p[9] = sqrt(1-y)*y, [x,0,2*pi,20], [y,-1,4,80], 'mode=cylindrical' # ding-dong surface
    #p[10] = (a+b*cos(x))*cos(y), (a+b*cos(x))*sin(y), b*sin(x), [x,0,pi*2,8], [y,0,pi*2,16]
    #p[11] = cos(y), sin(y), y/10.0 # (see issue 273)
    #p[12] = 2, 'mode=polar; color=x,y,t'
    #p[13] = 1, 'mode=spherical; color=0.5+x/2.0,0.5+y/2.0,0.5+z/2.0; style=solid'

    #cProfile.run("p.append(1, 'mode=polar')", 'plot.profile2')
    #cProfile.run("p.append(x**2+y**2)", 'plot.profile2')
    #s = Stats('plot.profile2')
    #s.sort_stats('cumulative').print_stats(20)
