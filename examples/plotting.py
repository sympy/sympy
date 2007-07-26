"""
Plotting Example

Note: In Python < 2.5, you will need the ctypes library
to use plotting. It is included with Python 2.5 and later.
"""
import sys
sys.path.append("..")

from time import sleep
from sympy import Symbol, sin, cos, log, pi, sqrt
from sympy.modules.plotting import Plot

if __name__ == "__main__":
    x, y = Symbol('x'), Symbol('y')

    p = Plot(width  = 300, height = 300,
             bounding_box = True, #grid = 'xy',
             wireframe = False,
             antialiasing = True,
             ortho = False)
             #ortho = True)

    #p[1]  =  1/x, [x,-3,3]
    #p[2]  =  x**2, [x,-3,3]
    #p[3]  =  x**2 + y**2, [x,-1,1,10], [y,-1,1,10]
    #p[4]  = -x**2 - y**2, [x,-1,1,10], [y,-1,1,10]
    #p[5]  = sin(x),x/10.0,cos(x), [x,-pi*4,pi*4,100], 'mode=parametric'
    #p[6]  = x*y**3-y*x**3, [x,-1.5,1.5,20], [y,-1.5,1.5,20]
    #p[7]  =  1 - x**2 + y**2, [x,-1,1], [y,-1,1]
    #p[8]  = -1 + x**2 - y**2, [x,-1,1], [y,-1,1]
    #p[9]  = 1.0, [x,0,2*pi,60], 'mode=polar'
    #p[10] = 1, [x,0,2*pi], [y,-2,2,20], 'mode=polar'
    #p[11] = -1/(x**2 + y**2), [x,-1,1], [y,-1,1]
    #p[12] = sin(x*y), [x,-6.282,6.282,50], [y,-6.282,6.282,50]
    #p[13] = 3, [x,0,2*pi,20], [y,0,pi,20], 'mode=spherical'
    #p[14] = 1/y, [x,0,6.282], [y,-0.25,0.25,4], 'mode=polar'
    p[15] = 1.5, [x,0,2*pi,20], [y,0,pi,20], 'mode=spherical'
