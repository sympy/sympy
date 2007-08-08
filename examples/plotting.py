"""
Plotting Examples

Note: In Python < 2.5, you will need the ctypes library
to use plotting. It is included with Python 2.5 and later.
"""
import sys
sys.path.append("..")

from sympy import symbols
from sympy import Plot
from sympy import sin, cos, Pi, sqrt

from time import sleep, clock

if __name__ == "__main__":

    x,y,z = symbols('xyz')
    axes_options = 'colored=true; label_ticks=true; label_axes=true; overlay=true; stride=0.5'
    #axes_options = 'colored=false; overlay=false; stride=(1.0, 0.5, 0.5)'
    #axes_options = 'none'
    p = Plot(width=600, height=600, ortho=False, invert_mouse_zoom=False, axes=axes_options)

    examples = []
    def example_wrapper(f):
        examples.append(f)
        return f

    @example_wrapper
    def ding_dong_surface():
        p[1] = sqrt(1.0-y)*y, [x,0,2*Pi,40], [y,-1,4,100], 'mode=cylindrical; style=solid'

    @example_wrapper
    def mirrored_ellipsoids():
        p[2] = x**2+y**2, 'style=solid'
        p[3] =-x**2-y**2, 'style=solid'

    @example_wrapper
    def mirrored_saddles():
        p[5] = x**2-y**2, [20], [20]
        p[6] = y**2-x**2, [20], [20]

    @example_wrapper
    def polar_circle():
        p[7] = 1, 'mode=polar'

    @example_wrapper
    def polar_flower():
        p[8] = 1.6*sin(4*x), [160], 'mode=polar'
        p[8].color = z, x, y, (0.5,0.5,0.5), (0.8,0.8,0.8), (x,y,None,z) # z is used for t

    @example_wrapper
    def simple_cylinder():
        p[9] = 1, 'mode=cylindrical'

    @example_wrapper
    def cylindrical_hyperbola():
        ## (note that polar is an alias for cylindrical)
        p[10] = 1/y, 'mode=polar', [x], [y,-2,2,20]

    @example_wrapper
    def extruded_hyperbolas():
        p[11] = 1/x, [x,-10,10,100], [1], 'style=wireframe'
        p[12] = -1/x, [x,-10,10,100], [1], 'style=solid'

    @example_wrapper
    def torus():
        a,b = 1, 0.5 # radius, thickness
        p[13] = (a+b*cos(x))*cos(y), (a+b*cos(x))*sin(y), b*sin(x), [x,0,Pi*2,40], [y,0,Pi*2,40]

    @example_wrapper
    def warped_torus():
        a,b = 2, 1 # radius, thickness
        p[13] = (a+b*cos(x))*cos(y), (a+b*cos(x))*sin(y), b*sin(x)+0.5*sin(4*y), [x,0,Pi*2,40], [y,0,Pi*2,40]

    @example_wrapper
    def parametric_spiral():
        p[14] = cos(y), sin(y), y/10.0, [y,-4*Pi,4*Pi,100]

    @example_wrapper
    def str_and_repr_demo():
        p[17] = x**2
        p[18] = eval(repr(p[17]))
        p[18].color = 0.9, 0.4, 0.4
        p[17].visible = False
        print str(p[18])
        print repr(p[18])

    @example_wrapper
    def lambda_vs_sympy_evaluation():
        start = clock()
        p[4] = x**2+y**2, [100], [100], 'style=solid'
        p.wait_for_calculations()
        print "lambda-based calculation took %s seconds." % (clock()-start)

        start = clock()
        p[4] = x**2+y**2, [100], [100], 'style=solid; use_sympy_eval'
        p.wait_for_calculations()
        print "sympy substitution-based calculation took %s seconds." % (clock()-start)

    def help_str():
        s =  ("\nPlot p has been created. Useful commands: \n"
              "    p[1] = x**2, print p, p.clear(), help(p) \n\n"
              "You can also run an example (source in plotting.py):\n\n")
        for i in xrange(len(examples)):
            s += "(%i) %s\n" % (i, examples[i].__name__)
        s += "\n"
        s += "e.g. >>> example(2)\n"
        s += "     >>> ding_dong_surface()\n"
        return s

    def example(i):
        if callable(i):
            p.clear()
            i()
        elif i >= 0 and i < len(examples):
            p.clear()
            examples[i]()
        else: print "Not a valid example.\n"
        print p

    mirrored_saddles()
    #ding_dong_surface()
    print help_str()

    #def profile_plotting():
        #import cProfile
        #from pstats import Stats
        #cProfile.run("p.append(1, 'mode=polar')", 'plot.profile2')
        #cProfile.run("p.append(x**2+y**2)", 'plot.profile2')
        #s = Stats('plot.profile2')
        #s.sort_stats('cumulative').print_stats(20)
