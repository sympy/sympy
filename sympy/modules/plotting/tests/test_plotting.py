import sys
sys.path.append(".")

disabled = False
try:
    from ctypes import *
except:
    disabled = True

from sympy import Symbol, log, sin, cos
x,y = Symbol('x'), Symbol('y')

class TestPlotting:
    def __init__(self):
        global disabled
        self.disabled = disabled

    def test_import(self):
        from sympy.modules.plotting import Plot

    def test_plot_2d(self):
        from sympy.modules.plotting import Plot
        Plot(x, [x, -5, 5, 10], show=False)
        
    def test_plot_2d_discontinuous(self):
        from sympy.modules.plotting import Plot
        Plot(1/x, [x, -1, 1, 2], show=False)
    
    def test_plot_3d(self):
        from sympy.modules.plotting import Plot
        Plot(x*y, [x, -5, 5, 10], [y, -5, 5, 10], show=False)

    def test_plot_3d_discontinuous(self):
        from sympy.modules.plotting import Plot
        Plot(1/x, [x, -3, 3], [y, -1, 1, 1], show=False)

    def _test_plot_2d_polar(self):
        from sympy.modules.plotting import Plot
        Plot(log(x), [x,0,6.282], 'mode=polar', show=False)
        Plot(1/x, [x,-1,1,4], 'mode=polar', show=False)

    def test_plot_3d_cylinder(self):
        from sympy.modules.plotting import Plot
        Plot(1/y, [x,0,6.282], [y,-1,1,6], 'mode=polar', show=False)

    def test_plot_3d_spherical(self):
        from sympy.modules.plotting import Plot
        Plot(1, [x,0,6.282], [y,0,3.141], 'mode=spherical', show=False)

    def test_plot_2d_parametric(self):
        from sympy.modules.plotting import Plot
        Plot(sin(x), cos(x), [x, 0, 6.282], 'mode=parametric', show=False)

    def test_plot_3d_parametric(self):
        from sympy.modules.plotting import Plot
        Plot(sin(x), cos(x), x/5.0, [x, 0, 6.282], 'mode=parametric', show=False)

    def test_plot_grid(self):
        from sympy.modules.plotting import Plot
        Plot(x, [x, -5, 5, 10], grid='xy', show=False)
