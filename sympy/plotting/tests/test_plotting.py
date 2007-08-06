import sys
sys.path.append(".")

disabled = False
try:
    from ctypes import *
except:
    disabled = True

from sympy import symbols, log, sin, cos
x,y = symbols('xy')

class TestPlotting:
    def __init__(self):
        global disabled
        self.disabled = disabled

    def test_import(self):
        from sympy import Plot

    def test_plot_2d(self):
        from sympy import Plot
        p=Plot(x, [x, -5, 5, 4], visible=False)
        p.wait_for_calculations()
        
    def test_plot_2d_discontinuous(self):
        from sympy import Plot
        p=Plot(1/x, [x, -1, 1, 2], visible=False)
        p.wait_for_calculations()
    
    def test_plot_3d(self):
        from sympy import Plot
        p=Plot(x*y, [x, -5, 5, 5], [y, -5, 5, 5], visible=False)
        p.wait_for_calculations()

    def test_plot_3d_discontinuous(self):
        from sympy import Plot
        p=Plot(1/x, [x, -3, 3, 6], [y, -1, 1, 1], visible=False)
        p.wait_for_calculations()

    def test_plot_2d_polar(self):
        from sympy import Plot
        p=Plot(1/x, [x,-1,1,4], 'mode=polar', visible=False)
        p.wait_for_calculations()

    def test_plot_3d_cylinder(self):
        from sympy import Plot
        p=Plot(1/y, [x,0,6.282,4], [y,-1,1,4], 'mode=polar;style=solid', visible=False)
        p.wait_for_calculations()

    def test_plot_3d_spherical(self):
        from sympy import Plot
        p=Plot(1, [x,0,6.282,4], [y,0,3.141,4], 'mode=spherical;style=wireframe', visible=False)
        p.wait_for_calculations()

    def test_plot_2d_parametric(self):
        from sympy import Plot
        p=Plot(sin(x), cos(x), [x, 0, 6.282, 4], visible=False)
        p.wait_for_calculations()

    def test_plot_3d_parametric(self):
        from sympy import Plot
        p=Plot(sin(x), cos(x), x/5.0, [x, 0, 6.282, 4], visible=False)
        p.wait_for_calculations()

    def _test_plot_log(self):
        from sympy import Plot
        p=Plot(log(x), [x,0,6.282,4], 'mode=polar', visible=False)
        p.wait_for_calculations()
