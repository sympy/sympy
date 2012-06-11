from sympy import (plot, pi, sin, cos, Symbol, Integral, summation, sqrt, log,
oo, LambertW, I)
from tempfile import NamedTemporaryFile
import warnings

def tmp_file(name=''):
    return NamedTemporaryFile(suffix='.png').name

def plot_and_save(name):
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    ###
    # Examples from the 'introduction' notebook
    ###

    p = plot(x, show=False)
    p = plot(x*sin(x),x*cos(x), show=False)
    p.extend(p)
    p[0].line_color = lambda a : a
    p[1].line_color='b'
    p.title = 'Big title'
    p.xlabel = 'the x axis'
    p[1].label = 'straight line'
    p.legend = True
    p.aspect_ratio = (1,1)
    p.xlim = (-15,20)
    p.save(tmp_file('%s_basic_options_and_colors.png' % name))

    p.extend(plot(x+1, show=False))
    p.append(plot((x+3,),(x**2,), show=False)[1])
    p.save(tmp_file('%s_plot_extend_append.png' % name))

    p[2] = x**2, (x, -2, 3)
    p.save(tmp_file('%s_plot_setitem.png' % name))

    p = plot(sin(x),(x,-2*pi,4*pi), show=False)
    p.save(tmp_file('%s_line_explicit.png' % name))

    p = plot(sin(x),(-2*pi,4*pi), show=False)
    p.save(tmp_file('%s_line_implicit_var.png' % name))

    p = plot(sin(x), show=False)
    p.save(tmp_file('%s_line_default_range.png' % name))

    p = plot((x,), (x*sin(x), x*cos(x)), show=False)
    p.save(tmp_file('%s_two_curves.png' % name))

    p = plot(sin(x),cos(x),x, show=False)
    p.save(tmp_file('%s_3d_line.png' % name))

    p = plot(x*y, show=False)
    p.save(tmp_file('%s_surface.png' % name))

    p = plot(x*sin(z),x*cos(z),z, show=False)
    p.save(tmp_file('%s_parametric_surface.png' % name))

    ###
    # Examples from the 'colors' notebook
    ###

    p = plot(sin(x), show=False)
    p[0].line_color = lambda a : a
    p.save(tmp_file('%s_colors_line_arity1.png' % name))

    p[0].line_color = lambda a, b : b
    p.save(tmp_file('%s_colors_line_arity2.png' % name))

    p = plot(x*sin(x), x*cos(x), (0,10), show=False)
    p[0].line_color = lambda a : a
    p.save(tmp_file('%s_colors_param_line_arity1.png' % name))

    p[0].line_color = lambda a, b : a
    p.save(tmp_file('%s_colors_param_line_arity2a.png' % name))

    p[0].line_color = lambda a, b : b
    p.save(tmp_file('%s_colors_param_line_arity2b.png' % name))

    p = plot(sin(x)+0.1*sin(x)*cos(7*x),
             cos(x)+0.1*cos(x)*cos(7*x),
             0.1*sin(7*x),
             (0, 2*pi),
             show=False)
    p[0].line_color = lambda a : sin(4*a)
    p.save(tmp_file('%s_colors_3d_line_arity1.png' % name))
    p[0].line_color = lambda a, b : b
    p.save(tmp_file('%s_colors_3d_line_arity2.png' % name))
    p[0].line_color = lambda a, b, c : c
    p.save(tmp_file('%s_colors_3d_line_arity3.png' % name))

    p = plot(sin(x)*y, (0,6*pi), show=False)
    p[0].surface_color = lambda a : a
    p.save(tmp_file('%s_colors_surface_arity1.png' % name))
    p[0].surface_color = lambda a, b : b
    p.save(tmp_file('%s_colors_surface_arity2.png' % name))
    p[0].surface_color = lambda a, b, c : c
    p.save(tmp_file('%s_colors_surface_arity3a.png' % name))
    p[0].surface_color = lambda a, b, c : sqrt((a-3*pi)**2+b**2)
    p.save(tmp_file('%s_colors_surface_arity3b.png' % name))

    p = plot(x*cos(4*y), x*sin(4*y), y,
             (-1, 1), (-1, 1), show=False)
    p[0].surface_color = lambda a : a
    p.save(tmp_file('%s_colors_param_surf_arity1.png' % name))
    p[0].surface_color = lambda a, b : a*b
    p.save(tmp_file('%s_colors_param_surf_arity2.png' % name))
    p[0].surface_color = lambda a, b, c : sqrt(a**2+b**2+c**2)
    p.save(tmp_file('%s_colors_param_surf_arity3.png' % name))

    ###
    # Examples from the 'advanced' notebook
    ###

    i = Integral(log((sin(x)**2+1)*sqrt(x**2+1)),(x,0,y))
    p = plot(i,(1,5), show=False)
    p.save(tmp_file('%s_advanced_integral.png' % name))

    s = summation(1/x**y,(x,1,oo))
    p = plot(s, (2,10), show=False)
    p.save(tmp_file('%s_advanced_inf_sum.png' % name))

    p = plot(summation(1/x,(x,1,y)), (2,10), show=False)
    p[0].only_integers = True
    p[0].steps = True
    p.save(tmp_file('%s_advanced_fin_sum.png' % name))


    ###
    # Test expressions that can not be translated to np and generate complex
    # results.
    ###


    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        plot(sqrt(sqrt(-x)), show=False).save(tmp_file())
        assert len(w) == 1
        assert "Complex values as arguments to numpy functions encountered." == str(w[-1].message)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        plot(LambertW(x), show=False).save(tmp_file())
        assert len(w) > 10 #TODO trace where do the other warnings come from
        assert "Complex values as arguments to python math functions encountered." == str(w[-1].message)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        plot(sqrt(LambertW(x)), show=False).save(tmp_file())
        assert len(w) == 1
        assert "Complex values as arguments to python math functions encountered." == str(w[-1].message)
    #TODO this one does not work properly:
    plot(sin(x)+I*cos(x)).save(tmp_file())


    ###
    # Test all valid input args for plot()
    ###

    # 2D line with all possible inputs
    plot(x              ).save(tmp_file())
    plot(x, (x,        )).save(tmp_file())
    plot(x, (   -10, 10)).save(tmp_file())
    plot(x, (x, -10, 10)).save(tmp_file())

    # two 2D lines
    plot((x,             ), (x+1, (x, -10, 10))).save(tmp_file())
    plot((x, (x,        )), (x+1, (x, -10, 10))).save(tmp_file())
    plot((x, (   -10, 10)), (x+1, (x, -10, 10))).save(tmp_file())
    plot((x, (x, -10, 10)), (x+1, (x, -10, 10))).save(tmp_file())
    plot([x, x+1],             ).save(tmp_file())
    plot([x, x+1], (x,        )).save(tmp_file())
    plot([x, x+1], (   -10, 10)).save(tmp_file())
    plot([x, x+1], (x, -10, 10)).save(tmp_file())

    # 2D and 3D parametric lines
    plot(x, 2*x).save(tmp_file())
    plot(x, 2*x, sin(x)).save(tmp_file())

    # surface and parametric surface
    plot(x*y).save(tmp_file())
    plot(x*y, 1/x, 1/y).save(tmp_file())

    # some bizarre combinatrions
    ## two 3d parametric lines and a surface
    plot(([x, -x], x**2, sin(x)), (x*y,)).save(tmp_file())


def test_matplotlib():
    try:
        import matplotlib
        import numpy
        if matplotlib.__version__ < '1.2.0':
            return
        plot_and_save('test')
    except ImportError:
        pass
