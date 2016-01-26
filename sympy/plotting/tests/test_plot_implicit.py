import warnings
from sympy import (plot_implicit, cos, Symbol, symbols, Eq, sin, re, And, Or, exp, I,
                   tan, pi)
from sympy.plotting.plot import unset_show
from tempfile import NamedTemporaryFile
from sympy.utilities.pytest import skip
from sympy.external import import_module

#Set plots not to show
unset_show()


def tmp_file(name=''):
    return NamedTemporaryFile(suffix='.png').name

def plot_and_save(expr, *args, **kwargs):
    name = kwargs.pop('name', '')
    p = plot_implicit(expr, *args, **kwargs)
    p.save(tmp_file(name))
    # Close the plot to avoid a warning from matplotlib
    p._backend.close()

def plot_implicit_tests(name):
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')
    #implicit plot tests
    plot_and_save(Eq(y, cos(x)), (x, -5, 5), (y, -2, 2), name=name)
    plot_and_save(Eq(y**2, x**3 - x), (x, -5, 5),
            (y, -4, 4), name=name)
    plot_and_save(y > 1 / x, (x, -5, 5),
            (y, -2, 2), name=name)
    plot_and_save(y < 1 / tan(x), (x, -5, 5),
            (y, -2, 2), name=name)
    plot_and_save(y >= 2 * sin(x) * cos(x), (x, -5, 5),
            (y, -2, 2), name=name)
    plot_and_save(y <= x**2, (x, -3, 3),
            (y, -1, 5), name=name)

    #Test all input args for plot_implicit
    plot_and_save(Eq(y**2, x**3 - x))
    plot_and_save(Eq(y**2, x**3 - x), adaptive=False)
    plot_and_save(Eq(y**2, x**3 - x), adaptive=False, points=500)
    plot_and_save(y > x, (x, -5, 5))
    plot_and_save(And(y > exp(x), y > x + 2))
    plot_and_save(Or(y > x, y > -x))
    plot_and_save(x**2 - 1, (x, -5, 5))
    plot_and_save(x**2 - 1)
    plot_and_save(y > x, depth=-5)
    plot_and_save(y > x, depth=5)
    plot_and_save(y > cos(x), adaptive=False)
    plot_and_save(y < cos(x), adaptive=False)
    plot_and_save(And(y > cos(x), Or(y > x, Eq(y, x))))
    plot_and_save(y - cos(pi / x))

    #Test plots which cannot be rendered using the adaptive algorithm
    #TODO: catch the warning.
    plot_and_save(Eq(y, re(cos(x) + I*sin(x))), name=name)

    with warnings.catch_warnings(record=True) as w:
        plot_and_save(x**2 - 1, legend='An implicit plot')
        assert len(w) == 1
        assert issubclass(w[-1].category, UserWarning)
        assert 'No labeled objects found' in str(w[0].message)

def test_line_color():
    x, y = symbols('x, y')
    p = plot_implicit(x**2 + y**2 - 1, line_color="green", show=False)
    assert p._series[0].line_color == "green"
    p = plot_implicit(x**2 + y**2 - 1, line_color='r', show=False)
    assert p._series[0].line_color == "r"

def test_matplotlib():
    matplotlib = import_module('matplotlib', min_module_version='1.1.0', catch=(RuntimeError,))
    if matplotlib:
        plot_implicit_tests('test')
        test_line_color()
    else:
        skip("Matplotlib not the default backend")
