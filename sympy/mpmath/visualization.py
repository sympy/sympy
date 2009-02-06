"""
Plotting (requires matplotlib)
"""

from mptypes import mpc, inf, isnan, isinf, arange, complex_types
from functions import sqrt, arg

from colorsys import hsv_to_rgb, hls_to_rgb

plot_ignore = (ValueError, ArithmeticError, ZeroDivisionError)

def plot(f, xlim=[-5,5], ylim=None, points=200, file=None, dpi=None,
    singularities=[]):
    r"""
    Shows a simple 2D plot of a function `f(x)` or list of functions
    `[f_0(x), f_1(x), \ldots, f_n(x)]` over a given interval
    specified by *xlim*. Some examples::

        plot(lambda x: exp(x)*li(x), [1, 4])
        plot([cos, sin], [-4, 4])
        plot([fresnels, fresnelc], [-4, 4])
        plot([sqrt, cbrt], [-4, 4])
        plot(lambda t: zeta(0.5+t*j), [-20, 20])
        plot([floor, ceil, abs, sign], [-5, 5])

    Points where the function raises a numerical exception or
    returns an infinite value are removed from the graph.
    Singularities can also be excluded explicitly
    as follows (useful for removing erroneous vertical lines)::

        plot(cot, ylim=[-5, 5])   # bad
        plot(cot, ylim=[-5, 5], singularities=[-pi, 0, pi])  # good

    For parts where the function assumes complex values, the
    real part is plotted with dashes and the imaginary part
    is plotted with dots.

    NOTE: This function requires matplotlib (pylab).
    """
    import pylab
    pylab.clf()
    if not isinstance(f, (tuple, list)):
        f = [f]
    a, b = xlim
    colors = ['b', 'r', 'g', 'm', 'k']
    for n, func in enumerate(f):
        x = arange(a, b, (b-a)/float(points))
        segments = []
        segment = []
        in_complex = False
        for i in xrange(len(x)):
            try:
                if i != 0:
                    for sing in singularities:
                        if x[i-1] <= sing and x[i] >= sing:
                            raise ValueError
                v = func(x[i])
                if isnan(v) or abs(v) > 1e300:
                    raise ValueError
                if isinstance(v, complex_types):
                    re = float(v.real)
                    im = float(v.imag)
                    if not in_complex:
                        in_complex = True
                        segments.append(segment)
                        segment = []
                    segment.append((float(x[i]), re, im))
                else:
                    if in_complex:
                        in_complex = False
                        segments.append(segment)
                        segment = []
                    segment.append((float(x[i]), v))
            except plot_ignore:
                if segment:
                    segments.append(segment)
                segment = []
        if segment:
            segments.append(segment)
        for segment in segments:
            x = [s[0] for s in segment]
            y = [s[1] for s in segment]
            if not x:
                continue
            c = colors[n % len(colors)]
            if len(segment[0]) == 3:
                z = [s[2] for s in segment]
                pylab.plot(x, y, '--'+c, linewidth=3)
                pylab.plot(x, z, ':'+c, linewidth=3)
            else:
                pylab.plot(x, y, c, linewidth=3)
    pylab.xlim(xlim)
    if ylim:
        pylab.ylim(ylim)
    pylab.xlabel('x')
    pylab.ylabel('f(x)')
    pylab.grid(True)
    if file:
        pylab.savefig(file, dpi=dpi)
    else:
        pylab.show()

def default_color_function(z):
    if isinf(z):
        return (1.0, 1.0, 1.0)
    if isnan(z):
        return (0.5, 0.5, 0.5)
    pi = 3.1415926535898
    a = (float(arg(z)) + pi) / (2*pi)
    a = (a + 0.5) % 1.0
    b = 1.0 - float(1/(1.0+abs(z)**0.3))
    return hls_to_rgb(a, b, 0.8)

def cplot(f, re=[-5,5], im=[-5,5], points=2000, color=default_color_function,
    verbose=False, file=None, dpi=None):
    """
    Plots the given complex-valued function *f* over a rectangular part
    of the complex plane specified by the pairs of intervals *re* and *im*.
    For example::

        cplot(lambda z: z, [-2, 2], [-10, 10])
        cplot(exp)
        cplot(zeta, [0, 1], [0, 50])

    By default, the complex argument (phase) is shown as color (hue) and
    the magnitude is show as brightness. You can also supply a
    custom color function (*color*). This function should take a
    complex number as input and return an RGB 3-tuple containing
    floats in the range 0.0-1.0.

    To obtain a sharp image, the number of points may need to be
    increased to 100,000 or thereabout. Since evaluating the
    function that many times is likely to be slow, the 'verbose'
    option is useful to display progress.

    NOTE: This function requires matplotlib (pylab).
    """
    import pylab
    pylab.clf()
    rea, reb = re
    ima, imb = im
    dre = reb - rea
    dim = imb - ima
    M = int(sqrt(points*dre/dim)+1)
    N = int(sqrt(points*dim/dre)+1)
    x = pylab.linspace(rea, reb, M)
    y = pylab.linspace(ima, imb, N)
    # Note: we have to be careful to get the right rotation.
    # Test with these plots:
    #   cplot(lambda z: z if z.real < 0 else 0)
    #   cplot(lambda z: z if z.imag < 0 else 0)
    w = pylab.zeros((N, M, 3))
    for n in xrange(N):
        for m in xrange(M):
            z = mpc(x[m], y[n])
            try:
                v = color(f(z))
            except plot_ignore:
                v = (0.5, 0.5, 0.5)
            w[n,m] = v
        if verbose:
            print n, "of", N
    pylab.imshow(w, extent=(rea, reb, ima, imb), origin='lower')
    pylab.xlabel('Re(z)')
    pylab.ylabel('Im(z)')
    if file:
        pylab.savefig(file, dpi=dpi)
    else:
        pylab.show()
