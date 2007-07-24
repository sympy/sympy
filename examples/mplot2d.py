import sys
sys.path.append("..")

from sympy import Symbol, Basic
from sample import sample

def mplot2d(f, var, show=True):
    """
    Plot a 2d function using matplotlib/Tk.
    """

    import warnings
    warnings.filterwarnings("ignore", "Could not match \S")

    try:
        import pylab as p
    except ImportError:
        raise ImportError("Matplotlib is required to use mplot2d.")

    if not isinstance(f, (tuple, list)):
        f = [f,]

    for f_i in f:
        x, y = sample(f_i, var)
        p.plot(x, y)
    
    p.draw()
    if show:
        p.show()

if __name__ == "__main__":
    from sympy import sqrt, sin, log, pi
    x = Symbol('x')
    
    #mplot2d(log(x), (x, 0, 2, 100))
    #mplot2d([sin(x), -sin(x)], (x, float(-2*pi), float(2*pi), 50))
    mplot2d([sqrt(x), -sqrt(x), sqrt(-x), -sqrt(-x)], (x, -40.0, 40.0, 80))
