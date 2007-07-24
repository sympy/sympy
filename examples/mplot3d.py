import sys
sys.path.append("..")

from sympy import Symbol, Basic
from sample import sample

def mplot3d(f, var1, var2, show=True):
    """
    Plot a 3d function using matplotlib/Tk.
    """

    import warnings
    warnings.filterwarnings("ignore", "Could not match \S")

    try:
        import pylab as p
        import matplotlib.axes3d as p3
    except ImportError:
        raise ImportError("Matplotlib is required to use mplot3d.")

    x, y, z = sample(f, var1, var2)

    fig = p.figure()
    ax = p3.Axes3D(fig)

    #ax.plot_surface(x,y,z) #seems to be a bug in matplotlib
    ax.plot_wireframe(x,y,z)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    if show:
        p.show()

if __name__ == "__main__":
    from sympy import sin
    x = Symbol('x')
    y = Symbol('y')

    mplot3d(x**2-y**2, (x, -10.0, 10.0, 20), (y, -10.0, 10.0, 20))
    #mplot3d(x**2+y**2, (x, -10.0, 10.0, 20), (y, -10.0, 10.0, 20))
    #mplot3d(sin(x)+sin(y), (x, -3.14, 3.14, 10), (y, -3.14, 3.14, 10))
