"""
Limited tests of the visualization module. Right now it just makes
sure that passing custom Axes works.

"""

# This test either prints something to the terminal or displays a plot,
# depending on whether matplotlib is installed or not.  Neither is ideal
# for a test, so let's just skip this entirely.

disabled = True
from sympy.mpmath import mp, fp

def test_axes():
    try:
        import matplotlib
        version = matplotlib.__version__.split("-")[0]
        version = version.split(".")[:2]
        if [int(_) for _ in version] < [0,99]:
            raise ImportError
        import pylab
    except ImportError:
        print("\nSkipping test (pylab not available or too old version)\n")
        return
    fig = pylab.figure()
    axes = fig.add_subplot(111)
    for ctx in [mp, fp]:
        ctx.plot(lambda x: x**2, [0, 3], axes=axes)
        assert axes.get_xlabel() == 'x'
        assert axes.get_ylabel() == 'f(x)'

    fig = pylab.figure()
    axes = fig.add_subplot(111)
    for ctx in [mp, fp]:
        ctx.cplot(lambda z: z, [-2, 2], [-10, 10], axes=axes)
    assert axes.get_xlabel() == 'Re(z)'
    assert axes.get_ylabel() == 'Im(z)'
