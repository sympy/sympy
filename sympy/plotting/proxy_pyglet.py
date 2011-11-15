from warnings import warn
from sympy.core.compatibility import SymPyDeprecationWarning
from pygletplot import PygletPlot

def Plot(*args, **kwargs):
    """ A temporary proxy for an interface under deprecation.

    This proxy is the one imported by `from sympy import *`.

    The Plot class will change in future versions of sympy to use the new
    plotting module. That new plotting module is already used by the
    plot() function (lowercase). To write code compatible with future versions
    of sympy use that function (plot() lowercase). Or if you want to use the
    old plotting module just import it directly:
    `from sympy.plotting.pygletplot import PygletPlot`

    To use Plot from the new plotting module do:
    `from sympy.plotting.plot import Plot`

    In future version of sympy you will also be able to use
    `from sympy.plotting import Plot` but in the current version this will
    import this proxy object. It's done for backward compatibility.

    The old plotting module is not deprecated. Only the location will
    change. The new location is sympy.plotting.pygletplot.
    """
    warn('This interface will change in future versions of sympy.'
         ' As a precatuion use the plot() function (lowercase).'
         ' See the docstring for details.',
         SymPyDeprecationWarning)

    return PygletPlot(*args, **kwargs)
