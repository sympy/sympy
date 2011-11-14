def Plot(*args, **kwargs):
    """ A temporary proxy for an interface under deprecation.

    The Plot class will change in future versions of sympy to use the new
    plotting module. That new plotting module is used also by the
    plot() function (lowercase). To write code compatible with future versions
    of sympy use that function (plot() lowercase). Or if you want to use the
    old plotting module just import it directly:
    from sympy.plotting.pygletplot import Plot

    The old plotting module is not deprecated. Only the location will
    change. The new location is sympy.plotting.pygletplot.

    The current location sympy.plotting will be used for the new plotting
    module (currently in sympy.plotting.newplot).
    """
    from warnings import warn
    warn('This interface will change in future versions of sympy.'
         ' As a precatuion use the plot() function (lowercase).')

    from pygletplot import Plot
    return Plot(*args, **kwargs)
