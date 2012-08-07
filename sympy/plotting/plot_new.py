def plot_line(*args, **kwargs):
    """
    Plots a function of a single variable.

    The plotting uses an adaptive algorithm which samples recursively to
    accurately plot the plot. The adaptive algorithm uses a random point near
    the midpoint of two points that has to be further sampled. Hence the same
    plots can appear slightly different.

    Usage
    =====
    plot_line((expr1, range), (expr2, range), ...)
    plot_line(expr, range, ...)

    Arguments
    =========
    ``expr`` : Expression representing the function of single variable
    ``range``: (x, 0 , 5), A 3 - tuple denoting the range of the free variable.

    If the ranges is not specified, then a default range of (-10, 10) is used.

    Named Parameters
    ===============
    ``adaptive``: Boolean. The default value is set to True. Set adaptive to False and
    specify ``nb_of_points`` if uniform sampling is required.

    ``depth``: int Recursion depth of the adaptive algorithm. A depth of value ``n``
    samples a maximum of `2^{n}` points.

    ``nb_of_points``: int. Used when the ``adaptive`` is set to False. The function
    is uniformly sampled at ``nb_of_point`` number of points.

    ``title`` : str. Title of the plot. It is set to the latex representation of
    the expression, if the plot has only one expression.

    ``xlabel`` : str. Label for the x - axis.

    ``line_color``: List of floats corresponding to the line color of the
    expressions. The length of the List should match the number of expressions.

    ``xscale``: {'linear', 'log'} Sets the scaling of the x - axis.

    ``yscale``: {'linear', 'log'} Sets the scaling if the y - axis.

    ``axis_center``: tuple of two floats denoting the coordinates of the center or
    {'center', 'auto'}

    ``xlim`` : tuple of two floats, denoting the x - axis limits.

    ``ylim`` : tuple of two floats, denoting the y - axis limits.

    Examples
    ========
    >>> plot_line(x**2, (x, -5, 5))

    Multiple plots
    >>> plot_line((x**2, (x, -6, 6)), (x, (x, -5, 5)))

    Multiple plots with same range.
    >>> plot_line(x**2, x, (x -5, 5))

    No adaptive sampling.
    >>> plot_line(x**2, adaptive = False, nb_of_points = 400)

    """


def plot_parametric(*args, **kwargs):
    """
    Plots a 2D parametric plot.

    The plotting uses an adaptive algorithm which samples recursively to
    accurately plot the plot. The adaptive algorithm uses a random point near
    the midpoint of two points that has to be further sampled. Hence the same
    plots can appear slightly different.

    Usage
    =====
    plot_parametric((expr_x, expr_y, range), ...)
    plot_parametric(expr_x, expr_y, range)

    Arguments
    =========
    ``expr_x`` : Expression representing the function along x.
    ``expr_y`` : Expression representing the function along y.
    ``range``: (u, 0 , 5), A 3 - tuple denoting the range of the parameter
    variable.
    If the range is not specified, then a default range of (-10, 10) is used.

    Named Parameters
    ===============
    ``adaptive``: Boolean. The default value is set to True. Set adaptive to
    False and specify ``nb_of_points`` if uniform sampling is required.

    ``depth``: int Recursion depth of the adaptive algorithm. A depth of
    value ``n`` samples a maximum of `2^{n}` points.

    ``nb_of_points``: int. Used when the ``adaptive`` is set to False. The
    function is uniformly sampled at ``nb_of_point`` number of points.

    ``title`` : str. Title of the plot.

    ``xlabel`` : str. Label for the x - axis.

    ``line_color``: List of floats corresponding to the line color of the
    expressions. The length of the List should match the number of expressions.

    ``xscale``: {'linear', 'log'} Sets the scaling of the x - axis.

    ``yscale``: {'linear', 'log'} Sets the scaling if the y - axis.

    ``axis_center``: tuple of two floats denoting the coordinates of the center or
    {'center', 'auto'}

    ``xlim`` : tuple of two floats, denoting the x - axis limits.

    ``ylim`` : tuple of two floats, denoting the y - axis limits.

    Examples
    ========

    >>> plot_parametric(cos(u), sin(u), (u, -5, 5))

    Multiple parametric plots.
    >>> plot_parametric((cos(u), sin(u), (u, -5, 5)), (cos(u), u, (u, -5, 5)))

    """


def plot3D_parametric(*args, **kwargs):
    """
    Plots a 3D parametric plot.

    Usage
    =====
    Single plot
    plot3D_parametric(expr_x, expr_y, expr_z, range, kwargs)

    Multiple plots
    plot3D_parametric((expr_x, expr_y, expr_z, range), ...)

    Arguments
    =========
    ``expr_x`` : Expression representing the function along x.
    ``expr_y`` : Expression representing the function along y.
    ``expr_z`` : Expression representing the function along z.
    ``range``: (u, 0 , 5), A 3 - tuple denoting the range of the parameter
    variable.
    If the range is not specified, then a default range of (-10, 10) is used.

    Named Parameters
    ===============
    ``nb_of_points``: int. Used when the ``adaptive`` is set to False. The
    function is uniformly sampled at ``nb_of_point`` number of points.

    ``title`` : str. Title of the plot.

    ``xlabel`` : str. Label for the x - axis.

    ``line_color``: List of floats corresponding to the line color of the
    expressions. The length of the List should match the number of expressions.

    ``xscale``: {'linear', 'log'} Sets the scaling of the x - axis.

    ``yscale``: {'linear', 'log'} Sets the scaling if the y - axis.

    ``axis_center``: tuple of two floats denoting the coordinates of the center or
    {'center', 'auto'}

    ``xlim`` : tuple of two floats, denoting the x - axis limits.

    ``ylim`` : tuple of two floats, denoting the y - axis limits.

    Examples
    ========
    >>> plot_parametric3D(cos(u), sin(u), u, (u, -5, 5))

    Multiple plots with default range.

    >>> plot_parametric3D((cos(u), sin(u), u), (sin(u), u**2, u))

    """

def plot3D(*args, **kwargs):
    """
    Plots a 3D surface plot.

    Usage
    =====
    Single plot
    plot3D(expr, range_x, range_y)

    Multiple plots.
    plot3D((expr, range_x, range_y), ...)

    Arguments
    =========
    ``expr`` : Expression representing the function along x.
    ``range_x``: (x, 0 , 5), A 3 - tuple denoting the range of the x
    variable.
    ``range_y``: (y, 0 , 5), A 3 - tuple denoting the range of the y
    variable.
    If the ranges is not specified, then a default range of (-10, 10) is used.

    Named Parameters
    ===============
    ``nb_of_points``: int. The function is uniformly sampled at ``nb_of_point``
    number of points.

    ``title`` : str. Title of the plot.

    ``xlabel`` : str. Label for the x - axis.

    ``line_color``: List of floats corresponding to the line color of the
    expressions. The length of the List should match the number of expressions.

    ``xscale``: {'linear', 'log'} Sets the scaling of the x - axis.

    ``yscale``: {'linear', 'log'} Sets the scaling if the y - axis.

    ``axis_center``: tuple of two floats denoting the coordinates of the center or
    {'center', 'auto'}

    ``xlim`` : tuple of two floats, denoting the x - axis limits.

    ``ylim`` : tuple of two floats, denoting the y - axis limits.

    Examples
    ========
    >>> plot_3D(x*y, (x, -5, 5), (y, -5, 5))

    Multiple plots
    >>> plot_3D((x**2 + y**2, (x, -5, 5), (y, -5, 5), (x*y))

    """


def plot3D_surface(*args, **kwargs):
    """
    Plots a 3D parametric surface plot.

    Usage
    =====
    Single plot.
    plot3D_surface(expr_x, expr_y, expr_z, range_x, range_y)

    Multiple plots.
    plot_parametric_surface((expr_x, expr_y, expr_z, range_x, range_y), ...)
    If the ranges are not specified, then a default range of (-10, 10) is used.

    Arguments
    =========
    ``expr`` : Expression representing the function along x.
   ``range_x``: (x, 0 , 5),  A 3 - tuple denoting the range of the x
    variable.
    ``range_y``: (y, 0 , 5),  A 3 - tuple denoting the range of the y
    variable.

    Named Parameters
    ===============
    ``nb_of_points``: int. The function is uniformly sampled at ``nb_of_point``
    number of points.

    ``title`` : str. Title of the plot.

    ``xlabel`` : str. Label for the x - axis.

    ``line_color``: List of floats corresponding to the line color of the
    expressions. The length of the List should match the number of expressions.

    ``xscale``: {'linear', 'log'} Sets the scaling of the x - axis.

    ``yscale``: {'linear', 'log'} Sets the scaling if the y - axis.

    ``axis_center``: tuple of two floats denoting the coordinates of the center or
    {'center', 'auto'}

    ``xlim`` : tuple of two floats, denoting the x - axis limits.

    ``ylim`` : tuple of two floats, denoting the y - axis limits.

    Examples
    ========
    >>> plot3D_surface(cos(u + v), sin(u - v), u - v, (u, -5, 5), (v, -5, 5))
    """
