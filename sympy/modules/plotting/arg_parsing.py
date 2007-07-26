from sympy import Basic, Symbol

from plot_function import PlotFunction

def parse_function(f):
    try:
        return Basic.sympify(f)
    except:
        return None

def parse_functions(args):
    t = []
    while len(args) > 0:
        f = parse_function(args[0])
        if f is None: break
        t.append(f)
        args.pop(0)
    return t

def parse_interval(i):
    if not isinstance(i, (tuple, list)) or len(i) not in [3, 4]:
        return None
    var, v_min, v_max, v_steps = None, None, None, None
    try:
        if not isinstance(i[0], Symbol):
            raise ValueError( "First element of interval must be a Symbol representing a dependent variable." )

        v_min = float(Basic.sympify(i[1]).evalf())
        v_max = float(Basic.sympify(i[2]).evalf())

        if len(i) == 4 and i[3] is not None:
             v_steps = int(i[3])
    except:
        return None

    return (i[0], v_min, v_max, v_steps)

def parse_intervals(args):
    t = []
    while len(args) > 0:
        i = parse_interval(args[0])
        if i is None:
            break
        t.append(i)
        args.pop(0)
    return t

def parse_option_string(s):
    if not isinstance(s, str):
        return None
    options = {}
    for token in s.split(';'):
        pieces = token.split('=')
        if len(pieces) == 1:
            option, value = pieces[0], ""
        elif len(pieces) == 2:
            option, value = pieces
        else:
            raise ValueError("Plot option string '%s' is malformed." % (s))
        options[option.strip()] = value.strip()
    return options

def parse_option_strings(args):
    t = {'mode':'cartesian'}
    while len(args) > 0:
        a = parse_option_string(args[0])
        if a is None: break
        t = dict(t, **a)
        args.pop(0)
    return t

def parse_plot_args(*args):
    """
    Regexp:    (function+, interval*, option_string*)+
    Examples:   x**2, -x
                x**2, [x, -5, 5]
                sin(r), [r, -pi/Real(2), pi/Real(2)], "mode=polar"
                x**2, [x, -5, 5], sin(r), "mode=polar"
    """
    args = list(args)
    while len(args) > 0:
        functions = parse_functions(args)
        if len(functions) == 0:
            raise ValueError( "Invalid argument list format." )
        intervals = parse_intervals(args)
        options = parse_option_strings(args)

        if len(intervals) == 0 and len(options) == 0 and len(args) > 0:
            raise ValueError( "Could not interpret token '%s'." % (str(args[0])) )

        if options['mode'] == 'parametric':
            yield PlotFunction(functions, intervals, options)
        else:
            for f in functions:
                yield PlotFunction(f, intervals, options)

def parse_function_args(*args):
    """
    Regexp:     function, interval*, option_string*
    Examples:   x**2, [x, -5, 5]
                x**2 - y**2, [x, -1, 1], [y, -1, 1]

    Note:       There can be multiple functions if option
                'type=parametric' is specified
    """
    args = list(args)
    if len(args) == 0:
        raise ValueError( "No arguments given." )
    functions = parse_functions(args)
    if len(functions) == 0:
        raise ValueError( "Invalid argument list format." )
    intervals = parse_intervals(args)
    options = parse_option_strings(args)

    if len(args) > 0:
        raise ValueError( "Could not interpret token '%s'." % (str(args[0])) )
    if options['mode'] == 'parametric':
        return PlotFunction(functions, intervals, options)
    elif len(functions) == 1:
        return PlotFunction(functions[0], intervals, options)
    else:
        raise Exception("Only one function was expected in this context.")
