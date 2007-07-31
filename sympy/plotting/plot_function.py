from sympy import Basic, Symbol

from plot_object import PlotObject
from plot_curve import PlotCurve
from plot_surface import PlotSurface
from plot_interval import PlotInterval
from plot_modes import get_plot_mode, fill_intervals, fill_i_vars

from time import clock

class PlotFunction(PlotObject):
    """
    """
    def __init__(self, *args, **kwargs):
        args, kwargs = extract_options(args, kwargs)
        self.visible = 'True' == kwargs.pop('visible', 'True')
        self.mode_str = kwargs.pop('mode', '')

        d_vars, intervals = interpret_args(args)
        i_vars, interval_bindings = find_independent_vars(d_vars, intervals)
        i_var_c, d_var_c = len(i_vars), len(d_vars)

        if isinstance(self.mode_str, str):
            mode, self.mode_str = get_plot_mode(self.mode_str, d_var_c, i_var_c)
            i_vars = fill_i_vars(self.mode_str, i_vars)
            intervals = fill_intervals(self.mode_str, i_vars, intervals)
        elif callable(self.mode_str):
            self.mode_str, mode = self.mode_str.__name__, self.mode_str
            if len(intervals) != len(i_vars):
                raise ValueError("Intervals must be provided for each independent " +
                                 "variable when using a custom plot mode.")
        else: raise ValueError("'mode' named argument must be a string or a callable.")

        try:
            f = mode(*(d_vars+i_vars))
        except Exception, e:
            raise ValueError( mode_init_error % (str(e)) )

        plot_object_choice = {1: PlotCurve, 2: PlotSurface}
        try: cls = plot_object_choice[len(intervals)]
        except: raise ValueError("Cannot plot in %i dimensions." % (dimensions))

        self.plot_object = cls(f, *intervals, **kwargs)

    def set_style(self, rs):
        self.plot_object.set_style(rs)

    def get_calc_state(self):
        return self.plot_object.get_calc_state()

    def draw(self):
        self.plot_object.draw()

interval_wrong_order = "Interval %s was given before any function(s)."
interpret_error = "Could not interpret %s as a function or interval."
mode_init_error = "Mode initialization failed. Inner exception: %s"

def interpret_args(args):
    functions, intervals = [], []
    for a in args:
        i = PlotInterval.try_parse(a)
        if i is not None:
            if len(functions) == 0:
                raise ValueError(interval_wrong_order % (str(i)))
            else:
                intervals.append(i)
        else:
            try:
                f = Basic.sympify(a)
                functions.append(f)
            except:
                raise ValueError(interpret_error % str(a))
    return functions, intervals

def find_independent_vars(functions, intervals):
    """
    i_vars collects the independent variables first in the
    order specified by any given intervals, followed by any
    remaining Symbols found in functions.
    
    interval_bindings maps each variable to an interval, or
    None if no interval given.
    """
    i_vars, interval_bindings = [], {}
    for i in intervals:
        if i.v not in interval_bindings:
            interval_bindings[i.v] = i
            i_vars.append(i.v)
    for f in functions:
        for a in f.atoms(type=Symbol):
            if a not in interval_bindings:
                interval_bindings[a] = None
                i_vars.append(a)
    return i_vars, interval_bindings

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

def extract_options(args, kwargs):
    nkwargs, nargs = {}, []
    for a in args:
        if isinstance(a, str):
            nkwargs = dict(nkwargs, **parse_option_string(a))
        else:
            nargs.append(a)
    nkwargs = dict(nkwargs, **kwargs)
    return nargs, nkwargs
