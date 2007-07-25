from sympy import Basic, Symbol
from plot_object import PlotObject

class PlotFunction(PlotObject):

    f = None
    intervals = []
    options = {}

    x_min, x_max = 0.0, 0.0
    y_min, y_max = 0.0, 0.0
    z_min, z_max = 0.0, 0.0

    def __new__(cls, f, intervals, options):
        subcls = PlotFunctionRegistry.get(options['mode'])
        if subcls is None:
            error_str = "Plot mode '%s' is not supported (or not registered)."
            raise Exception(error_str % (options['mode']))
        o = subcls(f, intervals, options)
        if 'visible' in options and options['visible'] == 'false':
            o.visible = False
            del options['visible']
        return o

    def __init__(self, *args):
        raise Exception("PlotFunction is not meant to be directly instantiated.")

    def __repr__(self):
        f_str = ""
        if isinstance(self.f, (tuple, list)):
            f_str = ", ".join( [str(f) for f in self.f] )
        else:
            f_str = str(self.f)
        i_str = ", ".join([ "[" + ( ",".join( [str(k) for k in i] ) )
                          + "]" for i in self.intervals ])
        o_str = options_str(self.options)

        a_list = [f_str]
        if i_str != "": a_list.append(i_str)
        if o_str != "": a_list.append(o_str)

        return ", ".join( a_list )

    def __str__(self):
        return self.__repr__()

class PlotFunctionRegistry(object):
    _r = {}

    @staticmethod
    def get(i):
        if i not in PlotFunctionRegistry._r:
            return None
        return PlotFunctionRegistry._r[i]

    @staticmethod
    def register(i, plot_type):
        assert isinstance(i, str)
        assert isinstance(plot_type, type)
        PlotFunctionRegistry._r[i] = plot_type

    @staticmethod
    def remove(i):
        del PlotFunctionRegistry._r[i]

    @staticmethod
    def list():
        return list(iter(PlotFunctionRegistry._r))

def get_vars(f):
    assert isinstance(f, Basic)
    return f.atoms(type=Symbol)

def count_vars(f):
    return len(get_vars(f))

def options_str(d):
    o_list = list("%s=%s" % (k, d[k]) for k in d)
    o_str = ";".join(o_list)
    return "'" + o_str + "'"

def vrange(a_min, a_max, a_steps):
    """
    Helper function which returns an array containing a_step+1
    elements ranging from a_min to a_max.
    """
    a_delta = (a_max-a_min)/float(a_steps)
    return list(a_min+a_delta*i for i in range(a_steps+1))

def interpolate(a_min, a_max, a_ratio):
    return a_min + a_ratio * (a_max - a_min)

def rinterpolate(a_min, a_max, a_value):
    a_range = a_max-a_min
    if a_range == 0:
        a_range = 1.0
    return (a_value - a_min) / float(a_range)

def interpolate_color(color1, color2, ratio):
    return [interpolate(color1[i], color2[i], ratio) for i in range(3)]

def fsubs(f, a, _a, b=None, _b=None):
    """
    Returns a float z = f(a) or z = f(a,b).
    """
    try:
        if b is None:
            return float(f.subs(a, _a)) # seems the best choice right now
            #return f.subs(a, _a).evalf(precision=10) # not a float val
            #return float(f.subs(a, _a).evalf(precision=10))
        else:
            return float(f.subs(a, _a).subs(b, _b))
            #return f.subs(a, _a).subs(b, _b).evalf(precision=10)
            #return float(f.subs(a, _a).subs(b, _b).evalf(precision=10))
    except:
        return None
