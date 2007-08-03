from sympy import Basic, Symbol, symbols
from math import sin, cos, tan, asin, acos, atan, log, sqrt
from util import interpolate
intrp = interpolate

default_color_schemes = {'rainbow': "intrp(0.45,0.95,z),intrp(0.45,0.9,y),intrp(0.5,0.98,x)",
                         'zfade': "intrp(0.2,0.8,z),0.3,intrp(0.8,0.2,z)"}

class ColorScheme(object):

    uniform = False

    def __new__(self, *args, **kwargs):
        e = Exception(("Could not parse color scheme arguments %s.")
                        % (', '.join(str(a) for a in args)))
        cls = None
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, (list,tuple)):
                args = arg
            elif isinstance(arg, str):
                if arg in default_color_schemes:
                    args = [default_color_schemes[arg]]
                cls = ColorSchemeLambda
        if not cls:
            if isinstance(args, (list,tuple)):
                if len(args) == 3:
                    a,b,c = args
                    if isinstance(a, str) and isinstance(b, str) and isinstance(c, str):
                        args = ["%s, %s, %s" % (a,b,c)]
                        cls = ColorSchemeLambda
                    else:
                        try:
                            fl = float
                            a,b,c = fl(a), fl(b), fl(c)
                            cls = ColorSchemeFlat
                        except:
                            args = tuple(list(args) + [symbols('xyz'),])
                if len(args) == 4:
                    a,b,c,xyzuv = args
                    if isinstance(xyzuv, Symbol):
                        xyzuv = tuple([xyzuv])
                    try:
                        assert len(xyzuv) <= 5 and len(xyzuv) > 0

                        if len(xyzuv) == 1:
                            xyzuv = (Symbol('unbound1'), Symbol('unbound2'),
                                     Symbol('unbound3'), xyzuv[0],
                                     Symbol('unbound4'))
                        elif len(xyzuv) == 2:
                            xyzuv = (Symbol('unbound1'), Symbol('unbound2'),
                                     Symbol('unbound3'), xyzuv[0], xyzuv[1])
                        elif len(xyzuv) == 3:
                            xyzuv = (xyzuv[0], xyzuv[1], xyzuv[2],
                                     Symbol('unbound1'), Symbol('unbound2'))
                        elif len(xyzuv) == 4:
                            xyzuv = (xyzuv[0], xyzuv[1], xyzuv[2],
                                     xyzuv[3], Symbol('unbound1'))
    
                        sy = Basic.sympify
                        a,b,c = sy(a),sy(b),sy(c)
                        for k in xyzuv:
                            assert isinstance(k, Symbol)
                        for k in [a,b,c]:
                            for t in k.atoms(type=Symbol):
                                assert t in xyzuv
                        args = a,b,c, xyzuv
                        cls = ColorSchemeSympy

                    except: raise e

            else: raise e

        if cls == None:
            raise e

        obj = object.__new__(cls)
        obj.args = args
        return obj
    
    def __repr__(self):
        args = list(self.args[::])
        if len(args) == 4 and isinstance(args[3], (list, tuple)):
            args[3] = tuple(a for a in args[3] if str(a).find('unbound') < 0)
        return ",".join(str(a) for a in args)
    
    def __str__(self):
        return repr(self)

class ColorSchemeFlat(ColorScheme):

    uniform = True

    def __init__(self, *a):
        rgb = self.args
        assert len(rgb) == 3
        rgb = list(rgb)
        for i in range(3):
            rgb[i] = float(rgb[i])
            assert rgb[i] >= 0.0
            assert rgb[i] <= 1.0
        self.rgb = rgb

    def __call__(self, x,y,z,u,v):
        return self.rgb

class ColorSchemeLambda(ColorScheme):

    def __init__(self, *a):
        evaluator_str = "lambda x,y,z,u,v: [%s]" % (self.args[0])
        self.evaluator = eval(evaluator_str)

    def __call__(self, x,y,z,u,v):
        return self.evaluator(x,y,z,u,v)

class ColorSchemeSympy(ColorScheme):

    def __init__(self, *a):
        r,g,b,self.xyzuv = self.args
        self.functions = [r,g,b]

    def __call__(self, x,y,z,u,v):
        d = dict(zip(self.xyzuv, (x,y,z,u,v)))
        return list(f.subs_dict(d) for f in self.functions)

if __name__ == "__main__":
    x,y,z,u,v,t = symbols('xyzuvt')

    s = [ ColorScheme('z','y','x'),
          ColorScheme(0.5,0.6,0.7),
          ColorScheme('0.5,0.6,0.7'),
          ColorScheme(u,u,u, [u]),
          ColorScheme(u,v,0, (u,v)),
          ColorScheme(z,y,x, (x,y,z)),
          ColorScheme(z,y*v,x*u, (x,y,z,u,v)),
          ColorScheme('z, y, x')  ]

    for c in s:
        print "%s - %s" % (str(c), c(1,2,3,10,10))

