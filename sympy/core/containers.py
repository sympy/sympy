from basic import Basic

class Tuple(Basic):
    """
    Wrapper around the builtin tuple object

    The Tuple is a subclass of Basic, so that it works well in the
    Sympy framework.  The wrapped tuple is available as self.args, but
    you can also access elements or slices with [:] syntax.

    >>> from sympy import symbols
    >>> from sympy.physics.secondquant import Tuple
    >>> a, b, c, d = symbols('a b c d')
    >>> Tuple(a, b, c)[1:]
    Tuple(b, c)
    >>> Tuple(a, b, c).subs(a, d)
    Tuple(d, b, c)

    """

    def __getitem__(self,i):
        if isinstance(i,slice):
            indices = i.indices(len(self))
            return Tuple(*[self.args[i] for i in range(*indices)])
        return self.args[i]

    def __len__(self):
        return len(self.args)

    def __contains__(self,item):
        return item in self.args


def tuple_wrapper(method):
    """
    Decorator that makes any tuple in arguments into Tuple
    """
    def wrap_tuples(*args, **kw_args):
        newargs=[]
        for arg in args:
            if type(arg) is tuple:
                newargs.append(Tuple(*arg))
            else:
                newargs.append(arg)
        return method(*newargs, **kw_args)
    return wrap_tuples
