
from basic import Basic, Composite
from methods import ArithMeths

class Add(ArithMeths, Composite, dict):

    def __new__(cls, *args, **options):
        """
        Use dict.__new__(Add) to create a mutable instance
        but after you are finished then call ._froze() method.
        """
        obj = dict.__new__(cls)
        [obj.update(a) for a in args]
        obj._froze()
        return obj

    def __init__(self, *args, **options):
        pass

    def update(self, a):
        if isinstance(a, dict):
            for k,v in a.items():
                try:
                    self[k] += v
                except KeyError:
                    self[k] = v
        else:
            a = Basic.sympify(a)
            try:
                self[a] += 1
            except KeyError:
                self[a] = 1

    def _froze(self):
        """
        Make instance immutable.
        """
        self.update = self.raise_immutable
        #todo: overwrite __setitem__, pop, etc..
        #todo: delete _froze method.

    # do we need _unfroze method? It could be dangerous.

    @staticmethod
    def raise_immutable(*args, **kws):
        raise TypeError('Add instance is immutable')

    def torepr(self):
        return '%s(%s)' % (self.__class__.__name__, dict(self))
