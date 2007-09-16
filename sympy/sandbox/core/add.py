
from basic import Composite
from methods import ArithMeths

class Add(ArithMeths, Composite):

    @classmethod
    def canonize(cls, *args, **kwds):
        args = tuple(map(cls.sympify, args))
        return (args,), kwds
