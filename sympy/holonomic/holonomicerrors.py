""" Common Exceptions for `holonomic` module. """

from __future__ import print_function, division

class BaseHolonomicError(Exception):

    def new(self, *args):
        raise NotImplementedError("abstract base class")

class NotPowerSeriesError(BaseHolonomicError):

    def __init__(self, holonomic, x0):
        self.holonomic = holonomic
        self.x0 = x0

    def __str__(self):
        s = 'A Power Series doesnot exists for '
        s += str(self.holonomic)
        s += ' about %s.' %self.x0
        return s

class NotHolonomicError(BaseHolonomicError):

    def __init__(self, m):
        self.m = m

    def __str__(self):
        return self.m
