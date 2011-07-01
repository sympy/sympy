from sympy import exp, sqrt, pi, S, Dummy
from crv import SingleContinuousPSpace

class NormalPSpace(SingleContinuousPSpace):
    _count = 0
    _name = 'x'
    def __new__(cls, mean, var, symbol = None):

        x = symbol or cls.create_symbol()
        pdf = exp(-(x-mean)**2 / (2*var)) / (sqrt(2*pi*var))
        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        obj.mean = mean
        obj.variance = var
        return obj

def Normal(mean, variance, symbol=None):
    return NormalPSpace(mean, variance, symbol).value
