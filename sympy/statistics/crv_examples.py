from sympy import exp, sqrt, pi, S, Dummy
from crv import SingleContinuousPSpace

class NormalPSpace(SingleContinuousPSpace):
    _count = 0
    _name = 'x'
    def __new__(cls, mean, std, symbol = None):

        x = symbol or cls.create_symbol()
        pdf = exp(-(x-mean)**2 / (2*std**2)) / (sqrt(2*pi)*std)
        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        obj.mean = mean
        obj.std = std
        obj.variance = std**2
        return obj

def Normal(mean, std, symbol=None):
    return NormalPSpace(mean, std, symbol).value
