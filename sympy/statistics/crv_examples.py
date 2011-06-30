from sympy import exp, sqrt, pi, S, Dummy
from crv import SingleContinuousPSpace

class NormalPSpace(SingleContinuousPSpace):
    def __new__(cls, mean, var, symbol = None):

        x = symbol or Dummy('x', real=True, finite=True)
        pdf = exp(-(x-mean)**2 / (2*var)) / (sqrt(2*pi*var))
        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        obj.mean = mean
        obj.variance = var
        return obj
