
from basic import Basic

class ArithMeths:

    def __pos__(self):
        return self
    def __add__(self, other):
        return Basic.Add(self, other)
