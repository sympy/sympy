from utils import memoizer_immutable_args
from basic import Basic, sympify
from number import Rational


class Integer(Rational, int):

    is_integer = True

    @memoizer_immutable_args('Integer.__new__')
    def __new__(cls, p):
        obj = int.__new__(cls, p)
        if p==0:
            obj.is_zero = True            
        elif p==1:
            obj.is_one = True
        elif p==2:
            obj.is_two = True
        return obj

    @property
    def p(self): return int(self)

    @property
    def q(self): return 1

    # relational methods

    def compare(self, other):
        if self is other: return 0
        c = cmp(self.__class__, other.__class__)
        if c: return c
        return int.__cmp__(self, int(other))

    # converter methods

    # int has __int__, __float__

    def evalf(self):
        return Basic.Float(int(self))

    # mathematical properties

    @property
    def is_even(self):
        return int.__mod__(self,2)==0

    @property
    def is_odd(self):
        return int.__mod__(self,2)==1

    @property
    def is_positive(self):
        return int(self)>0

    @property
    def is_negative(self):
        return int(self)<0
