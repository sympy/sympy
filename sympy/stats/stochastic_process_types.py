from sympy import Indexed, Mul, Eq, sympify, S, Lambda
from sympy.stats.frv import SingleFinitePSpace
from sympy.stats.frv_types import BernoulliDistribution
from sympy.stats.joint_rv import (JointDistributionHandmade, JointPSpace,
    StochasticProcess)
from sympy.functions.elementary.piecewise import Piecewise

class BernoulliProcess(StochasticProcess):
    def __init__(self, name, P, succ=1, fail=0):
        self.name = sympify(name)
        self.P = sympify(P)
        self.succ = sympify(succ)
        self.fail = sympify(fail)

    def __getitem__(self, key):
        return SingleFinitePSpace(Indexed(self.name, key), BernoulliDistribution(
            self.P, self.succ, self.fail)).value

    def joint_dist(self, *keys):
        if len(keys) == 1:
            return BernoulliDistribution(self.P, self.succ, self.fail)
        syms = [Indexed(self.name, i) for i in range(len(keys))]
        pdf = Lambda(syms, Mul(*[Piecewise((self.P, Eq(x, self.succ)),
            (S.One - self.P, Eq(x, self.fail)), (0, True)) for x in syms]))
        return JointDistributionHandmade(pdf)
