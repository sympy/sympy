# File: SineIntegrationRules.py
#
# This file is a Python code generated from the publicly available rule-based
# integrator (Rubi) system created by Albert Rich available at
# http://www.apmaths.uwo.ca/~arich/. This file is under the modified BSD
# #TODO: license (see the LICENSE file for more details).

from sympy import simplify, log
from sympy.core.symbol import Symbol, symbols
from sympy.core.symbol import Wild

def rubi_integrate(self, x):
    """
    rubi_integrate decides which category of rule does the
    given integrand (self) falls in and returns to that rule.
    """
    if self.is_Mul:
        if self.args[0].has(x) == False and len(self.args) == 2:
            # condition to check if the first arg is a constant,
            # returns in an is_Pow form so that it matches better.
            return (self.args[0]*int_alg_1_1(self.args[1], x))
    if self.is_Pow:
        return (int_alg_1_1(self, x))

def int_alg_1_1(self, x):
    # rule int_alg_1_1 indicates Algebraic > Linear Product > (a + b*x)**m category
    a = Wild('a', exclude = [x])
    b = Wild('b', exclude = [x])
    m = Wild('m', exclude = [x])
    simplify(self)
    type = self.match_between(((b*x)**m), (a + b*x)**m)
    b = type[b]
    m = type[m]
    a = type[a]
    if (m == -1):
        return log(a + b*x)/b
    elif (m != -1):
        return (a + b*x)**(m + 1)/(b*(m + 1))
    else:
        return Piecewise((log(a + b*x)/b, (m == -1)),((a + b*x)**(m+1)/(b*(m+1)), True))
