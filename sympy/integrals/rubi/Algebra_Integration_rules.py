# File: SineIntegrationRules.py
#
# This file is a Python code generated from the publicly available rule-based
# integrator (Rubi) system created by Albert Rich available at
# http://www.apmaths.uwo.ca/~arich/. This file is under the modified BSD
# #TODO: license (see the LICENSE file for more details).

from sympy import simplify
from sympy.core.symbol import Symbol, symbols
from sympy.core.symbol import Wild

def int_alg_1_1(self, x):
    # rule Algebraic_1_1 indicates Algebraic > Linear Product > (a + b*x)**m category
    a = Wild('a', exclude = [x])
    b = Wild('b', exclude = [x])
    m = Wild('m', exclude = [x])
    simplify(self)
    type = self.match_between(((b*x)**m), (a + b*x)**m)
    b = type[b]
    m = type[m]
    if type[a] == None:
        a = 0
    else:
        a = type[a]
    if (m == -1):
        return log(a + b*x)/b
    elif (m != -1):
        return (a + b*x)**(m+1)/(b*(m+1))
    else:
        return Piecewise((log(a + b*x)/b, (m == -1)),((a + b*x)**(m+1)/(b*(m+1)), True))
