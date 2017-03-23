# File: SineIntegrationRules.py
#
# This file is a Python code generated from the publicly available rule-based
# integrator (Rubi) system created by Albert Rich available at
# http://www.apmaths.uwo.ca/~arich/. This file is under the modified BSD
# #TODO: license (see the LICENSE file for more details).

from sympy import simplify

def Algebraic_1_1(self)
    # rule Algebraic_1_1 indicates Algebraic > Linear Product > (a + b*x)**m category
    a = wild('a', exclude = [x])
    b = wild('b', exclude = [x])
    m = wild('m', exclude = [x])
    simplify(self)
    type_1 = self.matches((b*x)**m)
        if type:
            a = 0
        if (m == -1):
            return log(a + b*x)/b
        elif (m != -1):
            Return (a + b*x)**(m+1)/(b*(m+1))
        else
            return Piecewise((log(a + b*x)/b, (m == -1)),((a + b*x)**(m+1)/(b*(m+1)), True))
