#######################
# Simplification code #
#######################
import sympy.logpy.train
import logpy
from logpy.variables import variables
from logpy import var
from logpy.assoccomm import eq_assoccomm as eqac
from logpy.goals import goalify
from sympy import assuming, ask

# TODO - make sure this doesn't happen prematurely
asko = goalify(ask)

def refine_one(expr, *assumptions, **kwargs):
    reduces = kwargs['reduces']
    vars = kwargs['vars']
    with assuming(*assumptions):
        with variables(*vars):
            source, target, condition = var(), var(), var()
            result = logpy.run(1, target, (reduces, source, target, condition),
                                          (eqac, source, expr),
                                          (asko, condition, True))
    return result[0] if result else expr
