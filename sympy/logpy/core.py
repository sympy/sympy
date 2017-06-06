#######################
# Simplification code #
#######################
from sympy import assuming, ask
try:
    import sympy.logpy.train
    from logpy import variables, var, run, goalify
    from logpy.assoccomm import eq_assoccomm as eqac
    logpy_exists = True
    # TODO - make sure this doesn't happen prematurely
    asko = goalify(ask)
except ImportError:
    logpy_exists = False


def refine_one(expr, *assumptions, **kwargs):
    if not logpy_exists:
        raise ImportError("Please install LogPy")
    reduces = kwargs['reduces']
    vars = kwargs['vars']
    with assuming(*assumptions):
        with variables(*vars):
            source, target, condition = var(), var(), var()
            result = run(1, target, (reduces, source, target, condition),
                                    (eqac, source, expr),
                                    (asko, condition, True))
    return result[0] if result else expr
