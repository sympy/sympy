from sympy.utilities.decorator import doctest_depends_on

rules = []
rules_applied = []

@doctest_depends_on(modules=('matchpy',))
def get_rubi_object():
    '''
    Returns rubi ManyToOneReplacer by adding all rules from different modules.

    Uncomment the lines to add integration capabilities of that module.

    Currently, there are parsing issues with special_function,
    derivative and miscellaneous_integration. Hence they are commented.
    '''
    from sympy.integrals.rubi.rules.integrand_simplification import integrand_simplification
    from sympy.integrals.rubi.rules.linear_products import linear_products
    from sympy.integrals.rubi.rules.quadratic_products import quadratic_products
    from sympy.integrals.rubi.rules.binomial_products import binomial_products
    from sympy.integrals.rubi.rules.trinomial_products import trinomial_products
    from sympy.integrals.rubi.rules.miscellaneous_algebraic import miscellaneous_algebraic
    from sympy.integrals.rubi.rules.exponential import exponential
    from sympy.integrals.rubi.rules.logarithms import logarithms
    from sympy.integrals.rubi.rules.sine import sine
    from sympy.integrals.rubi.rules.tangent import tangent
    from sympy.integrals.rubi.rules.secant import secant
    from sympy.integrals.rubi.rules.miscellaneous_trig import miscellaneous_trig
    from sympy.integrals.rubi.rules.inverse_trig import inverse_trig
    from sympy.integrals.rubi.rules.hyperbolic import hyperbolic
    from sympy.integrals.rubi.rules.inverse_hyperbolic import inverse_hyperbolic
    from sympy.integrals.rubi.rules.special_functions import special_functions
    #from sympy.integrals.rubi.rules.derivative import derivative
    #from sympy.integrals.rubi.rules.piecewise_linear import piecewise_linear
    from sympy.integrals.rubi.rules.miscellaneous_integration import miscellaneous_integration

    from matchpy import ManyToOneReplacer

    from multiprocessing import cpu_count
    from concurrent.futures import ProcessPoolExecutor as PoolExecutor

    rule_functions = [
        integrand_simplification,
        linear_products,
        quadratic_products,
        binomial_products,
        trinomial_products,
        miscellaneous_algebraic,
        exponential,
        logarithms,
        special_functions,
        sine,
        tangent,
        secant,
        miscellaneous_trig,
        inverse_trig,
        hyperbolic,
        inverse_hyperbolic,
        # piecewise_linear,
        miscellaneous_integration
    ]

    # Number of processors is equal to mp.cpu_count
    # pool = Pool(cpu_count() - 2) 

    def get_integration_rule(integration_rule_func):
        global rules
        global rules_applied
        rules += integration_rule_func(rules_applied)
    
    # pool.map(get_integration_rule, rule_functions)
    with PoolExecutor(max_workers=cpu_count() - 2) as executor:
        executor.map(get_integration_rule, rule_functions)

    rubi = ManyToOneReplacer(*rules)
    return rubi, rules_applied, rules_applied

if __name__ == '__main__':
    get_rubi_object()