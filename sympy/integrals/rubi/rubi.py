from sympy.external import import_module
matchpy = import_module("matchpy")
from sympy.utilities.decorator import doctest_depends_on
from sympy.core import Integer, Float
import inspect, re

if matchpy:
    from matchpy import (Operation, CommutativeOperation, AssociativeOperation,
        ManyToOneReplacer, OneIdentityOperation, CustomConstraint)
    from matchpy.expressions.functions import register_operation_iterator, register_operation_factory
    from sympy import Pow, Add, Integral, Basic, Mul, S
    from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf,
        exp, log, gamma, acosh, asinh, atanh, acoth, acsch, asech, cosh, sinh,
        tanh, coth, sech, csch, atan, acsc, asin, acot, acos, asec, fresnels,
        fresnelc, erfc, erfi)

    Operation.register(Integral)
    register_operation_iterator(Integral, lambda a: (a._args[0],) + a._args[1], lambda a: len((a._args[0],) + a._args[1]))

    Operation.register(Pow)
    OneIdentityOperation.register(Pow)
    register_operation_iterator(Pow, lambda a: a._args, lambda a: len(a._args))

    Operation.register(Add)
    OneIdentityOperation.register(Add)
    CommutativeOperation.register(Add)
    AssociativeOperation.register(Add)
    register_operation_iterator(Add, lambda a: a._args, lambda a: len(a._args))

    Operation.register(Mul)
    OneIdentityOperation.register(Mul)
    CommutativeOperation.register(Mul)
    AssociativeOperation.register(Mul)
    register_operation_iterator(Mul, lambda a: a._args, lambda a: len(a._args))

    Operation.register(exp)
    register_operation_iterator(exp, lambda a: a._args, lambda a: len(a._args))

    Operation.register(log)
    register_operation_iterator(log, lambda a: a._args, lambda a: len(a._args))

    Operation.register(gamma)
    register_operation_iterator(gamma, lambda a: a._args, lambda a: len(a._args))

    Operation.register(fresnels)
    register_operation_iterator(fresnels, lambda a: a._args, lambda a: len(a._args))

    Operation.register(fresnelc)
    register_operation_iterator(fresnelc, lambda a: a._args, lambda a: len(a._args))

    Operation.register(erfc)
    register_operation_iterator(erfc, lambda a: a._args, lambda a: len(a._args))

    Operation.register(erfi)
    register_operation_iterator(erfi, lambda a: a._args, lambda a: len(a._args))

    Operation.register(sin)
    register_operation_iterator(sin, lambda a: a._args, lambda a: len(a._args))

    Operation.register(cos)
    register_operation_iterator(cos, lambda a: a._args, lambda a: len(a._args))

    Operation.register(tan)
    register_operation_iterator(tan, lambda a: a._args, lambda a: len(a._args))

    Operation.register(cot)
    register_operation_iterator(cot, lambda a: a._args, lambda a: len(a._args))

    Operation.register(csc)
    register_operation_iterator(csc, lambda a: a._args, lambda a: len(a._args))

    Operation.register(sec)
    register_operation_iterator(sec, lambda a: a._args, lambda a: len(a._args))

    Operation.register(sinh)
    register_operation_iterator(sinh, lambda a: a._args, lambda a: len(a._args))

    Operation.register(cosh)
    register_operation_iterator(cosh, lambda a: a._args, lambda a: len(a._args))

    Operation.register(tanh)
    register_operation_iterator(tanh, lambda a: a._args, lambda a: len(a._args))

    Operation.register(coth)
    register_operation_iterator(coth, lambda a: a._args, lambda a: len(a._args))

    Operation.register(csch)
    register_operation_iterator(csch, lambda a: a._args, lambda a: len(a._args))

    Operation.register(sech)
    register_operation_iterator(sech, lambda a: a._args, lambda a: len(a._args))

    Operation.register(asin)
    register_operation_iterator(asin, lambda a: a._args, lambda a: len(a._args))

    Operation.register(acos)
    register_operation_iterator(acos, lambda a: a._args, lambda a: len(a._args))

    Operation.register(atan)
    register_operation_iterator(atan, lambda a: a._args, lambda a: len(a._args))

    Operation.register(acot)
    register_operation_iterator(acot, lambda a: a._args, lambda a: len(a._args))

    Operation.register(acsc)
    register_operation_iterator(acsc, lambda a: a._args, lambda a: len(a._args))

    Operation.register(asec)
    register_operation_iterator(asec, lambda a: a._args, lambda a: len(a._args))

    Operation.register(asinh)
    register_operation_iterator(asinh, lambda a: a._args, lambda a: len(a._args))

    Operation.register(acosh)
    register_operation_iterator(acosh, lambda a: a._args, lambda a: len(a._args))

    Operation.register(atanh)
    register_operation_iterator(atanh, lambda a: a._args, lambda a: len(a._args))

    Operation.register(acoth)
    register_operation_iterator(acoth, lambda a: a._args, lambda a: len(a._args))

    Operation.register(acsch)
    register_operation_iterator(acsch, lambda a: a._args, lambda a: len(a._args))

    Operation.register(asech)
    register_operation_iterator(asech, lambda a: a._args, lambda a: len(a._args))

    def sympy_op_factory(old_operation, new_operands, variable_name):
         return type(old_operation)(*new_operands)

    register_operation_factory(Basic, sympy_op_factory)

    @doctest_depends_on(modules=('matchpy',))
    def rubi_object():
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
        #from sympy.integrals.rubi.rules.special_function import special_function
        #from sympy.integrals.rubi.rules.derivative import derivative
        from sympy.integrals.rubi.rules.piecewise_linear import piecewise_linear
        #from sympy.integrals.rubi.rules.miscellaneous_integration import miscellaneous_integration

        rubi = ManyToOneReplacer()
        #rubi = integrand_simplification(rubi)
        rubi = linear_products(rubi)
        rubi = quadratic_products(rubi)
        rubi = binomial_products(rubi)
        #rubi = trinomial_products(rubi)
        #rubi = miscellaneous_algebraic(rubi)
        #rubi = exponential(rubi)
        #rubi = logarithms(rubi)
        #rubi = sine(rubi)
        #rubi = tangent(rubi)
        #rubi = secant(rubi)
        #rubi = miscellaneous_trig(rubi)
        #rubi = inverse_trig(rubi)
        #rubi = hyperbolic(rubi)
        #rubi = inverse_hyperbolic(rubi)
        #rubi = piecewise_linear(rubi)
        #rubi = miscellaneous_integration(rubi)

        return rubi

    rubi = rubi_object()

@doctest_depends_on(modules=('matchpy',))
def rubi_integrate(expr, var, showsteps=False):
    '''
    Rule based algorithm for integration. Integrates the expression by applying
    transformation rules to the expression.

    Returns `Integrate` if an expression cannot be integrated.

    Parameters
    ==========
    expr : integrand expression
    var : variable of integration

    Returns Integral object if unable to integrate.
    '''
    if isinstance(expr, (int, Integer)) or isinstance(expr, (float, Float)):
        return S(expr)*var

    result = rubi.replace(Integral(expr, var))

    return result


@doctest_depends_on(modules=('matchpy',))
def get_matching_rule_definition(expr, var):
    '''
    Prints the list or rules which match to `expr`.

    Parameters
    ==========
    expr : integrand expression
    var : variable of integration
    '''
    matcher = rubi.matcher
    miter = matcher.match(Integral(expr, var))
    for fun, e in miter:
        print("Rule matching: ")
        print(inspect.getsourcefile(fun))
        code, lineno = inspect.getsourcelines(fun)
        print("On line: ", lineno)
        print("\n".join(code))
        print("Pattern matching: ")
        pattno = int(re.match(r"^\s*rule(\d+)", code[0]).group(1))
        print(matcher.patterns[pattno-1])
        print(e)
        print()
