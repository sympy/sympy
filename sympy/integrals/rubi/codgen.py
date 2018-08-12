import re
import os

from sympy.core.expr import Basic
from sympy.integrals.rubi.rubi import rubi_object
from sympy.integrals.rubi.rubi import *
from sympy import *
from sympy.core.singleton import Singleton

from matchpy.utils import get_short_lambda_source
from matchpy.matching.code_generation import CodeGenerator
import re
class RubiCodeGenerator(CodeGenerator):
    def final_label(self, pattern_index, subst_name):
        label = self._matcher.patterns[pattern_index][1]
        if label is not None:
            return label
        elif label is None:
            return super().final_label(pattern_index, subst_name)
    def constraint_repr(self, constraint):
        if isinstance(constraint, CustomConstraint) and isinstance(constraint.constraint, type(lambda: 0)):
            src = get_short_lambda_source(constraint.constraint)
            mapping = {k: v for v, k in constraint._variables.items() }
            params = constraint._variables.keys()
            pstr = r'\b({})\b'.format('|'.join(map(re.escape, params)))
            new_src = re.sub(pstr, lambda m: 'subst{}[{!r}]'.format(self._substs, constraint._variables[m[0]]), src)
            return new_src, False
        return super().constraint_repr(constraint)
    def expr(self, expr):
        if isinstance(type(expr), Singleton):
            return 'S({!r})'.format(expr)
        return repr(expr)

    def get_args(self, operation, operation_type):
        if issubclass(operation_type, Integral):
            return '({0}._args[0],) + {0}._args[1]'.format(operation)
        if issubclass(operation_type, Basic):
            return '{}._args'.format(operation)
        return super().get_args(operation, operation_type)

GENERATED_TEMPLATE = '''
# -*- coding: utf-8 -*-

from sympy import *
from matchpy import *
from sympy.integrals.rubi.utility_function import *
from sympy.integrals.rubi.constraints import *
# from sympy.integrals.rubi.symbol import *
{}
{}
'''.strip()

def generate_code(rubi, output_file):
    generator = RubiCodeGenerator(rubi)
    global_code, code = generator.generate_code()
    code = GENERATED_TEMPLATE.format(global_code, code)

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(code)

def generate():
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
    # from sympy.integrals.rubi.rules.derivative import derivative
    # from sympy.integrals.rubi.rules.piecewise_linear import piecewise_linear
    from sympy.integrals.rubi.rules.miscellaneous_integration import miscellaneous_integration

    from matchpy import ManyToOneMatcher
    rules_applied = []

    rubi = ManyToOneMatcher()
    rubi = integrand_simplification(rules_applied, rubi, True)[0]
    rubi = linear_products(rules_applied, rubi, True)[0]
    rubi = quadratic_products(rules_applied, rubi, True)[0]
    generate_code(rubi, 'generated_1.py')

    # rubi = ManyToOneMatcher()
    # rubi = binomial_products(rules_applied, rubi, True)[0]
    # rubi = trinomial_products(rules_applied, rubi, True)[0]
    # generate_code(rubi, 'generated_2.py')

    # rubi = ManyToOneMatcher()
    # rubi = miscellaneous_algebraic(rules_applied, rubi, True)[0]
    # generate_code(rubi, 'generated_3.py')
