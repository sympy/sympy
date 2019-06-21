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

if not os.path.exists('generated.py'):
    r = rubi_object()
    rubi = r[0]
    rules =  r[2]
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

    generator = RubiCodeGenerator(rubi)
    global_code, code = generator.generate_code()
    code = GENERATED_TEMPLATE.format(global_code, code)

    with open('generated.py', 'w', encoding='utf-8') as f:
        f.write(code)

from generated import match_root
x = symbols('x')
ex = Integral(x , x)
for r, p in match_root(ex):
    print(r, p)

ex1 = Integral(x**2, x)
for r, p in match_root(ex1):
    print(r, p)