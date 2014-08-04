from sympy.core import pi, oo, symbols, Function, Rational, Integer, GoldenRatio, EulerGamma, Catalan, Lambda, Dummy, Eq
from sympy.functions import Piecewise, sin, cos, Abs, exp, ceiling, sqrt, gamma
from sympy.utilities.pytest import raises
from sympy.printing.rcode import RCodePrinter
from sympy.utilities.lambdify import implemented_function
from sympy.tensor import IndexedBase, Idx
from sympy.printing.fcode import fcode, FCodePrinter

# import test
from sympy import rcode

i, m = symbols('i m', integer=True, cls=Dummy)
x = IndexedBase('x')
y = IndexedBase('y')
i = Idx(i, m)
print((x[i].base.label))
print((x[i].indices))
print(m.dummy_index)
print(i.label.dummy_index)
code = rcode(x[i], assign_to=y[i])
print(code)
expected = (
        'for (i_%(icount)i in 1:m_%(mcount)i){\n'
    '   y[i_%(icount)i] = x[i_%(icount)i];\n'
    '}'
) % {'icount': i.label.dummy_index, 'mcount': m.dummy_index}
print(expected)
