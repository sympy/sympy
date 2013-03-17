from sympy.printing.theanocode import theano_code
from sympy import Symbol
import theano

def test_symbol():
    xs = Symbol('x')
    xt = theano_code(xs)
    assert isinstance(xt, theano.scalar.ScalarVariable)
    assert xt.name == xs.name
