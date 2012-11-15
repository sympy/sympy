from sympy.computations.core import Computation, unique, CompositeComputation
from sympy.rules.rl import flatten, unpack

a,b,c,d,e,f,g,h = 'abcdefgh'

class TComp(Computation):
    """ Test Computation class """
    def __init__(self, op, inputs, outputs):
        self.op = op
        self.inputs = tuple(inputs)
        self.outputs = tuple(outputs)

    def __str__(self):
        ins  = "["+', '.join(self.inputs) +"]"
        outs = "["+', '.join(self.outputs)+"]"
        return "%s -> %s -> %s"%(ins, str(self.op), outs)


def test_flatten():
    MM = TComp('minmax', (a, b), (d, e))
    A =  TComp('foo', (d,), (f,))
    B =  TComp('bar', (a, f), (g, h))
    C =  CompositeComputation(MM, A, B)
    C2 = MM+A+B
    assert len(C2.computations) == 2
    assert len(flatten(C2).computations) == 3
    assert C.inputs == C2.inputs
    assert C.outputs == C2.outputs

def test_unpack():

    A =  TComp('foo', (d,), (f,))
    C =  CompositeComputation(A)
    assert unpack(C) == A
