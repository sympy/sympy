from sympy.computations.core import Computation, unique, CompositeComputation

a,b,c,d,e,f,g,h = 'abcdefgh'

def test_Computation():
    C = Computation()
    assert hasattr(C, 'inputs')
    assert hasattr(C, 'outputs')
    assert hasattr(C, 'edges')

def test_unique():
    assert tuple(unique((1, 3, 1, 2))) == (1, 3, 2)

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

class TComposite(CompositeComputation):
    """ Test CompositeComputation class """
    def __init__(self, *computations):
        self.computations = computations

    def __str__(self):
        return "[[" + ", ".join(map(str, self.computations)) + "]]"

def test_testcomp():
    A = TComp('add', (a, b, c), (d,))
    assert A.inputs == (a, b, c)
    assert tuple(A.edges()) == ((a, A), (b, A), (c, A), (A, d))

def test_composite():
    A = TComp('add', (a, b, c), (d,))
    M = TComp('mul', (d, e), (f,))
    C = TComposite(A, M)
    assert tuple(C.inputs)  == (a, b, c, e)
    assert tuple(C.outputs) == (f,)
    assert tuple(C.edges()) == ((a, A), (b, A), (c, A), (A, d), (d, M), (e, M),
            (M, f))

def test_composite_dag():
    A = TComp('add', (a, b, c), (d,))
    M = TComp('mul', (d, e), (f,))
    C = TComposite(A, M)

    assert C.dag_io() == {A: set([M]), M: set()}
    assert C.dag_oi() == {M: set([A]), A: set()}

def test_toposort():
    A = TComp('add', (a, b, c), (d,))
    M = TComp('mul', (d, e), (f,))
    C = TComposite(A, M)

    assert tuple(C.toposort()) == (A, M)

def test_multi_out():
    MM = TComp('minmax', (a, b), (d, e))
    A =  TComp('foo', (d,), (f,))
    B =  TComp('bar', (a, f), (g, h))
    C =  TComposite(MM, A, B)
    assert C.inputs == (a, b)
    assert C.outputs == (e, g, h)
    assert tuple(C.toposort()) == (MM, A, B)
