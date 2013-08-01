"""

qasm.py - Functions to parse a set of qasm commands into a Sympy Circuit.

Examples taken from Chuang's page: http://www.media.mit.edu/quanta/qasm2circ/

Todo:
* Put in subscripts??
* Figure out the def boxes

The code returns a circuit and an associated list of labels.
>>> from sympy.physics.quantum.qasm import Qasm
>>> q = Qasm('qubit q0','qubit q1','h q0','cnot q0,q1')
>>> q.get_circuit()
CNOT(1,0)*H(1)

>>> q = Qasm('qubit q0','qubit q1','cnot q0,q1','cnot q1,q0','cnot q0,q1')
>>> q.get_circuit()
CNOT(1,0)*CNOT(0,1)*CNOT(1,0)
"""
from sympy.physics.quantum.gate import H, CNOT, X, Z, CGate, SWAP,S,T
from sympy.physics.quantum.circuitplot import Mz

def prod(c):
    def mul(a,b): return a*b # Can't import operator.mul b/c operator module in directory
    return reduce(mul,c,1)

def flip_index(i,n):
    """Reorder qubit indices from largest to smallest.
    >>> from sympy.physics.quantum.qasm import flip_index
    >>> flip_index(0,2)
    1
    >>> flip_index(1,2)
    0
    """
    return n-i-1

def isblank(line):
    """Returns True if the line is empty/blank/all whitespace.
    >>> from sympy.physics.quantum.qasm import isblank
    >>> isblank('   ')
    True
    >>> isblank(' _ ')
    False
    """
    return len(line.split()) == 0

def trim(line):
    """Remove everything following comment # characters in line.
    >>> from sympy.physics.quantum.qasm import trim
    >>> trim('nothing happens here')
    'nothing happens here'
    >>> trim('something #happens here')
    'something '
    """
    if not '#' in line: return line
    return line.split('#')[0]

def get_index(target,labels):
    """Get qubit labels from the rest of the line,
    and return their indices, properly flipped.
    >>> from sympy.physics.quantum.qasm import get_index
    >>> get_index('q0',['q0','q1'])
    1
    >>> get_index('q1',['q0','q1'])
    0
    """
    nq = len(labels)
    return flip_index(labels.index(target),nq)

def get_indices(targets,labels):
    return [get_index(t,labels) for t in targets]

def nonblank(args):
    for line in args:
        line = trim(line)
        if isblank(line): continue
        yield line
    return

def fullsplit(line):
    words = line.split()
    rest = ' '.join(words[1:])
    return fixcommand(words[0]),rest.split(',')

def fixcommand(c):
    """Fix Qasm command names.
    Remove all of forbidden characters from command c, and
    replace 'def' with 'qdef'.
    """
    forbidden_characters = ['-']
    for char in forbidden_characters:
        c = c.replace(char,'')
    if c == 'def':
        return 'qdef'
    return c

def stripquotes(s):
    """Replace explicit quotes in a string.
    >>> from sympy.physics.quantum.qasm import stripquotes
    >>> stripquotes("'S'") == 'S'
    True
    >>> stripquotes('"S"') == 'S'
    True
    >>> stripquotes('S') == 'S'
    True
    """
    s = s.replace('"','') # Remove second set of quotes?
    s = s.replace("'",'')
    return s

class Qasm(object):
    """
    >>> from sympy.physics.quantum.qasm import Qasm
    >>> q = Qasm('qubit q0','qubit q1','h q0','cnot q0,q1')
    >>> q.get_circuit()
    CNOT(1,0)*H(1)
    >>> q = Qasm('qubit q0','qubit q1','cnot q0,q1','cnot q1,q0','cnot q0,q1')
    >>> q.get_circuit()
    CNOT(1,0)*CNOT(0,1)*CNOT(1,0)
    """
    def __init__(self,*args,**kwargs):
        self.defs = {}
        self.circuit = []
        self.labels = []
        self.add(*args)
        self.kwargs = kwargs

    def add(self,*lines):
        for line in nonblank(lines):
            command,rest = fullsplit(line)
            if hasattr(self,command):
                function = getattr(self,command)
                function(*rest)
            elif self.defs.get(command):
                function = self.defs.get(command)
                fi = self.index(rest[0])
                self.circuit.append(function(fi))
            else:
                print "Function %s not defined. Skipping"

    def get_circuit(self): return prod(reversed(self.circuit))
    def get_labels(self): return list(reversed(self.labels))

    def plot(self):
        from sympy.physics.quantum.circuitplot import CircuitPlot
        circuit,labels = self.get_circuit(), self.get_labels()
        CircuitPlot(circuit,len(labels),labels=labels)

    def qubit(self,arg): self.labels.append(arg)
    def indices(self,args): return get_indices(args,self.labels)
    def index(self,arg): return get_index(arg,self.labels)
    def nop(self,*args): pass
    def x(self,arg): self.circuit.append(X(self.index(arg)))
    def z(self,arg): self.circuit.append(Z(self.index(arg)))
    def h(self,arg): self.circuit.append(H(self.index(arg)))
    def s(self,arg): self.circuit.append(S(self.index(arg)))
    def t(self,arg): self.circuit.append(T(self.index(arg)))
    def measure(self,arg): self.circuit.append(Mz(self.index(arg)))

    def cnot(self,a1,a2):self.circuit.append(CNOT(*self.indices([a1,a2])))
    def swap(self,a1,a2):self.circuit.append(SWAP(*self.indices([a1,a2])))
    def cphase(self,a1,a2):self.circuit.append(CPhase(*self.indices([a1,a2])))

    def cx(self,a1,a2):
        fi,fj = self.indices([a1,a2])
        self.circuit.append(CGate(fi,X(fj)))
    def cz(self,a1,a2):
        fi,fj = self.indices([a1,a2])
        self.circuit.append(CGate(fi,Z(fj)))

    def qdef(self,name,nq,symbol):
        from sympy.physics.quantum.circuitplot import CreateOneQubitGate
        nq = int(nq)
        command = fixcommand(name)
        symbol = stripquotes(symbol)
        if nq > 1:
            print "Def for nq>1 not defined ",nq
        else:
            self.defs[command] = CreateOneQubitGate(symbol)

if __name__ == '__main__':
    import doctest; doctest.testmod()
