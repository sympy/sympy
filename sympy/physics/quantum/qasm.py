"""

qasm.py - Functions to parse a set of qasm commands into a Sympy Circuit.

Examples taken from Chuang's page: http://www.media.mit.edu/quanta/qasm2circ/

Todo:
* Put in subscripts??
* Figure out the def boxes

The code returns a circuit and an associated list of labels.
from sympy.physics.quantum.qasm import qasm
>>> qasm('qubit q0','qubit q1','h q0','cnot q0,q1')
(CNOT(1,0)*H(1), ['q1', 'q0'])

>>> qasm('qubit q0','qubit q1','qubit q2','h q1','cnot q1,q2','cnot q0,q1',\
         'h q0','nop q1','measure q0','measure q1','c-x q1,q2','c-z q0,q2')
(C((2),Z(0))*C((1),X(0))*Mz(1)*Mz(2)*H(2)*CNOT(2,1)*CNOT(1,0)*H(1), ['q2', 'q1', 'q0'])

>>> qasm('qubit q0','qubit q1','cnot q0,q1','cnot q1,q0','cnot q0,q1')
(CNOT(1,0)*CNOT(0,1)*CNOT(1,0), ['q1', 'q0'])
"""
from sympy.physics.quantum.gate import *
from sympy.physics.quantum.circuitplot import Mz

def prod(c):
    def mul(a,b): return a*b
    return reduce(mul,c,1)

def flip_index(i,n):
    """Reorder qubit indices from largest to smallest.
    from sympy.physics.quantum.qasm import flip_index
    >>> flip_index(0,2)
    1
    >>> flip_index(1,2)
    0
    """
    return n-i-1

def isblank(line):
    """Returns True if the line is empty/blank/all whitespace.
    from sympy.physics.quantum.qasm import isblank
    >>> isblank('   ')
    True
    >>> isblank(' _ ')
    False
    """
    return len(line.split()) == 0

def trim(line):
    """Remove everything following comment # characters in line.
    from sympy.physics.quantum.qasm import trim
    >>> trim('nothing happens here')
    'nothing happens here'
    >>> trim('something #happens here')
    'something '
    """
    if not '#' in line: return line
    return line.split('#')[0]

def get_indices(rest,labels):
    """Get qubit labels from the rest of the line,
    and return their indices, properly flipped.
    from sympy.physics.quantum.qasm import get_indices
    >>> get_indices('q0',['q0','q1'])
    1
    >>> get_indices('q1',['q0','q1'])
    0
    """
    nq = len(labels)
    targets = rest.split(',')
    indices = [labels.index(target) for target in targets]
    if len(indices) == 1: return flip_index(indices[0],nq)
    return [flip_index(i,nq) for i in indices]

def qasm(*args,**kwargs):
    circuit = []
    labels = []
    commands = ['qubit','h','cnot','c-x','c-z','nop','measure']
    two_qubit_commands = ['cnot','c-x','c-z']
    for line in args:
        line = trim(line)
        if isblank(line): continue
        words = line.split()
        command = words[0]
        rest = ' '.join(words[1:])
        if command not in commands:
            print "Skipping unknown/unparsed command: ",command
        if command == 'qubit':
            labels.append(words[1])
        elif command == 'h':
            fi = get_indices(rest,labels)
            circuit.append(H(fi))
        elif command == 'cnot':
            fi,fj = get_indices(rest,labels)
            circuit.append(CNOT(fi,fj))
        elif command == 'c-x':
            fi,fj = get_indices(rest,labels)
            circuit.append(CGate(fi,X(fj)))
        elif command == 'c-z':
            fi,fj = get_indices(rest,labels)
            circuit.append(CGate(fi,Z(fj)))
        elif command == 'nop':
            pass
        elif command == 'measure':
            fi = get_indices(rest,labels)
            circuit.append(Mz(fi))
    return prod(reversed(circuit)),list(reversed(labels))

if __name__ == '__main__':
    import doctest; doctest.testmod()
