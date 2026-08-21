"""

qasm.py - Functions to parse a set of qasm commands into a SymPy Circuit.

Examples taken from Chuang's page: https://web.archive.org/web/20220120121541/https://www.media.mit.edu/quanta/qasm2circ/

The code returns a circuit and an associated list of labels.

>>> from sympy.physics.quantum.qasm import Qasm
>>> q = Qasm('qubit q0', 'qubit q1', 'h q0', 'cnot q0,q1')
>>> q.get_circuit()
CNOT(1,0)*H(1)

>>> q = Qasm('qubit q0', 'qubit q1', 'cnot q0,q1', 'cnot q1,q0', 'cnot q0,q1')
>>> q.get_circuit()
CNOT(1,0)*CNOT(0,1)*CNOT(1,0)
"""
from __future__ import annotations

__all__ = [
    'Qasm',
    ]

from math import prod
from typing import TYPE_CHECKING

from sympy.core.singleton import S as _S
from sympy.physics.quantum.gate import H, CNOT, X, Z, CGate, CGateS, SWAP, S, T,CPHASE
from sympy.physics.quantum.circuitplot import Mz

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator

    from sympy.core.expr import Expr

def read_qasm(lines: str) -> Qasm:
    return Qasm(*lines.splitlines())

def read_qasm_file(filename: str) -> Qasm:
    return Qasm(*open(filename).readlines())

def flip_index(i: int, n: int) -> int:
    """Reorder qubit indices from largest to smallest.

    >>> from sympy.physics.quantum.qasm import flip_index
    >>> flip_index(0, 2)
    1
    >>> flip_index(1, 2)
    0
    """
    return n-i-1

def trim(line: str) -> str:
    """Remove everything following comment # characters in line.

    >>> from sympy.physics.quantum.qasm import trim
    >>> trim('nothing happens here')
    'nothing happens here'
    >>> trim('something #happens here')
    'something '
    """
    if '#' not in line:
        return line
    return line.split('#')[0]

def get_index(target: str, labels: list[str]) -> int:
    """Get qubit labels from the rest of the line,and return indices

    >>> from sympy.physics.quantum.qasm import get_index
    >>> get_index('q0', ['q0', 'q1'])
    1
    >>> get_index('q1', ['q0', 'q1'])
    0
    """
    nq = len(labels)
    return flip_index(labels.index(target), nq)

def get_indices(targets: list[str], labels: list[str]) -> list[int]:
    return [get_index(t, labels) for t in targets]

def nonblank(args: tuple[str, ...]) -> Iterator[str]:
    for line in args:
        line = trim(line)
        if line.isspace():
            continue
        yield line
    return

def fullsplit(line: str) -> tuple[str, list[str]]:
    words = line.split()
    rest = ' '.join(words[1:])
    return fixcommand(words[0]), [s.strip() for s in rest.split(',')]

def fixcommand(c: str) -> str:
    """Fix Qasm command names.

    Remove all of forbidden characters from command c, and
    replace 'def' with 'qdef'.
    """
    forbidden_characters = ['-']
    c = c.lower()
    for char in forbidden_characters:
        c = c.replace(char, '')
    if c == 'def':
        return 'qdef'
    return c

def stripquotes(s: str) -> str:
    """Replace explicit quotes in a string.

    >>> from sympy.physics.quantum.qasm import stripquotes
    >>> stripquotes("'S'") == 'S'
    True
    >>> stripquotes('"S"') == 'S'
    True
    >>> stripquotes('S') == 'S'
    True
    """
    s = s.replace('"', '') # Remove second set of quotes?
    s = s.replace("'", '')
    return s

class Qasm:
    """Class to form objects from Qasm lines

    >>> from sympy.physics.quantum.qasm import Qasm
    >>> q = Qasm('qubit q0', 'qubit q1', 'h q0', 'cnot q0,q1')
    >>> q.get_circuit()
    CNOT(1,0)*H(1)
    >>> q = Qasm('qubit q0', 'qubit q1', 'cnot q0,q1', 'cnot q1,q0', 'cnot q0,q1')
    >>> q.get_circuit()
    CNOT(1,0)*CNOT(0,1)*CNOT(1,0)
    """
    def __init__(self, *args: str, **kwargs: str) -> None:
        self.defs: dict[str, Callable[..., Expr]] = {}
        self.circuit: list[Expr] = []
        self.labels: list[str] = []
        self.inits: dict[str, str] = {}
        self.add(*args)
        self.kwargs = kwargs

    def add(self, *lines: str) -> None:
        for line in nonblank(lines):
            command, rest = fullsplit(line)
            if command in self.defs: #defs come first, since you can override built-in
                function = self.defs[command]
                indices = self.indices(rest)
                if len(indices) == 1:
                    self.circuit.append(function(indices[0]))
                else:
                    self.circuit.append(function(indices[:-1], indices[-1]))
            elif hasattr(self, command):
                function = getattr(self, command)
                function(*rest)
            else:
                print("Function %s not defined. Skipping" % command)

    def get_circuit(self) -> Expr:
        return prod(reversed(self.circuit), start=_S.One)

    def get_labels(self) -> list[str]:
        return list(reversed(self.labels))

    def plot(self) -> None:
        from sympy.physics.quantum.circuitplot import CircuitPlot
        circuit, labels = self.get_circuit(), self.get_labels()
        CircuitPlot(circuit, len(labels), labels=labels, inits=self.inits)

    def qubit(self, arg: str, init: str | None = None) -> None:
        self.labels.append(arg)
        if init: self.inits[arg] = init

    def indices(self, args: list[str]) -> list[int]:
        return get_indices(args, self.labels)

    def index(self, arg: str) -> int:
        return get_index(arg, self.labels)

    def nop(self, *args: str) -> None:
        pass

    def x(self, arg: str) -> None:
        self.circuit.append(X(self.index(arg)))

    def z(self, arg: str) -> None:
        self.circuit.append(Z(self.index(arg)))

    def h(self, arg: str) -> None:
        self.circuit.append(H(self.index(arg)))

    def s(self, arg: str) -> None:
        self.circuit.append(S(self.index(arg)))

    def t(self, arg: str) -> None:
        self.circuit.append(T(self.index(arg)))

    def measure(self, arg: str) -> None:
        self.circuit.append(Mz(self.index(arg)))

    def cnot(self, a1: str, a2: str) -> None:
        self.circuit.append(CNOT(*self.indices([a1, a2])))

    def swap(self, a1: str, a2: str) -> None:
        self.circuit.append(SWAP(*self.indices([a1, a2])))

    def cphase(self, a1: str, a2: str) -> None:
        self.circuit.append(CPHASE(*self.indices([a1, a2])))

    def toffoli(self, a1: str, a2: str, a3: str) -> None:
        i1, i2, i3 = self.indices([a1, a2, a3])
        self.circuit.append(CGateS((i1, i2), X(i3)))

    def cx(self, a1: str, a2: str) -> None:
        fi, fj = self.indices([a1, a2])
        self.circuit.append(CGate(fi, X(fj)))

    def cz(self, a1: str, a2: str) -> None:
        fi, fj = self.indices([a1, a2])
        self.circuit.append(CGate(fi, Z(fj)))

    def defbox(self, *args: str) -> None:
        print("defbox not supported yet. Skipping: ", args)

    def qdef(self, name: str, ncontrols: str, symbol: str) -> None:
        from sympy.physics.quantum.circuitplot import CreateOneQubitGate, CreateCGate
        ncontrols = int(ncontrols)
        command = fixcommand(name)
        symbol = stripquotes(symbol)
        if ncontrols > 0:
            self.defs[command] = CreateCGate(symbol)
        else:
            self.defs[command] = CreateOneQubitGate(symbol)
