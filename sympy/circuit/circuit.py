# TODO : Store the final solution such that Voltage at node i can be accessed by V[i]
#        and Current through element between node u and v as I[(u, v)]
# TODO : Device implementation of mutually Coupled-Inductors, Controlled V-I sources
#        and corresponding changes in submatrices definition.
# TODO : Parsing netlist with 5th and 6th optional parameter for eg: ac/dc for V-I sources
#        and initial conditions for capacitors and inductors. Currently support only
#        for DC sources.
# TODO : Check the C matrix for opamp.
# TODO : Extend all the device classes to support future analyses.
#        Currently they just do book-keeping. And document them.
# TODO : Method for printing netlist.
# TODO : Limit line length to 80.
# TODO : Remove commands used for testing once done with them.
# TODO : Print method for all classes

__all__ = ['Circuit', 'Resistor', 'Capacitor', 'Inductor', 'VoltageSource', 'CurrentSource',
           'OpAmp', 'solve', 'parse_netlist', '_g_matrix', '_b_matrix', '_c_matrix',
           '_d_matrix', '_e_matrix', '_j_matrix', '_v_matrix', '_i_matrix', '_A_matrix',
           '_z_matrix']

from sympy import Symbol, pprint, Matrix
from sympy.matrices.dense import zeros


class Circuit:
    """
    Class for making a circuit object. Circuit can be instantiated in the
    following two ways --

    1. Pass a netlist file as the input to parse_netlist function.
    2. Instantiate by passing a list containing all the elements.

    """
    def __init__(self, elements):

        self.rlc_elements = []
        self.v_sources = []
        self.i_sources = []
        self.opamps = []
        self.node_count = 0
        for e in elements:
            if isinstance(e, Resistor) or isinstance(e, Inductor) or isinstance(e, Capacitor):
                self.rlc_elements.append(e)

            if isinstance(e, VoltageSource):
                self.v_sources.append(e)

            if isinstance(e, CurrentSource):
                self.i_sources.append(e)

            if isinstance(e, OpAmp):
                self.opamps.append(e)

        self.G = Matrix([])
        self.B = Matrix([])
        self.C = Matrix([])
        self.D = Matrix([])
        self.E = Matrix([])
        self.I = Matrix([])
        self.V = Matrix([])
        self.J = Matrix([])
        self.A = Matrix([])
        self.Z = Matrix([])
        self.x = Matrix([])
        self.Voltage = []
        self.Current = {}


def solve(circuit):
    """
    Solves the circuit for unknown node Voltages and unknown Currents flowing through
    Voltage Sources. For i independent nodes and j Voltage sources in the circuit, it
    solved circuit with first i elements of x Matrix as Voltages at
    i nodes and next j elements as Currents through j Voltage sources.

    Refer :: http://www.swarthmore.edu/NatSci/echeeve1/Ref/mna/MNA_All.html

    >>> from sympy.circuit import Circuit, VoltageSource, Inductor, Capacitor, Resistor, solve
    >>> elements = []
    >>> elements.append(VoltageSource('V1', '1', '0', '10'))
    >>> elements.append(Inductor('L1', '1', '2', '0.1'))
    >>> elements.append(Capacitor('C1', '2', '3', '0.0001'))
    >>> elements.append(Resistor('R1', '3', '0', '1000'))
    >>> my_cir = Circuit(elements)
    >>> solution = solve(my_cir)
    >>> solution.G
    [ 1/(L1*s),       -1/(L1*s),           0]
    [-1/(L1*s), C1*s + 1/(L1*s),       -C1*s]
    [        0,           -C1*s, C1*s + 1/R1]

    >>> solution.x
    [                                         V1]
    [V1*(C1*R1*s + 1)/(C1*L1*s**2 + C1*R1*s + 1)]
    [          C1*R1*V1*s/(C1*s*(L1*s + R1) + 1)]
    [            C1*V1*s/(-C1*s*(L1*s + R1) - 1)]

    """
    for element in circuit.rlc_elements:
        circuit.node_count = max([int(element.node1), int(element.node2), circuit.node_count])
    for element in circuit.v_sources:
        circuit.node_count = max([int(element.node1), int(element.node2), circuit.node_count])
    for element in circuit.i_sources:
        circuit.node_count = max([int(element.node1), int(element.node2), circuit.node_count])
    for element in circuit.opamps:
        circuit.node_count = max([int(element.node1), int(element.node2), int(element.node3), int(element.node4), circuit.node_count])
    circuit.node_count = circuit.node_count + 1

    circuit.G = _g_matrix(circuit.rlc_elements, circuit.node_count)
    circuit.B = _b_matrix(circuit.v_sources, circuit.opamps, circuit.node_count)
    circuit.C = _c_matrix(circuit.v_sources, circuit.opamps, circuit.node_count)
    circuit.D = _d_matrix(circuit.v_sources, circuit.opamps)
    circuit.E = _e_matrix(circuit.v_sources)
    circuit.J = _j_matrix(circuit.v_sources, circuit.opamps, circuit.node_count)
    circuit.V = _v_matrix(circuit.v_sources, circuit.node_count)
    circuit.I = _i_matrix(circuit.i_sources, circuit.node_count)

    circuit.A = _A_matrix(circuit.G, circuit.B, circuit.C, circuit.D)
    circuit.Z = _z_matrix(circuit.I, circuit.E)

    circuit.x = circuit.A.inv()*circuit.Z
    circuit.x.simplify()

    circuit.Voltage.append(0)
    for i in range(circuit.node_count - 1):
        circuit.Voltage.append(circuit.x[i, 0])
    for i in range(circuit.node_count, len(circuit.x)):
        circuit.Current[circuit.Z[i]] = circuit.x[i]
    for i in circuit.rlc_elements:
        if int(i.node1) > int(i.node2)
            circuit.Current[i.symbol] = (circuit.Voltage[int(i.node1)] - circuit.Voltage[int(i.node2)])/i.symbol
        if int(i.node2) > int(i.node1)
            circuit.Current[i.symbol] = (circuit.Voltage[int(i.node2)] - circuit.Voltage[int(i.node1)])/i.symbol
    return circuit



class Resistor:

    def __init__(self, name, node1, node2, value):
        self.symbol = Symbol(name)
        self.node1 = node1
        self.node2 = node2
        self.value = float(value)


class Inductor:

    def __init__(self, name, node1, node2, value):
        self.symbol = Symbol(name)
        self.node1 = node1
        self.node2 = node2
        self.value = float(value)


class Capacitor:

    def __init__(self, name, node1, node2, value):
        self.symbol = Symbol(name)
        self.node1 = node1
        self.node2 = node2
        self.value = float(value)

class VoltageSource:

    def __init__(self, name, node1, node2, value):
        self.symbol = Symbol(name)
        self.node1 = node1
        self.node2 = node2
        self.value = float(value)


class CurrentSource:

    def __init__(self, name, node1, node2, value):
        self.symbol = Symbol(name)
        self.node1 = node1
        self.node2 = node2
        self.value = float(value)


class OpAmp:

    def __init__(self, name, node1, node2, node3, node4, gain):
        self.symbol = Symbol(name)
        self.node1 = node1
        self.node2 = node2
        self.node3 = node3
        self.node4 = node4
        self.gain = float(gain)

def parse_netlist(netlist_file):
    """
    Parses a netlist file that describse a circuit and returns a list
    containing all the circuit elements.

    For more info on Netlist file :
    http://www.allaboutcircuits.com/vol_5/chpt_7/8.html
    """
    f = open(netlist_file)

    elements = []

    for line in f:

        if line[0] != 'O':
            name, node1, node2, value = line.split()

            if name[0] == 'R':
                element = Resistor(name, node1, node2, value)
                elements.append(element)

            if name[0] == 'L':
                element = Inductor(name, node1, node2, value)
                elements.append(element)

            if name[0] == 'C':
                element = Capacitor(name, node1, node2, value)
                elements.append(element)

            if name[0] == 'V':
                element = VoltageSource(name, node1, node2, value)
                elements.append(element)

            if name[0] == 'I':
                element = CurrentSource(name, node1, node2, value)
                elements.append(element)

        else:
            name, node1, node2, node3, node4, gain = line.split()
            element = OpAmp(name, node1, node2, node3, node4, gain)
            elements.append(element)

    f.close()
    return Circuit(elements)


# s = sigma + j*w --> The symbol used in Laplace Domain.
s = Symbol('s')


def _g_matrix(elements, node_count):
    """
    The Conductance matrix

    For a circuit with n nodes and m independent voltages :

    The G matrix is an nxn matrix formed in two steps --
    1. Each element in the diagonal matrix is equal to the sum of the
    conductance(one over the resistance) of each element connected to the
    corresponding node.  So the first diagonal element is the sum of
    conductances connected to node 1, the second diagonal element is the
    sum of conductances connected to node 2, and so on.

    2. The off diagonal elements are the negative conductance of the element
    connected to the pair of corresponding node.  Therefore a resistor between
    nodes 1 and 2 goes into the G matrix at location (1,2) and locations (2,1).
    """

    G = zeros(node_count - 1, node_count - 1)

    for e in elements:
        n1 = int(e.node1)
        n2 = int(e.node2)

        if isinstance(e, Resistor):
            g_elem = 1/e.symbol
        if isinstance(e, Inductor):
            g_elem = 1/(s*e.symbol)
        if isinstance(e, Capacitor):
            g_elem = s*e.symbol

        if n1 !=0 and n2 != 0:
            G[n1 - 1, n2 - 1] = G[n1 - 1, n2 - 1] - g_elem
            G[n2 - 1, n1 - 1] = G[n2 - 1, n1 - 1] - g_elem
        if n1 != 0:
            G[n1 - 1, n1 - 1] = G[n1 - 1, n1 - 1] + g_elem

        if  n2 != 0:
            G[n2 - 1, n2 - 1] = G[n2 - 1, n2 - 1] + g_elem
    return G


def _b_matrix(v_sources, opamps, node_count):
    """
    The B matrix is an nxm matrix with only 0, 1 and -1 elements. Each
    location in the matrix corresponds to a particular voltage source
    (first dimension) or a node (second dimension).  If the positive
    terminal of the ith voltage source is connected to node k, then the
    element (i,k) in the B matrix is a 1.  If the negative terminal of
    the ith voltage source is connected to node k, then the element (i,k)
    in the B matrix is a -1.  Otherwise, elements of the B matrix are zero.
    The output port of an OpAmp is also considered as a voltage source.
    """

    v_count = len(v_sources)
    o_count = len(opamps)

    if v_count + o_count == 0:
        return

    B = zeros(node_count - 1, v_count + o_count)

    for i in range(v_count):
        for j in range(node_count-1):
            if int(v_sources[i].node1) == j + 1:
                B[j, i] = 1
            elif int(v_sources[i].node2) == j + 1:
                B[j, i] = -1

    for i in range(o_count):
        for j in range(node_count - 1):
            if int(opamps[i].node1) == j + 1:
                B[j, i + v_count] = 1
            elif int(opamps[i].node2) == j + 1:
                B[j, i + v_count] = -1
            else:
                B[j, i + v_count] = 0

    return B


def _c_matrix(v_sources, opamps, node_count):
    """
    The C matrix is an nxm matrix with only 0, 1 and -1 elements. Each
    location in the matrix corresponds to a particular node (first dimension)
    or voltage source (second dimension).  For each indendent voltage source,
    if the positive terminal of the ith voltage source is connected to node k,
    then the element (k,i) in the C matrix is a 1; if the negative terminal
    of the ith voltage source is connected to node k, then the element (k,i)
    in the C matrix is a -1.  For each op-amp let the positive input  terminal
    be at node k and negative terminal at node j;  the corresponding (ith) row
    of the C matrix has a 1 at location corresponding to the positive terminal
    (k,i), and a -1 at the location corresponding to the negative terminal
    (j,i). Rest of the elements of the C matrix are zero.
    """

    v_count = len(v_sources)
    o_count = len(opamps)

    if v_count + o_count == 0:
        return

    C = zeros(v_count + o_count, node_count - 1)

    for i in range(v_count):
        for j in range(node_count - 1):
            if int(v_sources[i].node1) == j + 1:
                C[i, j] = 1
            elif int(v_sources[i].node2) == j + 1:
                C[i, j] = -1

    for i in range(o_count):
        for j in range(node_count - 1):
            if int(opamps[i].node1) == j + 1:
                C[i + v_count, j] = 1
            elif int(opamps[i].node2) == j + 1:
                C[i + v_count, j] = -1
            elif int(opamps[i].node3) == j + 1:
                C[i + v_count, j] = 1
            elif int(opamps[i].node4) == j + 1:
                C[i + v_count, j] = -1
            else:
                C[i + v_count, j] = 0

    return C


def _d_matrix(v_sources, opamps):
    """
    The D matrix is an mxm matrix that is composed entirely of zeros
    if the circuit contains only independent voltage sources.
    It can be non-zero if dependent sources are considered.
    """

    #TODO : Revamp this function when Controlled Voltage sources will be
    #       supported.

    v_count = len(v_sources)
    o_count = len(opamps)

    itm_count = v_count + o_count
    D = zeros(itm_count, itm_count)

    return D


def _e_matrix(v_sources):
    """
    The e matrix is an mx1 matrix with each element of the
    matrix equal in value to the corresponding independent voltage source.
    """

    v_count = len(v_sources)

    E = zeros(v_count, 1)

    for i in range(v_count):
        E[i] = v_sources[i].symbol

    return E


def _i_matrix(i_sources, node_count):
    """
    The i matrix is an nx1 matrix with each element of the matrix
    corresponding to a particular node.  The value of each element
    of i is determined by the sum of current sources into the
    corresponding node.  If there are no current sources connected
    to the node, the value is zero.
    """

    I = zeros(node_count - 1, 1)

    for q in range(node_count - 1):
        for i_src in i_sources:
            if int(i_src.node1) == q + 1:
                I[q, 0] = I[q, 0] - i_src.symbol
            elif int(i_src.node2) == q + 1:
                I[q, 0] = I[q, 0] + i_src.symbol

    return I


def _v_matrix(v_sources, node_count):
    """
    The v matrix is an nx1 matrix formed of the node voltages. Each element
    in v corresponds to the voltage at the equivalent node in the circuit
    (there is no entry for ground -- node 0).
    """

    V = zeros(node_count - 1, 1)

    for itm_no in range(node_count - 1):
        V[itm_no] = Symbol('V_' + str(itm_no + 1))

    return V


def _j_matrix(v_sources, opamps, node_count):
    """
    The j matrix is an mx1 matrix, with one entry for the current
    through each voltage source.
    """

    v_count = len(v_sources)
    o_count = len(opamps)

    if v_count + o_count == 0:
        return

    J = zeros(v_count + o_count, 1)

    for i in range(v_count):
        J[i] = Symbol('I_' + str(v_sources[i].symbol))

    for i in range(o_count):
        J[i + v_count] = Symbol('I_' + str(opamps[i].symbol))

    return J


def _A_matrix(G, B, C, D):
    """
    The A matrix:
    -- is (n+m)x(n+m) in size, and consists only of known quantities.
    -- the nxn part of the matrix in the upper left has only passive elements;
    elements connected to ground appear only on the diagonal;
    elements not connected to ground are both on the diagonal and off-diagonal
    terms.
    -- the rest of the A matrix (not included in the nxn upper left part)
    contains only 1, -1 and 0 (other values are possible if there are dependent
    current and voltage sources.
    """

    g_rows = G.shape[0]
    g_cols = G.shape[1]

    c_rows = C.shape[0]
    b_cols = B.shape[1]

    row_count = g_rows + c_rows
    col_count = g_cols + b_cols

    A = zeros(row_count, col_count)

    A[:g_rows, :g_cols] = G
    A[g_rows:, :g_cols] = C
    A[:g_rows, g_cols:] = B
    A[g_rows:, g_cols:] = D

    return A


def _z_matrix(I, E):
    """
    The z matrix
    -- is an (n+m)x1 vector that holds only known quantities.
    -- the top n elements are either zero or the sum and difference of
    independent current sources in the circuit.
    -- the bottom m elements represent the m independent voltage sources
    in the circuit.
    """

    row_count = I.shape[0] + E.shape[0]
    col_count = I.shape[1]

    Z = zeros(row_count, col_count)

    Z[:I.shape[0], :] = I
    Z[I.shape[0]:, :] = E

    return Z
