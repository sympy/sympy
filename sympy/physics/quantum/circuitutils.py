"""Primitive circuit operations on quantum circuits.

TODO:
* Add wrappers around the seq-based functions to use Mul
"""

from sympy import Wild, Integer, Tuple
from sympy.physics.quantum.gate import Gate

__all__ = [
    'kmp_table',
    'find_subcircuit_with_seq',
    'replace_subcircuit_with_seq',
    'conv2_symbolic_qubits_with_seq'
]

def kmp_table(word):
    """Build the 'partial match' table of the
       Knuth-Morris-Pratt algorithm.

    Note: This is applicable to strings or quantum circuits.
    """

    # Current position in subcircuit
    pos = 2
    # Beginning position of candidate substring that
    # may reappear later in word
    cnd = 0
    # The 'partial match' table that helps one determine
    # the next location to start substring search
    table = list()
    table.append(-1)
    table.append(0)

    while pos < len(word):
        if word[pos-1] == word[cnd]:
            cnd = cnd + 1
            table.append(cnd)
            pos = pos + 1
        elif cnd > 0:
            cnd = table[cnd]
        else:
            table.append(0)
            pos = pos + 1

    return table

def find_subcircuit_with_seq(circuit, subcircuit, start=0, end=0):
    """Finds the subcircuit in circuit, if it exists.

    If the subcircuit exists, the index of the start of
    the subcircuit in circuit is returned; otherwise,
    -1 is returned.  The algorithm that is implemented
    is the Knuth-Morris-Pratt algorithm.

    Parameters
    ==========
    circuit : tuple, Gate
        A tuple of Gates representing a quantum circuit
    subcircuit : tuple, Gate
        A tuple of Gates to check if circuit contains
    start : int
        The location to start looking for subcircuit.
        If start is the same or past end, -1 is returned.
    end : int
        The last place to look for a subcircuit.  If end
        is less than 1 (one), then the length of circuit
        is taken to be end.
    """

    if len(subcircuit) == 0 or len(subcircuit) > len(circuit):
        return -1

    if end < 1:
        end = len(circuit)

    # Location in circuit
    pos = start
    # Location in the subcircuit
    index = 0
    # 'Partial match' table
    table = kmp_table(subcircuit)

    while (pos + index) < end:
        if subcircuit[index] == circuit[pos + index]:
            index = index + 1
        else:
            pos = pos + index - table[index]
            index = table[index] if table[index] > -1 else 0

        if index == len(subcircuit):
            return pos

    return -1

def replace_subcircuit_with_seq(circuit, subcircuit, replace=(), pos=0):
    """Replaces a subcircuit with another subcircuit in circuit,
    if it exists.

    If multiple instances of subcircuit exists, the
    first instance is replaced.  A location to check may
    be optionally given.  If subcircuit can't be found,
    circuit is returned.

    Parameters
    ==========
    circuit : tuple, Gate
        A quantum circuit represented by a tuple of Gates
    subcircuit : tuple, Gate
        The circuit to be replaced
    replace : tuple, Gate
        The replacement circuit
    pos : int
        The location to start search and replace
        subcircuit, if it exists.  This may be used
        if it is known beforehand that multiple
        instances exist, and it is desirable to
        replace a specific instance.  If a negative number
        is given, pos will be defaulted to 0.
    """

    if pos < 0:
        pos = 0

    # Look for the subcircuit starting at pos
    loc = find_subcircuit_with_seq(circuit, subcircuit, start=pos)

    # If subcircuit was found
    if loc > -1:
        # Get the gates to the left of subcircuit
        left = circuit[0:loc]
        # Get the gates to the right of subcircuit
        right = circuit[loc + len(subcircuit):len(circuit)]
        # Recombine the left and right side gates into a circuit
        circuit = left + replace + right

    return circuit

def next_symbolic_index(cur):
    """Returns the "next" symbolic index."""

    # Don't really like this approach of converting Wild
    # to string.  Advise finding to way to deal with it.
    symb_str = repr(cur)
    next_ndx = Integer(symb_str[1:len(symb_str)-1]) + 1
    return Wild('i%r' % next_ndx)

def conv2_symbolic_qubits_with_seq(*seq, **kargs):
    """Returns the circuit with symbolic indices and the
    dictionary mapping symbolic indices to real indices.

    The mapping is 1 to 1 and onto (bijective).

    Parameters
    ==========
    seq : tuple, Gate/Integer/tuple
        A tuple of Gate, Integer, or tuple objects
    kargs : dict
        A dictionary of optional arguments.  It is only
        use to pass in an existing mapping of symbolic
        indices to real indices and a starting symbolic
        index.  The type of the starting symbolic index
        must be a Wild.

    All symbolic indices have the format 'i#', where # is
    some number >= 0.
    """

    cur_ndx = Wild('i-1')
    # keys are symbolic indices; values are real indices
    ndx_map = {}

    def create_inverse_map(symb_to_real_map):
        rev_items = lambda item: tuple([item[1], item[0]])
        return dict(map(rev_items, symb_to_real_map.items()))

    if 'start' in kargs:
        if not isinstance(kargs['start'], Wild):
            msg = 'Expected Wild for starting index, got %r.' % kargs['start']
            raise TypeError(msg)
        cur_ndx = kargs['start']

    if 'index_map' in kargs:
        if not isinstance(kargs['index_map'], dict):
            msg = ('Expected dict for existing map, got ' +
                   '%r.' % kargs['index_map'])
            raise TypeError(msg)
        ndx_map = kargs['index_map']

    # keys are real indices; keys are symbolic indices
    inv_map = create_inverse_map(ndx_map)

    sym_seq = ()
    for item in seq:
        # Nested items, so recurse
        if isinstance(item, Gate):
            sym_item, new_map, cur_ndx = conv2_symbolic_qubits_with_seq(
                                                 *item.args,
                                                 index_map=ndx_map,
                                                 start=cur_ndx)
            ndx_map.update(new_map)
            inv_map = create_inverse_map(ndx_map)

        elif isinstance(item, tuple) or isinstance(item, Tuple):
            sym_item, new_map, cur_ndx = conv2_symbolic_qubits_with_seq(
                                                 *item,
                                                 index_map=ndx_map,
                                                 start=cur_ndx)
            ndx_map.update(new_map)
            inv_map = create_inverse_map(ndx_map)

        elif item in inv_map:
            sym_item = inv_map[item]

        else:
            cur_ndx = next_symbolic_index(cur_ndx)
            ndx_map[cur_ndx] = item
            inv_map[item] = cur_ndx
            sym_item = cur_ndx

        if isinstance(item, Gate):
            sym_item = item.__class__(*sym_item)

        sym_seq = sym_seq + (sym_item,)

    return sym_seq, ndx_map, cur_ndx
