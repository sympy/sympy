"""Primitive circuit operations on quantum circuits."""

__all__ = [
    'kmp_table',
    'find_subcircuit',
    'remove_subcircuit'
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

def find_subcircuit(circuit, subcircuit, start=0, end=0):
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

def remove_subcircuit(circuit, subcircuit, pos=0):
    """Removes subcircuit from circuit, if it exists.

    If multiple instances of subcircuit exists, the
    first instance is removed.  A location to check may
    be optionally given.  If subcircuit can't be found,
    circuit is returned.

    Parameters
    ==========
    circuit : tuple, Gate
        A quantum circuit represented by a tuple of Gates
    subcircuit : tuple, Gate
        A quantum circuit to remove from circuit
    pos : int
        The location to remove subcircuit, if it exists.
        This may be used if it is known beforehand that
        multiple instances exist, and it is desirable
        to remove a specific instance.  If a negative number
        is given, pos will be defaulted to 0.
    """

    if pos < 0:
        pos = 0

    # Look for the subcircuit starting at pos
    loc = find_subcircuit(circuit, subcircuit, start=pos)

    # If subcircuit was found
    if loc > -1:
        # Get the gates to the left of subcircuit
        left = circuit[0:loc]
        # Get the gates to the right of subcircuit
        right = circuit[loc + len(subcircuit):len(circuit)]
        # Recombine the left and right side gates into a circuit
        circuit = left + right

    return circuit

