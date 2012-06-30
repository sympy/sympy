from sympy.ntheory import isprime, sieve

def _check_cycles_alt_sym(perm):
    """
    Checks for cycles of prime length p with n/2 < p < n-2.

    Helper function for the function is_alt_sym.

    See Also
    ========
    is_alt_sym
    """
    n = perm.size
    af = perm.array_form
    current_len = 0
    total_len = 0
    used = set()
    for i in xrange(n//2):
        if not i in used and i < n//2 - total_len:
            current_len = 1
            used.add(i)
            j = i
            while(af[j] != i):
                current_len += 1
                j = af[j]
                used.add(j)
            total_len += current_len
            if current_len > n//2 and current_len < n-2 and isprime(current_len):
                return True
    return False

def _distribute_gens_by_base(base, gens):
    k = len(base)
    S = []
    for i in xrange(k):
        S.append([])
    num_gens = len(gens)
    # for each generator, find the index of the
    # smallest (fixing the largest number of points)
    # basic stabilizer it belongs to
    stab_index = [0]*num_gens
    for i in xrange(num_gens):
        j = 0
        while j < k and gens[i](base[j]) == base[j]:
            j += 1
        stab_index[i] = j
    # distribute generators according to basic stabilizers
    for i in xrange(num_gens):
        index = stab_index[i]
        for j in xrange(index+1):
            S[j].append(gens[i])
    return S

def _orbits_transversals_from_bsgs(base, gens_distributed, transversals_only=False):
    from sympy.combinatorics.perm_groups import PermutationGroup
    k = len(base)
    transversals = [None]*k
    if transversals_only is False:
        basic_orbits = [None]*k
    for i in xrange(k):
        group = PermutationGroup(gens_distributed[i])
        transversals[i] = dict(group.orbit_transversal(base[i], pairs=True))
        if transversals_only is False:
            basic_orbits[i] = transversals[i].keys()
    if transversals_only:
        return transversals
    else:
        return basic_orbits, transversals

def _handle_precomputed_bsgs(base, strong_gens, transversals, basic_orbits, gens_distributed):
    if gens_distributed is None:
        gens_distributed = _distribute_gens_by_base(base, strong_gens)
    if transversals is None:
        if basic_orbits is None:
            basic_orbits, transversals = _orbits_transversals_from_bsgs(base, gens_distributed)
        else:
            transversals = _orbits_transversals_from_bsgs(base, gens_distributed, transversals_only=True)
    else:
        if basic_orbits is None:
            base_len = len(base)
            basic_orbits = [None]*base_len
            for i in xrange(base_len):
                basic_orbits[i] = transversals[i].keys()
    return transversals, basic_orbits, gens_distributed
