import itertools
from sympy.combinatorics.homomorphisms import homomorphism, GroupHomomorphism
from sympy.combinatorics.fp_groups import FpGroup
from sympy.combinatorics.perm_groups import PermutationGroup
from sympy import S

def group_isomorphism(G, H):
    '''
    Compute the isomorphism between 2 given groups.

    Arguments:
        G (a finite `FpGroup` or a `PermutationGroup`) -- First group
        H (a finite `FpGroup` or a `PermutationGroup`) -- Second group

    Returns:
    If the groups are not isomorphic, `None` is returned.
    If the groups are isomorphic, the corresponding `GroupHomomorphism` object is returned.

    Summary:
    Uses the approach suggested by Robert Tarjan to compute the isomorphism between two groups.
    First, the set of generators are mapped with the elements of `H` and we check if the mapping induces isomorphism.

    '''

    if isinstance(G, FpGroup):
        if G.order() == S.Infinity:
            raise NotImplementedError("Isomorphism methods are not implemented for infinite groups.")
        G = G._to_perm_group()[0]
    if isinstance(H, FpGroup):
        if H.order() == S.Infinity:
            raise NotImplementedError("Isomorphism methods are not implemented for infinite groups.")
        H = H._to_perm_group()[0]

    if not isinstance(G, PermutationGroup):
        raise TypeError("The group must be a PermutationGroup or a finite FpGroup")
    if not isinstance(H, PermutationGroup):
        raise TypeError("The group must be a PermutationGroup or a finite FpGroup")

    if (len(G) != len(H)) and (G.is_abelian != H.is_abelian):
        return None

    # Match the generators of `G` with the subsets of H
    gens = G.generators
    for subset in itertools.permutations(H, len(gens)):
        # Create a map
        func_map = dict(zip(gens, subset))
        counterexample = False
        # Loop through the mapped elements and try to extend the mapping
        # or to find a counter example.
        while not counterexample:
            extended_map = {}
            for g, h in itertools.product(func_map, func_map):
                if g*h in func_map:
                    if func_map[g]*func_map[h] != func_map[g*h]:
                        counterexample = True
                        break
                else:
                    extended_map[g*h] = func_map[g]*func_map[h]

            # Break when all the elements are mapped.
            if len(func_map) == len(G):
                break

            # Remove the duplicate elements in extended_map
            # and merge the extended_map with the map.
            image_len = len(set(func_map.values()) + set(extended_map.values()))
            if image_len != len(func_map) + len(extended_map):
                counterexample = True
                break
            func_map.update(extended_map)

        if not counterexample:
            # Trivial homomorphism is computed.
            return homomorphism(G, H, gens)

    return None

def is_isomorphism(G, H):
    '''
    Check if the given groups are isomorphic.

    Arguments:
        G (a finite `FpGroup` or a `PermutationGroup`) -- First group
        H (a finite `FpGroup` or a `PermutationGroup`) -- Second group

    Returns:
    `True` if the groups are isomorphic
    `False` if the groups are not isomorphic
    '''
    return bool(group_isomorphism(G, H))
