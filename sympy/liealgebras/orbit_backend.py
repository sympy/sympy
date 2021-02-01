import warnings
from sympy.external import import_module
from sympy.utilities.iterables import multiset_permutations
from sympy import Rational, Matrix, zeros, ones
np = import_module('numpy')

def _cast_to_type(ary, outtype, backend):
    if outtype not in ['sympy', 'numpy']:
        raise ValueError(f"outtype must be 'sympy' or 'numpy' (got '{outtype}')")

    if outtype == "sympy" and backend == "numpy":
        results = []
        for mat in ary:
            results.append(
                Matrix([[Rational(x) for x in mat]]).as_immutable()
            )
        return results
    elif outtype == "numpy" and backend == "sympy":
        return np.array(ary)
    else:
        return ary

#
#   Orbit methods
#

def _orbit_backend_sympy(algebra, weight, stabilizer):
    """Returns the orbit stabilizer theorem calculation using
    sympy objects in an easy to follow, but unperformant, algorithm.
    The docstring on `cartan_base.StandardCartan.orbit` explains the args."""
    # Generating the rotation matrices from simple roots
    # to avoid recalculating them each loop
    reflection_matrices = algebra._reflection_matrices()

    if stabilizer is not None:
        reflection_matrices = [reflection_matrices[i] for i in stabilizer]

    # Fill up master list of reflected roots by
    # continually operating on them until all weights in orbit are found
    master_list = [weight.as_immutable()]
    master_list_hash = set()
    while True:
        # we use this because we can't change master_list during
        # iteration (python rule)
        ref_list = []
        for w in master_list:
            # if we've seen w before, we've also see all its reflections
            if w in master_list_hash:
                continue
            for refl in reflection_matrices:
                reflected = w * refl
                if reflected not in ref_list and reflected not in master_list:
                    ref_list.append(reflected)
                    master_list_hash.add(w)

        # No new reflections have been found
        if len(ref_list) == 0:
            break

        master_list += ref_list
    return master_list

def _orbit_backend_numpy(algebra, weight, stabilizer, dtype):
    """Runs the orbit stabilizer theorem using an effecient numpy api empowered
    algorithm. The docstring on `cartan_base.StandardCartan.orbit` explains the args.
    """
    reflection_matrices = np.array(algebra._reflection_matrices(), dtype=dtype, order='c')
    # if stabilizers are passed, select relfection matricies corresponding to those
    # simple root indexes
    stab_refl_mats = reflection_matrices
    if stabilizer:
        stabilizer = np.array(stabilizer, dtype=int)
        stab_refl_mats = reflection_matrices[stabilizer]

    num_pos = algebra.roots() // 2
    # used to rotate basis
    omega_inv = np.array(algebra.omega_matrix().pinv(), dtype=dtype, order="c")
    weight = np.array(weight, dtype=dtype, order='c').squeeze()
    # rotate to dominant weyl chamber

    dom = _rotate_to_dominant(reflection_matrices, omega_inv, weight)

    # generate full orbit based on dominant weight
    return _full_orbit(stab_refl_mats, dom, num_pos)


def orbit(algebra, weight, stabilizer=None, dtype=int, backend="sympy", outtype="sympy"):
    """Returns the orbit stabilizer theorem algorithm with
    the option of using sympy or numpy. For larger groups
    such as TypeA-TypeD rank 10 or any of TypeE groups,
    numpy is preferred. For published documentation, see
    `sympy.liealgebra.cartan_base.StandardCartan.orbit`"""

    if backend == "sympy":
        orbit_ = _orbit_backend_sympy(algebra, weight, stabilizer)
    elif backend == "numpy":
        if np is None:
            warnings.warn("Numpy is not installed, falling back to sympy")
            return orbit(algebra, weight, stabilizer, dtype=int, backend="sympy", outtype="sympy")

        orbit_ = _orbit_backend_numpy(algebra, weight, stabilizer, dtype)
    else:
        raise ValueError(f"backend must be 'sympy' or 'numpy' (got '{backend}')")

    return _cast_to_type(orbit_, outtype, backend)

def _reflect(refl_mats, weights):
    """Returns all unique reflections of weights by
    map applying reflection matrices over the weights
    """
    # Nested for loop that is flatten to (N, dim)
    dim = weights.shape[-1]
    refs = np.array([refl_mats.dot(w.squeeze()) for w in weights]).reshape(-1,dim)
    comb = np.concatenate((refs, weights),axis=0)
    return np.unique(comb, axis=0)


def _rotate_to_dominant(refl_mats, omega_inv, weight, typeA=False):
    """Returns the dominant weight by rotating
    repeatedly across weyl chambers until all
    elements are positive.
    """

    if typeA:
        orbits = np.array(list(multiset_permutations(weight)))
        dom = _find_dom(orbits, omega_inv)
        if len(dom) == 0:
            # Check for lack of indexes here
            raise IndexError("Cannot find dominant weight in given orbit")
        return dom

    orbits = np.expand_dims(weight,axis=0)
    while True:
        orbits = _reflect(refl_mats, orbits)
        dom = _find_dom(orbits, omega_inv)
        if len(dom) > 0:
            return dom.squeeze()

def _find_dom(orbits, omega_inv):
    """Returns the dominant weight in list of weights, if exists, else
    returns an empty array"""
    for x in orbits:
        if np.all(np.dot(x, omega_inv) >= 0):
            return x
    return np.array([], dtype=orbits.dtype)

def _full_orbit(refl_mats, weight, num_pos):
    """Returns the full orbit of a weight by
    operating on the weight by the reflection matrices
    a number of times equal to the number of positive
    roots in the algebra."""
    orbit = np.expand_dims(weight, axis=0)
    for _ in range(num_pos):
        orbit = _reflect(refl_mats, orbit)
    return orbit

#
#   Root system methods
#

def _roots_system_backend_sympy(algebra):
    """For published documentation, see
    `sympy.liealgebra.cartan_base.StandardCartan.root_system`"""
    s_r = algebra.simple_roots()
    rank = algebra.rank

    orbits = set()
    for i in s_r:
        for r in algebra.orbit(i):
            orbits.add(r)

    zero_roots = [zeros(1, rank)] * rank

    orbits = [i * algebra.cocartan_matrix().T for i in orbits] + zero_roots


    # rotate back to the orthogonal basis for consistency
    omega_matrix = algebra.omega_matrix()
    # sort roots by their weight level, then general positive number order (2nd one is )
    sorbits = sorted(orbits, key=lambda x: (-algebra.root_level(x, "alpha"), tuple(x)))

    return [x * omega_matrix for x in sorbits]

def _roots_system_backend_numpy(algebra, dtype):
    """Numpy implemented rootsystem algorithm.
    For published documentation, see
    `sympy.liealgebra.cartan_base.StandardCartan.root_system`"""
    sr = np.array(algebra.simple_roots(), dtype=dtype, order='c')
    rank = algebra.rank

    # used to rotate to omega basis from orthogonal
    rot = np.array(algebra.cocartan_matrix().T, dtype=float)

    # used to rotate to alpha basis from omega
    rot2 = algebra.cartan_matrix().pinv()
    rot2 = np.array(rot2 * ones(rot2.rows, 1))

    #used to rotate from omega to orthogonal
    rot3 = np.array(algebra.omega_matrix(), dtype=dtype, order='c')

    # Alias for orbit function with kwargs set.
    orbit_func = lambda x: orbit(algebra, x, dtype=dtype, backend="numpy", outtype="numpy")

    # Get all orbits for each simple root and take unique weights
    orbits = np.array([i for x in sr for i in orbit_func(x)])
    orbits = np.unique(orbits, axis=0)

    # add in zero roots and rotate to omega basis for sorting
    roots = np.concatenate((orbits.dot(rot), *[np.zeros((1, rank))] * rank))

    # sort in omega basis by alpha sum, then by order of indexes (last part is convention)
    roots = np.array(sorted(roots, key=lambda x: (-x.dot(rot2), tuple(x))))[::-1]

    # return in orthogonal basis
    return roots.dot(rot3)


def root_system(algebra, dtype=int, backend="sympy", outtype="sympy"):
    if backend == "sympy":
        roots = _roots_system_backend_sympy(algebra)
    elif backend == "numpy":
        if np is None:
            warnings.warn("Numpy is not installed, falling back to sympy")
            return root_system(algebra, dtype, "sympy", "sympy")
        roots = _roots_system_backend_numpy(algebra, dtype)
    else:
        raise ValueError(f"backend must be 'sympy' or 'numpy' (got '{backend}')")

    return _cast_to_type(roots, outtype, backend)
