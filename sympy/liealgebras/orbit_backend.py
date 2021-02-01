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

def _orbit_backend(algebra, weight, stabilizer, dtype=None, backend="numpy"):
    """Runs the orbit stabilizer theorem using an effecient numpy api empowered
    algorithm. The docstring on `cartan_base.StandardCartan.orbit` explains the args.
    """
    name = str(algebra.series) + str(algebra.rank)

    if backend == "numpy" and dtype is None:
        raise ValueError("orbit needs a dtype if numpy backend is used")

    reflection_matrices = algebra._reflection_matrices()
    # used to rotate basis
    omega_inv = algebra.omega_matrix().pinv()
    stab_refl_mats = reflection_matrices
    # used to rotate to alpha basis from orthogonal
    rot = algebra.cartan_matrix().pinv()
    rot = (omega_inv * rot)


    # if stabilizers are passed, select relfection matricies corresponding to those
    # simple root indexes
    if stabilizer:
        stab_refl_mats = [reflection_matrices[i] for i in stabilizer]

    if backend == "numpy":
        weight = np.array(weight, dtype=dtype, order='c').squeeze()
        omega_inv = np.array(omega_inv, dtype=dtype, order="c")
        reflection_matrices = np.array(reflection_matrices, dtype=dtype, order='c')
        stab_refl_mats = np.array(stab_refl_mats, dtype=dtype, order='c')
        rot = np.array(rot, dtype=dtype, order='c')

    # rotate to dominant weyl chamber
    dom = _rotate_to_dominant(name, reflection_matrices, omega_inv, weight, backend=backend)
    if backend == "numpy":
        dom = np.array([dom]).squeeze()

    # generate full orbit based on dominant weight
    fullorbit = _full_orbit(name,stab_refl_mats, dom, algebra.roots() // 2, backend=backend)

    weight_level = lambda x: sum(x * rot)
    convention = sorted(fullorbit, key=lambda x: (tuple(x),weight_level(x)))
    return convention[-1:] + convention[:-1]


_cache = {}
def _reflect(refl_mats, weights, backend="numpy", key=None):
    """Returns all unique reflections of weights by
    map applying reflection matrices over the weights
    """
    if backend == "numpy":
        # Nested for loop that is flatten to (N, dim)
        dim = weights.shape[-1]
        refs = np.array([refl_mats.dot(w.squeeze()) for w in weights]).reshape(-1,dim)
        comb = np.concatenate((refs, weights),axis=0)
        return np.unique(comb, axis=0)
    else:
        results = set(weights)
        global _cache
        if key not in _cache:
            _cache[key]={}
        for w in weights:
            if _cache.get(w) is not None:
                reflections = _cache.get(w)
            else:
                reflections = [w * r for r in refl_mats]
                _cache[key][w] = reflections

            results = results.union(set(reflections))
        return list(results)

def _find_dom(orbits, omega_inv,backend="numpy"):
    """Returns the dominant weight in list of weights, if exists, else
    returns an empty array"""
    if backend=="numpy":
        for x in orbits:
            if np.all(np.dot(x, omega_inv) >= 0):
                return x.squeeze()
        return np.array([], dtype=orbits.dtype)
    else:
        for x in orbits:
            r = x * omega_inv
            if all ([i >= 0 for i in r.row(0)]):
                return x
        return []

def _rotate_to_dominant(name, refl_mats, omega_inv, weight, typeA=False, backend="numpy"):
    """Returns the dominant weight by rotating
    repeatedly across weyl chambers until all
    elements are positive.
    """

    if typeA and backend=="numpy":
        orbits = np.array(list(multiset_permutations(weight)))
        dom = _find_dom(orbits, omega_inv, backend)
        if len(dom) == 0:
            # Check for lack of indexes here
            raise IndexError("Cannot find dominant weight in given orbit")
        return dom
    if backend == "numpy":
        orbits = np.expand_dims(weight,axis=0)
    else:
        orbits = [weight]
    while True:
        orbits = _reflect(refl_mats, orbits, backend, key=name)
        dom = _find_dom(orbits, omega_inv, backend)
        if len(dom) > 0:
            return dom



def _full_orbit(name, refl_mats, weight, num_pos, backend="numpy"):
    """Returns the full orbit of a weight by
    operating on the weight by the reflection matrices
    a number of times equal to the number of positive
    roots in the algebra."""
    if backend == "numpy":
        orbit = np.expand_dims(weight,axis=0)
    else:
        orbit = [weight]
    for _ in range(num_pos):
        orbit = _reflect(refl_mats, orbit, backend=backend, key=name)
    return orbit


def orbit(algebra, weight, stabilizer=None, dtype=int, backend="sympy", outtype="sympy"):
    """Returns the orbit stabilizer theorem algorithm with
    the option of using sympy or numpy. For larger groups
    such as TypeA-TypeD rank 10 or any of TypeE groups,
    numpy is preferred. For published documentation, see
    `sympy.liealgebra.cartan_base.StandardCartan.orbit`"""

    if backend == "sympy":
        orbit_ = _orbit_backend(algebra, weight.as_immutable(), stabilizer, backend="sympy")
    elif backend == "numpy":
        if np is None:
            warnings.warn("Numpy is not installed, falling back to sympy")
            return orbit(algebra, weight, stabilizer, dtype=int, backend="sympy", outtype="sympy")

        orbit_ = _orbit_backend(algebra, weight, stabilizer, dtype)
    else:
        raise ValueError(f"backend must be 'sympy' or 'numpy' (got '{backend}')")

    return _cast_to_type(orbit_, outtype, backend)

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
