"""
Several methods to simplify expressions involving unit objects.
"""
from functools import reduce
from collections.abc import Iterable
from typing import Any, Dict, Optional, Union

from sympy import default_sort_key
from sympy.core.add import Add
from sympy.core.containers import Tuple
from sympy.core.expr import Expr
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.sorting import ordered
from sympy.core.sympify import sympify
from sympy.core.function import Function
from sympy.matrices.exceptions import NonInvertibleMatrixError
from sympy.physics.units.dimensions import Dimension, DimensionSystem
from sympy.physics.units.prefixes import Prefix
from sympy.physics.units.quantities import Quantity
from sympy.physics.units.unitsystem import UnitSystem
from sympy.utilities.iterables import sift
from sympy.physics.units.definitions import dimension_definitions


def _get_conversion_matrix_for_expr(expr, target_units, unit_system):
    from sympy.matrices.dense import Matrix

    dimension_system = unit_system.get_dimension_system()

    expr_dim = Dimension(unit_system.get_dimensional_expr(expr))
    dim_dependencies = dimension_system.get_dimensional_dependencies(expr_dim, mark_dimensionless=True)
    target_dims = [Dimension(unit_system.get_dimensional_expr(x)) for x in target_units]
    canon_dim_units = [i for x in target_dims for i in dimension_system.get_dimensional_dependencies(x, mark_dimensionless=True)]
    canon_expr_units = set(dim_dependencies)

    if not canon_expr_units.issubset(set(canon_dim_units)):
        return None

    seen = set()
    canon_dim_units = [i for i in canon_dim_units if not (i in seen or seen.add(i))]

    camat = Matrix([[dimension_system.get_dimensional_dependencies(i, mark_dimensionless=True).get(j, 0) for i in target_dims] for j in canon_dim_units])
    exprmat = Matrix([dim_dependencies.get(k, 0) for k in canon_dim_units])

    try:
        res_exponents = camat.solve(exprmat)
    except NonInvertibleMatrixError:
        return None

    return res_exponents


def convert_to(expr, target_units, unit_system="SI"):
    """
    Convert ``expr`` to the same expression with all of its units and quantities
    represented as factors of ``target_units``, whenever the dimension is compatible.

    ``target_units`` may be a single unit/quantity, or a collection of
    units/quantities.

    Examples
    ========

    >>> from sympy.physics.units import speed_of_light, meter, gram, second, day
    >>> from sympy.physics.units import mile, newton, kilogram, atomic_mass_constant
    >>> from sympy.physics.units import kilometer, centimeter
    >>> from sympy.physics.units import gravitational_constant, hbar
    >>> from sympy.physics.units import convert_to
    >>> convert_to(mile, kilometer)
    25146*kilometer/15625
    >>> convert_to(mile, kilometer).n()
    1.609344*kilometer
    >>> convert_to(speed_of_light, meter/second)
    299792458*meter/second
    >>> convert_to(day, second)
    86400*second
    >>> 3*newton
    3*newton
    >>> convert_to(3*newton, kilogram*meter/second**2)
    3*kilogram*meter/second**2
    >>> convert_to(atomic_mass_constant, gram)
    1.660539060e-24*gram

    Conversion to multiple units:

    >>> convert_to(speed_of_light, [meter, second])
    299792458*meter/second
    >>> convert_to(3*newton, [centimeter, gram, second])
    300000*centimeter*gram/second**2

    Conversion to Planck units:

    >>> convert_to(atomic_mass_constant, [gravitational_constant, speed_of_light, hbar]).n()
    7.62963087839509e-20*hbar**0.5*speed_of_light**0.5/gravitational_constant**0.5

    """
    from sympy.physics.units import UnitSystem
    unit_system = UnitSystem.get_unit_system(unit_system)

    if not isinstance(target_units, (Iterable, Tuple)):
        target_units = [target_units]

    def handle_Adds(expr):
        return Add.fromiter(convert_to(i, target_units, unit_system)
            for i in expr.args)

    if isinstance(expr, Add):
        return handle_Adds(expr)
    elif isinstance(expr, Pow) and isinstance(expr.base, Add):
        return handle_Adds(expr.base) ** expr.exp

    expr = sympify(expr)
    target_units = sympify(target_units)

    if isinstance(expr, Function):
        expr = expr.together()

    if not isinstance(expr, Quantity) and expr.has(Quantity):
        expr = expr.replace(lambda x: isinstance(x, Quantity),
            lambda x: x.convert_to(target_units, unit_system))

    def get_total_scale_factor(expr):
        if isinstance(expr, Mul):
            return reduce(lambda x, y: x * y,
                [get_total_scale_factor(i) for i in expr.args])
        elif isinstance(expr, Pow):
            return get_total_scale_factor(expr.base) ** expr.exp
        elif isinstance(expr, Quantity):
            return unit_system.get_quantity_scale_factor(expr)
        return expr

    depmat = _get_conversion_matrix_for_expr(expr, target_units, unit_system)
    if depmat is None:
        return expr

    expr_scale_factor = get_total_scale_factor(expr)
    return expr_scale_factor * Mul.fromiter(
        (1/get_total_scale_factor(u)*u)**p for u, p in
        zip(target_units, depmat))


def quantity_simplify(expr, across_dimensions: bool=False, unit_system=None):
    """Return an equivalent expression in which prefixes are replaced
    with numerical values and all units of a given dimension are the
    unified in a canonical manner by default. `across_dimensions` allows
    for units of different dimensions to be simplified together.

    `unit_system` must be specified if `across_dimensions` is True.

    Examples
    ========

    >>> from sympy.physics.units.util import quantity_simplify
    >>> from sympy.physics.units.prefixes import kilo
    >>> from sympy.physics.units import foot, inch, joule, coulomb
    >>> quantity_simplify(kilo*foot*inch)
    250*foot**2/3
    >>> quantity_simplify(foot - 6*inch)
    foot/2
    >>> quantity_simplify(5*joule/coulomb, across_dimensions=True, unit_system="SI")
    5*volt
    """

    if expr.is_Atom or not expr.has(Prefix, Quantity):
        return expr

    # replace all prefixes with numerical values
    p = expr.atoms(Prefix)
    expr = expr.xreplace({p: p.scale_factor for p in p})

    # replace all quantities of given dimension with a canonical
    # quantity, chosen from those in the expression
    d = sift(expr.atoms(Quantity), lambda i: i.dimension)
    for k in d:
        if len(d[k]) == 1:
            continue
        v = list(ordered(d[k]))
        ref = v[0]/v[0].scale_factor
        expr = expr.xreplace({vi: ref*vi.scale_factor for vi in v[1:]})

    if across_dimensions:
        # combine quantities of different dimensions into a single
        # quantity that is equivalent to the original expression

        if unit_system is None:
            raise ValueError("unit_system must be specified if across_dimensions is True")

        unit_system = UnitSystem.get_unit_system(unit_system)
        dimension_system: DimensionSystem = unit_system.get_dimension_system()
        dim_expr = unit_system.get_dimensional_expr(expr)
        dim_deps = dimension_system.get_dimensional_dependencies(dim_expr, mark_dimensionless=True)

        target_dimension: Optional[Dimension] = None
        for ds_dim, ds_dim_deps in dimension_system.dimensional_dependencies.items():
            if ds_dim_deps == dim_deps:
                target_dimension = ds_dim
                break

        if target_dimension is not None:
            target_unit = unit_system.derived_units.get(target_dimension)
        elif target_dimension is None:
            # There is no single regular dimension in the unit system that matches
            # the dimensions of the expression, so we need to build one out of 2
            # or more existing dimensions.
            # We'll do this by brute force combining the units in the unit system, flipping them around,
            # and finding the first that matches the expression's dimensionality.

            # TODO: To be true to the spirit of "simplification", we need to find all possible combinations
            # and then find the one that is "shortest".
            target_dimension = __find_compound_dimension(dimension_system, dim_deps)
            # Dimensions are special. When you do math with them, you get a Dimension object, not a normal sympy expression. The sympy expression is actually in the Dimension's symbol: `name`.

            # replace each dimension in target_dimension with derived unit from the unit system
            dim_map = { symbol: unit_system.derived_units.get(vars(dimension_definitions)[symbol.name]) for symbol in target_dimension.free_symbols }
            print(dim_map)
            target_unit = target_dimension.name.subs(dim_map)
            assert not isinstance(target_unit, Dimension)

            # find Adds that we can safely add together because they have the same dimensions
            adds = expr.atoms(Add)
            for add in adds:
                add_dim_deps = dimension_system.get_dimensional_dependencies(unit_system.get_dimensional_expr(add))
                if all(is_dimensionally_equivalent(add_dim_deps, arg, unit_system) for arg in add.args):
                    units = [unit for unit in unit_system.get_units_non_prefixed() if is_dimensionally_equivalent(unit, add_dim_deps, unit_system)]
                    if units:
                        new_add = add.xreplace({u: 1 for u in add.atoms(Quantity)}) * units[0]
                        expr = expr.xreplace({add: new_add})
                        print("REPLACED", add, "with", new_add)

        if target_unit:
            expr = convert_to(expr, target_unit, unit_system)

    return expr


def __find_compound_dimension(dimension_system: DimensionSystem, dimensional_dependencies) -> Dimension:
    # This can and should be switched out for something more intelligent.

    all_dimensions = dimension_system.base_dims + dimension_system.derived_dims
    # test all reciprocals of all dimensions
    for dim in all_dimensions:
        new_dim = 1 / dim
        new_dim_deps = dimension_system.get_dimensional_dependencies(new_dim)
        if new_dim_deps == dimensional_dependencies:
            print("FOUND COMPATIBLE: ", new_dim)
            return new_dim

    # test all 2 dim combinations of all dimensions
    for dA in all_dimensions:
        for dB in all_dimensions:
            if dA == dB:
                continue
            new_dim = dA * dB
            new_dim_deps = dimension_system.get_dimensional_dependencies(new_dim)
            if new_dim_deps == dimensional_dependencies:
                print("FOUND COMPATIBLE: ", new_dim)
                return new_dim

            new_dim = dA / dB
            new_dim_deps = dimension_system.get_dimensional_dependencies(new_dim)
            if new_dim_deps == dimensional_dependencies:
                print("FOUND COMPATIBLE: ", new_dim)
                return new_dim

            new_dim = dB / dA
            new_dim_deps = dimension_system.get_dimensional_dependencies(new_dim)
            if new_dim_deps == dimensional_dependencies:
                print("FOUND COMPATIBLE: ", new_dim)
                return new_dim


def check_dimensions(expr, unit_system="SI"):
    """Return expr if units in addends have the same
    base dimensions, else raise a ValueError."""
    # the case of adding a number to a dimensional quantity
    # is ignored for the sake of SymPy core routines, so this
    # function will raise an error now if such an addend is
    # found.
    # Also, when doing substitutions, multiplicative constants
    # might be introduced, so remove those now

    from sympy.physics.units import UnitSystem
    unit_system = UnitSystem.get_unit_system(unit_system)

    def addDict(dict1, dict2):
        """Merge dictionaries by adding values of common keys and
        removing keys with value of 0."""
        dict3 = {**dict1, **dict2}
        for key, value in dict3.items():
            if key in dict1 and key in dict2:
                   dict3[key] = value + dict1[key]
        return {key:val for key, val in dict3.items() if val != 0}

    adds = expr.atoms(Add)
    DIM_OF = unit_system.get_dimension_system().get_dimensional_dependencies
    for a in adds:
        deset = set()
        for ai in a.args:
            if ai.is_number:
                deset.add(())
                continue
            dims = []
            skip = False
            dimdict = {}
            for i in Mul.make_args(ai):
                if i.has(Quantity):
                    i = Dimension(unit_system.get_dimensional_expr(i))
                if i.has(Dimension):
                    dimdict = addDict(dimdict, DIM_OF(i))
                elif i.free_symbols:
                    skip = True
                    break
            dims.extend(dimdict.items())
            if not skip:
                deset.add(tuple(sorted(dims, key=default_sort_key)))
                if len(deset) > 1:
                    raise ValueError(
                        "addends have incompatible dimensions: {}".format(deset))

    # clear multiplicative constants on Dimensions which may be
    # left after substitution
    reps = {}
    for m in expr.atoms(Mul):
        if any(isinstance(i, Dimension) for i in m.args):
            reps[m] = m.func(*[
                i for i in m.args if not i.is_number])

    return expr.xreplace(reps)


def is_dimensionally_equivalent(A: Union[Expr, Dimension, Dict[Expr, int]], B: Union[Expr, Dimension, Dict[Expr, int]], unit_system):
    """Return True if A and B are dimensionally equivalent."""

    unit_system = UnitSystem.get_unit_system(unit_system)
    dimension_system: DimensionSystem = unit_system.get_dimension_system()

    def _get_dim_deps(expr: Union[Expr, Dimension, Dict[Expr, int]]) -> Dict[Expr, int]:
        if isinstance(expr, dict):
            return expr
        if isinstance(expr, Dimension):
            return dimension_system.get_dimensional_dependencies(expr)
        dim_expr = unit_system.get_dimensional_expr(expr)
        return dimension_system.get_dimensional_dependencies(dim_expr)

    dim_depsA = _get_dim_deps(A)
    dim_depsB = _get_dim_deps(B)

    return dim_depsA == dim_depsB
