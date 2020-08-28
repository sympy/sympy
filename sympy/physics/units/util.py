"""
Several methods to simplify expressions involving unit objects.
"""

from __future__ import division

from sympy import Add, Mul, Pow, Tuple, sympify
from sympy.core.compatibility import reduce, Iterable, ordered
from sympy.physics.units.dimensions import Dimension
from sympy.physics.units.prefixes import Prefix
from sympy.physics.units.quantities import Quantity
from sympy.utilities.iterables import sift


def _get_conversion_matrix_for_expr(expr, target_units, unit_system):
    from sympy import Matrix

    dimension_system = unit_system.get_dimension_system()

    expr_dim = Dimension(unit_system.get_dimensional_expr(expr))
    dim_dependencies = dimension_system.get_dimensional_dependencies(expr_dim, mark_dimensionless=True)
    target_dims = [Dimension(unit_system.get_dimensional_expr(x)) for x in target_units]
    canon_dim_units = [i for x in target_dims for i in dimension_system.get_dimensional_dependencies(x, mark_dimensionless=True)]
    canon_expr_units = {i for i in dim_dependencies}

    if not canon_expr_units.issubset(set(canon_dim_units)):
        return None

    seen = set([])
    canon_dim_units = [i for i in canon_dim_units if not (i in seen or seen.add(i))]

    camat = Matrix([[dimension_system.get_dimensional_dependencies(i, mark_dimensionless=True).get(j, 0) for i in target_dims] for j in canon_dim_units])
    exprmat = Matrix([dim_dependencies.get(k, 0) for k in canon_dim_units])

    res_exponents = camat.solve_least_squares(exprmat, method=None)
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
    7.62963085040767e-20*gravitational_constant**(-0.5)*hbar**0.5*speed_of_light**0.5

    """
    from sympy.physics.units import UnitSystem
    unit_system = UnitSystem.get_unit_system(unit_system)

    if not isinstance(target_units, (Iterable, Tuple)):
        target_units = [target_units]

    if isinstance(expr, Add):
        return Add.fromiter(convert_to(i, target_units, unit_system) for i in expr.args)

    expr = sympify(expr)

    if not isinstance(expr, Quantity) and expr.has(Quantity):
        expr = expr.replace(lambda x: isinstance(x, Quantity), lambda x: x.convert_to(target_units, unit_system))

    def get_total_scale_factor(expr):
        if isinstance(expr, Mul):
            return reduce(lambda x, y: x * y, [get_total_scale_factor(i) for i in expr.args])
        elif isinstance(expr, Pow):
            return get_total_scale_factor(expr.base) ** expr.exp
        elif isinstance(expr, Quantity):
            return unit_system.get_quantity_scale_factor(expr)
        return expr

    depmat = _get_conversion_matrix_for_expr(expr, target_units, unit_system)
    if depmat is None:
        return expr

    expr_scale_factor = get_total_scale_factor(expr)
    return expr_scale_factor * Mul.fromiter((1/get_total_scale_factor(u) * u) ** p for u, p in zip(target_units, depmat))


def quantity_simplify(expr):
    """Return an equivalent expression in which prefixes are replaced
    with numerical values and all units of a given dimension are the
    unified in a canonical manner.

    Examples
    ========

    >>> from sympy.physics.units.util import quantity_simplify
    >>> from sympy.physics.units.prefixes import kilo
    >>> from sympy.physics.units import foot, inch
    >>> quantity_simplify(kilo*foot*inch)
    250*foot**2/3
    >>> quantity_simplify(foot - 6*inch)
    foot/2
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

    return expr


def check_dimensions(expr, unit_system="SI"):
    """Return expr if there are not unitless values added to
    dimensional quantities, else raise a ValueError."""
    # the case of adding a number to a dimensional quantity
    # is ignored for the sake of SymPy core routines, so this
    # function will raise an error now if such an addend is
    # found.
    # Also, when doing substitutions, multiplicative constants
    # might be introduced, so remove those now

    from sympy.physics.units import UnitSystem
    unit_system = UnitSystem.get_unit_system(unit_system)

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
            for i in Mul.make_args(ai):
                if i.has(Quantity):
                    i = Dimension(unit_system.get_dimensional_expr(i))
                if i.has(Dimension):
                    dims.extend(DIM_OF(i).items())
                elif i.free_symbols:
                    skip = True
                    break
            if not skip:
                deset.add(tuple(sorted(dims)))
                if len(deset) > 1:
                    raise ValueError(
                        "addends have incompatible dimensions")

    # clear multiplicative constants on Dimensions which may be
    # left after substitution
    reps = {}
    for m in expr.atoms(Mul):
        if any(isinstance(i, Dimension) for i in m.args):
            reps[m] = m.func(*[
                i for i in m.args if not i.is_number])

    return expr.xreplace(reps)
