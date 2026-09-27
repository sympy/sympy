r"""
Conversion of expressions and equations between unit systems.

Explanation
===========

The symbols of a formula represent physical quantities, but the relations
among the quantities depend on the unit system. The function
:func:`convert_unit_system` rewrites a formula valid in a unit system into the
corresponding formula of another unit system. This is different from
:func:`~sympy.physics.units.util.convert_to`, which expresses a quantity with
other units of the same unit system.

The conversion is based on the dimension systems of the two unit systems. The
unit system with more independent dimensions is called *full*, the other one
*reduced*. For example, the SI is full with respect to Gaussian units, as in
Gaussian units the current is not a base dimension.

A quantity is represented by `X_f` in the full unit system and by `X_r` in
the reduced one. They are related by

.. math::
    X_f = K X_r, \qquad K = \prod_i \left( \frac{C_i}{c_i} \right)^{p_i}
    \prod_j G_j^{q_j}

where

- `C_i` are physical constants of the full unit system that are pure numbers
  in the reduced one, `c_i` being their values,
- `G_j` are physical constants that have the same dimension in both unit
  systems,
- the exponents are the solution of a linear system, which states that the
  dimension of `K` in the full unit system is the ratio of the dimensions of
  `X_f` and `X_r`.

Passing from the full to the reduced unit system, the constants `C_i` are
replaced by their values. In the opposite direction the symbols are divided
by `K`.
"""
from __future__ import annotations

from collections import defaultdict

from sympy.core.add import Add
from sympy.core.basic import Basic
from sympy.core.function import AppliedUndef, Derivative, UndefinedFunction
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.relational import Relational
from sympy.core.singleton import S
from sympy.core.sorting import default_sort_key
from sympy.core.symbol import Dummy, Symbol
from sympy.core.sympify import sympify
from sympy.functions.elementary.complexes import Abs, conjugate, im, re
from sympy.functions.elementary.miscellaneous import Max, Min
from sympy.functions.elementary.piecewise import Piecewise
from sympy.integrals.integrals import Integral
from sympy.matrices.dense import Matrix
from sympy.matrices.matrixbase import MatrixBase
from sympy.physics.units.definitions.unit_definitions import speed_of_light
from sympy.physics.units.dimensions import Dimension
from sympy.physics.units.quantities import PhysicalConstant, Quantity
from sympy.physics.units.unitsystem import UnitSystem


def _map_dependencies(dependencies, dimension_system):
    """
    Express the dimensional dependencies in terms of the base dimensions of
    ``dimension_system``.

    Parameters
    ==========

    dependencies : dict
        Dimensions and their exponents. The dimensions need not be base
        dimensions of ``dimension_system``.
    dimension_system : DimensionSystem

    Returns
    =======

    dict
        Base dimensions of ``dimension_system`` and their exponents.
        Dimensions unknown to ``dimension_system`` are left unchanged.

    Examples
    ========

    >>> from sympy.physics.units import current, length
    >>> from sympy.physics.units.systems.cgs import dimsys_cgs
    >>> from sympy.physics.units.unit_system_conversion import _map_dependencies
    >>> _map_dependencies({current: 2, length: -3}, dimsys_cgs)
    {Dimension(mass): 1, Dimension(time): -4}

    """
    result = defaultdict(int)
    for dim, exponent in dependencies.items():
        for base_dim, base_exponent in dimension_system.get_dimensional_dependencies(dim).items():
            result[base_dim] += exponent*base_exponent
    return {dim: exponent for dim, exponent in result.items() if exponent != 0}


def _is_contained(inner, outer):
    """
    Check that the base dimensions of the dimension system ``inner`` are
    independent dimensions in the dimension system ``outer``.

    Explanation
    ===========

    Every base dimension of ``inner`` is expressed with the base dimensions of
    ``outer``, the result is expressed again with the base dimensions of
    ``inner``. If ``outer`` does not merge dimensions that are distinguished
    by ``inner``, the initial dimension is recovered.

    Examples
    ========

    >>> from sympy.physics.units.systems.si import dimsys_SI
    >>> from sympy.physics.units.systems.cgs import dimsys_cgs
    >>> from sympy.physics.units.unit_system_conversion import _is_contained
    >>> _is_contained(dimsys_cgs, dimsys_SI)
    True
    >>> _is_contained(dimsys_SI, dimsys_cgs)
    False

    """
    for dim in inner.base_dims:
        dependencies = outer.get_dimensional_dependencies(dim)
        if _map_dependencies(dependencies, inner) != {dim: 1}:
            return False
    return True


def _solve_exponents(vectors, target):
    """
    Find the exponents ``p`` such that the sum of ``p[i]*vectors[i]`` is equal
    to ``target``.

    Parameters
    ==========

    vectors : list of dict
        Dimensional dependencies of some quantities, they have to be linearly
        independent.
    target : dict
        Dimensional dependencies to be obtained as a product of powers of the
        quantities.

    Returns
    =======

    list or None
        The exponents of the quantities, or ``None`` if there is no solution.

    Examples
    ========

    >>> from sympy.physics.units import length, mass, time
    >>> from sympy.physics.units.unit_system_conversion import _solve_exponents
    >>> velocity = {length: 1, time: -1}
    >>> action = {mass: 1, length: 2, time: -1}
    >>> _solve_exponents([velocity, action], {mass: 1, length: 3, time: -2})
    [1, 1]
    >>> _solve_exponents([velocity, action], {mass: 1}) is None
    True

    """
    if not target:
        return [S.Zero]*len(vectors)
    if not vectors:
        return None
    dims = sorted(set(target).union(*vectors), key=default_sort_key)
    matrix = Matrix([[vector.get(dim, 0) for vector in vectors] for dim in dims])
    rhs = Matrix([target.get(dim, 0) for dim in dims])
    if matrix.row_join(rhs).rank() != matrix.rank():
        return None
    solution, _ = matrix.gauss_jordan_solve(rhs)
    return list(solution)


def _independent(vectors):
    """
    Check that the dimensional dependencies contained in the list ``vectors``
    are linearly independent.
    """
    dims = sorted(set().union(*vectors), key=default_sort_key)
    matrix = Matrix([[vector.get(dim, 0) for vector in vectors] for dim in dims])
    return matrix.rank() == len(vectors)


def _is_defined(dimension, quantity):
    """
    Check that ``dimension`` is the dimension given to ``quantity`` by a unit
    system. If the unit system does not know the quantity, it returns a
    dimension with the name of the quantity.
    """
    return dimension.name != quantity.name


def _simplest(factors):
    """
    Return the factor with the lowest powers of the constants.

    Explanation
    ===========

    The addends of a sum may have different conversion factors. One of them
    is collected, the remaining part of the other ones is left inside the
    sum. The simplest factor is chosen in order to get the form which is
    commonly used, e.g. `E + v B/c` instead of `(c E + v B)/c`.
    """
    def weight(factor):
        return sum(abs(exponent) for base, exponent in factor.as_powers_dict().items()
                   if isinstance(base, Quantity))
    return min(factors, key=weight)


def _distribute(expr, factor):
    """
    Multiply ``expr`` by ``factor``. If ``expr`` is a sum, every addend is
    multiplied, so that the constants contained in ``factor`` can cancel.
    """
    if factor == 1:
        return expr
    if isinstance(expr, Add):
        return Add(*[arg*factor for arg in expr.args])
    return expr*factor


class _UnitSystemConverter:
    """
    Conversion of expressions between two unit systems, see
    :func:`convert_unit_system`.

    Explanation
    ===========

    The unit system with the larger number of independent dimensions is
    called *full*, the other one *reduced*.

    A quantity of given dimension is represented in the two unit systems by
    `X_f` and `X_r`, with `X_f = K X_r`. The factor `K` is a product of powers
    of physical constants of the full unit system, its exponents are found by
    comparing the dimensions of the quantity in the two unit systems. The
    constants contained in `K` are either dimensionless in the reduced unit
    system, or they have the same dimension in both unit systems.

    The expression is visited recursively by :meth:`_split`, which returns the
    conversion factor of every subexpression separated from the converted
    subexpression. In this way the factors of the terms of products and sums
    are combined before they are multiplied by the expression.

    Parameters
    ==========

    source : UnitSystem
        The unit system of the expressions to be converted.
    target : UnitSystem
        The unit system of the results.
    dimensions : dict
        The dimensions of symbols and undefined functions.
    constants : list or None
        The physical constants contained in the conversion factors. If
        ``None``, they are determined by :meth:`_set_constants`.

    Attributes
    ==========

    full, reduced : UnitSystem
        The two unit systems, sorted by number of independent dimensions.
    contracting : bool
        Whether the conversion is from the full to the reduced unit system.
    unit_constants : list
        The constants that are pure numbers in the reduced unit system.
    common_constants : list
        The constants with the same dimension in both unit systems.
    """

    def __init__(self, source, target, dimensions, constants):
        self.source = source
        self.target = target
        self.dimsys_source = source.get_dimension_system()
        self.dimsys_target = target.get_dimension_system()

        if _is_contained(self.dimsys_target, self.dimsys_source):
            self.contracting = True
            self.full, self.reduced = source, target
        elif _is_contained(self.dimsys_source, self.dimsys_target):
            self.contracting = False
            self.full, self.reduced = target, source
        else:
            raise NotImplementedError(
                "cannot convert from %s to %s: neither of the unit systems "
                "contains the base dimensions of the other one" % (source, target))
        self.dimsys_full = self.full.get_dimension_system()
        self.dimsys_reduced = self.reduced.get_dimension_system()

        self.dimensions = {}
        self.function_dimensions = {}
        for key, dim in dimensions.items():
            dim = sympify(dim)
            if dim == 1:
                dim = Dimension(1)
            if not isinstance(dim, Dimension):
                raise TypeError("expected a dimension for %s, got %s" % (key, dim))
            if isinstance(key, UndefinedFunction):
                self.function_dimensions[key] = dim
            else:
                self.dimensions[sympify(key)] = dim

        self._set_constants(constants)
        self._factors = {}

    def _quantity_dependencies(self, quantity):
        """
        Dimensional dependencies of the quantity in the full and in the
        reduced unit system.

        If one of the unit systems does not define the dimension of the
        quantity, the dimension defined by the other one is used.
        """
        dim_full = self.full.get_quantity_dimension(quantity)
        dim_reduced = self.reduced.get_quantity_dimension(quantity)
        if not _is_defined(dim_reduced, quantity):
            dim_reduced = dim_full
        elif not _is_defined(dim_full, quantity):
            dim_full = dim_reduced
        return (self.dimsys_full.get_dimensional_dependencies(dim_full),
                self.dimsys_reduced.get_dimensional_dependencies(dim_reduced))

    def _mismatch(self, deps_full, deps_reduced):
        """
        Dimensional dependencies in the full unit system of the factor
        relating the representations of a quantity in the two unit systems.

        Explanation
        ===========

        The dimensional dependencies in the reduced unit system are expressed
        with the base dimensions of the full unit system and subtracted from
        the ones in the full unit system. The result is empty if the quantity
        has the same meaning in both unit systems.

        Dimensions that are not known to one of the unit systems are
        expressed with the definition of the other one. An error is raised if
        the reduced unit system does not define a dimension that depends on
        the base dimensions it does not have.

        Parameters
        ==========

        deps_full : dict
            Dimensional dependencies of the quantity in the full unit system.
        deps_reduced : dict
            Dimensional dependencies of the quantity in the reduced unit
            system.
        """
        for dim in deps_reduced:
            if dim in self.dimsys_reduced.base_dims:
                continue
            for base_dim in self.dimsys_full.get_dimensional_dependencies(dim):
                if self.dimsys_reduced.get_dimensional_dependencies(base_dim) != {base_dim: 1}:
                    raise ValueError(
                        "dimension %s is not defined in the unit system %s" % (
                            dim.name, self.reduced))
        mismatch = defaultdict(int)
        for dim, exponent in deps_full.items():
            if dim in self.dimsys_full.base_dims:
                mismatch[dim] += exponent
                continue
            dependencies = self.dimsys_reduced.get_dimensional_dependencies(dim)
            for base_dim, base_exponent in _map_dependencies(dependencies, self.dimsys_full).items():
                mismatch[base_dim] += exponent*base_exponent
        for dim, exponent in _map_dependencies(deps_reduced, self.dimsys_full).items():
            mismatch[dim] -= exponent
        return {dim: exponent for dim, exponent in mismatch.items() if exponent != 0}

    def _known_constants(self):
        """
        Physical constants whose dimension is defined by the reduced unit
        system.
        """
        quantities = set(self.reduced._quantity_dimension_map)
        quantities.update(self.dimsys_reduced._quantity_dimension_map)
        return [q for q in quantities if isinstance(q, PhysicalConstant)]

    def _has_scale_factor(self, quantity):
        """
        Check that the scale factor of the quantity is defined by the reduced
        unit system.
        """
        return (quantity in self.reduced._quantity_scale_factors or
                quantity in self.dimsys_reduced._quantity_scale_factors)

    def _set_constants(self, constants):
        """
        Set the physical constants contained in the conversion factors.

        Explanation
        ===========

        If ``constants`` is ``None``, the constants are the ones that are
        dimensionless in the reduced unit system and not in the full one,
        followed by the speed of light. The constants equal to one are
        preferred. The constants whose dimension is a combination of the
        dimensions of the previous ones are skipped.

        The constants given by the user are required to be independent. They
        have to be either dimensionless in the reduced unit system, or have
        the same dimension in both unit systems.
        """
        if constants is None:
            unit_constants = []
            for constant in self._known_constants():
                deps_full, deps_reduced = self._quantity_dependencies(constant)
                if deps_full and not deps_reduced and self._has_scale_factor(constant):
                    unit_constants.append(constant)
            unit_constants.sort(key=lambda q: (
                self.reduced.get_quantity_scale_factor(q) != 1, default_sort_key(q)))
            constants = unit_constants + [speed_of_light]
            skip_dependent = True
        else:
            constants = list(constants)
            skip_dependent = False

        # Constants that are pure numbers in the reduced unit system:
        self.unit_constants = []
        # Constants with the same dimension in both unit systems:
        self.common_constants = []
        vectors = []
        for constant in constants:
            deps_full, deps_reduced = self._quantity_dependencies(constant)
            if not deps_full:
                if skip_dependent:
                    continue
                raise ValueError("%s is dimensionless in the unit system %s" % (constant, self.full))
            if not _independent(vectors + [deps_full]):
                if skip_dependent:
                    continue
                raise ValueError("the dimensions of the constants are not independent")
            if not deps_reduced:
                self.unit_constants.append(constant)
            elif not self._mismatch(deps_full, deps_reduced):
                self.common_constants.append(constant)
            elif skip_dependent:
                continue
            else:
                raise ValueError(
                    "%s should either be dimensionless in the unit system %s, "
                    "or have the same dimension in both unit systems" % (constant, self.reduced))
            vectors.append(deps_full)

        self.constants = self.unit_constants + self.common_constants
        self._vectors = [self._quantity_dependencies(c)[0] for c in self.constants]
        self._common_vectors_reduced = [self._quantity_dependencies(c)[1] for c in self.common_constants]

    def _factor_from_dependencies(self, deps_full, deps_reduced):
        """
        Factor multiplying a quantity with the given dimensional dependencies
        when the unit system is changed.

        Explanation
        ===========

        The exponents of the constants are determined by the mismatch of the
        dimensions, see :meth:`_mismatch`. Converting to the reduced unit
        system, the constants that are pure numbers are dropped, as the
        symbols are supposed to absorb them. Converting to the full unit
        system, the factor is the inverse of the product of the constants,
        divided by their values in the reduced unit system.
        """
        mismatch = self._mismatch(deps_full, deps_reduced)
        if not mismatch:
            return S.One
        key = tuple(sorted(mismatch.items(), key=default_sort_key))
        if key in self._factors:
            return self._factors[key]
        exponents = _solve_exponents(self._vectors, mismatch)
        if exponents is None:
            raise ValueError(
                "cannot relate the dimension %s of the unit system %s to the "
                "dimension %s of the unit system %s using the constants %s" % (
                    Dimension._from_dimensional_dependencies(deps_full), self.full,
                    Dimension._from_dimensional_dependencies(deps_reduced), self.reduced,
                    self.constants))
        if self.contracting:
            factor = Mul(*[c**p for c, p in zip(self.constants, exponents) if c in self.common_constants])
        else:
            factor = Mul(*[
                (self.reduced.get_quantity_scale_factor(c)/c)**p if c in self.unit_constants else c**(-p)
                for c, p in zip(self.constants, exponents)])
        self._factors[key] = factor
        return factor

    def factor(self, dimension):
        """
        Factor multiplying a symbol of given dimension when the unit system
        is changed.

        Examples
        ========

        >>> from sympy.physics.units import charge, magnetic_density, energy
        >>> from sympy.physics.units.systems.si import SI
        >>> from sympy.physics.units.systems.cgs import cgs_gauss
        >>> from sympy.physics.units.unit_system_conversion import _UnitSystemConverter
        >>> converter = _UnitSystemConverter(SI, cgs_gauss, {}, None)
        >>> converter.factor(magnetic_density)
        1/speed_of_light
        >>> converter.factor(energy)
        1
        >>> converter = _UnitSystemConverter(cgs_gauss, SI, {}, None)
        >>> converter.factor(charge)
        sqrt(coulomb_constant)
        >>> converter.factor(magnetic_density)
        speed_of_light/sqrt(coulomb_constant)

        """
        deps_full = self.dimsys_full.get_dimensional_dependencies(dimension)
        deps_reduced = self.dimsys_reduced.get_dimensional_dependencies(dimension)
        return self._factor_from_dependencies(deps_full, deps_reduced)

    def _declared_dimension(self, expr):
        """
        Dimension given by the user to ``expr``, or ``None``.
        """
        if expr in self.dimensions:
            return self.dimensions[expr]
        if isinstance(expr, AppliedUndef) and expr.func in self.function_dimensions:
            return self.function_dimensions[expr.func]
        return None

    def _constant_value(self, expr):
        r"""
        Value of a physical constant that is a multiple of the common
        constants in the reduced unit system.

        Explanation
        ===========

        Some constants of the full unit system are fixed by the conventions
        of the reduced unit system, for example in Gaussian units the vacuum
        permittivity is `1/(4 \pi)` and the magnetic constant is
        `4 \pi/c^2`. They are replaced by their value. ``None`` is returned
        for the other constants, which are handled like symbols.
        """
        if not isinstance(expr, PhysicalConstant) or not self._has_scale_factor(expr):
            return None
        deps_full, deps_reduced = self._quantity_dependencies(expr)
        if not self._mismatch(deps_full, deps_reduced):
            return None
        exponents = _solve_exponents(self._common_vectors_reduced, deps_reduced)
        if exponents is None:
            return None
        value = self.reduced.get_quantity_scale_factor(expr)
        for constant, exponent in zip(self.common_constants, exponents):
            value *= (constant/self.reduced.get_quantity_scale_factor(constant))**exponent
        return value

    def _split(self, expr):
        """
        Return the conversion factor of the quantity represented by ``expr``
        and the expression of the quantity in the target unit system.

        Explanation
        ===========

        The product of the two returned values is the converted expression.
        The rules are:

        - symbols and functions get the factor of their dimension, the
          arguments of the functions are the quantities of the target unit
          system;
        - the factor of a product is the product of the factors;
        - the simplest factor of the terms of a sum is collected, see
          :func:`_simplest`;
        - both sides of an equation are divided by the factor of the left
          hand side;
        - derivatives are divided by the factors of the variables, integrals
          are multiplied by them;
        - the arguments of the other functions are converted.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.units import charge, length
        >>> from sympy.physics.units.systems.si import SI
        >>> from sympy.physics.units.systems.cgs import cgs_gauss
        >>> from sympy.physics.units.unit_system_conversion import _UnitSystemConverter
        >>> q, r = symbols("q r")
        >>> converter = _UnitSystemConverter(cgs_gauss, SI, {q: charge, r: length}, None)
        >>> converter._split(q**2/r)
        (coulomb_constant, q**2/r)

        """
        dim = self._declared_dimension(expr)
        if dim is not None:
            if isinstance(expr, AppliedUndef):
                expr = expr.func(*[self._split(arg)[1] for arg in expr.args])
            return self.factor(dim), expr
        if isinstance(expr, Quantity):
            value = self._constant_value(expr)
            if value is not None:
                return S.One, value
            return self._factor_from_dependencies(*self._quantity_dependencies(expr)), expr
        if isinstance(expr, AppliedUndef):
            raise ValueError("the dimension of %s has not been specified" % expr)
        if not isinstance(expr, Basic) or not expr.args:
            return S.One, expr
        if isinstance(expr, Mul):
            factors, args = zip(*[self._split(arg) for arg in expr.args])
            return Mul(*factors), Mul(*args)
        if isinstance(expr, Pow) and expr.exp.is_number:
            factor, base = self._split(expr.base)
            return factor**expr.exp, base**expr.exp
        if isinstance(expr, (Add, Min, Max)):
            factors, args = zip(*[self._split(arg) for arg in expr.args])
            reference = _simplest(factors)
            return reference, expr.func(*[
                _distribute(arg, factor/reference) for factor, arg in zip(factors, args)])
        if isinstance(expr, Piecewise):
            factors, args = zip(*[self._split(arg) for arg, _ in expr.args])
            reference = _simplest(factors)
            return reference, Piecewise(*[
                (_distribute(arg, factor/reference), self._convert(cond))
                for factor, arg, (_, cond) in zip(factors, args, expr.args)])
        if isinstance(expr, (Abs, conjugate, re, im)):
            factor, arg = self._split(expr.args[0])
            return factor, expr.func(arg)
        if isinstance(expr, Relational):
            lhs_factor, lhs = self._split(expr.lhs)
            rhs_factor, rhs = self._split(expr.rhs)
            if lhs.is_number:
                return S.One, expr.func(lhs, rhs)
            return S.One, expr.func(lhs, _distribute(rhs, rhs_factor/lhs_factor))
        if isinstance(expr, Derivative):
            factor, function = self._split(expr.expr)
            for variable, count in expr.variable_count:
                factor /= self._split(variable)[0]**count
            return factor, Derivative(function, *expr.variable_count)
        if isinstance(expr, Integral):
            factor, function = self._split(expr.function)
            limits = []
            for variable, *bounds in expr.limits:
                variable_factor = self._split(variable)[0]
                factor *= variable_factor
                limits.append((variable, *[
                    _distribute(self._convert(bound), 1/variable_factor) for bound in bounds]))
            return factor, Integral(function, *limits)
        return S.One, expr.func(*[self._convert(arg) for arg in expr.args])

    def _convert(self, expr):
        """
        Convert ``expr`` to the target unit system.
        """
        factor, expr = self._split(expr)
        return _distribute(expr, factor)

    def convert(self, expr, dimension=None):
        """
        Convert ``expr`` to the target unit system, after checking that the
        dimensions of all its symbols are known.

        Parameters
        ==========

        expr : Expr, Relational
        dimension : Dimension, optional
            The dimension of the quantity represented by ``expr``. If given,
            the result is divided by the conversion factor of this dimension.
            It is ignored if ``expr`` is an equation.
        """
        declared = {key: Dummy() for key in self.dimensions if not isinstance(key, Symbol)}
        missing = expr.xreplace(declared).free_symbols - set(self.dimensions) - set(declared.values())
        if missing:
            raise ValueError(
                "the dimension of the following symbols has not been specified: %s" % (
                    ", ".join(str(i) for i in sorted(missing, key=default_sort_key))))
        factor, result = self._split(expr)
        if dimension is not None and not isinstance(expr, Relational):
            factor /= self.factor(dimension)
        return _distribute(result, factor)


def _apply(func, expr):
    """
    Apply ``func`` to ``expr``, or to its elements if it is a matrix, a list
    or a tuple.
    """
    if isinstance(expr, MatrixBase):
        return expr.applyfunc(func)
    if isinstance(expr, (list, tuple)):
        return type(expr)(_apply(func, i) for i in expr)
    return func(sympify(expr))


def convert_unit_system(expr, dimensions, source, target, dimension=None, constants=None):
    r"""
    Convert an expression or an equation from the form it has in the unit
    system ``source`` to the form it has in the unit system ``target``.

    Explanation
    ===========

    The form of physical laws depends on the unit system. For example, in the
    SI the force between two charges at distance `r` is

    .. math::
        F = \frac{1}{4 \pi \epsilon_0} \frac{q_1 q_2}{r^2}

    while in Gaussian units

    .. math::
        F = \frac{q_1 q_2}{r^2}

    The reason is that the unit systems do not have the same base
    dimensions. The SI has the current among its base dimensions, while in
    Gaussian units the current is a derived dimension, expressed in terms of
    length, mass and time. Moving to a unit system with fewer base dimensions
    turns some physical constants into pure numbers. The other way round,
    dimensionful constants have to be restored.

    Both cases are supported, provided that the dimension of every symbol in
    the expression is known.

    Parameters
    ==========

    expr : Expr, Relational, Matrix, list
        The expression to convert. Physical constants have to be given as
        quantities, e.g. ``speed_of_light``.
    dimensions : dict
        The dimension of every symbol and undefined function of ``expr``. Use
        ``1`` for dimensionless symbols.
    source : UnitSystem, str
        The unit system ``expr`` is written in.
    target : UnitSystem, str
        The unit system to convert to.
    dimension : Dimension, optional
        The dimension of the quantity represented by ``expr``. It is not
        needed for equations. If missing, it is assumed that the
        definition of the quantity is the same in both unit systems.
    constants : list, optional
        The physical constants used to express the conversion factors. The
        default choice is given by the constants that are dimensionless in
        one unit system and not in the other one, and the speed of light.

    Examples
    ========

    >>> from sympy import symbols, pi, Eq, Function
    >>> from sympy.physics.units import convert_unit_system
    >>> from sympy.physics.units import force, charge, length, velocity
    >>> from sympy.physics.units import vacuum_permittivity
    >>> from sympy.physics.units.systems.si import SI
    >>> from sympy.physics.units.systems.cgs import cgs_gauss
    >>> F, q1, q2, r = symbols("F q1 q2 r")
    >>> dims = {F: force, q1: charge, q2: charge, r: length}

    Coulomb's law from SI to Gaussian units:

    >>> eq = Eq(F, q1*q2/(4*pi*vacuum_permittivity*r**2))
    >>> convert_unit_system(eq, dims, SI, cgs_gauss)
    Eq(F, q1*q2/r**2)

    In the opposite direction the missing constant is restored. The Coulomb
    constant is used, as it is the one equal to one in Gaussian units:

    >>> convert_unit_system(Eq(F, q1*q2/r**2), dims, cgs_gauss, SI)
    Eq(F, coulomb_constant*q1*q2/r**2)

    A different choice of constants is possible:

    >>> from sympy.physics.units import speed_of_light
    >>> convert_unit_system(Eq(F, q1*q2/r**2), dims, cgs_gauss, SI,
    ...     constants=[vacuum_permittivity, speed_of_light])
    Eq(F, q1*q2/(4*vacuum_permittivity*pi*r**2))

    The Lorentz force and the Ampere-Maxwell law in one space dimension:

    >>> from sympy.physics.units import voltage, magnetic_density, current, time
    >>> from sympy.physics.units import magnetic_constant
    >>> q, v, x, t = symbols("q v x t")
    >>> E, B, J = symbols("E B J", cls=Function)
    >>> dims = {F: force, q: charge, v: velocity, x: length, t: time,
    ...     E: voltage/length, B: magnetic_density, J: current/length**2}
    >>> convert_unit_system(Eq(F, q*(E(x, t) + v*B(x, t))), dims, SI, cgs_gauss)
    Eq(F, q*(E(x, t) + v*B(x, t)/speed_of_light))
    >>> eq = Eq(B(x, t).diff(x), magnetic_constant*J(x, t)
    ...     + magnetic_constant*vacuum_permittivity*E(x, t).diff(t))
    >>> eq_gauss = convert_unit_system(eq, dims, SI, cgs_gauss)
    >>> eq_gauss
    Eq(Derivative(B(x, t), x), 4*pi*J(x, t)/speed_of_light + Derivative(E(x, t), t)/speed_of_light)
    >>> convert_unit_system(eq_gauss, dims, cgs_gauss, SI)
    Eq(Derivative(B(x, t), x), 4*coulomb_constant*pi*J(x, t)/speed_of_light**2 + Derivative(E(x, t), t)/speed_of_light**2)

    If the expression is not an equation, the dimension of the quantity it
    represents should be specified:

    >>> I = symbols("I")
    >>> dims = {I: current, r: length}
    >>> expr = magnetic_constant*I/(2*pi*r)
    >>> convert_unit_system(expr, dims, SI, cgs_gauss, dimension=magnetic_density)
    2*I/(speed_of_light*r)

    Notes
    =====

    The conversion factors are determined by dimensional analysis only. Some
    quantities differ by an additional numerical factor in some unit systems,
    for example the fields `\mathbf{D}` and `\mathbf{H}` in Gaussian units
    contain a factor `4 \pi` with respect to the ones of the SI. This has to
    be taken care of manually.

    See Also
    ========

    sympy.physics.units.util.convert_to

    """
    source = UnitSystem.get_unit_system(source)
    target = UnitSystem.get_unit_system(target)
    converter = _UnitSystemConverter(source, target, dimensions, constants)
    if dimension is not None:
        dimension = sympify(dimension)
        if dimension == 1:
            dimension = Dimension(1)
    return _apply(lambda i: converter.convert(i, dimension), expr)
