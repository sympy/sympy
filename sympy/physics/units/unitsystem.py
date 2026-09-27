"""
Unit system for physical quantities; include definition of constants.
"""
from __future__ import annotations

from sympy.core.add import Add
from sympy.core.function import (Application, Derivative, Function)
from sympy.core.mul import Mul
from sympy.core.numbers import Float
from sympy.core.power import Pow
from sympy.core.singleton import S
from sympy.core.sympify import sympify
from sympy.matrices.dense import Matrix
from sympy.physics.units.dimensions import _QuantityMapper

from .dimensions import Dimension, DimensionSystem
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sympy.core.expr import Expr
    from sympy.physics.units.quantities import Quantity


class UnitSystem(_QuantityMapper):
    """
    UnitSystem represents a coherent set of units.

    A unit system is basically a dimension system with notions of scales. Many
    of the methods are defined in the same way.

    It is much better if all base units have a symbol.

    Parameters
    ==========

    base_units : list
        The units of the base dimensions.
    units : list, optional
        The other units of the unit system.
    name : str, optional
        The name of the unit system, it is used by :meth:`get_unit_system`.
    descr : str, optional
        A description of the unit system.
    dimension_system : DimensionSystem
        The dimension system defining the relations among the dimensions.
    derived_units : dict, optional
        The units of the derived dimensions.
    defining_constants : dict, optional
        The physical constants that are pure numbers in the unit system, and
        their values.

    See Also
    ========

    contract
    """

    _unit_systems: dict[str, UnitSystem] = {}

    def __init__(self, base_units, units=(), name="", descr="", dimension_system=None, derived_units: dict[Dimension, Expr]={},
                 defining_constants: dict[Quantity, Expr]={}):

        UnitSystem._unit_systems[name] = self

        self.name = name
        self.descr = descr

        self._base_units = base_units
        self._dimension_system = dimension_system
        self._units = tuple(set(base_units) | set(units))
        self._base_units = tuple(base_units)
        self._derived_units = derived_units
        self._defining_constants = dict(defining_constants)
        self._parent: UnitSystem | None = None
        self._contraction: tuple | None = None

        super().__init__()

    def __str__(self):
        """
        Return the name of the system.

        If it does not exist, then it makes a list of symbols (or names) of
        the base dimensions.
        """

        if self.name != "":
            return self.name
        else:
            return "UnitSystem((%s))" % ", ".join(
                str(d) for d in self._base_units)

    def __repr__(self):
        return '<UnitSystem: %s>' % repr(self._base_units)

    def extend(self, base, units=(), name="", description="", dimension_system=None, derived_units: dict[Dimension, Expr]={}):
        """Extend the current system into a new one.

        Take the base and normal units of the current system to merge
        them to the base and normal units given in argument.
        If not provided, name and description are overridden by empty strings.
        """

        base = self._base_units + tuple(base)
        units = self._units + tuple(units)

        return UnitSystem(base, units, name, description, dimension_system, {**self._derived_units, **derived_units},
                          self._defining_constants)

    def contract(self, constants, base_units=(), units=(), name="", description="", derived_units: dict[Dimension, Expr]={}):
        r"""
        Create the unit system in which the given physical constants are
        pure numbers.

        Explanation
        ===========

        Every constant that becomes a pure number removes a base dimension.
        The dimensions and the scale factors of all quantities of the current
        unit system are expressed in terms of the remaining base dimensions.
        This is the way natural units are defined.

        Let `d_1, \ldots, d_n` be the base dimensions of the current unit
        system, and `C_1, \ldots, C_k` the constants. The dimensions of the
        constants, together with `n - k` further dimensions
        `b_1, \ldots, b_{n-k}`, are a basis of the space of the dimensions.
        The dimension of every quantity `Q` has a unique expression

        .. math::
            [Q] = \prod_i [C_i]^{p_i} \prod_j b_j^{q_j}

        In the new unit system the base dimensions are `b_j`, the dimension
        of the quantity is

        .. math::
            [Q] = \prod_j b_j^{q_j}

        and its scale factor is multiplied by

        .. math::
            \prod_i \left( \frac{c_i}{s_i} \right)^{p_i}

        where `s_i` is the scale factor of `C_i` in the current unit system
        and `c_i` is its value in the new one.

        The dimensions `b_j` are the dimensions of ``base_units``, followed by
        the base dimensions of the current unit system that are independent
        of the previous ones.

        The dimension system of the new unit system is created when it is
        used the first time.

        Parameters
        ==========

        constants : dict, list
            The physical constants and their values in the new unit system.
            If a list is given, the constants are set to one.
        base_units : list, optional
            Units whose dimensions are base dimensions of the new unit system.
            If they are not enough, the missing ones are picked from the base
            dimensions of the current unit system.
        units : list, optional
            Further units of the new unit system.
        name : str, optional
            The name of the new unit system.
        description : str, optional
            A description of the new unit system.
        derived_units : dict, optional
            The units of the derived dimensions of the new unit system.

        Returns
        =======

        UnitSystem
            The new unit system. Its defining constants are the ones of the
            current unit system and ``constants``.

        Raises
        ======

        ValueError
            If the dimensions of the constants and of the base units are not
            independent, or they are not defined in the current unit system.

        Examples
        ========

        >>> from sympy.physics.units import speed_of_light, second, meter
        >>> from sympy.physics.units import convert_to, time, length
        >>> from sympy.physics.units.systems import SI
        >>> unit_system = SI.contract([speed_of_light], [meter])
        >>> unit_system.defining_constants
        {speed_of_light: 1}
        >>> dimsys = unit_system.get_dimension_system()
        >>> dimsys.equivalent_dims(time, length)
        True
        >>> convert_to(second, meter, unit_system)
        299792458*meter
        >>> convert_to(speed_of_light, 1, unit_system)
        1

        The constants may have a value different from one. In reduced Planck
        units `8 \pi G = 1`:

        >>> from sympy import pi
        >>> from sympy.physics.units import gravitational_constant, hbar
        >>> from sympy.physics.units import boltzmann_constant
        >>> unit_system = SI.contract({speed_of_light: 1, hbar: 1,
        ...     boltzmann_constant: 1, gravitational_constant: 1/(8*pi)})
        >>> convert_to(gravitational_constant, 1, unit_system)
        1/(8*pi)

        See Also
        ========

        defining_constants
        sympy.physics.units.unit_system_conversion.convert_unit_system

        """
        if not isinstance(constants, dict):
            constants = dict.fromkeys(constants, S.One)
        constants = {constant: sympify(value) for constant, value in constants.items()}

        dimsys = self.get_dimension_system()
        old_base_dims = list(dimsys.base_dims)

        base_dims = []
        for unit in base_units:
            dimension = self.get_quantity_dimension(unit)
            if not dimension.name.is_Symbol:
                raise ValueError("the dimension of %s is not a single dimension" % unit)
            base_dims.append(dimension)

        columns = []
        for dimension in [*map(self.get_quantity_dimension, constants), *base_dims]:
            dependencies = dimsys.get_dimensional_dependencies(dimension)
            if not set(dependencies).issubset(old_base_dims):
                raise ValueError("%s is not defined in the unit system" % dimension)
            columns.append([dependencies.get(dim, 0) for dim in old_base_dims])
        number = len(columns)
        columns.extend(
            [int(i == j) for j in range(len(old_base_dims))] for i in range(len(old_base_dims)))
        _, pivots = Matrix(columns).T.rref()
        if pivots[:number] != tuple(range(number)):
            raise ValueError("the dimensions of constants and base units are not independent")
        base_dims.extend(old_base_dims[i - number] for i in pivots[number:])

        base_units = list(base_units)
        for unit in self._base_units:
            if unit not in base_units and self.get_quantity_dimension(unit) in base_dims:
                base_units.append(unit)
        unit_system = UnitSystem(
            base_units, self._units + tuple(constants) + tuple(units), name, description,
            None, derived_units, {**self._defining_constants, **constants})
        unit_system._parent = self
        unit_system._contraction = (constants, base_dims, [columns[i] for i in pivots])
        return unit_system

    def _contract_dimension_system(self, constants, base_dims, columns):
        """
        Dimension system of the unit system created by :meth:`contract`, it
        contains the dimensions and the scale factors of the quantities.

        Explanation
        ===========

        The columns of the matrix are the dimensions of the constants and the
        base dimensions of the new unit system, expressed with the base
        dimensions of the current one. The inverse matrix gives the exponents
        of constants and new base dimensions for every dimension of the
        current unit system. The exponents of the base dimensions are the new
        dimensional dependencies, the exponents of the constants determine
        the change of the scale factors.

        Only the quantities whose scale factor is set in the current unit
        system are considered. The scale factors of the other quantities are
        relative to them.

        Parameters
        ==========

        constants : dict
            The constants and their values.
        base_dims : list
            The base dimensions of the new unit system.
        columns : list
            The exponents of the base dimensions of the current unit system,
            for every constant and every dimension of ``base_dims``.
        """
        dimsys = self.get_dimension_system()
        old_base_dims = list(dimsys.base_dims)
        inverse = Matrix(columns).T.inv().tolist()
        number = len(constants)

        def exponents(dimension):
            dependencies = dimsys.get_dimensional_dependencies(dimension)
            if not set(dependencies).issubset(old_base_dims):
                return None
            vector = [dependencies.get(dim, 0) for dim in old_base_dims]
            return [sum(i*j for i, j in zip(row, vector) if j != 0) for row in inverse]

        dependencies = {}
        for dimension in dimsys.dimensional_dependencies:
            if dimension not in base_dims:
                dependencies[dimension] = {
                    dim: exponent for dim, exponent in zip(base_dims, exponents(dimension)[number:])
                    if exponent != 0}
        dimension_system = DimensionSystem(base_dims, dimensional_dependencies=dependencies)

        dimension_system._quantity_dimension_map.update(dimsys._quantity_dimension_map)
        dimension_system._quantity_dimension_map.update(self._quantity_dimension_map)
        ratios = [value/self.get_quantity_scale_factor(constant) for constant, value in constants.items()]
        factors: dict[Dimension, Expr | None] = {}
        for quantity in [*dimsys._quantity_scale_factors, *self._quantity_scale_factors]:
            dimension = self.get_quantity_dimension(quantity)
            if dimension not in factors:
                powers = exponents(dimension)
                factors[dimension] = None if powers is None else Mul(*[
                    ratio**power for ratio, power in zip(ratios, powers)])
            factor = factors[dimension]
            if factor is None:
                continue
            scale_factor = self.get_quantity_scale_factor(quantity)*factor
            if scale_factor.has(Float):
                scale_factor = scale_factor.evalf()
            dimension_system._quantity_scale_factors[quantity] = scale_factor
        dimension_system._quantity_scale_factors.update(constants)
        return dimension_system

    def get_dimension_system(self):
        """
        Return the dimension system of the unit system.
        """
        if self._dimension_system is None and self._contraction is not None:
            self._dimension_system = self._parent._contract_dimension_system(*self._contraction)
        return self._dimension_system

    def get_quantity_dimension(self, unit):
        qdm = self.get_dimension_system()._quantity_dimension_map
        if unit in qdm:
            return qdm[unit]
        return super().get_quantity_dimension(unit)

    def get_quantity_scale_factor(self, unit):
        qsfm = self.get_dimension_system()._quantity_scale_factors
        if unit in qsfm:
            return qsfm[unit]
        return super().get_quantity_scale_factor(unit)

    @staticmethod
    def get_unit_system(unit_system):
        if isinstance(unit_system, UnitSystem):
            return unit_system

        if unit_system not in UnitSystem._unit_systems:
            raise ValueError(
                "Unit system is not supported. Currently"
                "supported unit systems are {}".format(
                    ", ".join(sorted(UnitSystem._unit_systems))
                )
            )

        return UnitSystem._unit_systems[unit_system]

    @staticmethod
    def get_default_unit_system():
        return UnitSystem._unit_systems["SI"]

    @property
    def dim(self):
        """
        Give the dimension of the system.

        That is return the number of units forming the basis.
        """
        return len(self._base_units)

    @property
    def is_consistent(self):
        """
        Check if the underlying dimension system is consistent.
        """
        # test is performed in DimensionSystem
        return self.get_dimension_system().is_consistent

    @property
    def derived_units(self) -> dict[Dimension, Expr]:
        return self._derived_units

    @property
    def defining_constants(self) -> dict[Quantity, Expr]:
        """
        The physical constants that are pure numbers in the unit system,
        with their values.

        Examples
        ========

        >>> from sympy.physics.units.systems import SI, hartree_atomic_units
        >>> SI.defining_constants
        {}
        >>> hartree_atomic_units.defining_constants
        {coulomb_constant: 1, elementary_charge: 1, hbar: 1, electron_rest_mass: 1}

        See Also
        ========

        contract
        """
        return self._defining_constants

    def get_dimensional_expr(self, expr):
        from sympy.physics.units import Quantity
        if isinstance(expr, Mul):
            return Mul(*[self.get_dimensional_expr(i) for i in expr.args])
        elif isinstance(expr, Pow):
            return self.get_dimensional_expr(expr.base) ** expr.exp
        elif isinstance(expr, Add):
            return self.get_dimensional_expr(expr.args[0])
        elif isinstance(expr, Derivative):
            dim = self.get_dimensional_expr(expr.expr)
            for independent, count in expr.variable_count:
                dim /= self.get_dimensional_expr(independent)**count
            return dim
        elif isinstance(expr, Function):
            args = [self.get_dimensional_expr(arg) for arg in expr.args]
            if all(i == 1 for i in args):
                return S.One
            return expr.func(*args)
        elif isinstance(expr, Quantity):
            return self.get_quantity_dimension(expr).name
        return S.One

    def _collect_factor_and_dimension(self, expr):
        """
        Return tuple with scale factor expression and dimension expression.
        """
        from sympy.physics.units import Quantity
        if isinstance(expr, Quantity):
            return expr.scale_factor, expr.dimension
        elif isinstance(expr, Mul):
            factor = 1
            dimension = Dimension(1)
            for arg in expr.args:
                arg_factor, arg_dim = self._collect_factor_and_dimension(arg)
                factor *= arg_factor
                dimension *= arg_dim
            return factor, dimension
        elif isinstance(expr, Pow):
            factor, dim = self._collect_factor_and_dimension(expr.base)
            exp_factor, exp_dim = self._collect_factor_and_dimension(expr.exp)
            if self.get_dimension_system().is_dimensionless(exp_dim):
                exp_dim = 1
            return factor ** exp_factor, dim ** (exp_factor * exp_dim)
        elif isinstance(expr, Add):
            factor, dim = self._collect_factor_and_dimension(expr.args[0])
            for addend in expr.args[1:]:
                addend_factor, addend_dim = \
                    self._collect_factor_and_dimension(addend)
                if not self.get_dimension_system().equivalent_dims(dim, addend_dim):
                    raise ValueError(
                        'Dimension of "{}" is {}, '
                        'but it should be {}'.format(
                            addend, addend_dim, dim))
                factor += addend_factor
            return factor, dim
        elif isinstance(expr, Derivative):
            factor, dim = self._collect_factor_and_dimension(expr.args[0])
            for independent, count in expr.variable_count:
                ifactor, idim = self._collect_factor_and_dimension(independent)
                factor /= ifactor**count
                dim /= idim**count
            return factor, dim
        elif isinstance(expr, Application):
            # ``Application`` (rather than ``Function``) so that ``Min``/``Max``
            # and similar applied functions -- which are not ``Function``
            # subclasses -- are handled too, instead of falling through to the
            # ``else`` branch and being treated as opaque dimensionless objects.
            fds = [self._collect_factor_and_dimension(arg) for arg in expr.args]
            dims = [Dimension(1) if self.get_dimension_system().is_dimensionless(d[1])
                    else d[1] for d in fds]
            # ``_collect_factor_and_dimension`` must return a single
            # ``(factor, dimension)`` pair; the previous code splatted one
            # dimension per argument, which produced a malformed tuple (and a
            # crash when the result was unpacked) for multi-argument functions.
            # The arguments of a function must all carry the same dimension,
            # which is the dimension of the result.
            dim = dims[0] if dims else Dimension(1)
            for d in dims[1:]:
                if not self.get_dimension_system().equivalent_dims(dim, d):
                    raise ValueError(
                        'Dimension of "{}" is {}, but it should be {}'.format(
                            expr, d, dim))
            return expr.func(*(f[0] for f in fds)), dim
        elif isinstance(expr, Dimension):
            return S.One, expr
        else:
            return expr, Dimension(1)

    def get_units_non_prefixed(self) -> set[Quantity]:
        """
        Return the units of the system that do not have a prefix.
        """
        return set(filter(lambda u: not u.is_prefixed and not u.is_physical_constant, self._units))
