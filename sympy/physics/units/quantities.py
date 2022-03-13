"""
Physical quantities.
"""

from sympy.core.expr import AtomicExpr
from sympy.core.symbol import Symbol
from sympy.core.sympify import sympify
from sympy.physics.units.dimensions import _QuantityMapper
from sympy.physics.units.prefixes import Prefix
from sympy.utilities.exceptions import (sympy_deprecation_warning,
                                        SymPyDeprecationWarning,
                                        ignore_warnings)


class Quantity(AtomicExpr):
    """
    Physical quantity: can be a unit of measure, a constant or a generic quantity.
    """

    is_commutative = True
    is_real = True
    is_number = False
    is_nonzero = True
    _diff_wrt = True

    def __new__(cls, name, abbrev=None, dimension=None, scale_factor=None,
                latex_repr=None, pretty_unicode_repr=None,
                pretty_ascii_repr=None, mathml_presentation_repr=None,
                **assumptions):

        if not isinstance(name, Symbol):
            name = Symbol(name)

        # For Quantity(name, dim, scale, abbrev) to work like in the
        # old version of SymPy:
        if not isinstance(abbrev, str) and not \
                   isinstance(abbrev, Symbol):
            dimension, scale_factor, abbrev = abbrev, dimension, scale_factor

        if dimension is not None:
            sympy_deprecation_warning(
                """
                The 'dimension' argument to to Quantity() is deprecated.
                Instead use the unit_system.set_quantity_dimension() method.
                """,
                deprecated_since_version="1.3",
                active_deprecations_target="deprecated-quantity-dimension-scale-factor"
            )

        if scale_factor is not None:
            sympy_deprecation_warning(
                """
                The 'scale_factor' argument to to Quantity() is deprecated.
                Instead use the unit_system.set_quantity_scale_factors()
                method.
                """,
                deprecated_since_version="1.3",
                active_deprecations_target="deprecated-quantity-dimension-scale-factor"
            )

        if abbrev is None:
            abbrev = name
        elif isinstance(abbrev, str):
            abbrev = Symbol(abbrev)

        obj = AtomicExpr.__new__(cls, name, abbrev)
        obj._name = name
        obj._abbrev = abbrev
        obj._latex_repr = latex_repr
        obj._unicode_repr = pretty_unicode_repr
        obj._ascii_repr = pretty_ascii_repr
        obj._mathml_repr = mathml_presentation_repr

        if dimension is not None:
            # TODO: remove after deprecation:
            with ignore_warnings(SymPyDeprecationWarning):
                obj.set_dimension(dimension)

        if scale_factor is not None:
            # TODO: remove after deprecation:
            with ignore_warnings(SymPyDeprecationWarning):
                obj.set_scale_factor(scale_factor)
        return obj

    def set_dimension(self, dimension, unit_system="SI"):
        sympy_deprecation_warning(
            f"""
            Quantity.set_dimension() is deprecated. Use either
            unit_system.set_quantity_dimension() or
            {self}.set_global_dimension() instead.
            """,
            deprecated_since_version="1.5",
            active_deprecations_target="deprecated-quantity-methods",
        )
        from sympy.physics.units import UnitSystem
        unit_system = UnitSystem.get_unit_system(unit_system)
        unit_system.set_quantity_dimension(self, dimension)

    def set_scale_factor(self, scale_factor, unit_system="SI"):
        sympy_deprecation_warning(
            f"""
            Quantity.set_scale_factor() is deprecated. Use either
            unit_system.set_quantity_scale_factors() or
            {self}.set_global_relative_scale_factor() instead.
            """,
            deprecated_since_version="1.5",
            active_deprecations_target="deprecated-quantity-methods",
        )
        from sympy.physics.units import UnitSystem
        unit_system = UnitSystem.get_unit_system(unit_system)
        unit_system.set_quantity_scale_factor(self, scale_factor)

    def set_global_dimension(self, dimension):
        _QuantityMapper._quantity_dimension_global[self] = dimension

    def set_global_relative_scale_factor(self, scale_factor, reference_quantity):
        """
        Setting a scale factor that is valid across all unit system.
        """
        from sympy.physics.units import UnitSystem
        scale_factor = sympify(scale_factor)
        # replace all prefixes by their ratio to canonical units:
        scale_factor = scale_factor.replace(
            lambda x: isinstance(x, Prefix),
            lambda x: x.scale_factor
        )
        scale_factor = sympify(scale_factor)
        UnitSystem._quantity_scale_factors_global[self] = (scale_factor, reference_quantity)
        UnitSystem._quantity_dimensional_equivalence_map_global[self] = reference_quantity

    @property
    def name(self):
        return self._name

    @property
    def dimension(self):
        from sympy.physics.units import UnitSystem
        unit_system = UnitSystem.get_default_unit_system()
        return unit_system.get_quantity_dimension(self)

    @property
    def abbrev(self):
        """
        Symbol representing the unit name.

        Prepend the abbreviation with the prefix symbol if it is defines.
        """
        return self._abbrev

    @property
    def scale_factor(self):
        """
        Overall magnitude of the quantity as compared to the canonical units.
        """
        from sympy.physics.units import UnitSystem
        unit_system = UnitSystem.get_default_unit_system()
        return unit_system.get_quantity_scale_factor(self)

    def _eval_is_positive(self):
        return True

    def _eval_is_constant(self):
        return True

    def _eval_Abs(self):
        return self

    def _eval_subs(self, old, new):
        if isinstance(new, Quantity) and self != old:
            return self

    @staticmethod
    def get_dimensional_expr(expr, unit_system="SI"):
        sympy_deprecation_warning(
            """
            Quantity.get_dimensional_expr() is deprecated. It is now
            associated with UnitSystem objects. The dimensional relations
            depend on the unit system used. Use
            unit_system.get_dimensional_expr() instead.
            """,
            deprecated_since_version="1.5",
            active_deprecations_target="deprecated-quantity-methods",
        )
        from sympy.physics.units import UnitSystem
        unit_system = UnitSystem.get_unit_system(unit_system)
        return unit_system.get_dimensional_expr(expr)

    @staticmethod
    def _collect_factor_and_dimension(expr, unit_system="SI"):
        """Return tuple with scale factor expression and dimension expression."""
        sympy_deprecation_warning(
            """
            Quantity._collect_factor_and_dimension() is deprecated. This
            method has been moved to the UnitSystem class. Use
            unit_system._collect_factor_and_dimension(expr) instead.
            """,
            deprecated_since_version="1.5",
            active_deprecations_target="deprecated-quantity-methods",
        )
        from sympy.physics.units import UnitSystem
        unit_system = UnitSystem.get_unit_system(unit_system)
        return unit_system._collect_factor_and_dimension(expr)

    def _latex(self, printer):
        if self._latex_repr:
            return self._latex_repr
        else:
            return r'\text{{{}}}'.format(self.args[1] \
                          if len(self.args) >= 2 else self.args[0])

    def convert_to(self, other, unit_system="SI"):
        """
        Convert the quantity to another quantity of same dimensions.

        Examples
        ========

        >>> from sympy.physics.units import speed_of_light, meter, second
        >>> speed_of_light
        speed_of_light
        >>> speed_of_light.convert_to(meter/second)
        299792458*meter/second

        >>> from sympy.physics.units import liter
        >>> liter.convert_to(meter**3)
        meter**3/1000
        """
        from .util import convert_to
        return convert_to(self, other, unit_system)

    @property
    def free_symbols(self):
        """Return free symbols from quantity."""
        return set()
