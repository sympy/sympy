"""Computational algebraic field theory. """

__all__ = [
    'minpoly', 'minimal_polynomial', 'primitive_element',

    'field_isomorphism',

    'to_number_field', 'isolate',
]

# Before `sympy.polys.numberfields` was a package, it was a module, containing
# what is now defined in the `minpoly`, `isomorphism` and `numbers` modules.
# To maintain backwards compatibility, we have to import everything that was
# available in the old `sympy.polys.numberfields` module.
from .minpoly import (
    _choose_factor, _separate_sq, _minimal_polynomial_sq,
    _minpoly_op_algebraic_element, _invertx, _muly, _minpoly_pow,
    _minpoly_add, _minpoly_mul, _minpoly_sin, _minpoly_cos, _minpoly_tan,
    _minpoly_exp, _minpoly_rootof, _minpoly_compose,
    minimal_polynomial, _minpoly_groebner, minpoly, _switch_domain,
    _linsolve, primitive_element
)
from .isomorphism import (
    is_isomorphism_possible, field_isomorphism_pslq, field_isomorphism_factor,
    field_isomorphism,
)
from .numbers import (
    to_number_field, IntervalPrinter, isolate
)
