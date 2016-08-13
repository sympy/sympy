"""
The Holonomic module is intended to deal with Holonomic Functions along
with various operations on them like addition, multiplication, composition,
integration and differentiation. The module also implements various kinds of
conversions like converting a Holonomic Function to an other form and the
other way around.
"""

from .holonomic import (DifferentialOperator, HolonomicFunction, DifferentialOperators,
    from_hyper, from_meijerg, expr_to_holonomic)
from .recurrence import RecurrenceOperators, RecurrenceOperator, HolonomicSequence
