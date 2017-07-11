r"""
The :py:mod:`~sympy.holonomic` module is intended to deal with holonomic functions along
with various operations on them like addition, multiplication, composition,
integration and differentiation. The module also implements various kinds of
conversions such as converting holonomic functions to a different form and the
other way around.
"""

from .holonomic import DifferentialOperator, DifferentialOperators, \
    HolonomicFunction, expr_to_holonomic, from_hyper, from_meijerg
from .recurrence import HolonomicSequence, RecurrenceOperator, \
    RecurrenceOperators
