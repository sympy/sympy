"""Biomechanics extension for SymPy.

Includes biomechanics-related constructs which allows users to extend multibody
models created using `sympy.physics.mechanics` into biomechanical or
musculoskeletal models involding musculotendons and activation dynamics.

"""

from .characteristic import (
   fl_T_de_groote_2016,
)


__all__ = [
   # Musculotendon characteristic curve functions
   'fl_T_de_groote_2016',
]
