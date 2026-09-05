"""Cylindrical algebraic decomposition (CAD).

The projection phase is in :mod:`sympy.polys.cad.projection`.
"""
from __future__ import annotations

from .projection import (projection_sets, mccallum_projection,
    hong_projection, squarefree_basis)

__all__ = ['projection_sets', 'mccallum_projection', 'hong_projection',
    'squarefree_basis']
