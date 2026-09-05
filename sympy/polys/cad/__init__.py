"""Cylindrical algebraic decomposition (CAD).

The projection phase is in :mod:`sympy.polys.cad.projection` and the exact
real algebraic sample points used by the lifting phase are in
:mod:`sympy.polys.cad.samplepoints`.
"""
from __future__ import annotations

from .projection import (projection_sets, mccallum_projection,
    hong_projection, squarefree_basis)
from .samplepoints import SamplePoint

__all__ = ['projection_sets', 'mccallum_projection', 'hong_projection',
    'squarefree_basis', 'SamplePoint']
