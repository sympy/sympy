"""Cylindrical algebraic decomposition (CAD).

The projection phase is in :mod:`sympy.polys.cad.projection`, the exact real
algebraic sample points in :mod:`sympy.polys.cad.samplepoints` and the
lifting phase, with the main entry point
:func:`cylindrical_algebraic_decomposition`, in
:mod:`sympy.polys.cad.lifting`. Quantifier elimination and decision of
formulas built on top of it are in :mod:`sympy.polys.cad.qe`.
"""
from __future__ import annotations

from .projection import (projection_sets, mccallum_projection,
    hong_projection, squarefree_basis)
from .samplepoints import SamplePoint
from .lifting import (cylindrical_algebraic_decomposition, CAD, CADCell,
    NotWellOriented)
from .qe import quantifier_elimination, decide, sample_points, solution_set

__all__ = ['projection_sets', 'mccallum_projection', 'hong_projection',
    'squarefree_basis', 'SamplePoint', 'cylindrical_algebraic_decomposition',
    'CAD', 'CADCell', 'NotWellOriented', 'quantifier_elimination', 'decide',
    'sample_points', 'solution_set']
