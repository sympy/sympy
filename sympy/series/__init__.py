"""A module that handles series: find a limit, order the series etc.
"""
from .order import Order
from .limits import limit, Limit
from .gruntz import gruntz
from .series import series
from .approximants import approximants, pade_approximant
from .residues import residue
from .sequences import SeqPer, SeqFormula, sequence, SeqAdd, SeqMul
from .fourier import fourier_series
from .formal import fps
from .limitseq import difference_delta, limit_seq

from sympy.core.singleton import S
EmptySequence = S.EmptySequence

O = Order

__all__ = ['Order', 'O', 'limit', 'Limit', 'gruntz', 'series', 'approximants',
        'pade_approximant', 'residue', 'EmptySequence', 'SeqPer', 'SeqFormula',
        'sequence', 'SeqAdd', 'SeqMul', 'fourier_series', 'fps', 'difference_delta',
        'limit_seq'
        ]
