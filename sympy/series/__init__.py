"""A module that handles series: find a limit, order the series etc.
"""
from .approximants import approximants
from .formal import fps
from .fourier import fourier_series
from .gruntz import gruntz
from .limits import Limit, limit
from .limitseq import difference_delta, limit_seq
from .order import Order
from .residues import residue
from .sequences import EmptySequence, SeqAdd, SeqFormula, SeqMul, SeqPer, \
    sequence
from .series import series

O = Order

__all__ = ['Order', 'O', 'limit', 'Limit', 'gruntz', 'series', 'residue',
           'EmptySequence', 'SeqPer', 'SeqFormula', 'sequence',
           'SeqAdd', 'SeqMul', 'fourier_series', 'fps', 'difference_delta',
           'limit_seq']
