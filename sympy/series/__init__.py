"""A module that handles series: find a limit, order the series etc.
"""
from .order import Order
from .limits import limit, Limit
from .gruntz import gruntz
from .series import series
from .residues import residue
from .sequences import (EmptySequence, SeqPer, SeqFormula, sequence, SeqAdd,
                        SeqMul)
from .fourier import fourier_series
from .formal import fps
from .limitseq import difference_delta, limitseq

O = Order
dd = difference_delta

__all__ = ['Order', 'O', 'limit', 'Limit', 'gruntz', 'series', 'residue',
           'EmptySequence', 'SeqPer', 'SeqFormula', 'sequence',
           'SeqAdd', 'SeqMul', 'fourier_series', 'fps', 'difference_delta',
           'dd', 'limitseq']
