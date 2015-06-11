"""A module that handles series: find a limit, order the series etc.
"""
from .order import Order
from .limits import limit, Limit
from .gruntz import gruntz
from .series import series
from .residues import residue
from .sequences import (EmptySequence, SeqPer, SeqFormula, sequence, SeqAdd,
                        SeqMul)

O = Order

__all__ = ['Order', 'O', 'limit', 'Limit', 'gruntz', 'series', 'residue',
           'EmptySequence', 'SeqPer', 'SeqFormula', 'sequence',
           'SeqAdd', 'SeqMul']
