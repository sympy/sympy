"""A module that handles series: find a limit, order the series etc.
"""
from .order import Order
from .limits import limit, Limit
from .gruntz import gruntz
from .series import series
from .residues import residue
from .sequences import EmptySequence, SeqPer

O = Order
