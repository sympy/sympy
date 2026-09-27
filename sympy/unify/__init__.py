""" Unification in SymPy

See sympy.unify.core docstring for algorithmic details

See http://matthewrocklin.com/blog/work/2012/11/01/Unification/ for discussion
"""
from __future__ import annotations

from .usympy import unify, rebuild
from .rewrite import rewriterule

__all__ = [
    'unify', 'rebuild',

    'rewriterule',
]
