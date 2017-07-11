""" Rewrite Rules

DISCLAIMER: This module is experimental. The interface is subject to change.

A rule is a function that transforms one expression into another

    Rule :: Expr -> Expr

A strategy is a function that says how a rule should be applied to a syntax
tree. In general strategies take rules and produce a new rule

    Strategy :: [Rules], Other-stuff -> Rule

This allows developers to separate a mathematical transformation from the
algorithmic details of applying that transformation. The goal is to separate
the work of mathematical programming from algorithmic programming.

Submodules

strategies.rl         - some fundamental rules
strategies.core       - generic non-SymPy specific strategies
strategies.traverse   - strategies that traverse a SymPy tree
strategies.tools      - some conglomerate strategies that do depend on SymPy
"""

from . import branch, rl, traverse
from .core import chain, condition, debug, do_one, exhaust, minimize, \
    null_safe, tryit
from .rl import distribute, flatten, glom, rebuild, rm_id, sort, unpack
from .tools import canon, typed
from .util import new
