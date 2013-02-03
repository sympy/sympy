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

rules.rl         - some fundamental rules
rules.strat_pure - generic non-SymPy specific strategies
rules.traverse   - strategies that traverse a SymPy tree
rules.strat      - some conglomerate strategies that do depend on SymPy
"""

import rl
import traverse
from rl import rm_id, unpack, flatten, sort, glom, distribute, rebuild
from util import new
from strat import (canon, condition, debug, typed, chain, null_safe, do_one,
        exhaust)
import branch
