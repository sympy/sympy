from .core import exhaust, multiplex
from .traverse import top_down
from collections.abc import Generator
from typing import Any, Callable


def canon(*rules) -> Callable[..., Generator[Any, Any, None]]:
    """ Strategy for canonicalization

    Apply each branching rule in a top-down fashion through the tree.
    Multiplex through all branching rule traversals
    Keep doing this until there is no change.
    """
    return exhaust(multiplex(*map(top_down, rules)))
