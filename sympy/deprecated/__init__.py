"""This module contains deprecations that could not stay in their original
module for some reason.

Such reasons include:
- Original module had to be removed.
- Adding @deprecated to a declaration caused an import cycle.

Since no modules in SymPy ever depend on deprecated code, SymPy always imports
this last, after all other modules have been imported.
"""

__all__ = []

from .class_registry import C, ClassRegistry
__all__ += ["C", "ClassRegistry"]
