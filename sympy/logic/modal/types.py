"""
Layer 1: Modal Type System.
Replaces untyped predicate logic with a stratified type universe supporting second-order quantification.
"""
from __future__ import annotations

from typing import Any
from sympy.core.basic import Basic
from sympy.core.symbol import Symbol
from sympy.logic.boolalg import Boolean

class Type(Basic):
    """Base class for all types in the modal type universe."""
    pass

class Universe(Type):
    """
    A cumulative type universe Level (Type₀, Type₁, ...).
    Universe(0) is the base level (predicative). Impredicative quantification
    requires elevated universe levels.
    """
    def __new__(cls, level: int) -> 'Universe':
        from sympy.core.numbers import Integer
        if isinstance(level, Integer):
            level = int(level)
        if not isinstance(level, int) or level < 0:
            raise ValueError("Universe level must be a non-negative integer.")
        return super().__new__(cls, level) # type: ignore

    @property
    def level(self) -> int:
        return int(self.args[0]) # type: ignore

    def __str__(self) -> str:
        return f"Type_{self.level}"

    def __repr__(self) -> str:
        return f"Universe({self.level})"

class BoolType(Type):
    """The type of propositions/booleans."""
    def __new__(cls) -> 'BoolType':
        return super().__new__(cls)

    def __str__(self) -> str:
        return "Bool"

    def __repr__(self) -> str:
        return "BoolType()"

class FunctionType(Type):
    """The type of functions (e.g. from Universe(0) to BoolType)."""
    def __new__(cls, domain: Type, codomain: Type) -> 'FunctionType':
        if not isinstance(domain, Type) or not isinstance(codomain, Type):
            raise TypeError("domain and codomain must be of type Type")
        return super().__new__(cls, domain, codomain)

    @property
    def domain(self) -> Type:
        return self.args[0]

    @property
    def codomain(self) -> Type:
        return self.args[1]

    def __str__(self) -> str:
        return f"({self.domain} -> {self.codomain})"

    def __repr__(self) -> str:
        return f"FunctionType({repr(self.domain)}, {repr(self.codomain)})"


class TypedSymbol(Symbol):
    """A symbol with an associated type."""
    def __new__(cls, name: str, type: Type, **kwargs: Any) -> 'TypedSymbol':
        obj = super().__new__(cls, name, **kwargs)
        obj._type = type
        return obj

    @property
    def type(self) -> Type:
        return self._type


class PredicateVariable(TypedSymbol, Boolean):
    """
    A second-order variable ranging over predicates of a given type.
    """
    def __new__(cls, name: str, type: Type, **kwargs: Any) -> 'PredicateVariable':
        if not isinstance(type, FunctionType) and not isinstance(type, BoolType):
            raise TypeError("PredicateVariable type must be a FunctionType or BoolType")
        return super().__new__(cls, name, type, **kwargs)

    def __call__(self, *args: Any) -> 'ModalPredicate':
        # Apply the predicate to arguments
        return ModalPredicate(self.name, self.type, args)


class ModalPredicate(Boolean):
    """
    A typed predicate whose validity is relative to a Kripke frame.
    Created by applying a PredicateVariable to arguments.
    """
    def __new__(cls, name: str, type: Type, args_tuple: tuple[Any, ...]) -> 'ModalPredicate':
        obj = super().__new__(cls, *args_tuple)
        obj._name = name
        obj._type = type
        return obj

    @property
    def name(self) -> str:
        return self._name

    @property
    def type(self) -> Type:
        return self._type


class GuardedFixedPoint(Boolean):
    """
    A controlled fixed-point constructor permitting the diagonal constructions
    required by Löb's theorem, guarded to ensure well-foundedness.
    """
    def __new__(cls, operator: Any, name: str) -> 'GuardedFixedPoint':
        obj = super().__new__(cls, operator)
        obj._name = name
        return obj

    @property
    def operator(self) -> Any:
        return self.args[0]

    @property
    def name(self) -> str:
        return self._name
