"""Static type checks for subs(). mypy checks them; they never execute."""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing_extensions import assert_type
    from sympy import And, Basic, Expr, Symbol
    from sympy.logic.boolalg import Boolean

    x, y = Symbol("x"), Symbol("y")
    e = x**2

    d = {x: y + 1}
    assert_type(e.subs(d), Expr)
    assert_type(e.subs(x, 1, simultaneous=True), Expr)
    assert_type(e.subs({x: 1}, simultaneous=True), Expr)
    assert_type(e.subs(e.args[0], 1), Expr)  # args are typed Basic

    assert_type(e.subs(x, Basic()), Basic)
    b = {x: Basic()}
    assert_type(e.subs(b), Basic)

    assert_type(And(x, y).subs(Basic(), True), Boolean)
    assert_type(And(x, y).subs(x, True, simultaneous=True), Boolean)
