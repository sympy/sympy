from .decorators import deprecated
from .expr import Expr

@deprecated(useinstead="sympy.physics.quantum.trace.Tr",
    deprecated_since_version="1.10", issue=22330)
class Tr(Expr):
    def __new__(cls, *args):
        from sympy.physics.quantum.trace import Tr
        return Tr(*args)
