#sympy.vector.kind

from sympy.core.kind import Kind, _NumberKind
from sympy.core.mul import Mul

class VectorKind(Kind):
    """
    Kind for all vector objects in SymPy.

    Parameters
    ==========
    No parameters.

    Examples
    ========

    Any instance of Vector class has kind ``VectorKind``:

    >>> from sympy.vector.coordsysrect import CoordSys3D
    >>> Sys = CoordSys3D('Sys')
    >>> Sys.i.kind
    VectorKind

    Operations between instances of Vector keep also have the kind ``VectorKind``:

    >>> from sympy.core.add import Add
    >>> v1 = Sys.i * 2 + Sys.j * 3 + Sys.k * 4
    >>> v2 = Sys.i * Sys.x + Sys.j * Sys.y + Sys.k * Sys.z
    >>> v1.kind
    VectorKind
    >>> v2.kind
    VectorKind
    >>> Add(v1, v2).kind
    VectorKind

    Subclasses of Vector also have the kind ``VectorKind``, such as
    Cross, VectorAdd, VectorMul or VectorZero.

    See Also
    ========

    sympy.core.kind.Kind
    sympy.matrices.kind.MatrixKind

    """
    def __new__(cls):
        return super().__new__(cls)

    def __repr__(self):
        return "VectorKind"

@Mul._kind_dispatcher.register(_NumberKind, VectorKind)
def num_vec_mul(k1, k2):
    """
    The result of a multiplication between a number and a Vector should be of VectorKind.
    """
    return VectorKind()
