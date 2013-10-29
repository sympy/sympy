from __future__ import print_function, division

from sympy import simplify
from sympy.tensor.components import *

__all__ = ["christoffel", "riemann_tensor", "ricci_curvature_tensor", "ricci_scalar"]

def christoffel(indices, metric, coordinate_list):
    """Calculate the christoffel of some metric

    Examples
    ========
    >>> from sympy import symbols, sin, cos
    >>> r, theta, phi = symbols("r theta phi")
    >>> g = Metric([Index("_", "i"), Index("_", "j")], [[r**2, 0], [0, r**2 * sin(theta)**2]])
    >>> christoffel([Index("^", "i"), Index("_", "j"), Index("_", "k")], g, [theta, phi])
    Tensor[^i, _j, _k] = [[[0, 0], [0, -1.0*sin(theta)*cos(theta)]], [[0, 1.0*cos(theta)/sin(theta)], [1.0*cos(theta)/sin(theta), 0]]]
    >>> christoffel([Index("_", "i"), Index("_", "j"), Index("_", "k")], g, [theta, phi])
    ValueError: Indices not in the right form as a Christoffel symbols
    """
    if (len(indices) != 3
        or indices[0].is_covariant()
        or indices[1].is_contravariant()
        or indices[2].is_contravariant()):
        raise ValueError("Indices not in the right form as a Christoffel symbol")
    m,k,l = Index("_", "m"), Index("_", "k"), Index("_", "l")
    mu,iu = Index("^", "m"), Index("^", "i")
    pgpx = partial_diff(metric.lowerlower().change_indices([m, k]),
                        Tensor([l], coordinate_list))
    return (metric.upperupper().change_indices([iu, mu])*
            (pgpx + pgpx.change_indices([m,l,k]).switch_indices([m,k,l])
                  - pgpx.change_indices([k,l,m]).switch_indices([m,k,l]))
            *0.5).einstein_summation().change_indices(indices)

def riemann_tensor(indices, christoffel, coordinate_list):
    """Calculate the Riemann curvature tensor

    Examples
    ========
    >>> from sympy import symbols, sin, cos
    >>> r, theta, phi = symbols("r theta phi")
    >>> g = Metric([Index("_", "i"), Index("_", "j")], [[r**2, 0], [0, r**2 * sin(theta)**2]])
    >>> gamma = christoffel([Index("^", "i"), Index("_", "j"), Index("_", "k")], g, [theta, phi])
    >>> iemann_tensor([Index("^", "i"), Index("_", "j"), Index("_", "k"), Index("_", "l")], gamma, [theta,phi])
    Tensor[^i, _j, _k, _l] = [[[[0, 0], [0, 0]], [[0, 1.0*sin(theta)**2], [-1.0*sin(theta)**2, 0]]], [[[0, -1.0], [1.0, 0]], [[0, 0], [0, 0]]]]
    """
    if (len(indices) != 4
        or indices[0].is_covariant()
        or indices[1].is_contravariant()
        or indices[2].is_contravariant()
        or indices[3].is_contravariant()):
        raise ValueError("Indices not in the right form as a Riemann tensor")
    i,j,k,l = Index("^", "i"), Index("_", "j"), Index("_", "k"), Index("_", "l")
    du,dd = Index("^", "d"), Index("_", "d")
    return (partial_diff(christoffel.change_indices([i,l,j]), Tensor([k], coordinate_list)).switch_indices([i,j,k,l])
            - partial_diff(christoffel.change_indices([i,k,j]), Tensor([l], coordinate_list)).switch_indices([i,j,k,l])
            + (christoffel.change_indices([i,k,dd])*christoffel.change_indices([du,l,j])).einstein_summation().switch_indices([i,j,k,l])
            - (christoffel.change_indices([i,l,dd])*christoffel.change_indices([du,k,j])).einstein_summation().switch_indices([i,j,k,l])
            ).change_indices(indices)

def ricci_curvature_tensor(indices, riemann):
    """Calculate the Ricci curvature tensor

    Examples
    ========
    >>> from sympy import symbols, sin, cos
    >>> r, theta, phi = symbols("r theta phi")
    >>> g = Metric([Index("_", "i"), Index("_", "j")], [[r**2, 0], [0, r**2 * sin(theta)**2]])
    >>> gamma = christoffel([Index("^", "i"), Index("_", "j"), Index("_", "k")], g, [theta, phi])
    >>> rabcd = riemann_tensor([Index("^", "i"), Index("_", "j"), Index("_", "k"), Index("_", "l")], gamma, [theta,phi])
    >>> ricci_curvature_tensor([Index("_", "i"), Index("_", "j")], rabcd)
    Tensor[_i, _j] = [[1.0, 0], [0, 1.0*sin(theta)**2]]
    """
    if (len(indices) != 2
        or indices[0].is_contravariant()
        or indices[1].is_contravariant()):
        raise ValueError("Indices not in the right form as a Ricci curvature tensor")
    i,j = Index("_", "i"), Index("_", "j")
    du,dd = Index("^", "d"), Index("_", "d")
    return riemann.change_indices([du,i,dd,j]).einstein_summation().change_indices(indices)

def ricci_scalar(metric, ricci_tensor):
    """Calculate the Ricci scalar

    Examples
    ========
    >>> from sympy import symbols, sin, cos
    >>> r, theta, phi = symbols("r theta phi")
    >>> g = Metric([Index("_", "i"), Index("_", "j")], [[r**2, 0], [0, r**2 * sin(theta)**2]])
    >>> gamma = christoffel([Index("^", "i"), Index("_", "j"), Index("_", "k")], g, [theta, phi])
    >>> rabcd = riemann_tensor([Index("^", "i"), Index("_", "j"), Index("_", "k"), Index("_", "l")], gamma, [theta,phi])
    >>> rab = ricci_curvature_tensor([Index("_", "i"), Index("_", "j")], rabcd)
    >>> ricci_scalar(g, rab)
    2.0/r**2

    >>> t, r, rs, theta, phi = symbols("t r rs theta phi")
    >>> coordinates = [t, r, theta, phi]
    >>> g = Metric([Index("_", "i"), Index("_", "j")], [[(1-rs/r), 0, 0, 0],
                                                        [0, -1/(1-rs/r), 0, 0],
                                                        [0, 0, -r**2, 0],
                                                        [0, 0, 0, - r**2 * sin(theta)**2]])
    >>> gamma = christoffel([Index("^", "i"), Index("_", "j"), Index("_", "k")], g, coordinates)
    >>> rabcd = riemann_tensor([Index("^", "i"), Index("_", "j"), Index("_", "k"), Index("_", "l")], gamma, coordinates)
    >>> rab = ricci_curvature_tensor([Index("_", "i"), Index("_", "j")], rabcd)
    >>> ricci_scalar(g, rab)
    0
    """
    ju,ku,jd,kd = Index("^", "j"), Index("^", "k"), Index("_", "j"), Index("_", "k")
    return simplify((metric.upperupper().change_indices([ju,ku])*ricci_tensor.change_indices([kd,jd])).einstein_summation())

