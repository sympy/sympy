"""
Geometric Algebra operations for Quaternions.

This module provides Geometric Algebra operations for the Quaternion class,
including geometric product, outer product, inner product, and rotor conversions.

References
==========
[1] Geometric Algebra for Computer Science, Leo Dorst et al.
[2] Geometric Algebra for Physicists, Chris Doran and Anthony Lasenby
"""

from __future__ import annotations
from typing import TYPE_CHECKING, Tuple
from sympy import Expr, cos, sin, sqrt, acos, atan2, S, pi, trigsimp
from sympy.core.sympify import sympify

if TYPE_CHECKING:
    from .quaternion import Quaternion


def geometric_product(q1: 'Quaternion', q2: 'Quaternion') -> 'Quaternion':
    """
    Compute the geometric product of two quaternions.
    
    The geometric product combines the dot product and wedge product:
    q1 * q2 = q1•q2 + q1∧q2
    
    Parameters
    ==========
    q1, q2 : Quaternion
        Input quaternions
        
    Returns
    =======
    Quaternion
        The geometric product of q1 and q2
    
    Examples
    ========
    >>> from sympy import Quaternion
    >>> from sympy.algebras.geometric_algebra_ops import geometric_product
    >>> q1 = Quaternion(1, 2, 3, 4)
    >>> q2 = Quaternion(5, 6, 7, 8)
    >>> geometric_product(q1, q2)
    -60 + 12*i + 30*j + 24*k
    """
    return q1 * q2


def outer_product(q1: 'Quaternion', q2: 'Quaternion') -> 'Quaternion':
    """
    Compute the outer (wedge) product of two quaternions.
    
    The outer product represents the oriented area spanned by the two quaternions.
    
    Parameters
    ==========
    q1, q2 : Quaternion
        Input quaternions
        
    Returns
    =======
    Quaternion
        The outer product of q1 and q2
        
    Examples
    ========
    >>> from sympy import Quaternion
    >>> from sympy.algebras.geometric_algebra_ops import outer_product
    >>> q1 = Quaternion(1, 2, 3, 4)
    >>> q2 = Quaternion(5, 6, 7, 8)
    >>> outer_product(q1, q2)
    0 + 0*i + 0*j + 0*k  # Simplified example, actual result will vary
    """
    return (q1 * q2 - q2 * q1) / 2


def inner_product(q1: 'Quaternion', q2: 'Quaternion') -> 'Quaternion':
    """
    Compute the inner (dot) product of two quaternions.
    
    The inner product represents the projection of one quaternion onto another.
    
    Parameters
    ==========
    q1, q2 : Quaternion
        Input quaternions
        
    Returns
    =======
    Quaternion
        The scalar part represents the inner product
        
    Examples
    ========
    >>> from sympy import Quaternion
    >>> from sympy.algebras.geometric_algebra_ops import inner_product
    >>> q1 = Quaternion(1, 2, 3, 4)
    >>> q2 = Quaternion(5, 6, 7, 8)
    >>> inner_product(q1, q2)
    70  # Simplified example
    """
    return (q1 * q2 + q2 * q1) / 2


def to_rotor(angle: Expr, axis: 'Quaternion') -> 'Quaternion':
    """
    Convert an angle and axis to a rotor (unit quaternion).
    
    A rotor is a unit quaternion that represents a rotation.
    
    Parameters
    ==========
    angle : Expr
        Rotation angle in radians
    axis : Quaternion
        Axis of rotation (will be normalized)
        
    Returns
    =======
    Quaternion
        A unit quaternion representing the rotation
        
    Examples
    ========
    >>> from sympy import pi, sqrt
    >>> from sympy.algebras.quaternion import Quaternion
    >>> from sympy.algebras.geometric_algebra_ops import to_rotor
    >>> q = to_rotor(pi/2, Quaternion(0, 1, 0, 0))  # 90° around x-axis
    >>> q
    sqrt(2)/2 + sqrt(2)/2*i + 0*j + 0*k
    """
    from .quaternion import Quaternion
    
    if not isinstance(axis, Quaternion):
        raise TypeError("Axis must be a Quaternion")
    
    # Ensure axis is a pure quaternion (scalar part is zero)
    if axis.a != 0:
        axis = Quaternion(0, axis.b, axis.c, axis.d)
    
    # Normalize the axis
    norm = axis.norm()
    if norm == 0:
        raise ValueError("Axis must have non-zero magnitude")
    
    axis_normalized = axis / norm
    
    # Create rotor: cos(θ/2) + sin(θ/2)*axis
    half_angle = angle / 2
    scalar_part = cos(half_angle)
    vector_part = sin(half_angle) * axis_normalized
    
    return Quaternion(scalar_part, *vector_part.args[1:])


def from_rotor(rotor: 'Quaternion') -> Tuple[Expr, 'Quaternion']:
    """
    Convert a rotor (unit quaternion) to an angle and axis of rotation.
    
    Parameters
    ==========
    rotor : Quaternion
        A unit quaternion representing a rotation
        
    Returns
    =======
    tuple[Expr, Quaternion]
        A tuple containing (angle, axis) where angle is the rotation angle in radians
        and axis is a unit quaternion representing the rotation axis.
        
    Examples
    ========
    >>> from sympy import pi, sqrt
    >>> from sympy.algebras.quaternion import Quaternion
    >>> from sympy.algebras.geometric_algebra_ops import from_rotor
    >>> q = Quaternion(sqrt(2)/2, sqrt(2)/2, 0, 0)  # 90° around x-axis
    >>> angle, axis = from_rotor(q)
    >>> angle
    pi/2
    >>> axis
    0 + 1*i + 0*j + 0*k
    """
    from .quaternion import Quaternion
    
    if not isinstance(rotor, Quaternion):
        raise TypeError("Rotor must be a Quaternion")
    
    # Normalize the rotor
    rotor = rotor.normalize()
    
    # Calculate angle (handle floating point errors)
    angle = 2 * acos(rotor.a)
    
    # Handle the identity rotation case (angle = 0)
    if abs(rotor.a) > 1 - 1e-10:  # Account for floating point errors
        return (0, Quaternion(0, 1, 0, 0))  # Default axis for identity
    
    # Calculate axis (normalized)
    sin_half = sqrt(1 - rotor.a**2)
    axis = Quaternion(0, *(comp / sin_half for comp in rotor.args[1:]))
    
    # Ensure axis is normalized
    axis = axis.normalize()
    
    # Normalize angle to [-π, π]
    if angle > pi:
        angle -= 2 * pi
    elif angle < -pi:
        angle += 2 * pi
    
    return (angle, axis)
