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
from sympy import Expr, cos, sin, sqrt, acos, atan2, S, pi, trigsimp, Abs
from sympy.core.sympify import sympify

if TYPE_CHECKING:
    from .quaternion import Quaternion


def geometric_product(q1: 'Quaternion', q2: 'Quaternion') -> 'Quaternion':
    """
    Compute the geometric product of two quaternions.
    
    In quaternion algebra, the geometric product is simply the quaternion multiplication.
    For vectors in 3D space, the geometric product combines dot and wedge products,
    but for quaternions it's just the standard quaternion product.
    
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


def commutator_product(q1: 'Quaternion', q2: 'Quaternion') -> 'Quaternion':
    """
    Compute the commutator product (analogous to outer/wedge product) of two quaternions.
    
    The commutator [q1, q2] = q1*q2 - q2*q1 captures the non-commutative part
    of quaternion multiplication and is related to the cross product for pure quaternions.
    
    Parameters
    ==========
    q1, q2 : Quaternion
        Input quaternions
        
    Returns
    =======
    Quaternion
        The commutator of q1 and q2
        
    Examples
    ========
    >>> from sympy import Quaternion
    >>> from sympy.algebras.geometric_algebra_ops import commutator_product
    >>> q1 = Quaternion(0, 1, 0, 0)  # i
    >>> q2 = Quaternion(0, 0, 1, 0)  # j
    >>> commutator_product(q1, q2)  # Should give 2*k
    0 + 0*i + 0*j + 2*k
    """
    return q1 * q2 - q2 * q1


def anticommutator_product(q1: 'Quaternion', q2: 'Quaternion') -> 'Quaternion':
    """
    Compute the anticommutator product (analogous to inner/dot product) of two quaternions.
    
    The anticommutator {q1, q2} = q1*q2 + q2*q1 captures the commutative part
    of quaternion multiplication.
    
    Parameters
    ==========
    q1, q2 : Quaternion
        Input quaternions
        
    Returns
    =======
    Quaternion
        The anticommutator of q1 and q2
        
    Examples
    ========
    >>> from sympy import Quaternion
    >>> from sympy.algebras.geometric_algebra_ops import anticommutator_product
    >>> q1 = Quaternion(1, 2, 3, 4)
    >>> q2 = Quaternion(5, 6, 7, 8)
    >>> anticommutator_product(q1, q2)
    -120 + 24*i + 60*j + 48*k
    """
    return q1 * q2 + q2 * q1


def dot_product(q1: 'Quaternion', q2: 'Quaternion') -> Expr:
    """
    Compute the dot product of two quaternions (treating them as 4D vectors).
    
    Parameters
    ==========
    q1, q2 : Quaternion
        Input quaternions
        
    Returns
    =======
    Expr
        The scalar dot product
        
    Examples
    ========
    >>> from sympy import Quaternion
    >>> from sympy.algebras.geometric_algebra_ops import dot_product
    >>> q1 = Quaternion(1, 2, 3, 4)
    >>> q2 = Quaternion(5, 6, 7, 8)
    >>> dot_product(q1, q2)
    70
    """
    return q1.a*q2.a + q1.b*q2.b + q1.c*q2.c + q1.d*q2.d


def to_rotor(angle: Expr, axis: 'Quaternion') -> 'Quaternion':
    """
    Convert an angle and axis to a rotor (unit quaternion for rotation).
    
    A rotor represents a rotation by angle θ around the given axis.
    The formula is: R = cos(θ/2) + sin(θ/2) * (normalized axis vector)
    
    Parameters
    ==========
    angle : Expr
        Rotation angle in radians
    axis : Quaternion
        Axis of rotation (should be a pure quaternion, will be normalized)
        
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
    
    # Extract vector part (ignore scalar part for axis)
    axis_vector = Quaternion(0, axis.b, axis.c, axis.d)
    
    # Calculate magnitude of vector part
    axis_magnitude = sqrt(axis.b**2 + axis.c**2 + axis.d**2)
    
    if axis_magnitude == 0:
        raise ValueError("Axis must have non-zero vector component")
    
    # Normalize the axis vector
    axis_normalized = axis_vector / axis_magnitude
    
    # Create rotor: cos(θ/2) + sin(θ/2)*normalized_axis
    half_angle = angle / 2
    scalar_part = cos(half_angle)
    vector_part = sin(half_angle) * axis_normalized
    
    return Quaternion(scalar_part, vector_part.b, vector_part.c, vector_part.d)


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
        and axis is a pure quaternion representing the normalized rotation axis.
        
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
    
    # Normalize the rotor to ensure it's a unit quaternion
    rotor = rotor.normalize()
    
    # Calculate the magnitude of the vector part
    vector_magnitude = sqrt(rotor.b**2 + rotor.c**2 + rotor.d**2)
    
    # Handle the identity rotation case (pure scalar quaternion)
    if vector_magnitude == 0:
        return (S.Zero, Quaternion(0, 1, 0, 0))  # Default to x-axis, zero angle
    
    # Calculate angle: θ = 2*atan2(|v|, w) where v is vector part, w is scalar part
    angle = 2 * atan2(vector_magnitude, Abs(rotor.a))
    
    # Calculate normalized axis
    axis = Quaternion(0, rotor.b/vector_magnitude, rotor.c/vector_magnitude, rotor.d/vector_magnitude)
    
    # Handle the case where the scalar part is negative (equivalent rotation with opposite angle/axis)
    if rotor.a < 0:
        angle = 2*pi - angle
        axis = -axis
    
    # Normalize angle to [0, 2π) or [-π, π] as preferred
    # Here we'll use [0, 2π)
    angle = angle % (2*pi)
    
    return (angle, axis)


def quaternion_exp(q: 'Quaternion') -> 'Quaternion':
    """
    Compute the exponential of a quaternion.
    
    For a quaternion q = w + v (where v is the vector part),
    exp(q) = exp(w) * (cos(|v|) + (v/|v|)*sin(|v|))
    
    Parameters
    ==========
    q : Quaternion
        Input quaternion
        
    Returns
    =======
    Quaternion
        The exponential of the quaternion
    """
    from .quaternion import Quaternion
    from sympy import exp
    
    if not isinstance(q, Quaternion):
        raise TypeError("Input must be a Quaternion")
    
    # Extract scalar and vector parts
    w = q.a
    v_mag = sqrt(q.b**2 + q.c**2 + q.d**2)
    
    # Handle the case where vector part is zero
    if v_mag == 0:
        return Quaternion(exp(w), 0, 0, 0)
    
    # Calculate exponential
    exp_w = exp(w)
    cos_v = cos(v_mag)
    sin_v = sin(v_mag)
    
    # Normalize vector part
    v_normalized = Quaternion(0, q.b/v_mag, q.c/v_mag, q.d/v_mag)
    
    return exp_w * (Quaternion(cos_v, 0, 0, 0) + sin_v * v_normalized)


def quaternion_log(q: 'Quaternion') -> 'Quaternion':
    """
    Compute the natural logarithm of a quaternion.
    
    For a unit quaternion q = cos(θ) + sin(θ)*v (where v is a unit vector),
    log(q) = θ*v
    
    For a general quaternion q = |q| * (cos(θ) + sin(θ)*v),
    log(q) = log(|q|) + θ*v
    
    Parameters
    ==========
    q : Quaternion
        Input quaternion (should be non-zero)
        
    Returns
    =======
    Quaternion
        The natural logarithm of the quaternion
    """
    from .quaternion import Quaternion
    from sympy import log, atan2
    
    if not isinstance(q, Quaternion):
        raise TypeError("Input must be a Quaternion")
    
    q_norm = q.norm()
    if q_norm == 0:
        raise ValueError("Cannot take logarithm of zero quaternion")
    
    # Normalize the quaternion
    q_normalized = q / q_norm
    
    # Extract vector part magnitude
    v_mag = sqrt(q_normalized.b**2 + q_normalized.c**2 + q_normalized.d**2)
    
    # Handle scalar quaternion case
    if v_mag == 0:
        return Quaternion(log(q_norm), 0, 0, 0)
    
    # Calculate the angle
    theta = atan2(v_mag, q_normalized.a)
    
    # Calculate logarithm
    log_norm = log(q_norm)
    vector_part = theta * Quaternion(0, q_normalized.b/v_mag, q_normalized.c/v_mag, q_normalized.d/v_mag)
    
    return Quaternion(log_norm, vector_part.b, vector_part.c, vector_part.d)
