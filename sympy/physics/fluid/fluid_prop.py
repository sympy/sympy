"""
Module to handle fluid properties.
"""
from sympy.core import S, Symbol, diff

def density(mass, volume):
    """
    Density: The density, or more precisely, the volumetric mass density, of a substance is its mass per unit volume
    ================================================================================================================
    Input Paramters :
     m - Mass of the Fluid
     v - Volume of the Fluid
    Output Parameter :
     d - Density of the Fluid
    ================================================================================================================
    """
    return mass/volume


def specific_weight(mass, volume , gravity):
    """
    Specific Weight: The specific weight (also known as the unit weight) is the weight per unit volume of a material.
    ================================================================================================================
    Input Paramters :
     m - Mass of the Fluid
     v - Volume of the Fluid
     g - Acceleration due to Gravity
    Output Parameter :
     w - Specific Weight
    ================================================================================================================
    """
    return (mass*gravity)/volume
