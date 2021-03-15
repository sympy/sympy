"""
This module has all functions related to inference of optic waves.

**Contains**

* Young's double-slit experiment

"""

__all__ = ['fringe_spacing']


def fringe_spacing(wavelength, slits_distance, screen_distance):
    """
    This function provides the fringe spacing
    when the following parameters are provided.

    ===========================================

    Parameters:

    wavelength : Wavelength of light
    slits_distance : Distance between two slits
    screen_distance : Distance between slits and screen

    """
    return ((wavelength*slits_distance)/screen_distance)


def wavelength(fringe_spacing, slits_distance, screen_distance):
    """
    This function provides the wavelength of light
    when the following parameters are provided.

    ===========================================

    Parameters:

    fringe_spacing : Fringe spacing
    slits_distance : Distance between two slits
    screen_distance : Distance between slits and screen

    """
    return ((fringe_spacing*screen_distance)/slits_distance)
