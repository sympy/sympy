from sympy.vector import gradient, divergence,curl

def laplacian(funct):
    """
    

    Parameters
    ----------
    funct : CoordSys3D type
        The function that we are finding the laplacian of

    Returns
    -------
    CoordSys3D type
    This returns the laplacian of the function, which is the divergence of the gradient of the function.
        

    """
    return divergence(gradient(funct))

def hasNoExtrema(funct):
    """
    

    Parameters
    ----------
    funct : CoordSys3D type
        The function that we are finding if it has extrema

    Returns
    -------
    CoordSys3D type
    This returns True if the laplacian of the function is 0 (If the laplacian of the function is 0, then the gradient field of the function has no sinks nor sources, which means there are no maximums or minimums)
        

    """
    if laplacian(funct) == 0:
        return True
    return False
