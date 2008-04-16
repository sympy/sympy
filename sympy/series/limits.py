from gruntz import gruntz            

def limit(e, z, z0, dir="+"):
    """
    Compute the limit of e(z) at the point z0.

    z0 can be any expression, including oo and -oo.

    For dir="+" (default) it calculates the limit from the right
    (z->z0+) and for dir="-" the limit from the left (z->z0-). For infinite z0
    (oo or -oo), the dir argument doesn't matter.
    """
    return gruntz(e, z, z0, dir)
