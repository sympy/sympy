
def flatten(cc):
    """ Flatten T(a, b, T(c, d), T2(e)) to T(a, b, c, d, T2(e)) """
    cls = cc.__class__
    computations = []
    for c in cc.computations:
        if c.__class__ == cls:
            computations.extend(c.computations)
        else:
            computations.append(c)
    return cls(*computations)

def unpack(cc):
    """ Rule to unpack singleton computations """
    if len(cc.computations) == 1:
        return cc.computations[0]
    else:
        return cc
