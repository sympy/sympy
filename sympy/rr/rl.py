# Generic rules for SymPy

def rmid(isid):
    def ident_remove(expr):
        """ Remove identities """
        ids = map(isid, expr.args)
        if sum(ids) == 0:           # No identities. Common case
            return expr
        elif sum(ids) != len(ids):  # there is at least one non-identity
            return expr.__class__(*[arg for arg, x in zip(expr.args, ids)
                                        if  not x])
        else:
            first_id = (arg for arg, x in zip(expr.args, ids) if x).next()
            return expr.__class__(first_id)

    return ident_remove

def frequencies(coll):
    counts = {}
    for elem in coll:
        counts[elem] = counts.get(elem, 0) + 1
    return counts

def glom(mkglom):
    def conglomerate(expr):
        """ Conglomerate together identical args x + x -> 2x """
        freqs = frequencies(expr.args)
        return expr.__class__(*[arg if freqs[arg] == 1
                                    else mkglom(freqs[arg], arg)
                                    for arg in freqs])
    return conglomerate
