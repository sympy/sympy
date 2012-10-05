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
