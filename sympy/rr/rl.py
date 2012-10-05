# Generic rules for SymPy

def rmid(isid):
    def ident_remove(expr):
        """ Remove identities """
        ids = map(isid, expr.args)
        if sum(ids) > 0 and len(expr.args) > 1:
            return expr.__class__(*[arg for arg, x in zip(expr.args, ids) if x])
        else:
            return expr
    return ident_remove
