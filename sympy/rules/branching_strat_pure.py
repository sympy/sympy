
def exhaust(brule):
    """ Apply a rule repeatedly until it has no effect """
    def exhaust_brl(expr):
        seen = set([expr])
        for nexpr in brule(expr):
            if nexpr not in seen:
                seen.add(nexpr)
                for nnexpr in exhaust_brl(nexpr):
                    yield nnexpr
        if seen == set([expr]):
            yield expr
    return exhaust_brl

