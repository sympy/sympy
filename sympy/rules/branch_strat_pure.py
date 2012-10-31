
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

def debug(brule, file=None):
    """ Print out before and after expressions each time rule is used """
    if file is None:
        from sys import stdout
        file = stdout
    def debug_brl(expr):
        for result in brule(expr):
            if result != expr:
                file.write("Rule: %s\n"%brule.func_name)
                file.write("In: %s\nOut: %s\n\n"%(expr, result))
            yield result
    return debug_brl

def multiplex(*brules):
    """ Multiplex many branching rules into one """
    def multiplex_brl(expr):
        seen = set([])
        for brl in brules:
            for nexpr in brl(expr):
                if nexpr not in seen:
                    seen.add(nexpr)
                    yield nexpr
    return multiplex_brl

def condition(cond, brule):
    """ Only apply rule if condition is true """
    def conditioned_brl(expr):
        if cond(expr):
            for x in brule(expr): yield x
        else:
            pass
    return conditioned_brl

def notempty(brule):
    def notempty_brl(expr):
        yielded = False
        for nexpr in brule(expr):
            yielded = True
            yield nexpr
        if not yielded:
            yield expr
    return notempty_brl
