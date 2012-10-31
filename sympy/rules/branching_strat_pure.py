
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

