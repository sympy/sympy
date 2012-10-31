from branch_strat_pure import exhaust, multiplex
from branch_traverse import top_down

def canon(*rules):
    """ Strategy for canonicalization

    Apply each rule in a bottom_up fashion through the tree.
    Do each one in turn.
    Keep doing this until there is no change.
    """
    return exhaust(multiplex(*map(top_down, rules)))
