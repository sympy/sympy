from strat_pure import *
from traverse import bottom_up, top_down

def canon(*rules):
    """ Strategy for canonicalization

    Apply each rule in a bottom_up fashion through the tree.
    Do each one in turn.
    Keep doing this until there is no change.
    """
    return exhaust(chain(*map(bottom_up, rules)))
