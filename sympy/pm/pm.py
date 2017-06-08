from sympy.utilities.iterables import dict_merge
from sympy.core.compatibility import _nodes


def rubi_match(subject, pattern, var=None):
    '''
    Commutative pattern matcher for Rubi Integration.

    :param var: variable of integration
    '''
    return next(_rubi_match(subject, pattern, var, sig={}))


def _rubi_match(s, p, var, head=None, sig={}):
    '''
    Helper function for `rubi_match`.

    :param head: head of pattern. Eg: Add, Mul, log, etc
    '''
    if p.is_Wild:
        yield dict_merge(sig, {p: s})
    elif p == var and s == var:
        yield sig
    elif p == var:
        yield None
    else:
        p_partition = p._partition_args(var)
        s_partition = s._partition_args(var)

        if p_partition.keys() != s_partition.keys():
            sig = None
        else:
            if len(p_partition['Constant']) > 0:
                if p_partition['Constant'][0].is_Wild:
                    sig[p_partition['Constant'][0]] = s_partition['Constant'][0]

            keys = p_partition.keys()

            for h in keys:
                # Match subexpressions having same head
                if len(s_partition[h]) != len(p_partition[h]):
                    sig = None
                    break
                else:
                    #arrange expression in ascending order
                    s_partition[h].sort(key=lambda l: _nodes(l))
                    p_partition[h].sort(key=lambda l: _nodes(l))
                    for i in s_partition[h]:
                        matched = None # expression in subject which is matched
                        for j in p_partition[h]:
                            for res in _rubi_match(i, j, var, h, sig):
                                if res != None:
                                    sig = dict_merge(sig, res)
                                    matched = j
                                    break
                            if matched != None: # if matched, break the loop
                                break
                        if matched == None:
                            sig = None
                            break
                        else:
                            p_partition[h].remove(matched) # if matched, remove from search space

            yield sig
