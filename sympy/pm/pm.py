def rubi_match(subject, pattern, var=None):
    '''
    Commutative pattern matcher for Rubi Integration.

    :param var: variable of integration
    '''
    return next(_rubi_match(subject, pattern, var, sig={}))

def _rubi_match(s, p, var, head=None, sig={}):
    '''
    Helper function for `rubi_match`.
    '''

    if p.is_Wild:
        yield {**sig, **{p: s}}
    elif p == var and s == var:
        yield sig
    elif p == var:
        yield None
    else:
        p_partition = p._partition_args(var)
        s_partition = s._partition_args(var)

        if p_partition.keys() != s_partition.keys():
            yield None

        else:
            if len(p_partition['Constant']) > 0:
                if p_partition['Constant'][0].is_Wild:
                    sig[p_partition['Constant'][0]] = s_partition['Constant'][0]

            keys = p_partition.keys()

            for k in keys:
                match = None
                matched = None
                if len(s_partition[k]) != len(p_partition[k]):
                    yield None
                    break
                else:
                    for i in s_partition[k]:
                        for j in p_partition[k]:
                            for res in _rubi_match(i, j, var, k, sig):
                                match = res
                                if match != None:
                                    matched = j # expression which is matched
                                    break
                            if match != None: # if matched, break the loop
                                break
                        if match == None:
                            yield None
                            break
                        else:
                            p_partition[k].remove(matched) # if matched, remove from search space
                            sig = {**sig, **match}
                    if match == None and len(s_partition[k]) > 0:
                        yield None
                        break

            yield sig
