import itertools as it

def interleave(seqs, pass_exceptions=()):
    iters = it.imap(iter, seqs)
    while iters:
        newiters = []
        for itr in iters:
            try:
                yield next(itr)
                newiters.append(itr)
            except (StopIteration,) + tuple(pass_exceptions):
                pass
        iters = newiters

def unique(seq):
    seen = set()
    for item in seq:
        if item not in seen:
            seen.add(item)
            yield item

def fnmap(fns, val):
    return (fn(val) for fn in fns)

