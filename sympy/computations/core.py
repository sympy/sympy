import itertools

def unique(seq):
    seen = set()
    for item in seq:
        if item not in seen:
            seen.add(item)
            yield item

class Computation(object):

    inputs  = None
    outputs = None

    def edges(self):
        inedges  = ((i, self) for i in self.inputs)
        outedges = ((self, o) for o in self.outputs)
        return itertools.chain(inedges, outedges)

class CompositeComputation(Computation):

    @property
    def inputs(self):
        return unique(itertools.chain(*[c.inputs for c in self.computations]))

    @property
    def outputs(self):
        return unique(itertools.chain(*[c.outputs for c in self.computations]))

    @property
    def edges(self):
        return itertools.chain(*[c.edges() for c in self.computations])
