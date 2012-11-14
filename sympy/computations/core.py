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

    def _input_outputs(self):
        allin = tuple(unique(itertools.chain(
                        *[c.inputs  for c in self.computations])))
        allout = tuple(unique(itertools.chain(
                        *[c.outputs for c in self.computations])))
        inputs  = tuple(i for i in allin  if i not in allout)
        outputs = tuple(o for o in allout if o not in allin)
        return inputs, outputs

    @property
    def inputs(self):
        return self._input_outputs()[0]

    @property
    def outputs(self):
        return self._input_outputs()[1]

    def edges(self):
        return itertools.chain(*[c.edges() for c in self.computations])
