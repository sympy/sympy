import itertools

def unique(seq):
    seen = set()
    for item in seq:
        if item not in seen:
            seen.add(item)
            yield item

def intersect(a, b):
    return not not set(a).intersection(set(b))

class Computation(object):

    inputs  = None
    outputs = None

    def edges(self):
        inedges  = ((i, self) for i in self.inputs)
        outedges = ((self, o) for o in self.outputs)
        return itertools.chain(inedges, outedges)

    @property
    def variables(self):
        return itertools.chain(self.inputs, self.outputs)

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

    @property
    def variables(self):
        return unique(itertools.chain(
                        *[c.variables for c in self.computations]))

    def edges(self):
        return itertools.chain(*[c.edges() for c in self.computations])

    def dag_io(self):
        """ Return a dag of computations from inputs to outputs

        returns {A: {Bs}} such that A must occur before each of the Bs
        """
        return {A: set([B for B in self.computations
                          if intersect(A.outputs, B.inputs)])
                    for A in self.computations}

    def dag_oi(self):
        """ Return a dag of computations from outputs to inputs

        returns {A: {Bs}} such that A requires each of the Bs before it runs
        """
        return {A: set([B for B in self.computations
                          if intersect(A.inputs, B.outputs)])
                    for A in self.computations}

    def toposort(self):
        from sympy.utilities.iterables import _toposort
        return _toposort(self.dag_io())
