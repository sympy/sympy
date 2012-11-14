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
    """ An interface for a Computation

    Computations have inputs and outputs
    """

    inputs  = None
    outputs = None

    def edges(self):
        inedges  = ((i, self) for i in self.inputs)
        outedges = ((self, o) for o in self.outputs)
        return itertools.chain(inedges, outedges)

    @property
    def variables(self):
        return itertools.chain(self.inputs, self.outputs)

    def __add__(self, other):
        return CompositeComputation(self, other)

    def __str__(self):
        ins  = "["+', '.join(self.inputs) +"]"
        outs = "["+', '.join(self.outputs)+"]"
        return "%s -> %s -> %s"%(ins, str(self.__class__.__name__), outs)

class CompositeComputation(Computation):
    """ A computation composed of other computations """

    def __init__(self, *computations):
        self.computations = computations

    def _input_outputs(self):
        """ Find the inputs and outputs of the complete computation """
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

    def __str__(self):
        return "[[" + ", ".join(map(str, self.computations)) + "]]"

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
        """ Order computations in an executable order """
        from sympy.utilities.iterables import _toposort
        return _toposort(self.dag_io())
