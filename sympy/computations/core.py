class Computation(object):

    def inputs(self):
        raise NotImplementedError()

    def outputs(self):
        raise NotImplementedError()

    def edges(self):
        import itertools
        inedges  = ((i, self) for i in self.inputs())
        outedges = ((self, o) for o in self.outputs())
        return itertools.chain(inedges, outedges)
