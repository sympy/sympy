from sympy.core import Basic, Expr
from sympy.printing.pretty.stringpict import prettyForm
from sympy.printing.pretty.pretty_symbology import pretty_symbol

class U(Basic):
    """Up decorator for indices."""
    pass

class D(Basic):
    """Down decorator for indices."""
    pass

class MetricTangentSpace(Basic):
    """A tangent space with a metric.

    It is usable both for abstract Penrose index notation that is independent
    of coordinates as well as for coordinate dependent representation.
    """
    #__init__(self, name, metric_matrix):

    # XXX args aliases that can not be implemented in __init__ because
    # the rewrite rules do not call __init__
    @property
    def name(self):
        return self.args[0]

    @property
    def metric_matrix(self):
        return self.args[1]

    def metric(self, *indices):
        return MetricTensorComponent(self, *indices)

class TensorComponent(Expr):
    """
    It is usable both for abstract Penrose index notation that is independent
    of coordinates as well as for coordinate dependent representation.
    """
    #__init__(self, name, tangent_space, *indices):

    # XXX args aliases that can not be implemented in __init__ because
    # the rewrite rules do not call __init__
    @property
    def name(self):
        return self.args[0]

    @property
    def indices(self):
        return self.args[2:]

    def _pretty(self, printer, *args):
        tensor = prettyForm(pretty_symbol(self.name))
        up = down = prettyForm(' '*tensor.width())
        for i in self.indices:
            i_pretty = prettyForm(pretty_symbol(i.args[0].name))
            if isinstance(i, U):
                up = prettyForm(*up.right(i_pretty))
                down = prettyForm(*down.right(' '*i_pretty.width()))
            else:
                down = prettyForm(*down.right(i_pretty))
                up = prettyForm(*up.right(' '*i_pretty.width()))
            tensor = prettyForm(*tensor.right(' '*i_pretty.width()))
        return prettyForm(*prettyForm.stack(up, tensor, down))

class MetricTensorComponent(TensorComponent):
    #__init__(self, tangent_space, *indices):

    # XXX args aliases that can not be implemented in __init__ because
    # the rewrite rules do not call __init__
    name = 'g' # for printing

    @property
    def indices(self):
        return self.args[1:]
