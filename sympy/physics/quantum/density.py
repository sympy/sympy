from sympy import Tuple, Add, Matrix, log, expand
from sympy.printing.pretty.stringpict import prettyForm
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.operator import HermitianOperator, OuterProduct
from sympy.physics.quantum.represent import represent
from matrixutils import numpy_ndarray, scipy_sparse_matrix, to_numpy

class Density(HermitianOperator):
    """Density operator for representing mixed states."""

    def states(self):
        return Tuple(*[arg[0] for arg in self.args])

    def probs(self):
        return Tuple(*[arg[1] for arg in self.args])

    def get_state(self, index):
        state = self.args[index][0]
        return state

    def get_prob(self, index):
        prob = self.args[index][1]
        return prob

    def operate_on(self, op):
        new_args = [(op*state, prob) for (state, prob) in self.args]
        return Density(*new_args)

    def doit(self, **hints):
        terms = []
        for (state, prob) in self.args:
            terms.append(prob*state*Dagger(state))
        return Add(*terms)

    def _represent(self, **options):
        return represent(self.doit(), **options)

    def _print_operator_name_latex(self, printer, *args):
        return printer._print(r'\rho', *args)

    def _print_operator_name_pretty(self, printer, *args):
        return prettyForm(u"\u03C1")


def entropy(density):
    """Compute the entropy of a density matrix.
    
    This computes -Tr(density*ln(density)) using the eigenvalue decomposition
    of density, which is given as either a Density instance or a matrix
    (numpy.ndarray, sympy.Matrix or scipy.sparse).
    """
    if isinstance(density, Density):
        density = represent(density, format='numpy')
        
    if isinstance(density, scipy_sparse_matrix):
        density = to_numpy(density)

    if isinstance(density, Matrix):
        eigvals = density.eigenvals().keys()
        return expand(-sum(e*log(e) for e in eigvals))

    elif isinstance(density, numpy_ndarray):
        import numpy as np
        eigvals = np.linalg.eigvals(density)
        return -np.sum(eigvals*np.log(eigvals))

    else:
        raise ValueError("numpy.ndarray, scipy.sparse or sympy matrix expected")

