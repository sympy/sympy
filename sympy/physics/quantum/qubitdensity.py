"""Qubit specific things for density matrices.

WARNING: This module may not work.  It used to work, but we have refactored
the density matrix classes and this needs to be updated as the density
matrix work is completed.
"""

def matrix_to_density(mat):
    """
       Works by finding the eigenvectors and eigenvalues of the matrix.
       We know we can decompose rho by doing:
           sum(EigenVal*|Eigenvect><Eigenvect|)
    """
    eigen = mat.eigenvects()
    args = [[matrix_to_qubit(Matrix([vector,])),x[0]] for x in eigen for vector in x[2] if x[0] != 0]
    return Density(*args)

def reduced_density(density, qubit, **options):
    """Compute the reduced density matrix by tracing out one qubit."""

    def find_index_that_is_projected(j, k, qubit):
         bit_mask = 2**qubit - 1
         return ((j >> qubit) << (1 + qubit)) + (j & bit_mask) + (k << qubit)

    old_density = represent(density, **options)
    old_size = old_density.cols
    new_size = old_size/2
    new_density = Matrix().zeros(new_size)
    for i in xrange(new_size):
        for j in xrange(new_size):
            for k in xrange(2):
                col = find_index_that_is_projected(j,k,qubit)
                row = find_index_that_is_projected(i,k,qubit)
                new_density[i,j] += old_density[row, col]          
    return Matrix(new_density)
    
def entropy_of_entanglement(density, qubit, **options):
    from sympy.physics.quantum.density import entropy
    rho = reduced_density(density, qubit, **options)
    return entropy(rho)

def test_matrix_to_density():
    assert matrix_to_density(represent(Density([Qubit('00'),1]), nqubits=2)) ==\
     Density([Qubit('00'),1])
    
    den = Density([Qubit('00'),1], [(Qubit('00')+Qubit('11'))/sqrt(2),1])
    assert represent(matrix_to_density(represent(den), nqubits=2))\
      == represent(den, nqubits=2)