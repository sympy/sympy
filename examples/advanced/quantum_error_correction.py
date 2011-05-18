"""
Notes to self:
Pauli class in sympy/physics/paulialgebra.py does not use matrices for
their representation.  An operation such as Pauli(1)*Qubit(1) will
not be evaluated with apply_operators.  

For the matrix reps, look at file sympy/physics/matrices.py

Based on PHD thesis by Daniel Gottesman
It seems the generators are a set of operators used to diagonose an error,
i.e. they are used for a syndrome measurement.  Generators anticommunte with
errors, so if an error occurred, the error syndrome measurement would 
measure an eigenvalue of -1.  If an error did not occur, an error syndrome
measurement yields an eigenvalue of 1.  The generators themselves are
a set of errors since they use the error basis.  
"""

from sympy import Matrix, I
from sympy.physics.matrices import msigma
from sympy.physics.quantum.qubit import Qubit, qubit_to_matrix, matrix_to_qubit
from sympy.physics.quantum.tensorproduct import TensorProduct

# Matrix operations evaluate themselves, no need for apply_operators
# Playing around with 5 qubit stabilizer code

# Number of generators are (n-k), where I want to encode k bit with n bits
# Use <matrix instance>.eye(n) to return identity matrix

# Generators commute with each other
# The size of the stabilizers set is 2^(n-k), and each operator that fixes 
# an error is some unordered product of the generators, where each 
# operators may consist of a product of 1 to (n-k) generators
# In other words, given (n-k) generators, the operators are such that
# (n-k) C j where j = 1, ... , (n-k), (nCr is the choose operation)
# 2^(n-k) = sum from i = 0 to i = (n-k): (n-k) C i

# At least for now
# Create a function to create the different combinations of numbers
# Use a bitmap.  
# Ex. 0101 means combine generators 1 and 3 to form operator
# Ex. 1100 means combine generators 3 and 4 to form operator
# The function returns an array of tuples, these tuples could be 
# indices into an array of generators to form the operator

# This allows one to create the codeword, which is a sum of the generators
# operating on a n bit 0 (5 bits, M|00000>)

# Returns a list of tuples.  The tuples are the generators to combine as
# a product to form an operator in the stabilizer
def generator_combinations(numGenerators):
    # Start with an empty list
    combinationList = []
    totalOperators = pow(2, numGenerators)
  
    largestMask = 0x01 << numGenerators

    for i in range(totalOperators):
        bitmap = i
        mask = 0x01
        currentGenerator = 1
        currentCombo = ()

        while currentGenerator <= numGenerators:
            # If map was nonzero, then include this generator in the combo
            if mask & bitmap:
                currentCombo = currentCombo + (currentGenerator,)
            mask = mask << 1
            currentGenerator = currentGenerator + 1

        # Add combo to the list
        list.append(combinationList, currentCombo)

    return combinationList
        
identity = Matrix( ((1, 0), (0, 1)))
sigma_x = msigma(1)
sigma_y = msigma(2)
sigma_z = msigma(3)

m1 = TensorProduct(sigma_x, sigma_z, sigma_z, sigma_x, identity)
m2 = TensorProduct(identity, sigma_x, sigma_z, sigma_z, sigma_x)
m3 = TensorProduct(sigma_x, identity, sigma_x, sigma_z, sigma_z)
m4 = TensorProduct(sigma_z, sigma_x, identity, sigma_x, sigma_z)

generators = [m1, m2, m3, m4]

# Test the function generator_combinations
print generator_combinations(4)
# Ok it seems to work

combos = generator_combinations(4)

stabilizer_set = []

identity_5 = TensorProduct(identity, identity, identity, identity, identity)
for a_combo in combos:
    a_stabilizer = identity_5

    for which_gen in a_combo:
        a_stabilizer = a_stabilizer * generators[which_gen - 1]

    list.append(stabilizer_set, a_stabilizer)
#print 'Stabilizer set: ', stabilizer_set

zero_codeword = map((lambda op: 
                         matrix_to_qubit(
                             op * qubit_to_matrix(Qubit('00000')))),
                    stabilizer_set)

print ''
print '5 qubit codeword for 0:'
print zero_codeword

# Possible to generate new codes by "pasting" codes together

# Concept still confusing to me: how do you find the generators

print 'End of quantum error correction example.'
