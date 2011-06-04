"""
Notes to self:
Based on PHD thesis by Daniel Gottesman
It seems the generators are a set of operators used to diagonose an error,
i.e. they are used for a syndrome measurement.  Generators anticommunte with
errors, so if an error occurred, the error syndrome measurement would 
measure an eigenvalue of -1.  If an error did not occur, an error syndrome
measurement yields an eigenvalue of 1.  The generators themselves are
a set of errors since they use the error basis.  
"""

from sympy.physics.quantum.qubit import Qubit, IntQubit
from sympy.physics.quantum.applyops import apply_operators
from sympy.physics.quantum.gate import X, Y, Z

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

m1_5 = X(3)*Z(2)*Z(1)*X(0)
m2_5 = X(4)*Z(3)*Z(2)*X(1)
m3_5 = Z(4)*Z(3)*X(2)*X(0)
m4_5 = Z(4)*X(3)*X(1)*Z(0)

# Difference between Mermin book and PHD thesis
# Mermin book has Z in place of X and X in place of Z operators
'''
m1 = Z(4)*X(3)*X(2)*Z(1); m1
m2 = Z(3)*X(2)*X(1)*Z(0); m2
m3 = Z(4)*Z(2)*X(1)*X(0); m3
m4 = X(4)*Z(3)*Z(1)*X(0); m4
'''

generators = [m1_5, m2_5, m3_5, m4_5]

# Test the function generator_combinations
#print generator_combinations(4)
# Ok it seems to work

combos = generator_combinations(4)

stabilizer_set = []

for a_combo in combos:
    a_stabilizer = 1

    for which_gen in a_combo:
        a_stabilizer = a_stabilizer * generators[which_gen - 1]

    list.append(stabilizer_set, a_stabilizer)
#print 'Stabilizer set: ', stabilizer_set

zero_codeword = map((lambda op: apply_operators(op*Qubit('00000'))),
                    stabilizer_set)

# IntQubit parameters: # to be represented, # of qubits to represent number
zero_codeword_intqubit = map((lambda op: apply_operators(op*IntQubit(0, 5))),
                             stabilizer_set)

qubit_5_not = X(4)*X(3)*X(2)*X(1)*X(0)

one_codeword = map((lambda basis: apply_operators(qubit_5_not*basis)),
                   zero_codeword)

one_codeword_intqubit = map((lambda basis: apply_operators(qubit_5_not*basis)), 
                            zero_codeword_intqubit)

#one_codeword_intqubit = [IntQubit(one_codeword[0]),
# IntQubit(one_codeword[1])]

print ''
print '5 qubit codeword for 0:'
print zero_codeword
print ''
print zero_codeword_intqubit
print ''
print '5 qubit codeword for 1:'
print one_codeword
print ''
print one_codeword_intqubit
print ''

# 7 qubit codeword generators
m1_7 = X(3)*X(2)*X(1)*X(0)
m2_7 = X(5)*X(4)*X(1)*X(0)
m3_7 = X(6)*X(4)*X(2)*X(0)
m4_7 = Z(3)*Z(2)*Z(1)*Z(0)
m5_7 = Z(5)*Z(4)*Z(1)*Z(0)
m6_7 = Z(6)*Z(4)*Z(2)*Z(0)

# m4_7, m5_7, m6_7 are used to detect phase changes.  They have no effect
# on codewords with no errors in them, so applying m4_7, m5_7, and m6_7 will
# not produce a change in the codeword; differences are observed only when an
# error occurs.
generators_7 = [m1_7, m2_7, m3_7]

combos_7 = generator_combinations(3)

stabilizer_set_7 = []

for a_combo in combos_7:
    a_stabilizer = 1

    for which_gen in a_combo:
        a_stabilizer = a_stabilizer * generators_7[which_gen - 1]

    list.append(stabilizer_set_7, a_stabilizer)

zero_codeword_7 = map((lambda op: apply_operators(op*Qubit('0000000'))),
                      stabilizer_set_7)
qubit_7_not = X(6)*X(5)*X(4)*X(3)*X(2)*X(1)*X(0)
one_codeword_7 = map((lambda basis: apply_operators(qubit_7_not*basis)),
                     zero_codeword_7)

print '7 qubit codeword example:'
print '|0> = ', zero_codeword_7
print '|1> = ', one_codeword_7

print ''

# 3 qubit to 8 qubit codeword
m1_3_8 = X(7)X(6)X(5)X(4)X(3)X(2)X(1)X(0)
m2_3_8 = Z(7)Z(6)Z(5)Z(4)Z(3)Z(2)Z(1)Z(0)
m3_3_8 = Z(7)Y(6)Z(5)Y(4)X(3)X(1)
m4_3_8 = Y(7)Z(6)X(5)Y(3)Z(2)X(1)
m5_3_8 = Y(7)Z(5)X(4)Z(3)X(2)Y(1)

# For a 3 qubit codeword, there are 8 codewardes (2^8) to encode

# Possible to generate new codes by "pasting" codes together

# Concept still confusing to me: how do you find the generators

print 'End of quantum error correction example.'
