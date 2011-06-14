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

# Generators commute with each other.
# The size of the stabilizer set is 2^(n-k), and each operator that fixes
# an error is some unordered product of the generators, where each
# operator may consist of a product of 1 to (n-k) generators.
# In other words, given (n-k) generators, the operators are such that
# (n-k) C j where j = 1, ... , (n-k), (nCr is the choose operation)
# 2^(n-k) = sum from i = 0 to i = (n-k): (n-k) C i

# A function to create different combinations of numbers (nCr).
# Used to determine which generators to combine to form an operator.
# Ex. 0101 means combine generators 1 and 3 to form an operator
# Ex. 1100 means combine generators 3 and 4 to form an operator
# The function returns a list of tuples.  These tuples are
# indices into an list/array of generators.
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

# 1 qubit to 5 qubit codeword
m1_5 = X(4)*Z(3)*Z(2)*X(1)
m2_5 = X(3)*Z(2)*Z(1)*X(0)
m3_5 = X(4)*X(2)*Z(1)*Z(0)
m4_5 = Z(4)*X(3)*X(1)*Z(0)

generators = [m1_5, m2_5, m3_5, m4_5]

combos = generator_combinations(generators.__len__())

stabilizer_op_set = []

for a_combo in combos:
    a_stabilizer = 1

    for which_gen in a_combo:
        a_stabilizer = a_stabilizer * generators[which_gen - 1]

    list.append(stabilizer_op_set, a_stabilizer)

zero_codeword = map((lambda op: apply_operators(op*Qubit('00000'))),
                    stabilizer_op_set)

zero_codeword_intqubit = map((lambda op: apply_operators(op*IntQubit(0, 5))),
                             stabilizer_op_set)

qubit_5_not = X(4)*X(3)*X(2)*X(1)*X(0)

one_codeword = map((lambda basis: apply_operators(qubit_5_not*basis)),
                   zero_codeword)

one_codeword_intqubit = map((lambda basis: apply_operators(qubit_5_not*basis)),
                            zero_codeword_intqubit)

print ''
print '5 qubit codeword example:'
print IntQubit(0), ' = ', reduce((lambda x, y: x + y), zero_codeword)
print ''
print IntQubit(0), ' = ', reduce((lambda x, y: x + y), zero_codeword_intqubit)
print ''
print IntQubit(1), ' = ', reduce((lambda x, y: x + y), one_codeword)
print ''
print IntQubit(1), ' = ', reduce((lambda x, y: x + y), one_codeword_intqubit)
print ''

# 1 qubit to 7 qubit codeword generators
m1_7 = X(6)*X(5)*X(4)*X(3)
m2_7 = X(6)*X(5)*X(2)*X(1)
m3_7 = X(6)*X(4)*X(2)*X(0)
m4_7 = Z(6)*Z(5)*Z(4)*Z(3)
m5_7 = Z(6)*Z(5)*Z(2)*Z(1)
m6_7 = Z(6)*Z(4)*Z(2)*Z(0)

# m4_7, m5_7, m6_7 are used to detect phase changes.  They have no effect
# on codewords with no errors in them, so applying m4_7, m5_7, and m6_7 will
# not produce a change in the codeword.
generators_7 = [m1_7, m2_7, m3_7]

combos_7 = generator_combinations(generators_7.__len__())

stabilizer_op_set_7 = []

for a_combo in combos_7:
    a_stabilizer = 1

    for which_gen in a_combo:
        a_stabilizer = a_stabilizer * generators_7[which_gen - 1]

    list.append(stabilizer_op_set_7, a_stabilizer)

zero_codeword_7 = map((lambda op: apply_operators(op*Qubit('0000000'))),
                      stabilizer_op_set_7)

qubit_7_not = X(6)*X(5)*X(4)*X(3)*X(2)*X(1)*X(0)

one_codeword_7 = map((lambda basis: apply_operators(qubit_7_not*basis)),
                     zero_codeword_7)

print '7 qubit codeword example:'
print IntQubit(0),' = ', reduce((lambda x, y: x + y), zero_codeword_7)
print ''
print IntQubit(1),' = ', reduce((lambda x, y: x + y), one_codeword_7)
print ''

# 3 to 8 qubit codeword
m1_3_8 = X(7)*X(6)*X(5)*X(4)*X(3)*X(2)*X(1)*X(0)
m2_3_8 = Z(7)*Z(6)*Z(5)*Z(4)*Z(3)*Z(2)*Z(1)*Z(0)
m3_3_8 = X(6)*X(4)*Y(3)*Z(2)*Y(1)*Z(0)
m4_3_8 = X(6)*Z(5)*Y(4)*X(2)*Z(1)*Y(0)
m5_3_8 = Y(6)*X(5)*Z(4)*X(3)*Z(2)*Y(0)

# These operators, together with the generators, generate the
# normalizer.  The stabilizers are a subset of the normalizer.
x1_3_8 = X(7)*X(6)*Z(2)*Z(0)
x2_3_8 = X(7)*X(5)*Z(4)*Z(1)
x3_3_8 = X(7)*Z(4)*X(3)*Z(2)
z1_3_8 = Z(6)*Z(4)*Z(2)*Z(0)
z2_3_8 = Z(5)*Z(4)*Z(1)*Z(0)
z3_3_8 = Z(3)*Z(2)*Z(1)*Z(0)

# For a 3 qubit codeword, there are 8 codewords (2^3) to encode
generators_3_8 = [m1_3_8, m2_3_8, m3_3_8, m4_3_8, m5_3_8]

combos_3_8 = generator_combinations(generators_3_8.__len__())

stabilizer_op_set_3_8 = []

for a_combo in combos_3_8:
    a_stabilizer = 1

    for which_gen in a_combo:
        a_stabilizer = a_stabilizer * generators_3_8[which_gen - 1]

    list.append(stabilizer_op_set_3_8, a_stabilizer)

zero_codeword_3_8 = map((lambda op: apply_operators(op*Qubit('00000000'))),
                        stabilizer_op_set_3_8)

all_codewords_3_8 = []
for number in range(pow(2, 3)):
    a_codeword = zero_codeword_3_8

    if number & 0x01:
        a_codeword = map((lambda st: apply_operators(x3_3_8 * st)),
                         a_codeword)
    if number & 0x02:
        a_codeword = map((lambda st: apply_operators(x2_3_8 * st)),
                         a_codeword)
    if number & 0x04:
        a_codeword = map((lambda st: apply_operators(x1_3_8 * st)),
                         a_codeword)

    list.append(all_codewords_3_8, (number, a_codeword))

print '3 to 8 qubit codeword example:'
for (a_number, a_codeword) in all_codewords_3_8:
    print IntQubit(a_number, 3),' = ', reduce((lambda x, y: x + y), a_codeword)
    print ''

print ''

# For a [4, 2, 2] code, 2 qubit words are converted to 4 qubit words.
# The 2 stabilizers are one with all X operators and one with all Z operators.
m1_2_4 = X(3)*X(2)*X(1)*X(0)
m2_2_4 = Z(3)*Z(2)*Z(1)*Z(0)

# The operators below help generate the normalizer set.
# According to page 88 of Daniel Gottesman PHD thesis, the normalizer for
# the [4, 2, 2] codeword can be corrected such that each normalizer has an
# even number of X, Y, and Z operators because it is has an even n (n = 4).
x1_2_4 = X(3)*X(2)
x2_2_4 = X(3)*X(1)
z1_2_4 = Z(2)*Z(0)
z2_2_4 = Z(1)*Z(0)

generators_2_4 = [m1_2_4, m2_2_4]

combos_2_4 = generator_combinations(generators_2_4.__len__())

stabilizer_op_set_2_4 = []

for a_combo in combos_2_4:
    a_stabilizer = 1

    for which_gen in a_combo:
        a_stabilizer = a_stabilizer * generators_2_4[which_gen - 1]

    list.append(stabilizer_op_set_2_4, a_stabilizer)

zero_codeword_2_4 = map((lambda op: apply_operators(op*Qubit('0000'))),
                        stabilizer_op_set_2_4)

all_codewords_2_4 = []
for number in range(pow(2, 2)):
    a_codeword = zero_codeword_2_4

    if number & 0x02:
        a_codeword = map((lambda st: apply_operators(x2_2_4 * st)),
                         a_codeword)
    if number & 0x01:
        a_codeword = map((lambda st: apply_operators(x1_2_4 * st)),
                         a_codeword)

    list.append(all_codewords_2_4, (number, a_codeword))

print '2 to 4 qubit codeword example:'
for (a_number, a_codeword) in all_codewords_2_4:
    print IntQubit(a_number, 2),' = ', reduce((lambda x, y: x + y), a_codeword)
    print ''

print ''

print 'End of quantum error correction example.'
