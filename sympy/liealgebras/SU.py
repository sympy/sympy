#copy of file su.py, name of which was changed

from sympy.matrices import zeros, ones, Matrix, randMatrix, matrices
from sympy import KroneckerDelta, sqrt, nsimplify
import random

def su_generator(n, index_of_generator):
    ind = index_of_generator

    def lambda1(i, j, mu, nu):
        return  KroneckerDelta(j,mu)*KroneckerDelta(i,nu)+KroneckerDelta(j,nu)*KroneckerDelta(i,mu)

    def lambda2(i, j, mu, nu):
        return (-1j)*(KroneckerDelta(j,mu)*KroneckerDelta(i,nu)-KroneckerDelta(j,nu)*KroneckerDelta(i,mu))

    j = ind%n
    i = int((ind-j)/n)+1
    j = j+1
    #generators from equation 4 [2]
    if (j < i):
        answer = Matrix([lambda1(i, j, 1, 1)])
        for nu in range(2, n+1):
            answer = answer.col_join(Matrix([lambda1(i, j, 1, nu)]))
        for mu in range(2, n+1):
            row = [lambda1(i, j, mu, 1)]
            for nu in range(2, n+1):
                row.append(lambda1(i, j, mu, nu))
            answer = answer.row_join(Matrix(row))
        return answer
        print(answer)
    if (j > i):
        answer = Matrix([lambda2(i, j, 1, 1)])
        for nu in range(2, n+1):
            answer = answer.col_join(Matrix([lambda2(i, j, 1, nu)]))
        for mu in range(2, n+1):
            row = [lambda2(i, j, mu, 1)]
            for nu in range(2, n+1):
                row.append(lambda2(i, j, mu, nu))
            answer = answer.row_join(Matrix(row))
        return answer
        print(answer)
    #generators from equation 5 [2]
    if (i == j):
        k = i+1
        coefficient = nsimplify(sqrt(2/(k*k-k)))
        row = [1*coefficient]
        for l in range(1, n):
            row.append(0)
        answer = Matrix(row)
        for l in range(1, n):
            row = []
            for m in range(0, n):
                if (l == m and l < k-1):
                    row.append(1*coefficient)
                elif (l == m and l == k-1):
                    row.append(-(k-1)*coefficient)
                else:
                    row.append(0)
            answer = answer.row_join(Matrix(row))
        return answer

#Example no. 1
#building of generators for su(3)
for i in range(0, 8):
    print(su_generator(3,i))

for i in range(0, 3):
    print(su_generator(2,i))

def su_generator_matrix(n):
    answer = Matrix([su_generator(n, 0)[0, 0]])
    row = []
    for k in range(0, n*n-1):
        row.append(su_generator(n, k)[0, 0])
    answer = Matrix(row)
    answer = answer.transpose()
    for i in range(0, n):
        for j in range(0, n):
            if (i+j != 0):
                row = []
                for k in range(0, n*n-1):
                    row.append(su_generator(n, k)[i, j])
                answer = answer.col_join(Matrix(row).transpose())
    return answer

print(su_generator_matrix(3))

#random su(n) element constructor
#can be used for testing is_su(n) function
def random_su_element(n):
    ans = zeros(n,n)
    for i in range(n*n-1):
        rand_coefficient = random.random()
        ans = ans+rand_coefficient*su_generator(n, i)
    return ans

#su(n) generator with coefficients for generators
def su_with_coefficients(n, list_of_coefficients):
    ans = zeros(n, n)
    for i in range(n*n-1):
        ans = ans+list_of_coefficients[i]*su_generator(n, i)
    return ans

#Example no.2
#Random element of su(n)
random_su_element(2)

#Example no.3
#su(2) with coefficient for every generation matrix
su_with_coefficients(2,[1,2,3])

#function which takes matrix of the group and gives answer, is it su(n) or not
def is_su(given_matrix, matrix_size):
    gen_matrix = su_generator_matrix(matrix_size)
    test_row = gen_matrix.row(-1)
    gen_matrix.row_del(-1)
    eq_answer = []
    for i in range(matrix_size):
        for j in range(matrix_size):
            eq_answer.append(given_matrix[i, j])
    equation_ans = Matrix(eq_answer[:matrix_size*matrix_size-1])
    equation_answer = equation_ans.col(0)
    solution = gen_matrix.LUsolve(equation_answer)
    calc_ans = test_row*solution
    if (calc_ans[0] == eq_answer[matrix_size*matrix_size-1]):
        return 0
    else:
        return 1

#Example no.4
#We will look if a random_su_element(3) is an element of SU(3)
is_su(random_su_element(3), 3)

#links:
#https://en.wikipedia.org/wiki/Special_unitary_group[1]
#https://arxiv.org/pdf/math-ph/0205016.pdf
