from sympy.matrices import zeros, ones, Matrix, randMatrix, matrices

#generator of SU(n) group
#now works only for SU(2), but since all the other functions can be
#used for every SU(n) group, i will have to write only a new generator
#here

def SU_generator(n, request):
    generator = []
    list_of_generators = []
    if (n == 2):
        generator.append(Matrix([[0, 1j], [1j, 0]]))
        generator.append(Matrix([[0, -1j], [1j, 0]]))
        generator.append(Matrix([[1j, 0], [0, -1j]]))
    gen_matrix = Matrix([[0, 0, 1j], [1j, -1j, 0], [1j, 1j, 0], [0, 0, -1j]])
    if (request == 'generator'):
        return generator
    if (request == 'gen_matrix'):
        return gen_matrix

#simple SU(n) group constructor
def SU(n):
    generator = []
    generator = SU_generator(n, 'generator')
    ans = zeros(n,n)
    for i in range(n*n-1):
        ans = ans+generator[i]
    return ans

#SU(n) group generator with coefficients for group generator
def SU_with_coefficients(n, list_of_coefficients):
    generator = []
    generator = SU_generator(n, 'generator')
    ans = zeros(n, n)
    for i in range(n*n-1):
        ans = ans+list_of_coefficients[i]*generator[i]
    return ans

#Example no.1
#Simple SU(2) group
SU(2)

#Example no.2
#SU(2) group with coefficient for every generation matrix
SU_with_coefficients(2,[1,2,3])

#function which takes matrix of the group and gives answer, is it SU(n)
#or not
def is_SU(given_matrix, matrix_size):
    gen_matrix = SU_generator(matrix_size,'gen_matrix')
    test_row = gen_matrix.row(-1)
    gen_matrix.row_del(-1)
    eq_answer = []
    for i in range(matrix_size):
        for j in range(matrix_size):
            eq_answer.append(given_matrix[i][j])
    equation_ans = Matrix(eq_answer[:3])
    equation_answer = equation_ans.col(0)
    solution = gen_matrix.LUsolve(equation_answer)
    calc_ans = test_row*solution
    if (calc_ans[0] == eq_answer[matrix_size*matrix_size-1]):
        return 0
    else:
        return 1

#Example no.3
#We will look if matrix[[0,1],[1,0]] is matrix of  SU(2)
is_SU([[0,1],[1,0]],2)

#links:
#https://en.wikipedia.org/wiki/Special_unitary_group
