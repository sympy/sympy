from sympy import Matrix
from sympy.matrices import eigen

a = Matrix([[1, 2, 3], [3, 2, 1], [3, 1, 2]])
print(eigen._eigenvals(a))
