import sympy as sp
from sympy.polys.subresultants_qq_zz import sylvester

def rearrange_rows(matrix):
    num_rows = matrix.shape[0]
    def leading_zeros_count(row):
        count = 0
        for elem in row:
            if elem == 0:
                count += 1
            else:
                break
        return count

    for i in range(num_rows):
        target_zeros = leading_zeros_count(matrix.row(i))
        for j in range(i + 1, num_rows):
            current_zeros = leading_zeros_count(matrix.row(j))
            if current_zeros == target_zeros:
                matrix.row_swap(j, (i + 2) if i + 2 < num_rows else j)
    return matrix

def routh_hurwitz_matrix(poly, var):
    real = sp.simplify((poly.subs(var, sp.I * var) + poly.subs(var, -sp.I * var)) / 2)
    imag = sp.simplify((poly.subs(var, -sp.I * var) - poly.subs(var, sp.I * var)) / 2 * sp.I)
    a = rearrange_rows(sylvester(real, imag, var))
    if a.shape[0] == 3:
        a.row_swap(0, 2)
    minors = [a.minor(i, i) for i in range(a.shape[0])]
    unique_inequalities = unique_inequalities = {sp.simplify(minor > 0) for minor in minors}
    return unique_inequalities
