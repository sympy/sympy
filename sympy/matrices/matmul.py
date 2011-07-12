from matexpr import MatrixExpr, ShapeError, matrixify, Identity
from sympy.core import Mul

class MatMul(MatrixExpr, Mul):

    def __new__(cls, *args):

        # Check that the shape of the args is consistent
        matrices = [arg for arg in args if arg.is_Matrix]

        for i in range(len(matrices)-1):
            A,B = matrices[i:i+2]
            if A.m != B.n:
                raise ShapeError("Matrices %s and %s are not aligned"%(A, B))

        expr = matrixify(Mul.__new__(cls, *args))
        nonmatrices = [arg for arg in expr.args if not arg.is_Matrix]
        matrices = [arg for arg in expr.args if arg.is_Matrix]

        old = []
        new = matrices
        while(new != old):
            passnext=False
            old = new
            new = []
            # Flattening logic
            for i, mat in enumerate(old):
                if passnext:
                    passnext=False
                    continue
                if (mat.is_Inverse): # Cancel inverses if easy
                    if(i+1<len(old) and mat.arg == old[i+1]):
                        new.append(Identity(mat.n)) # add Id
                        passnext = True
                    elif i>0 and mat.arg == old[i-1]:
                        new.pop() # remove i-1 from the new list
                        new.append(Identity(mat.n)) # add Id
                    else:
                        new.append(mat) # Pass through

                else:
                    new.append(mat) # Pass through
            if all(mat.is_Identity for mat in new):
                new = [new[0]]
            else:
                new = [mat for mat in new if not mat.is_Identity] # clear ident

        return matrixify(Mul.__new__(cls, *(nonmatrices+new)))

    @property
    def shape(self):
        matrices = [arg for arg in self.args if arg.is_Matrix]
        return (matrices[0].n, matrices[-1].m)

    def _check_shape(self):
        matrices = [arg for arg in self.args if arg.is_Matrix]
        for A, B in zip(matrices[:-1], matrices[1:]):
            if A.m != B.n:
                return False
        return all(mat._check_shape() for mat in matrices)

from inverse import Inverse
