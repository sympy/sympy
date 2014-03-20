from __future__ import print_function, division

from sympy import Basic, Expr, S, Q
from .matexpr import ShapeError
from  .matexpr import MatrixExpr

class Determinant(Expr):
    """Matrix Determinant

    Represents the determinant of a matrix expression.

    >>> from sympy import MatrixSymbol, Determinant, eye
    >>> A = MatrixSymbol('A', 3, 3)
    >>> Determinant(A)
    Determinant(A)

    >>> Determinant(eye(3)).doit()
    1
    """

    def __new__(cls, mat):
        if not mat.is_Matrix:
            raise TypeError("Input to Determinant, %s, not a matrix" % str(mat))

        if not mat.is_square:
            raise ShapeError("Det of a non-square matrix")

        n=mat.rows;

        matrix=[[0 for x in range(n)] for x in range(n)]

        for i in range(0,n):
            for j in range(0,n):
                matrix[i][j]=str(mat) +"["+str(i)+"]["+str(j)+"]";

        return Determinant.det_sub(matrix);

    @classmethod
    def det_sub(self,matrix):
        n=len(matrix);
        det_mat="";

        if(n==1):
            det_mat=matrix[0][0];
        elif(n==2):
            det_mat=matrix[0][0]+"*" + matrix[1][1]+" - " + matrix[0][1]+"*" + matrix[1][0];
        else:
            sign="+";

            for i in range(0,n):
                if(i==0):
                    det_mat=matrix[0][0]+"("+Determinant.det_sub(Determinant.getM(matrix,0,0))+")";
                else:
                    if(i%2==0):
                        sign="+";
                    else:
                        sign="-";
                    det_mat=det_mat+sign+matrix[0][i]+"("+Determinant.det_sub(Determinant.getM(matrix,0,i))+")";

        return det_mat;

    @classmethod
    def getM(self,matrix,p,q):
        n=len(matrix);

        matrixA=[[0 for x in range(n-1)] for x in range(n-1)]

        countR=0;
        countC=0;

        for i in range(1,n):
                countC=0;
                for j in range(0,n):
                    if(j==q):
                        continue;
                    matrixA[countR][countC]=matrix[i][j];
                    countC+=1;
                countR+=1;

        return matrixA;



    @property
    def arg(self):
        return self.args[0]

    def doit(self, expand=False):
        try:
            return self.arg._eval_determinant()
        except (AttributeError, NotImplementedError):
            return self

def det(matexpr):
    """ Matrix Determinant

    >>> from sympy import MatrixSymbol, det, eye
    >>> A = MatrixSymbol('A', 3, 3)
    >>> det(A)
    Determinant(A)

    >>> det(eye(3))
    1
    """

    return Determinant(matexpr).doit()
