#OrthoCopy.py

import sys
from sympy import expand,Mul,Add,Symbol,S,Expr,Wild,Pow,diff,trigsimp,\
                  simplify,Matrix

from GA import *
from GAsympy import *

def bilinear_expand(expr):
    """
    If a sympy 'Expr' is of the form:

    expr = expr_0+expr_1*a_1+...+expr_n*a_n

    where all the a_j are noncommuting symbols in basis then

    (expr_0,...,expr_n) and (1,a_1,...,a_n) are returned.  Note that
    expr_j*a_j does not have to be of that form, but rather can be any
    Mul with a_j as a factor (it doen not have to be a postmultiplier).
    expr_0 is the scalar part of the expression.
    """

    if expr.is_commutative: #commutative expr only contains expr_0
        return((expr,),(ONE,))

    expr = expand(expr)
    if isinstance(expr,Mul): #expr only contains one term
        x = (coefs,bases) = expr.args_cnc()
        coefs = Mul(*coefs)
        bases = Mul(*bases)
    elif isinstance(expr,Symbol): #term is Symbol
        coefs = ONE
        bases = expr
    elif isinstance(expr,Add): #expr has multiple terms
        coefs = []
        bases = []
        for arg in expr.args:
            term = arg.args_cnc()
            coefs.append(Mul(*term[0]))
            bases.append(Mul(*term[1]))

    """
    if not isinstance(coefs,list): #convert single coef to list
        coefs = [coefs]
    if not isinstance(bases,list): #convert single base to list
        bases = [bases]
    """
    coefs = tuple(coefs)
    bases = tuple(bases)
    return(coefs,bases)

class OrthoCopy(object):
    copy_indices = []
    bases        = {}
    bases_inv    = {}

    @staticmethod
    def build_base_maps(copy_num):
        if copy_num in OrthoCopy.copy_indices:
            return
        OrthoCopy.copy_indices.append(copy_num)

        for grade in MV.blades[1:]:
            for base in grade:
                copy_base_str = ''
                for index in MV.blade_to_index[base]:
                    copy_base_str += str(MV.basis_vectors[index])+'__('+str(copy_num)+')^'
                copy_base = Symbol(copy_base_str[:-1],commutative=False)
                OrthoCopy.bases[copy_base] = (base,copy_num,len(MV.blade_to_index[base]))
                OrthoCopy.bases_inv[(copy_num,base)] = copy_base
        return

    @staticmethod
    def mul_common_bases(base1,base2): #Multipy two bases in same subspace
        (MV_base1,copy_num1) = OrthoCopy.bases[base1][0:1]
        (MV_base2,copy_num2) = OrthoCopy.bases[base2][0:1]
        MV1 = MV(MV_base1)
        MV2 = MV(MV_base2)
        obj = (MV1*MV2).obj
        (coefs,bases) = linear_expand(obj)
        new_obj = 0
        for (coef,base) in zip(coefs,bases):
            new_obj += coef*OrthoCopy.bases_inv[(copy_num1,base]
        return(new_obj)

    def __init__(self,mv,copy_num=1):
        if copy_num == 0:
            self.obj = mv
        else:
            OrthoCopy.build_base_maps(copy_num)
            obj = mv.obj
            obj = expand(obj)
            (coefs,bases) = linear_expand(obj)
            self.obj = ZERO
            for (coef,base) in zip(coefs,bases):
                if base == MV.ONE:
                    self.obj += coef*base
                else:
                    self.obj += coef*OrthoCopy.bases_inv[(copy_num,base)]
            self.obj = self.obj.subs({MV.ONE:1})

    def __str__(self):
        obj = expand(self.obj)
        (coefs,bases) = bilinear_expand(obj)
        rezip = zip(bases,coefs)
        rezip.sort()
        obj_str = ''
        for (base,coef) in rezip:
            term_str = str(coef*base)
            if term_str[0] != '-':
                term_str = '+'+term_str
            obj_str += term_str
        if obj_str[0] == '+':
            obj_str = obj_str[1:]
        return(obj_str)

    def __mul__(self,OC):
        prod = expand(self.obj*OC.obj)
        (coefs,bases) = bilinear_expand(prod)
        new_prod = 0
        for (coef,base) in zip(coefs,bases):
            base_lst = base.args_cnc()[1]
            copy_lst = []
            grade_lst = []
            for base_base in base_lst:
                if base_base != []:
                    (copy_num,grade) = OrthoCopy.bases[base_base][1:]
                    copy_lst.append(copy_num)
                    grade_lst.append(grade)
            if copy_lst == sorted(copy_lst):
                new_prod += coef*base
            else:
                print 'unsorted copy_lst =',copy_lst
        return(OrthoCopy(new_prod,0))

    def __add__(self,OC):
        return(OrthoCopy(self.obj+OC.obj))

    def __sub__(self,OC):
        return(OrthoCopy(self.obj-OC.obj))

    def __radd__(self,a):
        return(OrthoCopy(a+self.obj)

    def __rsub__(self,a):
        return(OrthoCopy(a-self.obj)



(ex,ey) = MV.setup('e_x e_y','[1,1]')

b = MV('b','bivector')
s = 1+b

S1 = OrthoCopy(s,1)
S2 = OrthoCopy(s,2)
S3 = OrthoCopy(s,3)

P = S1*S2*S3
print 'P =',P
