#linear_transformations

import sys
import inspect
import types
import itertools
import metric
from sympy import collect, expand, symbols, Matrix, Symbol, Function, S, Add
from copy import copy
import printer
import mv
import lt

def aprint(a):
    out = ''
    for ai in a:
        out += str(ai)+','
    print '['+out[:-1]+']'
    return


class Lt(object):

    mat_fmt = False
    init_slots = {'ga': (None, 'Name of metric (geometric algebra)'),
                  'f': (False, 'True if Lt if function of coordinates.')}

    @staticmethod
    def setup(ga):
        #coords = [Symbol('mu_' + str(x)) for x in ga.coords]
        coords = ga.coords
        x = sum([coords[i] * ga.basis[i] for i in ga.n_range])
        return coords, x

    @staticmethod
    def format(mat_fmt=False):
        Lt.mat_fmt = mat_fmt
        return

    def __init__(self, *kargs, **kwargs):

        kwargs = metric.test_init_slots(Lt.init_slots, **kwargs)

        mat_rep = kargs[0]
        ga = kwargs['ga']
        self.fct_flg = kwargs['f']

        self.Ga = ga
        self.coords = ga.lt_coords
        self.X = ga.lt_x

        self.spinor = False
        if isinstance(mat_rep, dict):  # Dict input
            self.mv_dict = {}
            self.lt_dict = {}
            for key in mat_rep:
                self.lt_dict[key] = mat_rep[key]
        elif isinstance(mat_rep, list):  # List input
            self.lt_dict = {}
            self.mv_dict = {}
            if not isinstance(mat_rep[0], list):
                for (lt_i, base) in zip(mat_rep, self.Ga.basis):
                    self.lt_dict[base] = lt_i
            else:
                mat_rep = map(list, zip(*mat_rep))  # Transpose list of lists
                for (row, base1) in zip(mat_rep, self.Ga.basis):
                    tmp = 0
                    for (col, base2) in zip(row, self.Ga.basis):
                        tmp += col * base2
                    self.lt_dict[base1] = tmp
        elif isinstance(mat_rep, mv.Mv):  # Spinor input
            self.spinor = True
            self.R = mat_rep
            self.Rrev = mat_rep.rev()
        elif isinstance(mat_rep, str):  # String input
            root = mat_rep
            mat_rep = []
            for row in self.Ga.coords:
                row_str = ''
                for col in self.Ga.coords:
                    row_str += root + '_' + str(row) + str(col) + ','
                if self.fct_flg:
                    fct_names = row_str[:-1].split(',')
                    mat_rep.append([Function(fct_name)(*self.Ga.coords) for fct_name in fct_names])
                else:
                    mat_rep.append(list(symbols(row_str[:-1])))
            self.__init__(mat_rep, ga=self.Ga)
        else:  # F is a multivector function to be tested for linearity
            F = mat_rep
            a = mv.Mv('a', 'vector', ga=self.Ga)
            b = mv.Mv('b', 'vector', ga=self.Ga)
            if F(a + b) == F(a) + F(b):
                self.lt_dict = {}
                for base in self.Ga.basis:
                    self.lt_dict[base] = (F(mv.Mv(base, ga=self.Ga))).obj
            else:
                raise ValueError(str(mat_rep) + ' is not supported for Lt definition\n')

    def __call__(self, v, obj=False):

        if isinstance(v, mv.Mv) and self.Ga.name != v.Ga.name:
                raise ValueError('In Lt call Lt and argument refer to different vector spaces')

        if self.spinor:
            if not isinstance(v, mv.Mv):
                v = mv.Mv(v, ga=self.Ga)
            R_v_Rrev = self.R * v * self.Rrev
            if obj:
                return R_v_Rrev.obj
            else:
                return R_v_Rrev

        if isinstance(v, mv.Mv):
            if v.is_vector():
                lt_v = v.obj.xreplace(self.lt_dict)
                if obj:
                    return lt_v
                else:
                    return mv.Mv(lt_v, ga=self.Ga)
            else:
                mv_obj = v.obj
        else:
            mv_obj = mv.Mv(v, ga=self.Ga).obj

        if len(self.mv_dict) == 0:  # Build dict for linear transformation of multivector
            self.mv_dict = copy(self.lt_dict)
            for key in self.Ga.blades[2:]:
                for blade in key:
                    index = self.Ga.blades_to_indexes_dict[blade]
                    lt_blade = self(self.Ga.basis[index[0]], obj=True)
                    for i in index[1:]:
                        lt_blade = self.Ga.wedge(lt_blade, self(self.Ga.basis[i], obj=True))
                    self.mv_dict[blade] = lt_blade

        lt_v = mv_obj.xreplace(self.mv_dict)
        if obj:
            return lt_v
        else:
            return mv.Mv(lt_v, ga=self.Ga)

    def __add__(self, LT):

        if self.Ga.name != LT.Ga.name:
            raise ValueError("Attempting addition of Lt's from different geometric algebras")

        self_add_LT = copy(self.lt_dict)
        for key in LT.lt_dict.keys():
            if key in self_add_LT:
                self_add_LT[key] = collect(self_add_LT[key] + LT.lt_dict[key], self.Ga.basis)
            else:
                self_add_LT[key] = LT.lt_dict[key]
        return(Lt(self_add_LT, ga=self.Ga))

    def __sub__(self, LT):

        if self.Ga.name != LT.Ga.name:
            raise ValueError("Attempting subtraction of Lt's from different geometric algebras")

        self_add_LT = copy(self.lt_dict)
        for key in LT.lt_dict.keys():
            if key in self_add_LT:
                self_add_LT[key] = collect(self_add_LT[key] - LT.lt_dict[key], self.Ga.basis)
            else:
                self_add_LT[key] = -LT.lt_dict[key]
        return(Lt(self_add_LT, ga=self.Ga))

    def __mul__(self, LT):

        if isinstance(LT, Lt):

            if self.Ga.name != LT.Ga.name:
                raise ValueError("Attempting multiplication of Lt's from different geometric algebras")
            self_mul_LT = {}
            for base in LT.lt_dict:
                self_mul_LT[base] = self(LT(base, obj=True), obj=True)
            for key in self_mul_LT:
                self_mul_LT[key] = expand(self_mul_LT[key]).collect(self.Ga.basis)
            return(Lt(self_mul_LT, ga=self.Ga))
        else:
            self_mul_LT = {}
            for key in self.lt_dict:
                self_mul_LT[key] = LT * self.lt_dict[key]
            return(Lt(self_mul_LT, ga=self.Ga))

    def __rmul__(self, LT):

        if not isinstance(LT, Lt):
            self_mul_LT = {}
            for key in self.lt_dict:
                self_mul_LT[key] = LT * self.lt_dict[key]
            return(Lt(self_mul_LT, ga=self.Ga))
        else:
            raise TypeError('Cannot have LT as left argument in Lt __rmul__\n')

    def _repr_latex_(self):
        latex_str = printer.GaLatexPrinter.latex(self)
        if r'\begin{align*}' not in latex_str:
            latex_str = r'\begin{equation*} ' + latex_str + r' \end{equation*}'
        return latex_str

    def det(self):  # det(L) defined by L(I) = det(L)I

        lt_I = self(self.Ga.i, obj=True)
        det_lt_I = lt_I.subs(self.Ga.i.obj, S(1))
        return(det_lt_I)

    def tr(self):  # tr(L) defined by tr(L) = grad|L(x)

        connect_flg = self.Ga.connect_flg
        self.Ga.connect_flg = False

        F_x = mv.Mv(self(self.Ga.lt_x, obj=True), ga=self.Ga)
        tr_F = (self.Ga.grad | F_x).scalar()
        self.Ga.connect_flg = connect_flg
        return(tr_F)

    def adj(self):

        self_adj = []
        self.Ga.dot_mode = '|'
        for e_j in self.Ga.basis:
            s = S(0)
            for (e_i, er_i) in zip(self.Ga.basis, self.Ga.r_basis):
                s += er_i * self.Ga.dot(e_j, self(e_i, obj=True))
            if self.Ga.is_ortho:
                self_adj.append(expand(s))
            else:
                self_adj.append(expand(s) / self.Ga.inorm)
        return Lt(self_adj, ga=self.Ga)

    def Lt_str(self):

        if self.spinor:
            return 'R = ' + str(self.R)
        else:
            pre = 'Lt('
            s = ''
            for base in self.Ga.basis:
                if base in self.lt_dict:
                    s += pre + str(base) + ') = ' + str(mv.Mv(self.lt_dict[base], ga=self.Ga)) + '\n'
                else:
                    s += pre + str(base) + ') = 0\n'
            return s[:-1]

    def Lt_latex_str(self):

        if self.spinor:
            #return str(self.R)
            s = '\\left \\{ \\begin{array}{ll} '
            for base in self.Ga.basis:
                s += 'L \\left ( ' + str(base) + '\\right ) =& ' + str(self.R * mv.Mv(base, ga=self.Ga) * self.Rrev) + ' \\\\ '
            s = s[:-3] + ' \\end{array} \\right \\} \n'
            return s
        else:
            s = '\\left \\{ \\begin{array}{ll} '
            for base in self.Ga.basis:
                if base in self.lt_dict:
                    s += 'L \\left ( ' + str(base) + '\\right ) =& ' + str(mv.Mv(self.lt_dict[base], ga=self.Ga)) + ' \\\\ '
                else:
                    s += 'L \\left ( ' + str(base) + '\\right ) =& 0 \\\\ '
            s = s[:-3] + ' \\end{array} \\right \\} \n'
            return s

    def Fmt(self, fmt=1, title=None):

        self.fmt = fmt
        latex_str = printer.GaLatexPrinter.latex(self)

        if printer.GaLatexPrinter.ipy:
            if title is None:
                if r'\begin{align*}' not in latex_str:
                    latex_str = r'\begin{equation*} ' + latex_str + r' \end{equation*}'
            else:
                if r'\begin{align*}' not in latex_str:
                    latex_str = r'\begin{equation*} ' + title + ' = ' + latex_str + r' \end{equation*}'
                else:
                    latex_str = latex_str.replace(r'\begin{align*}', r'\begin{align*} ' + title)
                    latex_str = latex_str.replace('&', '=&', 1)
            latex_str = latex_str.replace(r'\bfrac', r'\frac')
            from IPython.core.display import display, Math
            display(Math(latex_str))
        else:
            if title is not None:
                print title + ' = ' + latex_str
            else:
                    print latex_str
        return

    def __str__(self):
        if printer.GaLatexPrinter.latex_flg:
            Printer = printer.GaLatexPrinter
        else:
            Printer = printer.GaPrinter
        return Printer().doprint(self)

    def matrix(self):

        if self.spinor:
            self.lt_dict = {}
            for base in self.Ga.basis:
                self.lt_dict[base] = self(base).simplify()
            self.spinor = False
            mat = self.matrix()
            self.spinor = True
            return(mat)
        else:
            mat_rep = []
            for base in self.Ga.basis:
                if base in self.lt_dict:
                    row = []
                    image = (self.lt_dict[base])
                    if isinstance(image, mv.Mv):
                        image = image.obj
                    (coefs, bases) = metric.linear_expand(image)
                    for base in self.Ga.basis:
                        try:
                            i = bases.index(base)
                            row.append(coefs[i])
                        except:
                            row.append(0)
                    mat_rep.append(row)
                else:
                    mat_rep.append(self.Ga.n * [0])
            return(Matrix(mat_rep).transpose())


class Mlt(object):
    """
    A multilinear transformation (mlt) is a multilinear multivector function of
    a list of vectors (*kargs) F(v_1,...,v_r) where for any argument slot
    j we have (a is a scalar and u_j a vector)
          F(v_1,...,a*v_j,...,v_r) = a*F(v_1,...,v_j,...,v_r)
          F(v_1,...,v_j+u_j,...,v_r) = F(v_1,...,v_j,...,v_r) + F(v_1,...,u_j,...,v_r).
    If F and G are two mlt's with the same number of argument slots then the sum is
          (F+G)F(v_1,...,v_r) = F(v_1,...,v_r) + G(v_1,...,v_r).
    If F and G two mlt's with r and s argument slots then their product is
          (F*G)(v_1,...,v_r,...,v_(r+s)) = F(v_1,...,v_r)*G(v_(r+1),...,v_(r+s)),
    where * is any of the multivector multiplicative operations.
    The derivative of a mlt with is defined as the directional derivative with respect
    with respect to the coordinate vector (we assume F is implicitely a function of the
    coordinates)
        F(v_1,...,v_r;v_(r+1)) = (v_(r+1)|grad)F(v_1,...,v_j,...,v_r).
    The contraction of a mlt between slots j and k is defined as the
    geometric derivative of F with respect to slot k and the inner geometric
    derivative with respect to slot j (this gives the standard tensor
    definition of contraction for the case that F is a scalar function)
        Contract(i,j,F) = grad_i|(grad_j*F(v_1,...,v_i,...,v_j,...,v_r))
                        = grad_j|(grad_i*F(v_1,...,v_i,...,v_j,...,v_r)).
    This returns a mlt with slot i and j removed.
    """

    @staticmethod
    def subs(Ga, anew):
        #  Generate coefficient substitution list for new Mlt slot
        #  vectors (arguments) where anew is a list of slot vectors
        #  to be substituted for the old slot vectors.
        #  This is used when one wishes to substitute specific vector
        #  values into the Mlt such as the basis/reciprocal basis vectors.
        sub_lst = []
        i = 0
        for a in anew:
            acoefs = a.get_coefs(1)
            sub_lst += zip(Ga.pdiffs[i], acoefs)
            i += 1
        return sub_lst

    @staticmethod
    def increment_slots(nargs, Ga):
        # Increment cache of available slots (vector variables) if needed for Mlt class
        n_a = len(Ga.a)
        if n_a < nargs:
            for i in range(n_a, nargs):
                #  New slot variable with coefficients a_{n_a}__k
                a = Ga.mv('a_' + str(i + 1), 'vector')
                #  Append new slot variable a_j
                Ga.a.append(a)
                #  Append slot variable coefficients a_j__k for purpose
                #  of differentiation
                coefs = a.get_coefs(1)
                Ga.pdiffs.append(coefs)
                Ga.acoefs += coefs
        return

    @staticmethod
    def extact_basis_indexes(Ga):
        base_indexes = []
        for base in Ga.basis:
            base_str = str(base)
            base_str = base_str.replace(r'\boldsymbol','')
            base_str = base_str.replace('{','')
            base_str = base_str.replace('}','')
            i = base_str.find('_') + 1
            if i == 0:
               base_indexes.append(base_str)
            else:
                if base_str[i] == '_':
                    i += 1
                base_indexes.append(base_str[i:])
        return base_indexes

    def Mlt_str(self):
        return str(self.fvalue)

    def Mlt_latex_str(self):
        if self.nargs <= 1:
            return printer.latex(self.fvalue)
        expr_lst = Mlt.expand_expr(self.fvalue,self.Ga)
        latex_str = '\\begin{align*} '
        first = True
        cnt = 1 #  Component count on line
        for term in expr_lst:
            coef_str = str(term[0])
            coef_latex = printer.latex(term[0])
            term_add_flg = isinstance(term[0],Add)
            if term_add_flg:
                coef_latex = r'\left ( ' + coef_latex + r'\right ) '
            if first:
                first = False
            else:
                if coef_str[0].strip() is not '-' or term_add_flg:
                    coef_latex = ' + ' + coef_latex
            for aij in term[1]:
                coef_latex += printer.latex(aij) + ' '
            if cnt == 1:
                latex_str += ' & ' + coef_latex
            else:
                latex_str += coef_latex
            if cnt%self.lcnt == 0:
                latex_str += '\\\\ '
                cnt = 1
            else:
                cnt += 1
        if self.lcnt == len(expr_lst) or self.lcnt == 1:
            latex_str = latex_str[:-3]
        latex_str = latex_str + ' \\end{align*} \n'
        return  latex_str

    def Fmt(self, lcnt=1, title=None):
        """
        Set format for printing of Tensors -

            lcnt = Number of components per line

        Usage for tensor T example is -

            T.fmt('2','T')

        output is

            print 'T = '+str(A)

        with two components per line.  Works for both standard printing and
        for latex.
        """
        self.lcnt = lcnt
        latex_str = printer.GaLatexPrinter.latex(self)
        self.lcnt = 1

        if printer.GaLatexPrinter.ipy:
            if title is None:
                if r'\begin{align*}' not in latex_str:
                    latex_str = r'\begin{equation*} ' + latex_str + r' \end{equation*}'
            else:
                if r'\begin{align*}' not in latex_str:
                    latex_str = r'\begin{equation*} ' + title + ' = ' + latex_str + r' \end{equation*}'
                else:
                    latex_str = latex_str.replace(r'\begin{align*}', r'\begin{align*} ' + title)
                    latex_str = latex_str.replace('&', '=&', 1)
            from IPython.core.display import display, Math
            display(Math(latex_str))
        else:
            if title is not None:
                print title + ' = ' + latex_str
            else:
                print latex_str
        return

    @staticmethod
    def expand_expr(expr,ga):
        lst_expr = []
        expr = expand(expr)
        for term in expr.args:
            coef = S(1)
            a_lst = []
            for factor in term.args:
                if factor in ga.acoefs:
                    a_lst.append(factor)
                else:
                    coef *= factor
            a_lst = tuple([x for x in a_lst if x in ga.acoefs])
            b_lst = tuple([ga.acoefs.index(x) for x in a_lst])
            lst_expr.append((coef,a_lst,b_lst))
        lst_expr = sorted(lst_expr, key=lambda x: x[2])
        new_lst_expr = []
        previous = (-1,)
        first = True
        for term in lst_expr:
            if previous == term[2]:
                coef += term[0]
                previous = term[2]
            else:
                if not first:
                    new_lst_expr.append((coef, a))
                else:
                    first = False
                coef = term[0]
                previous = term[2]
                a = term[1]
        new_lst_expr.append((coef, a))
        return new_lst_expr

    def __init__(self, f, Ga, nargs=None, fct=False):
        #  f is a function, a multivector, a string, or a component expression
        #  self.f is a function or None such as T | a_1 where T and a_1 are vectors
        #  self.fvalue is a component expression such as
        #  T_x*a_1__x+T_y*a_1__y+T_z*a_1__z for a rank 1 tensor in 3 space and all
        #  symbols are sympy real scalar symbols
        self.Ga = Ga
        self.lcnt = 1
        if isinstance(f, mv.Mv):
            if f.is_vector(): # f is vector T = f | a1
                self.nargs = 1
                Mlt.increment_slots(1, Ga)
                self.fvalue = (f | Ga.a[0]).obj
                self.f = None
            else: #  To be inplemented for f a general pure grade mulitvector
                self.nargs = nargs
                self.fvalue = f
                self.f = None
        elif isinstance(f, lt.Lt): #  f is linear transformation T = a1 | f(a2)
            self.nargs = 2
            Mlt.increment_slots(2, Ga)
            self.fvalue = (Ga.a[0] | f(Ga.a[1])).obj
            self.f = None
        elif isinstance(f, str) and nargs is not None:
            self.f = None
            self.nargs = nargs
            Mlt.increment_slots(nargs, Ga)
            if nargs > 1: #  General tensor of rank > 1
                t_indexes = nargs * [Mlt.extact_basis_indexes(self.Ga)]
                self.fvalue = 0
                for (t_index,a_prod) in zip(itertools.product(*t_indexes),
                                            itertools.product(*self.Ga.pdiffs)):
                    if fct: #  Tensor field
                        coef = Function(f+'_'+''.join(map(str,t_index)),real=True)(*self.Ga.coords)
                    else: #  Constant Tensor
                        coef = symbols(f+'_'+''.join(map(str,t_index)),real=True)
                    coef *= reduce(lambda x, y: x*y, a_prod)
                    self.fvalue += coef
            else: #  General tensor of rank = 1
                self.fvalue = 0
                for (t_index,a_prod) in zip(Mlt.extact_basis_indexes(self.Ga),self.Ga.pdiffs[0]):
                    if fct: #  Tensor field
                        coef = Function(f+'_'+''.join(map(str,t_index)),real=True)(*self.Ga.coords)
                    else: #  Constant Tensor
                        coef = symbols(f+'_'+''.join(map(str,t_index)),real=True)
                    self.fvalue += coef * a_prod
        else:
            if isinstance(f, types.FunctionType): #  Tensor defined by general multi-linear function
                args, kargs, kwargs, defaults = inspect.getargspec(f)
                self.nargs = len(args)
                self.f = f
                Mlt.increment_slots(self.nargs, Ga)
                self.fvalue = f(*tuple(Ga.a[0:self.nargs]))
            else: # Tensor defined by component expression
                self.f = None
                self.nargs = nargs
                Mlt.increment_slots(self.nargs, Ga)
                self.fvalue = f

    def __str__(self):
        if printer.GaLatexPrinter.latex_flg:
            Printer = printer.GaLatexPrinter
        else:
            Printer = printer.GaPrinter
        return Printer().doprint(self)

    def __call__(self, *kargs):
        if len(kargs) == 0:
            return self.fvalue
        if self.f is not None:
            return self.f(*kargs)
        else:
            sub_lst = []
            for (x, ai) in zip(kargs, self.Ga.pdiffs):
                for (r_base, aij) in zip(self.Ga.r_basis_mv, ai):
                    sub_lst.append((aij, (r_base | x).scalar()))
            return self.fvalue.subs(sub_lst,simultaneous=True)

    def __add__(self, X):
        if isinstance(Mlt, X):
            if self.nargs == X.nargs:
                return Mlt(self.fvalue + X.fvalue, self.Ga, self.nargs)
            else:
                raise ValueError('In Mlt add number of args not the same\n')
        else:
            raise TypeError('In Mlt add second argument not an Mkt\n')

    def __sub__(self, X):
        if isinstance(Mlt, X):
            if self.nargs == X.nargs:
                return Mlt(self.fvalue - X.fvalue, self.Ga, self.nargs)
            else:
                raise ValueError('In Mlt sub number of args not the same\n')
        else:
            raise TypeError('In Mlt sub second argument not an Mlt\n')

    def __mul__(self, X):
        if isinstance(X, Mlt):
            nargs = self.nargs + X.nargs
            Mlt.increment_slots(nargs, self.Ga)
            self_args = self.Ga.a[:self.nargs]
            X_args = X.Ga.a[self.nargs:nargs]
            value = expand(self(*self_args) * X(*X_args))
            return Mlt(value, self.Ga, nargs)
        else:
            return Mlt(X * self.fvalue, self.Ga, self.nargs)

    def __xor__(self, X):
        if isinstance(X, Mlt):
            nargs = self.nargs + X.nargs
            Mlt.increment_slots(nargs, self.Ga)
            value = self(*self.Ga.a[:self.nargs]) ^ X(X.Ga.a[self.nargs:nargs])
            return Mlt(value, self.Ga, nargs)
        else:
            return Mlt(X * self.fvalue, self.Ga, self.nargs)

    def __or__(self, X):
        if isinstance(X, Mlt):
            nargs = self.nargs + X.nargs
            Mlt.increment_slots(nargs, self.Ga)
            value = self(*self.Ga.a[:self.nargs]) | X(X.Ga.a[self.nargs:nargs])
            return Mlt(value, self.Ga, nargs)
        else:
            return Mlt(X * self.fvalue, self.Ga, self.nargs)

    def _repr_latex_(self):
        latex_str = printer.GaLatexPrinter.latex(self)
        if r'\begin{align*}' not in latex_str:
            latex_str = r'\begin{equation*} ' + latex_str + r' \end{equation*}'
        return latex_str

    def dd(self):
        Mlt.increment_slots(self.nargs + 1, self.Ga)
        dd_fvalue = (self.Ga.a[self.nargs] | self.Ga.grad) * self.fvalue
        return Mlt(dd_fvalue, self.Ga, self.nargs + 1)

    def pdiff(self, slot):
        # Take geometric derivative of mlt with respect to slot argument
        self.Ga.dslot = slot - 1
        return self.Ga.grad * self.Ga.mv(self.fvalue)

    @staticmethod
    def remove_slot(mv, slot, nargs, ga):
        if slot == nargs:
            return mv
        for islot in range(slot, nargs):
            mv = mv.subs(zip(ga.pdiffs[islot], ga.pdiffs[islot - 1]))
        return mv

    def contract(self, slot1, slot2):
        min_slot = min(slot1, slot2)
        max_slot = max(slot1, slot2)
        cnargs = self.nargs - 2
        self.Ga.dslot = min_slot - 1
        grad_self = self.Ga.grad * self.Ga.mv(self.fvalue)
        grad_self = Mlt.remove_slot(grad_self.obj, min_slot, self.nargs, self.Ga)
        self.Ga.dslot = max_slot - 2
        div_grad_self = self.Ga.grad | self.Ga.mv(grad_self)
        div_grad_self = Mlt.remove_slot(div_grad_self.obj, max_slot - 1, self.nargs - 1, self.Ga)
        return Mlt(div_grad_self, self.Ga, cnargs)

    def cderiv(self):
        Mlt.increment_slots(self.nargs + 1, self.Ga)
        agrad = self.Ga.a[self.nargs] | self.Ga.grad
        CD = Mlt((agrad * self.Ga.mv(self.fvalue)).obj, self.Ga, self.nargs + 1)
        if CD != 0:
            CD = CD.fvalue
        for i in range(self.nargs):
            args = self.Ga.a[:self.nargs]
            tmp = agrad * self.Ga.a[i]
            if tmp.obj != 0:
                args[i] = tmp
                CD = CD - self(*args)
        CD = Mlt(CD, self.Ga, self.nargs + 1)
        return CD

    def expand(self):
        self.fvalue = expand(self.fvalue)
        return self

    def comps(self):
        basis = self.Ga.mv()
        rank = self.nargs
        ndim = len(basis)
        i_indexes = rank*[range(ndim)]
        indexes = rank*[basis]
        i = 1
        output = ''
        for (e,i_index) in zip(itertools.product(*indexes),itertools.product(*i_indexes)):
            if i_index[-1] % ndim == 0: print ''
            output += str(i)+':'+str(i_index)+':'+str(self(*e)) + '\n'
            i += 1
        return output
