#mv.py

import sys
import itertools
import copy
from operator import itemgetter, mul, add
from numpy.linalg import matrix_rank
from sympy import Symbol, Function, S, expand, Add, Mul, Pow, Basic, \
    sin, cos, sinh, cosh, sqrt, trigsimp, \
    simplify, diff, Rational
from sympy import N as Nsympy
import metric
import printer

half = Rational(1, 2)

modules = \
"""
from sympy import symbols, sin, Function
from mv import Mv
from ga import Ga, half
from printer import Eprint, xdvi
from lt import Lt
"""

########################### Multivector Class ##########################


class Mv(object):
    """
    Wrapper class for multivector objects (self.obj) so that it is easy
    to overload operators (*,^,|,<,>) for the various multivector
    products and for printing.  Also provides an __init__ fuction to
    easily instanciate multivector objects.  Additionally, the functionality
    of the multivector derivative have been added via the special vector
    'grad' so that one can take the geometric derivative of a multivector
    function 'A' by applying 'grad' from the left, 'grad*A', or the
    right 'A*grad' for both the left and right derivatives.  The operator
    between the 'grad' and the 'A' can be any of the multivector product
    operators.

    If 'f' is a scalar function 'grad*f' is the usual gradient of a function.
    If 'A' is a vector function 'grad|f' is the divergence of 'A' and
    '-I*(grad^A)' is the curl of 'A' (I is the pseudo scalar for the geometric
    algebra)
    """

    ################### Multivector initialization #####################

    restore = False
    init_slots = {'f': (False, 'True if function of coordinates'),
                  'ga': (None, 'Geometric algebra to be used with multivectors'),
                  'coords': (None, 'Coordinates to be used with multivector function'),
                  'recp': (None, 'Normalization for reciprocal vector')}

    @staticmethod
    def setup(ga):
        """
        Set up constant mutilvectors reqired for multivector class for
        a given geometric algebra, 'ga'.
        """
        Mv.fmt = 1

        basis = [Mv(x, ga=ga) for x in ga.basis]
        I = Mv(ga.iobj, ga=ga)  # default pseudoscalar
        x = Mv('XxXx', 'vector', ga=ga)  # testing vectors
        # return default basis vectors and grad vector if coords defined
        return I, basis, x

    @staticmethod
    def get_Ga(name):
        return(Mv.ga[name])

    @staticmethod
    def Format(mode=1):
        Mv.fmt = mode
        return

    def characterise_Mv(self):
        if self.char_Mv:
            return
        obj = self.obj
        if obj.is_commutative:
            self.i_grade = 0
            self.is_blade_rep = True
            self.grades = [0]
            return
        if isinstance(obj, Add):
            args = obj.args
        else:
            args = [obj]
        grades = []
        self.is_blade_rep = True
        for term in args:
            if term.is_commutative:
                if 0 not in grades:
                    grades.append(0)
            else:
                c, nc = term.args_cnc(split_1=False)
                blade = nc[0]
                if blade in self.Ga.blades_lst:
                    grade = self.Ga.blades_to_grades_dict[blade]
                    if not grade in grades:
                        grades.append(grade)
                else:
                    self.char_Mv = True
                    self.is_blade_rep = False
                    self.i_grade = None
                    return
        if len(grades) == 1:
            self.i_grade = grades[0]
        else:
            self.i_grade = None
        self.grades = grades
        self.char_Mv = True
        return

    def make_grade(self, *kargs, **kwargs):
        root = kargs[0] + '__'
        grade = kargs[1]
        self.i_grade = grade

        if isinstance(kwargs['f'], bool) and not kwargs['f']:
            self.obj = sum([Symbol(root + super_script, real=True) * base
                            for (super_script, base) in zip(self.Ga.blade_super_scripts[grade], self.Ga.blades[grade])])
        else:
            if isinstance(kwargs['f'], bool):
                self.obj = sum([Function(root + super_script, real=True)(*self.Ga.coords) * base
                    for (super_script, base) in zip(self.Ga.blade_super_scripts[grade], self.Ga.blades[grade])])
            else:
                self.obj = sum([Function(root + super_script, real=True)(kwargs['f']) * base
                    for (super_script, base) in zip(self.Ga.blade_super_scripts[grade], self.Ga.blades[grade])])
        return

    def make_scalar(self, *kargs, **kwargs):
        if isinstance(kargs[0],str):
            if 'f' in kwargs and kwargs['f']:
                self.obj = Function(kargs[0])(*self.Ga.coords)
            else:
                self.obj = Symbol(kargs[0], real=True)
        else:
            self.obj = kargs[0]
        return

    def make_vector(self, *kargs, **kwargs):
        self.make_grade(*(kargs[0], 1), **kwargs)
        return

    def make_bivector(self, *kargs, **kwargs):
        self.make_grade(*(kargs[0], 2), **kwargs)
        return

    def make_pseudo_scalar(self, *kargs, **kwargs):
        self.make_grade(*(kargs[0], self.Ga.n), **kwargs)
        return

    def make_multivector(self, *kargs, **kwargs):
        self.make_scalar(kargs[0], **kwargs)
        tmp = self.obj
        for grade in self.Ga.n_range:
            self.make_grade(*(kargs[0], grade + 1), **kwargs)
            tmp += self.obj
        self.obj = tmp
        return

    def make_spinor(self, *kargs, **kwargs):
        self.make_scalar(kargs[0], **kwargs)
        tmp = self.obj
        for grade in self.Ga.n_range:
            if (grade + 1) % 2 == 0:
                self.make_grade(*(kargs[0], grade + 1), **kwargs)
                tmp += self.obj
        self.obj = tmp
        return

    def make_odd(self, *kargs, **kwargs):
        self.make_scalar(kargs[0], **kwargs)
        tmp = S(0)
        for grade in self.Ga.n_range:
            if (grade + 1) % 2 == 1:
                self.make_grade(*(kargs[0], grade + 1), **kwargs)
                tmp += self.obj
        self.obj = tmp
        return

    init_dict = {'scalar': make_scalar,
                 'vector': make_vector,
                 'bivector': make_bivector,
                 'pseudo': make_pseudo_scalar,
                 'mv': make_multivector,
                 'spinor': make_spinor,
                 'even': make_spinor,
                 'odd': make_odd,
                 'grade': make_grade}

    def __init__(self, *kargs, **kwargs):

        if 'ga' not in kwargs:
            raise ValueError("Geometric algebra key inplut 'ga' required")

        kwargs = metric.test_init_slots(Mv.init_slots, **kwargs)

        self.Ga = kwargs['ga']
        self.recp = kwargs['recp']  # Normalization for reciprocal vectors

        self.char_Mv = False
        self.i_grade = None  # if pure grade mv, grade value
        self.grades = None  # list of grades in mv
        self.is_blade_rep = True  # flag for blade representation
        self.blade_flg = None  # if is_blade is called flag is set
        self.versor_flg = None  # if is_versor is called flag is set
        self.coords = self.Ga.coords
        self.fmt = 1

        if len(kargs) == 0:  # default constructor 0
            self.obj = S(0)
            self.i_grade = 0
        elif len(kargs) == 1 and not isinstance(kargs[0], str):  # copy constructor
            x = kargs[0]
            if isinstance(x, Mv):
                self.obj = x.obj
                self.is_blade_rep = x.is_blade_rep
                self.i_grade = x.i_grade
            else:
                self.obj = x
                self.is_blade_rep = True
                self.characterise_Mv()
        else:
            if kargs[1] not in Mv.init_dict:
                raise ValueError('"' + kargs[1] + '" not an allowed multivector type.')
            mode = kargs[1]
            kargs = [kargs[0]] + list(kargs[2:])
            Mv.init_dict[mode](self, *kargs, **kwargs)
            self.characterise_Mv()

    ################# Multivector member functions #####################

    def base_rep(self):
        if self.is_blade_rep:
            self.obj = self.Ga.blade_to_base_rep(self.obj)
            self.is_blade_rep = False
            return self
        else:
            return self

    def blade_rep(self):
        if self.is_blade_rep:
            return self
        else:
            self.obj = self.Ga.base_to_blade_rep(self.obj)
            self.is_blade_rep = True
            return self

    def __eq__(self, A):
        if not isinstance(A, Mv):
            return False
        self.obj = expand(self.obj)
        A.obj = expand(A.obj)
        if self.obj == A.obj:
            return True
        else:
            return False

    def __neg__(self):
        return Mv(-self.obj, ga=self.Ga)

    def __add__(self, A):

        if (not isinstance(A, Mv)) and (not isinstance(A, Dop)):
            return Mv(self.obj + A, ga=self.Ga)

        if self.Ga.name != A.Ga.name:
            raise ValueError('In + operation Mv arguments are not from same geometric algebra')

        if isinstance(A, Dop):
            return A.Add(self.obj, A)

        if self.is_blade_rep == A.is_blade_rep:
            return Mv(self.obj + A.obj, ga=self.Ga)
        else:
            if self.is_blade_rep:
                A = A.blade_rep()
            else:
                self = self.blade_rep()
            return Mv(self.obj + A.obj, ga=self.Ga)

    def __radd__(self, A):
        return(self + A)

    def __add_ab__(self, A):  # self += A
        self.obj += A.obj
        self.char_Mv = False
        self.characterise_Mv()
        return(self)

    def __sub__(self, A):

        if (not isinstance(A, Mv)) and (not isinstance(A, Dop)):
            return Mv(self.obj - A, ga=self.Ga)

        if self.Ga.name != A.Ga.name:
            raise ValueError('In - operation Mv arguments are not from same geometric algebra')

        if isinstance(A, Dop):
            return A.Add(self.obj, A, mode='-')

        if self.is_blade_rep == A.is_blade_rep:
            return Mv(self.obj - A.obj, ga=self.Ga)
        else:
            if self.is_blade_rep:
                A = A.blade_rep()
            else:
                self = self.blade_rep()
            return Mv(self.obj - A.obj, ga=self.Ga)

    def __rsub__(self, A):
        return -self + A

    def __sub_ab__(self, A):  # self -= A
        self.obj -= A.obj
        self.char_Mv = False
        self.characterise_Mv()
        return(self)

    def __mul__(self, A):

        if (not isinstance(A, Mv)) and (not isinstance(A, Dop)):
            return Mv(expand(A * self.obj), ga=self.Ga)

        if self.Ga.name != A.Ga.name:
            raise ValueError('In * operation Mv arguments are not from same geometric algebra')

        if isinstance(A, Dop):
            return A.Mul(self, A, mode='*')

        if self.is_blade_rep and A.is_blade_rep:
            self = self.base_rep()
            A = A.base_rep()

            selfxA = Mv(self.Ga.mul(self.obj, A.obj), ga=self.Ga)
            selfxA.is_blade_rep = False
            selfxA = selfxA.blade_rep()

            self = self.blade_rep()
            A = A.blade_rep()
        elif self.is_blade_rep:
            self = self.base_rep()

            selfxA = Mv(self.Ga.mul(self.obj, A.obj), ga=self.Ga)
            selfxA.is_blade_rep = False
            selfxA = selfxA.blade_rep()

            self = self.blade_rep()
        elif A.is_blade_rep:
            A = A.base_rep()

            selfxA = Mv(self.Ga.mul(self.obj, A.obj), ga=self.Ga)
            selfxA.is_blade_rep = False
            selfxA = selfxA.blade_rep()

            A = A.blade_rep()
        else:
            selfxA = Mv(self.Ga.mul(self.obj, A.obj), ga=self.Ga)

        return selfxA

    def __rmul__(self, A):
            return Mv(expand(A * self.obj), ga=self.Ga)

    def __mul_ab__(self, A):  # self *= A
        self.obj *= A.obj
        self.char_Mv = False
        self.characterise_Mv()
        return(self)

    def __div__(self, A):
        self_div = Mv(self.obj, ga=self.Ga)
        self_div.obj /= A
        return(self_div)

    def __str__(self):
        if printer.GaLatexPrinter.latex_flg:
            Printer = printer.GaLatexPrinter
        else:
            Printer = printer.GaPrinter
        return Printer().doprint(self)

    def Mv_str(self):
        # str representation of multivector
        if self.i_grade == 0:
            return str(self.obj)
        self.obj = expand(self.obj)
        self.characterise_Mv()
        self.obj = metric.Simp.apply(self.obj)
        if self.is_blade_rep or self.Ga.is_ortho:
            base_keys = self.Ga.blades_lst
            grade_keys = self.Ga.blades_to_grades_dict
        else:
            base_keys = self.Ga.bases_lst
            grade_keys = self.Ga.bases_to_grades_dict
        if isinstance(self.obj, Add):  # collect coefficients of bases
            if self.obj.is_commutative:
                return self.obj
            args = self.obj.args
            terms = {}  # dictionary with base indexes as keys
            grade0 = S(0)
            for arg in args:
                c, nc = arg.args_cnc()
                if len(c) > 0:
                    c = reduce(mul, c)
                else:
                    c = S(1)
                if len(nc) > 0:
                    base = nc[0]
                    if base in base_keys:
                        index = base_keys.index(base)
                        if index in terms:
                            (c_tmp, base, g_keys) = terms[index]
                            terms[index] = (c_tmp + c, base, g_keys)
                        else:
                            terms[index] = (c, base, grade_keys[base])
                else:
                    grade0 += c
            if grade0 != S(0):
                terms[0] = (grade0, S(1), -1)
            terms = terms.items()

            sorted_terms = sorted(terms, key=itemgetter(0))  # sort via base indexes

            s = str(sorted_terms[0][1][0] * sorted_terms[0][1][1])
            if self.fmt == 3:
                s = ' ' + s + '\n'
            if self.fmt == 2:
                s = ' ' + s
            old_grade = sorted_terms[0][1][2]
            for (key, (c, base, grade)) in sorted_terms[1:]:
                term = str(c * base)
                if self.fmt == 2 and old_grade != grade:  # one grade per line
                    old_grade = grade
                    s += '\n'
                if term[0] == '-':
                    term = ' - ' + term[1:]
                else:
                    term = ' + ' + term
                if self.fmt == 3:  # one base per line
                    s += term + '\n'
                else:  # one multivector per line
                    s += term
            if s[-1] == '\n':
                s = s[:-1]
            return s
        else:
            return str(self.obj)

    def Mv_latex_str(self):
        self.first_line = True

        def append_plus(c_str):
            if self.first_line:
                self.first_line = False
                return c_str
            else:
                c_str = c_str.strip()
                if c_str[0] == '-':
                    return ' ' + c_str
                else:
                    return ' + ' + c_str

        # str representation of multivector
        self.obj = expand(self.obj)
        self.characterise_Mv()
        self.obj = metric.Simp.apply(self.obj)

        if self.is_blade_rep or self.Ga.is_ortho:
            base_keys = self.Ga.blades_lst
            grade_keys = self.Ga.blades_to_grades_dict
        else:
            base_keys = self.Ga.bases_lst
            grade_keys = self.Ga.bases_to_grades_dict
        if isinstance(self.obj, Add):
            args = self.obj.args
        else:
            args = [self.obj]
        terms = {}  # dictionary with base indexes as keys
        grade0 = S(0)
        for arg in args:
            c, nc = arg.args_cnc(split_1=False)
            if len(c) > 0:
                c = reduce(mul, c)
            else:
                c = S(1)
            if len(nc) > 0:
                base = nc[0]
                if base in base_keys:
                    index = base_keys.index(base)
                    if index in terms:
                        (c_tmp, base, g_keys) = terms[index]
                        terms[index] = (c_tmp + c, base, g_keys)
                    else:
                        terms[index] = (c, base, grade_keys[base])
            else:
                grade0 += c
        if grade0 != S(0):
            terms[-1] = (grade0, S(1), 0)
        terms = terms.items()

        sorted_terms = sorted(terms, key=itemgetter(0))  # sort via base indexes

        if len(sorted_terms) == 1 and sorted_terms[0][1][2] == 0:  # scalar
            return printer.latex(printer.coef_simplify(sorted_terms[0][1][0]))

        lines = []
        old_grade = -1
        s = ''
        for (index, (coef, base, grade)) in sorted_terms:
            coef = printer.coef_simplify(coef)
            #coef = simplify(coef)
            l_coef = printer.latex(coef)
            if l_coef == '1' and base != S(1):
                l_coef = ''
            if base == S(1):
                l_base = ''
            else:
                l_base = printer.latex(base)
            if isinstance(coef, Add):
                cb_str = '\\left ( ' + l_coef + '\\right ) ' + l_base
            else:
                cb_str = l_coef + ' ' + l_base
            if self.fmt == 3:  # One base per line
                lines.append(append_plus(cb_str))
            elif self.fmt == 2:  # One grade per line
                if grade != old_grade:
                    old_grade = grade
                    if not self.first_line:
                        lines.append(s)
                    s = append_plus(cb_str)
                else:
                    s += append_plus(cb_str)
            else:  # One multivector per line
                s += append_plus(cb_str)
        if self.fmt == 2:
            lines.append(s)
        if self.fmt >= 2:
            if len(lines) == 1:
                return lines[0]
            s = ' \\begin{align*} '
            for line in lines:
                s += ' & ' + line + ' \\\\ '
            s = s[:-3] + ' \\end{align*} \n'
        return s

    def __xor__(self, A):  # wedge (^) product

        if (not isinstance(A, Mv)) and (not isinstance(A, Dop)):
            return Mv(A * self.obj, ga=self.Ga)

        if self.Ga.name != A.Ga.name:
            raise ValueError('In ^ operation Mv arguments are not from same geometric algebra')

        if isinstance(A, Dop):
            return A.Mul(self, A, mode='^')

        self = self.blade_rep()
        A = A.blade_rep()
        self_W_A = self.Ga.wedge(self.obj, A.obj)
        self_W_A = Mv(self_W_A, ga=self.Ga)
        return self_W_A

    def __rxor__(self, A):  # wedge (^) product
        if not isinstance(A, Mv):
            return Mv(A * self.obj, ga=self.Ga)
        else:
            return A * self

    def __or__(self, A):  # dot (|) product
        if (not isinstance(A, Mv)) and (not isinstance(A, Dop)):
            return Mv(ga=self.Ga)

        if self.Ga.name != A.Ga.name:
            raise ValueError('In | operation Mv arguments are not from same geometric algebra')

        self.Ga.dot_mode = '|'

        if isinstance(A, Dop):
            return A.Mul(self, A, mode='|')

        self = self.blade_rep()
        A = A.blade_rep()
        self_dot_A = Mv(self.Ga.dot(self.obj, A.obj), ga=self.Ga)
        return self_dot_A

    def __ror__(self, A):  # dot (|) product
        if not isinstance(A, Mv):
            return Mv(ga=self.Ga)
        else:
            return A | self

    def __lt__(self, A):  # left contraction (<)

        if (not isinstance(A, Mv)) and (not isinstance(A, Dop)):  # sympy scalar
            return Mv(A * self.obj, ga=self.Ga)

        if self.Ga.name != A.Ga.name:
            raise ValueError('In < operation Mv arguments are not from same geometric algebra')

        self.Ga.dot_mode = '<'

        if isinstance(A, Dop):
            return A.Mul(self, A, mode='<')

        self = self.blade_rep()
        A = A.blade_rep()
        self_lc_A = Mv(self.Ga.dot(self.obj, A.obj), ga=self.Ga)
        return self_lc_A

    def __gt__(self, A):  # right contraction (>)


        if (not isinstance(A, Mv)) and (not isinstance(A, Dop)):  # sympy scalar
            return self.Ga.mv(A * self.scalar())

        if self.Ga.name != A.Ga.name:
            raise ValueError('In > operation Mv arguments are not from same geometric algebra')

        self.Ga.dot_mode = '>'

        if isinstance(A, Dop):
            return A.Mul(self, A, mode='>')

        self = self.blade_rep()
        A = A.blade_rep()
        self_rc_A = Mv(self.Ga.dot(self.obj, A.obj), ga=self.Ga)
        return self_rc_A

    def collect(self):
        # group coeffients of blades of multivector
        # so there is only one coefficient per grade
        self.obj = expand(self.obj)
        if self.is_blade_rep or Mv.Ga.is_ortho:
            c = self.Ga.blades_lst
        else:
            c = self.Ga.bases_lst
        self.obj = self.obj.collect(c)
        return self

    def is_scalar(self):
        grades = self.Ga.grades(self.obj)
        if len(grades) == 1 and grades[0] == 0:
            return True
        else:
            return False

    def is_vector(self):

        self = self.blade_rep()
        self.characterise_Mv()

        if self.i_grade is not None and self.i_grade == 1:
            return True
        else:
            #grades = self.Ga.grades(self.obj)
            if len(self.grades) == 1 and self.grades[0] == 1:
                return True
            else:
                return False

    def is_blade(self):  # True is self is blade, otherwise False
        # sets self.blade_flg and returns value
        if self.blade_flg is not None:
            return self.blade_flg
        else:
            if self.is_versor():
                if self.i_grade is not None:
                    self.blade_flg = True
                else:
                    self.blade_flage = False
            else:
                self.blade_flg = False
            return self.blade_flg

    def is_versor(self):  # Test for versor (geometric product of vectors)
        """
        This follows Leo Dorst's test for a versor.
        Leo Dorst, 'Geometric Algebra for Computer Science,' p.533
        Sets self.versor_flg and returns value
        """
        if self.versor_flg is not None:
            return self.versor_flg
        self.characterise_Mv()
        even = True
        odd = True
        # test to see if self is pure even or pure odd multivector
        for grade in self.grades:
            if grade % 2 == 0:
                odd = False
            else:
                even = False
        if even is False and odd is False:
            self.versor_flg = False
            return self.versor_flg

        # see if self*x*self.rev() returns a vector for x an arbitrary vector
        test = self * self.Ga.mv_x * self.rev()
        self.versor_flg = test.is_vector()
        return self.versor_flg

    def scalar(self):
        # return scalar part of multivector
        # as sympy expression
        return self.Ga.scalar_part(self.obj)

    def get_grade(self, r):
        # return r-th grade of multivector as
        # a multivector
        return Mv(self.Ga.get_grade(self.obj, r), ga=self.Ga)

    def get_coefs(self, grade):
        (coefs, bases) = metric.linear_expand(self.obj)
        bases_lst = self.Ga.blades_lst
        cb = zip(coefs, bases)
        cb = sorted(cb, key=lambda x: self.Ga.blades[grade].index(x[1]))
        (coefs, bases) = zip(*cb)
        return coefs

    def proj(self, bases_lst):
        bases_lst = [x.obj for x in bases_lst]
        (coefs, bases) = metric.linear_expand(self.obj)
        obj = 0
        for (coef, base) in zip(coefs, bases):
            if base in bases_lst:
                obj += coef * base
        return Mv(obj, ga=self.Ga)

    def dual(self):
        return self.Ga.i * self

    def even(self):
        # return even parts of multivector
        return Mv(self.Ga.even_odd(self.obj, True), ga=self.Ga)

    def odd(self):
        # return odd parts of multivector
        return Mv(self.Ga.even_odd(self.obj, False), ga=self.Ga)

    def rev(self):
        self = self.blade_rep()
        return Mv(self.Ga.reverse(self.obj), ga=self.Ga)

    def diff(self, coord):
        Dself = Mv(ga=self.Ga)
        if coord not in self.Ga.coords:
            if self.Ga.par_coords is None:
                Dself.obj = diff(self.obj, coord)
            else:
                Dself.obj = diff(self.obj, coord)
                for x_coord in self.Ga.coords:
                    f = self.Ga.par_coords[x_coord]
                    if f != S(0):
                        tmp1 = self.Ga.pDiff(self.obj, x_coord)
                        tmp2 = diff(f, coord)
                        Dself.obj += tmp1 * tmp2
            Dself.characterise_Mv()
            return Dself
        else:
            Dself.obj = self.Ga.pDiff(self.obj, coord)
            Dself.characterise_Mv()
            return Dself

    def pdiff(self, var):
        return Mv(self.Ga.pDiff(self.obj, var), ga=self.Ga)

    def Grad(self, coords, mode='*', left=True):
        """
        Returns various derivatives (*,^,|,<,>) of multivector functions
        with respect to arbitrary coordinates, 'coords'.  This would be
        used where you have a multivector function of both the basis
        coordinate set and and auxilliary coordinate set.  Consider for
        example a linear transformation in which the matrix coefficients
        depend upon the manifold coordinates, but the vector being
        transformed does not and you wish to take the divergence of the
        linear transformation with respect to the linear argument.
        """
        return Mv(self.Ga.Diff(self, mode, left, coords=coords), ga=self.Ga)

    def exp(self, hint='+'):  # Calculate exponential of multivector
        """
        Only works if square of multivector is a scalar.  If square is a
        number we can determine if square is > or < zero and hence if
        one should use trig or hyperbolic functions in expansion.  If
        square is not a number use 'hint' to determine which type of
        functions to use in expansion
        """
        self_sq = self * self
        if self_sq.is_scalar():
            sq = self_sq.obj
            if sq.is_number:
                if sq > S(0):
                    norm = sqrt(sq)
                    value = self.obj / norm
                    return Mv(cosh(norm) + sinh(norm) * value, ga=self.Ga)
                elif sq < S(0):
                    norm = sqrt(-sq)
                    value = self.obj / norm
                    return Mv(cos(norm) + sin(norm) * value, ga=self.Ga)
                else:
                    return Mv(S(1), 'scalar', ga=self.Ga)
            else:
                norm = metric.square_root_of_expr(sq)
                value = self.obj / norm
                if hint == '+':
                    return Mv(cosh(norm) + sinh(norm) * value, ga=self.Ga)
                else:
                    return Mv(cos(norm) + sin(norm) * value, ga=self.Ga)
        else:
            raise ValueError('"' + str(self) + '**2" is not a scalar in exp.')

    def set_coef(self, igrade, ibase, value):
        if self.blade_rep:
            base = self.Ga.blades[igrade][ibase]
        else:
            base = self.Ga.bases[igrade][ibase]
        (coefs, bases) = metric.linear_expand(self.obj)
        bases_lst = list(bases)  # python 2.5
        if base in bases:
            self.obj += (value - coefs[bases_lst.index(base)]) * base
        else:
            self.obj += value * base
        return

    def Fmt(self, fmt=1, title=None):
        """
        Set format for printing of multivectors -

            fmt = 1 - One multivector per line
            fmt = 2 - One grade per line
            fmt = 3 - one base per line

        Usage for multivector A example is -

            A.fmt('2','A')

        output is

            print 'A = '+str(A)

        with one grade per line.  Works for both standard printing and
        for latex.
        """
        self.fmt = fmt
        printer.GaLatexPrinter.fmt = self.fmt
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
            from IPython.core.display import display, Math
            display(Math(latex_str))
        else:
            if title is not None:
                print title + ' = ' + latex_str
            else:
                print latex_str
        return

    def _repr_latex_(self):
        latex_str = printer.GaLatexPrinter.latex(self)
        if r'\begin{align*}' not in latex_str:
            latex_str = r'\begin{equation*} ' + latex_str + r' \end{equation*}'
        return latex_str

    def norm2(self):
        reverse = self.rev()
        product = self * reverse
        if product.is_scalar():
            return product.func(lambda coefficient: abs(coefficient))
        else:
            raise TypeError('"(' + str(product) + ')**2" is not a scalar in norm2.')

    def norm(self):
        reverse = self.rev()
        product = self * reverse
        if product.is_scalar():
            return metric.square_root_of_expr(product.obj)
        else:
            raise TypeError('"(' + str(product) + ')" is not a scalar in norm.')

    def inv(self):
        reverse = self.rev()
        product = self * reverse
        if(product.is_scalar()):
            return reverse.func(lambda coefficient: coefficient / product.obj)
        else:
            raise TypeError('"(' + str(product) + ')" is not a scalar.')

    def func(self, fct):  # Apply function, fct, to each coefficient of multivector
        (coefs, bases) = metric.linear_expand(self.obj)
        s = S(0)
        for (coef, base) in zip(coefs, bases):
            s += fct(coef) * base
        fct_self = Mv(s, ga=self.Ga)
        fct_self.characterise_Mv()
        return fct_self

    def trigsimp(self):
        return self.func(trigsimp)

    def simplify(self, modes=simplify):
        (coefs, bases) = metric.linear_expand(self.obj)
        obj = S(0)
        if isinstance(modes, list) or isinstance(modes, tuple):
            for (coef, base) in zip(coefs, bases):
                for mode in modes:
                    coef = mode(coef)
                obj += coef * base
        else:
            for (coef, base) in zip(coefs, bases):
                obj += modes(coef) * base
        self.obj = obj
        return self

    def subs(self, d):
        # For each scalar coef of the multivector apply substitution argument d
        (coefs, bases) = metric.linear_expand(self.obj)
        obj = S(0)
        for (coef, base) in zip(coefs, bases):
            obj += coef.subs(d) * base
        self.obj = obj
        return self

    def expand(self):
        self.obj = expand(self.obj)
        return self

    def list(self):
        (coefs, bases) = metric.linear_expand(self.obj)
        indexes = []
        key_coefs = []
        for (coef, base) in zip(coefs, bases):
            if base in self.Ga.basis:
                index = self.Ga.basis.index(base)
                key_coefs.append((coef, index))
                indexes.append(index)

        for index in self.Ga.n_range:
            if index not in indexes:
                key_coefs.append((S(0), index))

        key_coefs = sorted(key_coefs, key=itemgetter(1))
        coefs = [x[0] for x in key_coefs]
        return coefs

###################### Differential Operator Class #####################


def mod_convert(l, mod):
    m = 0
    for i in reversed(l):
        m = mod * m + i
    return m


def pdop_sort(coef_lst, pdiff_lst):
    if len(pdiff_lst) == 1:
        return [(coef_lst[0], pdiff_lst[0])]

    mod = max(list(itertools.chain(*pdiff_lst))) + 1
    key_lst = [mod_convert(x, mod) for x in pdiff_lst]
    return [(x, y) for (z, x, y) in sorted(zip(key_lst, coef_lst, pdiff_lst))]


def pdiff_str(pdiffs, coords):
    s = 'D{'
    nz_cnt = 0
    for (pdiff, coord) in zip(pdiffs, coords):
        if pdiff > 0:
            nz_cnt += 1
            if pdiff == 1:
                s += str(coord) + ','
            else:
                s += str(coord) + '^' + str(pdiff) + ','
    if nz_cnt > 0:
        s = s[:-1] + '}'
    else:
        s = ''
    return s


def pdiff_latex_str(pdiffs, coords):
    if len(pdiffs) > 1:
        nz_cnt = sum(pdiffs)
    else:
        nz_cnt = pdiffs[0]
    if nz_cnt == 0:
        return ''
    if nz_cnt > 1:
        s = '\\bfrac{\partial^{' + str(nz_cnt) + '}}{'
    else:
        s = '\\bfrac{\partial}{'
    for (pdiff, coord) in zip(pdiffs, coords):
        if pdiff > 0:
            if pdiff == 1:
                s += '\\partial ' + printer.latex(coord) + ' '
            else:
                s += '\\partial ' + printer.latex(coord) + '^{' + str(pdiff) + '}'
    return s + '} '


class Dop(object):
    """
    Differential operator class for multivectors.  The operators are of
    the form

        D = D^{i_{1}...i_{n}}\partial_{i_{1}...i_{n}}

    where the D^{i_{1}...i_{n}} are multivector functions of the coordinates
    x_{1},...,x_{n} and \partial_{i_{1}...i_{n}} are partial derivative
    operators

        \partial_{i_{1}...i_{n}} =
             \partial^{i_{1}+...+i_{n}}/\partial{x_{1}^{i_{1}}}...\partial{x_{n}^{i_{n}}}.

    If * is any multivector multiplicative operation then the operator D
    operates on the multivector function F by the following definitions

        D*F = D^{i_{1}...i_{n}}*\partial_{i_{1}...i_{n}}F

    returns a multivector and

        F*D = F*D^{i_{1}...i_{n}}\partial_{i_{1}...i_{n}}

    returns a differential operator.  If the 'cmpflg' in the operator is
    set to 'True' the operation returns

        F*D = (\partial_{i_{1}...i_{n}}F)*D^{i_{1}...i_{n}}

    a multivector function.  For example the representation of the grad
    operator in 3d would be:

        D^{i_{1}...i_{n}} = [e__x,e__y,e__z]
        \partial_{i_{1}...i_{n}} = [(1,0,0),(0,1,0),(0,0,1)].

    See LaTeX documentation for definitions of operator algebraic
    operations +, -, *, ^, |, <, and >.
    """

    init_slots = {'ga': (None, 'Associated geometric algebra'),
                  'cmpflg': (False, 'Complement flag for Dop'),
                  'debug': (False, 'True to print out debugging information')}

    def __init__(self, *kargs, **kwargs):

        kwargs = metric.test_init_slots(Dop.init_slots, **kwargs)

        self.coefs = kargs[0]  # List of coeffients for Dop
        self.pdiffs = kargs[1]  # List of partial derivatives for Dop
        self.cmpflg = kwargs['cmpflg']  # Complement flag (default False)
        self.Ga = kwargs['ga']  # Associated geometric algebra
        self.fmt = 1

        if self.Ga is None:
            raise ValueError('In Dop ga must be defined to instanciate.')

    def remove_zero_coefs(self):

        new_coefs = []
        new_pdiffs = []
        for (coef, pd) in zip(self.coefs, self.pdiffs):
            if coef != S(0):
                new_coefs.append(coef)
                new_pdiffs.append(pd)
        self.coefs = new_coefs
        self.pdiffs = new_pdiffs
        return

    def consolidate_coefs(self):
        """
        Remove zero coefs and consolidate coefs with repeated pdiffs.
        """
        new_coefs = []
        new_pdiffs = []
        for (coef, pd) in zip(self.coefs, self.pdiffs):
            if coef != S(0):
                if pd in new_pdiffs:
                    index = new_pdiffs.index(pd)
                    new_coefs[index] += coef
                else:
                    new_coefs.append(coef)
                    new_pdiffs.append(pd)
        self.coefs = new_coefs
        self.pdiffs = new_pdiffs
        self.TSimplify()
        return

    @staticmethod
    def Add(dop1, dop2, mode='+'):

        if dop1.Ga.name != dop2.Ga.name:
            raise ValueError('In Dop + operation Dop arguments are not from same geometric algebra')

        if isinstance(dop1, Dop) and isinstance(dop2, Dop):
            sum = Dop(list(dop1.coefs), list(dop1.pdiffs), ga=dop1.Ga,
                      cmpflg=dop1.cmpflg)

            for (coef, pdiff) in zip(dop2.coefs, dop2.pdiffs):
                if pdiff in sum.pdiffs:
                    index = sum.pdiffs.index(pdiff)
                    if mode == '+':
                        sum.coefs[index] += coef
                    else:
                        sum.coefs[index] -= coef
                else:
                    sum.pdiffs.append(pdiff)
                    if mode == '+':
                        sum.coefs.append(coef)
                    else:
                        sum.coefs.append(-coef)

        elif isinstance(dop2, Mv):
            sum = Dop(list(dop1.coefs), list(dop1.pdiffs), ga=dop1.Ga,
                      cmpflg=dop1.cmpflg)

            if self.Ga.pd0 in sum.pdiffs:
                index0 = sum.pdiffs.index(self.Ga.pd0)
                if mode == '+':
                    sum.coefs[index0] += dop2.obj
                else:
                    sum.coefs[index0] += dop2.obj
            else:
                sum.pdiffs.append(self.Ga.pd0)
                if mode == '+':
                    sum.coefs.append(dop2.obj)
                else:
                    sum.coefs.append(-dop2.obj)
        else:
            dop2_coefs = [-x for x in dop2.coefs]
            sum = Dop(dop2_coefs, list(dop2.pdiffs), ga=dop2.Ga,
                      cmpflg=dop2.cmpflg)

            if self.Ga.pd0 in sum.pdiffs:
                index0 = sum.pdiffs.index(self.Ga.pd0)
                sum.coefs[index0] += dop1.obj
            else:
                sum.pdiffs.append(self.Ga.pd0)
                sum.coefs.append(dop1.obj)

        sum.remove_zero_coefs()
        return sum

    def __add__(self, dop):
        return Dop.Add(self, dop)

    def __neg__(self):

        ncoefs = [-x for x in self.coefs]

        neg = Dop(ncoefs, list(self.pdiffs), ga=self.Ga,
                  cmpflg=self.cmpflg)

        return neg

    def __sub__(self, dop):
        return Dop.Add(self, dop, mode='-')

    @ staticmethod
    def Mul(dopl, dopr, mode='*'):  # General multiplication of Dop's

        if dopl.Ga.name != dopr.Ga.name:
            raise ValueError('In Dop ' + mode + 'operation Dop arguments are not from same geometric algebra')

        if isinstance(dopr, Mv):
            if not dopl.cmpflg:  # dopr is multivector then dopl*dopr is multivector
                DmodeF = S(0)
                if dopl.Ga.dslot == -1:  # Grad with respect to coordinates
                    for (coef, pdiff) in zip(dopl.coefs, dopl.pdiffs):
                        pdopl = dopl.Ga.pDiff(dopr.obj, pdiff)
                        DmodeF += dopl.Ga.Mul(coef, pdopl, mode=mode)
                    DmodeF = Mv(DmodeF, ga=dopl.Ga).blade_rep()
                else:  # Grad with respect to Mlt *kargs slot
                    for (coef, pdiff) in zip(dopl.coefs, dopl.Ga.pdiffs[dopl.Ga.dslot]):
                        pdopl = dopl.Ga.pDiff(dopr.obj, pdiff)
                        DmodeF += dopl.Ga.Mul(coef, pdopl, mode=mode)
                    DmodeF = Mv(DmodeF, ga=dopl.Ga).blade_rep()
                    dopl.Ga.dslot = -1  # Reset to derivative with respect to coordinates
                return DmodeF
            else:  # dopr is multivector then dopl*dopr is differential operator
                DmodeF = S(0)
                for (coef, pdiff) in zip(dopl.coefs, dopl.pdiffs):
                    DmodeF += dopl.Ga.Mul(coef, dopl.Ga.pDiff(dopr.obj, pdiff), mode=mode)

                coefs = [dopl.Ga.Mul(dopr.obj, x, mode=mode) for x in dopl.coefs]
                pdiffs = list(dopl.pdiffs)

                if dopl.Ga.pd0 in pdiffs:
                    index = pdiffs.index(self.Ga.pd0)
                    coefs[index] += DmodeF
                else:
                    coefs.append(DmodeF)
                    pdiffs.append(dopl.Ga.pd0)

                DmodeF = Dop(coefs, pdiffs, ga=dopl.Ga, cmpflg=dopl.cmpflg)
                DmodeF.consolidate_coefs()
                return DmodeF

        if isinstance(dopl, Mv):
            if not dopr.cmpflg:  # dopl is multivector then dopl*dopr is differential operator
                new_coefs = [dopl.Ga.Mul(dopl.obj, x, mode=mode) for x in dopr.coefs]

                FmodeD = Dop(new_coefs, list(dopr.pdiffs), ga=dopr.Ga,
                             cmpflg=dopr.cmpflg)

                FmodeD.consolidate_coefs()
                return FmodeD
            else:  # dopl is multivector then dopl*dopr is multivector
                FmodeD = S(0)
                if dopr.Ga.dslot == -1:  # Grad with respect to coordinates
                    for (coef, pdiff) in zip(dopr.coefs, dopr.pdiffs):
                        pdopr = dopl.Ga.pDiff(dopl.obj, pdiff)
                        FmodeD += dopl.Ga.Mul(pdopr, coef, mode=mode)
                else:  # Grad with respect to Mlt *kargs slot
                    for (coef, pdiff) in zip(dopl.coefs, dopl.Ga.pdiffs[dopr.Ga.dslot]):
                        pdopl = dopl.Ga.pDiff(dopr.obj, pdiff)
                        DmodeF += dopl.Ga.Mul(coef, pdopl, mode=mode)
                    DmodeF = Mv(DmodeF, ga=dopl.Ga).blade_rep()
                    dopr.Ga.dslot = -1  # Reset to derivative with respect to coordinates
                FmodeD = Mv(FmodeD, ga=self.Ga).blade_rep()
                return FmodeD

        # dopl and dopr are both differential operators

        coefs = []
        pdiffs = []

        # See eq (72) in LaTeX docs

        for (l_coef, l_pdiff) in zip(dopl.coefs, dopl.pdiffs):
            for (r_coef, r_pdiff) in zip(dopr.coefs, dopr.pdiffs):
                coef1 = dopl.Ga.Mul(l_coef, dopl.Ga.pDiff(r_coef, l_pdiff), mode=mode)
                coef2 = dopl.Ga.Mul(l_coef, r_coef, mode=mode)
                pdiff = map(add, l_pdiff, r_pdiff)
                if r_pdiff in pdiffs:
                    index = pdiffs.index(r_pdiff)
                    coefs[index] += coef1
                else:
                    pdiffs.append(r_pdiff)
                    coefs.append(coef1)
                if pdiff in pdiffs:
                    index = pdiffs.index(pdiff)
                    coefs[index] += coef2
                else:
                    pdiffs.append(pdiff)
                    coefs.append(coef2)

        dopl_mode_dopr = Dop(coefs, pdiffs, ga=dopr.Ga,
                             cmpflg=dopl.cmpflg or dopr.cmpflg)

        dopl_mode_dopr.consolidate_coefs()

        return dopl_mode_dopr

    def TSimplify(self):
        new_coefs = []
        for coef in self.coefs:
            new_coefs.append(metric.Simp.apply(coef))
        self.coefs = new_coefs
        return

    def __mul__(self, dopr):  # * geometric product
        return Dop.Mul(self, dopr, mode='*')

    def __xor__(self, dopr):  # ^ outer product
        return Dop.Mul(self, dopr, mode='^')

    def __or__(self, dopr):  # | inner product
        return Dop.Mul(self, dopr, mode='|')

    def __lt__(self, dopr):  # < left contraction
        return Dop.Mul(self, dopr, mode='<')

    def __gt__(self, dopr):  # > right contraction
        return Dop.Mul(self, dopr, mode='>')

    def __str__(self):
        if printer.GaLatexPrinter.latex_flg:
            Printer = printer.GaLatexPrinter
        else:
            Printer = printer.GaPrinter

        return Printer().doprint(self)

    def _repr_latex_(self):
        if printer.GaLatexPrinter.latex_flg:
            Printer = printer.GaLatexPrinter
            return Printer().doprint(self)
        return

    def Dop_str(self):
        if len(self.coefs) == 0:
            return '0'

        dop_sorted = pdop_sort(self.coefs, self.pdiffs)

        s = ''
        for (coef, pdop) in dop_sorted:
            pd_str = pdiff_str(pdop, self.Ga.coords)
            if pd_str == '':
                coef_str = str(Mv(coef, ga=self.Ga))
            elif isinstance(coef, Mul) or \
                  isinstance(coef, Symbol) or \
                  isinstance(coef, Rational) or \
                  isinstance(coef, Function) or \
                  isinstance(coef, Pow):
                coef_str = str(Mv(coef, ga=self.Ga)) + '*'
            elif isinstance(coef, Add):
                coef_str = '(' + str(Mv(coef, ga=self.Ga)) + ')*'
            elif coef == S(1):
                coef_str = ''
            elif isinstance(coef, Rational):
                coef_str = str(coef) + '*'
            s += coef_str + printer.Eprint.Deriv(pd_str) + ' + '

        return s[:-3]

    def Dop_latex_str(self):
        if len(self.coefs) == 0:
            return ' 0 '

        dop_sorted = pdop_sort(self.coefs, self.pdiffs)

        s = ''
        first_flg = False
        pop = ' '
        for (coef, pdop) in dop_sorted:
            pd_str = pdiff_latex_str(pdop, self.Ga.coords)
            if pd_str == '':
                if coef == S(1):
                    if first_flg:
                        pop = ' + '
                    coef_str = '1'
                elif coef == -S(1):
                    pop = ' - '
                    coef_str = '1'
                elif isinstance(coef, Mul) or \
                      isinstance(coef, Rational) or \
                      isinstance(coef, Symbol) or \
                      isinstance(coef, Function) or \
                      isinstance(coef, Pow):
                    coef_str = str(Mv(coef, ga=self.Ga)) + ' '
                    if first_flg:
                        if coef_str[0] != '-':
                            pop = ' + '
                        else:
                            pop = ' '
                elif isinstance(coef, Add):
                    coef_str = str(Mv(coef, ga=self.Ga))
                    if first_flg:
                        if coef_str[0] != '-':
                            pop = ' + '
                        else:
                            pop = ' '
            else:
                if coef == S(1):
                    if first_flg:
                        pop = ' + '
                    coef_str = ''
                elif coef == -S(1):
                    pop = ' - '
                    coef_str = ''
                elif isinstance(coef, Mul) or \
                      isinstance(coef, Rational) or \
                      isinstance(coef, Symbol) or \
                      isinstance(coef, Function) or \
                      isinstance(coef, Pow):
                    coef_str = str(Mv(coef, ga=self.Ga)) + ' '
                    if first_flg:
                        if coef_str[0] != '-':
                            pop = ' + '
                        else:
                            pop = ' '
                elif isinstance(coef, Add):
                    if first_flg:
                        pop = ' + '
                    coef_str = '\\paren{' + str(Mv(coef, ga=self.Ga)) + '}'

            first_flg = True
            s += pop + coef_str + pd_str + ' '
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

    @staticmethod
    def basic(ga):
        r_basis = list(ga.r_basis)

        if not ga.is_ortho:
            r_basis = [x / ga.inorm for x in r_basis]
        if ga.norm:
            r_basis = [x / e_norm for (x, e_norm) in zip(r_basis, ga.e_norm)]

        ga.lgrad = Dop(r_basis, ga.pdx, ga=ga)
        ga.rgrad = Dop(r_basis, ga.pdx, ga=ga, cmpflg=true)
        return ga.lgrad, ga.rgrad

################################# Alan Macdonald's additions #########################


def Nga(x, prec=5):
    if isinstance(x, Mv):
        Px = Mv(x, ga=x.Ga)
        Px.obj = Nsympy(x.obj, prec)
        return(Px)
    else:
        return(Nsympy(x, prec))


def Com(A, B):  # Commutator
    return((A * B - B * A) / S(2))


def rank(M):  # Return rank of matrix M.
    return matrix_rank(M)


def printeigen(M):    # Print eigenvalues, multiplicities, eigenvectors of M.
    evects = M.eigenvects()
    for i in range(len(evects)):                   # i iterates over eigenvalues
        print('Eigenvalue =', evects[i][0], '  Multiplicity =', evects[i][1], ' Eigenvectors:')
        for j in range(len(evects[i][2])):         # j iterates over eigenvectors of a given eigenvalue
            result = '['
            for k in range(len(evects[i][2][j])):  # k iterates over coordinates of an eigenvector
                result += str(trigsimp(evects[i][2][j][k]).evalf(3))
                if k != len(evects[i][2][j]) - 1:
                    result += ', '
            result += '] '
            print(result)


def printGS(M, norm=False):  # Print Gram-Schmidt output.
    from sympy import GramSchmidt
    global N
    N = GramSchmidt(M, norm)
    result = '[ '
    for i in range(len(N)):
        result += '['
        for j in range(len(N[0])):
            result += str(trigsimp(N[i][j]).evalf(3))
            if j != len(N[0]) - 1:
                result += ', '
        result += '] '
        if j != len(N[0]) - 1:
            result += ' '
    result += ']'
    print(result)


def printrref(matrix, vars="xyzuvwrs"):   # Print rref of matrix with variables.
    rrefmatrix = matrix.rref()[0]
    rows, cols = rrefmatrix.shape
    if len(vars) < cols - 1:
        print('Not enough variables.')
        return
    for i in range(rows):
        result = ''
        for j in range(cols - 1):
            result += str(rrefmatrix[i, j]) + vars[j]
            if j != cols - 2:
                result += ' + '
        result += ' = ' + str(rrefmatrix[i, cols - 1])
        print(result)


def correlation(u, v, dec=3):  # Compute the correlation coefficient of vectors u and v.
    rows, cols = u.shape
    uave = 0
    vave = 0
    for i in range(rows):
        uave += u[i]
        vave += v[i]
    uave = uave / rows
    vave = vave / rows
    ulocal = u[:, :]  # Matrix copy
    vlocal = v[:, :]
    for i in range(rows):
        ulocal[i] -= uave
        vlocal[i] -= vave
    return ulocal.dot(vlocal) / (ulocal.norm() * vlocal.norm()). evalf(dec)
