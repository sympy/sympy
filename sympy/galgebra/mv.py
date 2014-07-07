#mv.py

import sys
sys.path.append('./')
import itertools
import copy
import numbers
import operator
from compiler.ast import flatten
from operator import itemgetter, mul, add
from numpy.linalg import matrix_rank
from sympy import Symbol, Function, S, expand, Add, Mul, Pow, Basic, \
    sin, cos, sinh, cosh, sqrt, trigsimp, \
    simplify, diff, Rational, Expr
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

    fmt = 1
    latex_flg = False
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
        Mv.latex_flg = True
        Mv.fmt = mode
        return

    @staticmethod
    def Mul(A, B, op):
        if not isinstance(A, Mv):
            A = B.Ga.mv(A)
        if not isinstance(B, Mv):
            B = A.Ga.mv(B)

        if op == '*':
            return A * B
        elif op == '^':
            return A ^ B
        elif op == '|':
            return A | B
        elif op == '<':
            return A < B
        elif op == '>':
            return A > B
        else:
            raise ValeError('Operation ' + op + 'not allowed in Mv.Mul!')
        return

    def characterise_Mv(self):
        if self.char_Mv:
            return
        obj = self.obj
        if isinstance(obj, numbers.Number):
            self.i_grade = 0
            self.is_blade_rep = True
            self.grades = [0]
            return
        if  obj.is_commutative:
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
                if isinstance(x, Expr):
                   self.obj = x
                else:
                    self.obj = S(x)
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
            if not self.is_scalar():
                return False
            if expand(self.obj) == expand(A):
                return True
            else:
                return False
        if self.is_blade_rep != A.is_blade_rep:
            self = self.blade_rep()
            A = A.blade_rep()
        coefs, bases = metric.linear_expand(self.obj)
        Acoefs, Abases = metric.linear_expand(A.obj)
        if len(bases) != len(Abases):
            return False
        if set(bases) != set(Abases):
            return False
        for base in bases:
            index = bases.index(base)
            indexA = Abases.index(base)
            if expand(coefs[index]) != expand(Acoefs[index]):
                return False
        return True

    def __neg__(self):
        return Mv(-self.obj, ga=self.Ga)

    def __add__(self, A):

        if (not isinstance(A, Mv)) and (not isinstance(A, Dop)):
            return Mv(self.obj + A, ga=self.Ga)

        if self.Ga.name != A.Ga.name:
            raise ValueError('In + operation Mv arguments are not from same geometric algebra')

        if isinstance(A, Dop):
            return Dop.Add(A, self)

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
            return Dop.Add(self, -A)

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
            return A.Mul(self, A, op='*')

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

    def __repr__(self):
        return str(self)

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
                terms[-1] = (grade0, S(1), -1)
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
            return A.Mul(self, A, op='^')

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
            return A.Mul(self, A, op='|')

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
            return A.Mul(self, A, op='<')

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
            return A.Mul(self, A, op='>')

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

    def components(self):
        (coefs, bases) = metric.linear_expand(self.obj)
        bases_lst = self.Ga.blades_lst
        cb = zip(coefs, bases)
        cb = sorted(cb, key=lambda x: self.Ga.blades_lst0.index(x[1]))
        terms = []
        for (coef, base) in cb:
            terms.append(self.Ga.mv(coef * base))
        return terms

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
        if Mv.latex_flg:
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
        else:
            printer.GaPrinter.fmt = self.fmt
            if title is not None:
                print title + ' = ' + str(self)
            else:
                print self
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

################ Scalar Partial Differential Operator Class ############

class Sdop(object):

    init_slots = {'ga': (None, 'Associated geometric algebra'),
                  'fmt': (1, '1 for normal formating')}

    ga = None
    str_mode = False

    @staticmethod
    def setGa(ga):
        Sdop.ga = ga
        Pdop.setGa(ga)
        return

    def TSimplify(self):
        new_terms = []
        for (coef, pdiff) in self.terms:
            new_terms.append((Simp.apply(coef), pdiff))
        self.terms = new_terms
        return

    @staticmethod
    def consolidate_coefs(sdop):
        """
        Remove zero coefs and consolidate coefs with repeated pdiffs.
        """
        if isinstance(sdop, Sdop):
            terms = sdop.terms
        else:
            terms = sdop

        new_coefs = []
        new_pdiffs = []
        for (coef, pd) in terms:
            if coef != S(0):
                if pd in new_pdiffs:
                    index = new_pdiffs.index(pd)
                    new_coefs[index] += coef
                else:
                    new_coefs.append(coef)
                    new_pdiffs.append(pd)
        new_terms = zip(new_coefs, new_pdiffs)

        if isinstance(sdop, Sdop):
            return Sdop(new_terms, ga=sdop.Ga)
        else:
            return new_terms

    def sort_terms(self):
        self.terms.sort(key=operator.itemgetter(1), cmp=Pdop.compare)
        return

    def Sdop_str(self):
        if len(self.terms) == 0:
            return '0'

        self.sort_terms()
        s = ''
        for (coef, pdop) in self.terms:
            pd_str = str(pdop)

            if coef == S(1):
                s += pd_str
            elif coef == S(-1):
                s += '-' + pd_str
            else:
                if isinstance(coef, Add):
                    s += '(' + str(coef) + ')*' + pd_str
                else:
                    s += str(coef) + '*' + pd_str
            s += ' + '

        s = s.replace('+ -','- ')
        s = s[:-3]
        if Sdop.str_mode:
            if len(self.terms) > 1 or isinstance(self.terms[0][0], Add):
                s = '(' + s + ')'
        return s

    def Sdop_latex_str(self):
        if len(self.terms) == 0:
            return '0'

        self.sort_terms()

        s = ''
        for (coef, pdop) in self.terms:
            pd_str = str(pdop)
            if coef == S(1):
                if pd_str == '':
                    s += '1'
                else:
                    s += pd_str
            elif coef == S(-1):
                if pd_str == '':
                    s += '-1'
                else:
                    s += '-' + pd_str
            else:
                s += str(coef) + ' ' + pd_str
            s += ' + '

        s = s.replace('+ -','- ')
        return s[:-3]

    def __str__(self):
        if printer.GaLatexPrinter.latex_flg:
            Printer = printer.GaLatexPrinter
        else:
            Printer = printer.GaPrinter

        return Printer().doprint(self)

    def __repr__(self):
        return str(self)

    def _repr_latex_(self):
        latex_str = printer.GaLatexPrinter.latex(self)
        if r'\begin{align*}' not in latex_str:
            latex_str = r'\begin{equation*} ' + latex_str + r' \end{equation*}'
        return latex_str

    def __init__(self, *kargs, **kwargs):

        kwargs = metric.test_init_slots(Sdop.init_slots, **kwargs)

        self.Ga = kwargs['ga']  # Associated geometric algebra (coords)
        self.fmt = kwargs['fmt']  # Output format

        if self.Ga is None:
            if Sdop.ga is None:
                raise ValueError('In Sdop.__init__ self.Ga must be defined.')
            else:
                self.Ga = Sdop.ga

        if len(kargs[0]) == 0:  # identity Dop
            self.terms = [(S(1), self.Ga.Pdop_identity)]
        else:
            if len(kargs) == 2:
                if len(kargs[0]) != len(kargs[1]):
                    raise ValueError('In Sdop.__init__ coefficent list and Pdop list must be same length.')
                self.terms = zip(kargs[0],kargs[1])
            elif len(kargs) == 1:
                self.terms = kargs[0]
            else:
                raise ValueError('In Sdop.__init__ length of kargs must be 1 or 2.')

    def __call__(self, arg):
        if isinstance(arg, Sdop):
            if self.Ga != arg.Ga:
                raise ValueError('In Sdop.__call__  self.Ga != arg.Ga.')
            terms = []
            for (coef, pdiff) in self.terms:
                new_terms = pdiff(arg.terms)
                new_terms = [ (coef * x[0], x[1]) for x in new_terms]
                terms += new_terms
            return Sdop(terms, ga=self.Ga)
        else:
            return sum([x[0] * x[1](arg) for x in self.terms])


    def __neg__(self):
        return Sdop([(-x[0], x[1]) for x in self.terms], ga=self.Ga)

    @staticmethod
    def Add(sdop1, sdop2):
        if isinstance(sdop1, Sdop) and isinstance(sdop1, Sdop):
            if sdop1.Ga != sdop2.Ga:
                raise ValueError('In Sdop.Add sdop1.Ga != sdop2.Ga.')
            coefs1, pdiffs1 = zip(*sdop1.terms)
            coefs2, pdiffs2 = zip(*sdop2.terms)

            pdiffs1 = list(pdiffs1)
            pdiffs2 = list(pdiffs2)

            pdiffs = pdiffs1 + [x for x in pdiffs2 if x not in pdiffs1]
            coefs = len(pdiffs) * [S(0)]

            for pdiff in pdiffs1:
                index = pdiffs.index(pdiff)
                coef = coefs1[pdiffs1.index(pdiff)]
                coefs[index] += coef

            for pdiff in pdiffs2:
                index = pdiffs.index(pdiff)
                coef = coefs2[pdiffs2.index(pdiff)]
                coefs[index] += coef

            sdop_sum = Sdop(coefs, pdiffs, ga=sdop1.Ga)
        elif isinstance(sdop1, Sdop):
            coefs, pdiffs = zip(*sdop1.terms)
            if sdop1.Ga.Pdop_identity in pdiffs:
                index = pdiffs.index(sdop1.Ga.Pdop_identity)
                coef[index] += sdop2
            else:
                coef.append(sdop2)
                pdiff.append(sdop1.Ga.Pdop_identity)
            return Sdop(coefs, pdiffs, ga=sdop1.Ga)
        else:
            coefs, pdiffs = zip(*sdop2.terms)
            if sdop2.Ga.Pdop_identity in pdiffs:
                index = pdiffs.index(sdop2.Ga.Pdop_identity)
                coef[index] += sdop1
            else:
                coef.append(sdop1)
                pdiff.append(sdop2.Ga.Pdop_identity)
            sdop_sum = Sdop(coefs, pdiffs, ga=sdop2.Ga)

        return Sdop.consolidate_coefs(sdop_sum)

    def __eq__(self, sdop):
        if isinstance(sdop, Sdop):
            if self.Ga != sdop.Ga:
                return False
            self = Sdop.consolidate_coefs(self)
            sdop = Sdop.consolidate_coefs(sdop)
            if len(self.terms) != len(sdop.terms):
                return False
            if set(self.terms) != set(sdop.terms):
                return False
            return True
        else:
            return False

    def __add__(self, sdop):
        return Sdop.Add(self, sdop)

    def __radd__(self, sdop):
        return Sdop(self, sdop)

    def __add_ab__(self, sdop):
        if isinstance(sdop, Sdop):
            if self.Ga != sdop.Ga:
                raise ValueError('In Sdop.__add_ab__ self.Ga != sdop.Ga.')

            coefs, pdiffs = zip(*self.terms)
            pdiffs = list(pdiffs)
            coefs = list(coefs)

            for (coef, pdiff) in sdop.terms:
                if pdiff in pdiffs:
                    index = pdiffs.index(pdiff)
                    coefs[index] += coef
                else:
                    pdiffs.append(pdiff)
                    coefs.append(coef)
            self.term = zip(coefs, pdiffs)
            self = Sdop.consolidate_coefs(self)
            return

        elif isinstance(sdop, tuple):
            self.term.append(sdop)
            self = Dfop.consolidate_coefs(self)
            return

        else:
            self.terms.append((sdop, self.Ga.Pdop_identity))
            self = Sdop.consolidate_coefs(self)
            return

    def __sub__(self, sdop):
        return Sdop.Add(self, -sdop)

    def __rsub__(self, sdop):
        return Sdop.Add(-self, sdop)

    def __mul__(sdopl, sdopr):
        if isinstance(sdopl, Sdop) and isinstance(sdopr, Sdop):
            if sdopl.Ga != sdopr.Ga:
                raise ValueError('In Sdop.__mul__ Sdop arguments are not from same geometric algebra')
            terms = []
            for (coef, pdiff) in sdopl.terms:
                Dsdopl = pdiff(sdopr.terms)  # list of terms
                Dsdopl = [(coef * x[0], x[1]) for x in Dsdopl]
                terms += Dsdopl
            product = Sdop(terms, ga=sdopl.Ga)
            return Sdop.consolidate_coefs(product)
        else:
            if not isinstance(sdopl, Sdop):  # sdopl is a scalar
                terms = [(sdopl * x[0], x[1]) for x in sdopr.terms]
                product = Sdop(terms, ga=sdopr.Ga)  # returns Sdop
                return Sdop.consolidate_coefs(product)
            else:  # sdopr is a scalar or a multivector
                return sum([x[0] * x[1](sdopr) for x in sdopl.terms])  # returns scalar

    def __rmul__(self,sdop):
        terms = [(sdop * x[0], x[1]) for x in self.terms]
        return Sdop(terms, ga=self.Ga)

#################### Partial Derivative Operator Class #################

class Pdop(object):
    """
    Partial derivative class for multivectors.  The partial derivatives
    are of the form

        \partial_{i_{1}...i_{n}} =
            \partial^{i_{1}+...+i_{n}}/\partial{x_{1}^{i_{1}}}...\partial{x_{n}^{i_{n}}}.

    If i_{j} = 0 then the partial derivative does not contain the x^{i_{j}}
    coordinate.

    The partial derivative is represented by a dictionary with coordinates
    for keys and key value are the number of times one differentiates with
    respect to the key.
    """

    ga = None

    init_slots = {'ga': (None, 'Associated geometric algebra')}

    @staticmethod
    def setGa(ga):
        Pdop.ga = ga
        return

    @staticmethod
    def compare(pdop1, pdop2):  # compare two Pdops
        if pdop1.order > pdop2.order:
            return 1
        if pdop1.order < pdop2.order:
            return -1

        keys1 = pdop1.pdiffs.keys()
        keys2 = pdop2.pdiffs.keys()
        lkeys1 = len(keys1)
        lkeys2 = len(keys2)

        if lkeys1 == lkeys2:
            s1 = ''.join([str(pdop1.Ga.coords.index(x)) for x in keys1])
            s2 = ''.join([str(pdop1.Ga.coords.index(x)) for x in keys2])
            if s1 < s2:
                return -1
            else:
                return 1
        else:
            if lkeys1 < lkeys2:
                return 1
            else:
                return -1

    def __eq__(self,A):
        if isinstance(A, Pdop) and self.Ga.name == A.Ga.name and self.pdiffs == A.pdiffs:
            return True
        else:
            if len(self.pdiffs) == 0 and A == S(1):
                return True
            return False

    def __init__(self, *kargs, **kwargs):

        kwargs = metric.test_init_slots(Pdop.init_slots, **kwargs)

        self.Ga = kwargs['ga']  # Associated geometric algebra
        self.fmt = 1
        self.order = 0

        if self.Ga is None:
            if Pdop.ga is None:
                raise ValueError('In Pdop.__init__ self.Ga must be defined.')
            else:
                self.Ga = Pdop.ga  # use geometric algebra of class Pdop

        if kargs[0] is None:  # Pdop is the identity (1)
            self.pdiffs = {}
        elif isinstance(kargs[0], dict):  # Pdop defined by dictionary
            self.pdiffs = kargs[0]
        else:  # Pdop defined by list of integers
            self.pdiffs = {}
            for x in kargs[0]:
                if x != 0:
                    ix = kargs.index[x]
                    self.pdiffs[self.Ga.coords[ix]] = x

        for x in self.pdiffs:  # self.order is total number of differentiations
            self.order += self.pdiffs[x]

    def factor(self):
        """
        If partial derivative operator self.order > 1 factor out first
        order differential operator.  Needed for application of partial
        derivative operator to product of sympy expression and partial
        differential operator.  For example if D = Pdop({x:3}) then

            (Pdop({x:2}),Pdop({x:1})) = D.factor()
        """
        if self.order == 1:
            return S(0), self
        else:
            x = self.pdiffs.keys()[0]
            self.order -= 1
            n = self.pdiffs[x]
            if n == 1:
                del self.pdiffs[x]
            else:
                self.pdiffs[x] -= 1
            return self, self.Ga.Pdiffs[x]

    def __call__(self, arg):

        if self.pdiffs == {}:
            return arg  # result is Pdop identity (1)

        if isinstance(arg, Pdop):  # arg is Pdop
            if self.Ga.name != arg.Ga.name:
                raise ValueError('In Pdop.__call__ arguments do not belong to same geometric algebra.')
            elif arg.pdiffs == {}:  # arg is one
                return S(0)  # derivative is zero
            else:  # arg is partial derivative
                pdiffs = copy.copy(arg.pdiffs)
                for key in self.pdiffs:
                    if key in pdiffs:
                        pdiffs[key] += self.pdiffs[key]
                    else:
                        pdiffs[key] = self.pdiffs[key]
            return Pdop(pdiffs,ga=self.Ga)  # result is Pdop

        elif isinstance(arg, Mv):  # arg is multivector
            for x in self.pdiffs:
                for i in range(self.pdiffs[x]):
                    arg = self.Ga.pDiff(arg, x)
            return arg  # result is multivector

        elif isinstance(arg, (Expr, Symbol, numbers.Number)):  # arg is sympy expression
            for x in self.pdiffs:
                arg = diff(arg,x,self.pdiffs[x])
            return arg  # derivative is sympy expression

        elif isinstance(arg,list):  # arg is list of tuples (coef, partial derivative)
            D = copy.deepcopy(self)
            terms = copy.deepcopy(arg)

            while True:
                D, D0 = D.factor()
                k = 0
                for term in terms:
                    dc = D0(term[0])
                    pd = D0(term[1])
                    tmp = []
                    if dc != 0:
                        tmp.append((dc,term[1]))
                    if pd != 0 :
                        tmp.append((term[0],pd))
                    terms[k] = tmp
                    k += 1
                terms = [i for o in terms for i in o]  # flatten list one level
                if D == 0:
                    break
            terms = Sdop.consolidate_coefs(terms)
            return terms  # result is list of tuples (coef, partial derivative)
        elif isinstance(arg, Sdop):  # arg is scalar differential operator
            if self.Ga != arg.Ga:
                raise ValueError('In Pdop.__call__ self.Ga != arg.Ga.')
            return self(arg.terms)  # result is list of tuples (coef, partial derivative)
        else:
            raise ValueError('In Pdop.__call__ type(arg) = ' + str(type(arg)) + ' not allowed.')

    def __mul__(self, pdop):  # functional product of self and arg (self*arg)
        return self(pdop)

    def __rmul__(self, pdop):  # functional product of arg and self (arg*self)
        if isinstance(pdop, Pdop):
            return pdop(self)
        return Sdop([(pdop, self)], ga=self.Ga)

    def Pdop_str(self):
        if self.order == 0:
            return 'D{}'
        s = 'D'
        for x in self.pdiffs:
            s += '{' + str(x) + '}'
            n = self.pdiffs[x]
            if n > 1:
                s += '^' + str(n)
        return s

    def Pdop_latex_str(self):
        if self.order == 0:
            return ''
        s = r'\bfrac{\partial'
        if self.order > 1:
            s += '^{' + str(self.order) + '}'
        s += '}{'
        keys = self.pdiffs.keys()
        keys.sort(key=(self.Ga.coords + keys).index)
        for key in keys:
            i = self.pdiffs[key]
            s += r'\partial ' + str(key)
            if i > 1:
                s += '^{' + str(i) + '}'
        s += '}'
        return s

    def __str__(self):
        if printer.GaLatexPrinter.latex_flg:
            Printer = printer.GaLatexPrinter
        else:
            Printer = printer.GaPrinter
        return Printer().doprint(self)

    def __repr__(self):
        return str(self)

################# Multivector Differential Operator Class ##############

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
                  'debug': (False, 'True to print out debugging information'),
                  'fmt': (1, '1 for normal formating')}

    ga = None

    @staticmethod
    def setGa(ga):  # set geometric algebra globally for all Dop's
        Dop.ga = ga
        Sdop.setGa(ga)
        return

    @staticmethod
    def flatten_one_level(lst):
        return [inner for outer in lst for inner in outer]
        """
        new_lst = []
        for x in lst:
            if isinstance(x[0], list):
                new_lst += x
            else:
                new_lst.append(x)
        return new_lst
        """

    def __init__(self, *kargs, **kwargs):

        kwargs = metric.test_init_slots(Dop.init_slots, **kwargs)

        self.cmpflg = kwargs['cmpflg']  # Complement flag (default False)
        self.Ga = kwargs['ga']  # Associated geometric algebra

        if self.Ga is None:
            if Dop.ga is None:
                raise ValueError('In Dop.__init__ self.Ga must be defined.')
            else:
                self.Ga = Dop.ga

        self.fmt = kwargs['fmt']  # Output format

        if len(kargs[0]) == 0:  # identity Dop
            self.terms = [(S(1),self.Ga.pdop_identity)]
        else:
            if len(kargs) == 2:
                if len(kargs[0]) != len(kargs[1]):
                    raise ValueError('In Dop.__init__ coefficent list and Pdop list must be same length.')
                self.terms = zip(kargs[0],kargs[1])
            elif len(kargs) == 1:
                if isinstance(kargs[0][0][0], Mv):  # Mv expansion [(Mv, Pdop)]
                    self.terms = kargs[0]
                elif isinstance(kargs[0][0][0], Sdop):  # Sdop expansion [(Sdop, Mv)]
                    coefs = []
                    pdiffs = []
                    for (sdop, mv) in kargs[0]:
                        for (coef, pdiff) in sdop.terms:
                            if pdiff in pdiffs:
                                index = pdiffs.index(pdiff)
                                coefs[index] += coef * mv
                            else:
                                pdiffs.append(pdiff)
                                coefs.append(coef * mv)
                    self.terms = zip(coefs, pdiffs)
                else:
                    raise ValueError('In Dop.__init__ kargs[0] form not allowed.')
            else:
                raise ValueError('In Dop.__init__ length of kargs must be 1 or 2.')

    """
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
    """

    def consolidate_coefs(self):
        """
        Remove zero coefs and consolidate coefs with repeated pdiffs.
        """
        new_coefs = []
        new_pdiffs = []
        for (coef, pd) in self.terms:
            if isinstance(coef, Mv) and coef.is_scalar():
                coef = coef.obj
            if coef != S(0):
                if pd in new_pdiffs:
                    index = new_pdiffs.index(pd)
                    new_coefs[index] += coef
                else:
                    new_coefs.append(coef)
                    new_pdiffs.append(pd)

        self.terms = zip(new_coefs, new_pdiffs)
        #self.TSimplify()
        return

    def blade_rep(self):
        N = len(self.blades)
        coefs = N * [[]]
        bases = N * [0]
        for term in self.terms:
            for (coef, base) in metric.linear_expand(self.terms[0].obj, mode=False):
                index = self.blades.index(base)
                coefs[index] = coef
                bases[index] = base

    @staticmethod
    def Add(dop1, dop2):

        if isinstance(dop1, Dop) and isinstance(dop2, Dop):
            if dop1.Ga.name != dop2.Ga.name:
                raise ValueError('In Dop.Add Dop arguments are not from same geometric algebra')

            if dop1.cmpflg != dop2.cmpflg:
                raise ValueError('In Dop.Add complement flags have different values.')

            coefs1, pdiffs1 = zip(*dop1.terms)
            coefs2, pdiffs2 = zip(*dop2.terms)

            pdiffs1 = list(pdiffs1)
            pdiffs2 = list(pdiffs2)

            pdiffs = pdiffs1 + [x for x in pdiffs2 if x not in pdiffs1]
            coefs = len(pdiffs) * [S(0)]

            for pdiff in pdiffs1:
                index = pdiffs.index(pdiff)
                coef = coefs1[pdiffs1.index(pdiff)]
                coefs[index] += coef

            for pdiff in pdiffs2:
                index = pdiffs.index(pdiff)
                coef = coefs2[pdiffs2.index(pdiff)]
                coefs[index] += coef

            return Dop(coefs, pdiffs, cmpflg=dop1.cmpflg, ga=dop1.Ga)
        else:
            if isinstance(dop1, Dop):  # dop1 is Dop
                if not isinstance(dop2, Mv):
                    dop2 = dop1.Ga.mv(dop2)
                dop2 = Dop([dop2], [dop1.Ga.Pdop_identity], cmpflg=dop1.cmpflg, ga=dop1.Ga)
            else:  # dop2 is Dop
                if not isinstance(dop1, Mv):
                    dop1 = dop2.Ga.mv(dop1)
                dop1 = Dop([dop1], [dop2.Ga.Pdop_identity], cmpflg=dop2.cmpflg, ga=dop2.Ga)
            return Dop.Add(dop1, dop2)

    def __add__(self, dop):
        return Dop.Add(self, dop)

    def __radd__(self, dop):
        return Dop.Add(dop, self)

    def __neg__(self):

        coefs, pdiffs = zip(*self.terms)

        coefs = [-x for x in coefs]

        neg = Dop(coefs, pdiffs, ga=self.Ga,
                  cmpflg=self.cmpflg)

        return neg

    def __sub__(self, dop):
        return Dop.Add(self, -dop)

    def __rsub__(self, dop):
        return Dop.Add(dop, -self)

    @staticmethod
    def Mul(dopl, dopr, op='*'):  # General multiplication of Dop's
        # cmpflg is True if the Dop operates on the left argument and
        # False if the Dop operates on the right argument

        if isinstance(dopl, Dop) and isinstance(dopr, Dop):
            if dopl.Ga != dopr.Ga:
                raise ValueError('In Dop.Mul Dop arguments are not from same geometric algebra')
            if dopl.cmpflg != dopr.cmpflg:
                raise ValueError('In Dop.Mul Dop arguments do not have same cmplfg')
            if not dopl.cmpflg:  # dopl and dopr operate on right argument
                terms = []
                for (coef, pdiff) in dopl.terms:
                    Ddopl = pdiff(dopr.terms)  # list of terms
                    Ddopl = [(Mv.Mul(coef, x[0], op=op), x[1]) for x in Ddopl]
                    terms += Ddopl
                product = Dop(terms, ga=dopl.Ga)
            else:  # dopl and dopr operate on left argument
                terms = []
                for (coef, pdiff) in dopr.terms:
                    Ddopr = pdiff(dopl.terms)  # list of terms
                    Ddopr = [(Mv.Mul(x[0], coef, op=op), x[1]) for x in Ddopr]
                    terms += Ddopr
                product = Dop(terms, ga=dopr.Ga, cmpflg=True)
        else:
            if not isinstance(dopl, Dop):  # dopl is a scalar or Mv and dopr is Dop
                if isinstance(dopl, Mv) and dopl.Ga != dopr.Ga:
                    raise ValueError('In Dop.Mul Dop arguments are not from same geometric algebra')
                else:
                    dopl = dopr.Ga.mv(dopl)

                if not dopr.cmpflg:  # dopr operates on right argument
                    terms = [(Mv.Mul(dopl, x[0], op=op), x[1]) for x in dopr.terms]
                    return Dop(terms, ga=dopr.Ga)  # returns Dop
                else:
                    product = sum([Mv.Mul(x[1](dopl), x[0], op=op) for x in dopr.terms])  # returns multivector
            else:  # dopr is a scalar or a multivector

                if isinstance(dopr, Mv) and dopl.Ga != dopr.Ga:
                    raise ValueError('In Dop.Mul Dop arguments are not from same geometric algebra')

                if not dopl.cmpflg:  # dopl operates on right argument
                    return sum([Mv.Mul(x[0], x[1](dopr), op=op) for x in dopl.terms])  # returns multivector
                else:
                    terms = [(Mv.Mul(x[0], dopr, op=op), x[1]) for x in dopl.terms]
                    product = Dop(terms, ga=dopl.Ga, cmpflg=True)  # returns Dop complement
        if isinstance(product, Dop):
            product.consolidate_coefs()
        return product

    def TSimplify(self):
        new_terms = []
        for (coef, pdiff) in self.terms:
            new_terms.append((metric.Simp.apply(coef), pdiff))
        self.terms = new_terms
        return

    def __mul__(self, dopr):  # * geometric product
        return Dop.Mul(self, dopr, op='*')

    def __rmul__(self, dopl):  # * geometric product
        return Dop.Mul(dopl, self, op='*')

    def __xor__(self, dopr):  # ^ outer product
        return Dop.Mul(self, dopr, op='^')

    def __rxor__(self, dopl):  # ^ outer product
        return Dop.Mul(dopl, self, op='^')

    def __or__(self, dopr):  # | inner product
        return Dop.Mul(self, dopr, op='|')

    def __ror__(self, dopl):  # | inner product
        return Dop.Mul(dopl, self, op='|')

    def __lt__(self, dopr):  # < left contraction
        return Dop.Mul(self, dopr, op='<')

    def __gt__(self, dopr):  # > right contraction
        return Dop.Mul(self, dopr, op='>')

    def __eq__(self, dop):
        if isinstance(dop, Dop):
            if self.Ga != dop.Ga:
                return False
            self = Sdop.consolidate_coefs(self)
            dop = Sdop.consolidate_coefs(dop)
            if len(self.terms) != len(dop.terms):
                return False
            if set(self.terms) != set(dop.terms):
                return False
            return True
        else:
            return False

    def __str__(self):
        if printer.GaLatexPrinter.latex_flg:
            Printer = printer.GaLatexPrinter
        else:
            Printer = printer.GaPrinter

        return Printer().doprint(self)

    def __repr__(self):
        return str(self)

    def _repr_latex_(self):
        latex_str = printer.GaLatexPrinter.latex(self)
        if r'\begin{align*}' not in latex_str:
            latex_str = r'\begin{equation*} ' + latex_str + r' \end{equation*}'
        return latex_str

    def is_scalar(self):
        for x in self.terms:
            if isinstance(x[0], Mv) and not x[0].is_scalar():
                return False
        return True

    def components(self):
        dop_lst = []
        for (sdop, base) in self.Dop_mv_expand():
            new_coefs = []
            new_pdiffs = []
            for (coef, pdiff) in sdop.terms:
                if pdiff in new_pdiffs:
                    index = new_pdiffs.index(pdiff)
                    new_coefs[index] += coef * base
                else:
                    new_pdiffs.append(pdiff)
                    new_coefs.append(coef * base)
            new_coefs = [Mv(x, ga=self.Ga) for x in new_coefs]
            terms = zip(new_coefs, new_pdiffs)
            dop_lst.append(Dop(terms, ga=self.Ga))
        return tuple(dop_lst)

    def Dop_mv_expand(self):
        coefs = []
        bases = []
        self.consolidate_coefs()
        for (coef, pdiff) in self.terms:
            if isinstance(coef, Mv) and not coef.is_scalar():
                mv_terms = metric.linear_expand(coef.obj, mode=False)
                for (mv_coef, mv_base) in mv_terms:
                    if mv_base in bases:
                        index = bases.index(mv_base)
                        coefs[index] += Sdop([(mv_coef, pdiff)], ga=self.Ga)
                    else:
                        bases.append(mv_base)
                        coefs.append(Sdop([(mv_coef, pdiff)], ga=self.Ga))
            else:
                if isinstance(coef, Mv):
                    mv_coef = coef.obj
                else:
                    mv_coef = coef
                if S(1) in bases:
                    index = bases.index(S(1))
                    coefs[index] += Sdop([(mv_coef, pdiff)], ga=self.Ga)
                else:
                    bases.append(S(1))
                    coefs.append(Sdop([(mv_coef, pdiff)], ga=self.Ga))
        terms = zip(coefs, bases)
        return sorted(terms, key=lambda x: self.Ga.blades_lst0.index(x[1]))

    def Dop_str(self):
        if len(self.terms) == 0:
            return ' 0 '

        mv_terms = self.Dop_mv_expand()
        s = ''

        for (sdop, base) in mv_terms:
            str_sdop = str(sdop)
            if base == S(1):
                s += str_sdop
            else:
                if len(sdop.terms) > 1:
                    if self.cmpflg:
                        s += '(' + str_sdop + ')*' + str(base)
                    else:
                        s += str(base) + '*(' + str_sdop + ')'
                else:
                    if str_sdop[0] == '-' and not isinstance(sdop.terms[0][0], Add):
                        if self.cmpflg:
                            s += str_sdop + '*' + str(base)
                        else:
                            s += '-' + str(base) + '*' + str_sdop[1:]
                    else:
                        if self.cmpflg:
                            s += str_dop + '*' + str(base)
                        else:
                            s += str(base) + '*' + str_sdop
            s += ' + '

        s = s.replace('+ -','-')
        return s[:-3]

    def Dop_latex_str(self):
        if len(self.terms) == 0:
            return ' 0 '

        self.consolidate_coefs()

        mv_terms = self.Dop_mv_expand()
        s = ''

        for (sdop, base) in mv_terms:
            str_sdop = str(sdop)
            if base == S(1):
                s += str_sdop
            else:
                if str_sdop == '1':
                    s += str(base)
                if str_sdop == '-1':
                    s += '-' + str(base)
                    if str_sdop[1:] != '1':
                        s += ' ' + str_sdop[1:]
                else:
                    if len(sdop.terms) > 1:
                        if self.cmpflg:
                            s += r'\left ( ' + str_sdop + r'\right ) ' + str(base)
                        else:
                            s += str(base) + ' ' + r'\left ( ' + str_sdop + r'\right ) '
                    else:
                        if str_sdop[0] == '-' and not isinstance(sdop.terms[0][0], Add):
                            if self.cmpflg:
                                s += str_sdop + str(base)
                            else:
                                s += '-' + str(base) + ' ' + str_sdop[1:]
                        else:
                            if self.cmpflg:
                                s += str_sdop + ' ' + str(base)
                            else:
                                s += str(base) + ' ' + str_sdop
            s += ' + '

        s = s.replace('+ -','-')
        Sdop.str_mode = False
        return s[:-3]

    def Fmt(self, fmt=1, title=None):

        self.fmt = fmt
        print 'Mv.latex_flg =', Mv.latex_flg
        if Mv.latex_flg:
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
