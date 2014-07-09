# ga.py

import operator
import copy
from sympy import diff, Rational, Symbol, S, Mul, Pow, Add, \
    collect, expand, eye, trigsimp, sin, cos, sinh, cosh, \
    symbols, sqrt, Abs, numbers
from collections import OrderedDict
#from sympy.core.compatibility import combinations
from itertools import combinations
import printer
import metric
import mv
import lt

half = Rational(1, 2)
one = S(1)
zero = S(0)


def is_bases_product(w):
    nc_w = w.args_cnc()
    nc = nc_w[1]
    return len(nc) == 2 or len(nc) == 1 and nc[0].is_Pow and nc[0].exp == 2


class auto_update_dict(dict):
    """
    auto_update_dict creats entries to a dictionary on the fly.  When the
    dictionary is called and the key used is not one of the existing keys
    The function self.f_update(key) is called to evaluate the key.  The
    result is then added to the dictionary so that self.f_update is not
    used to evaluate the same key again.

    The __init__ is used to input self.f_update for a given dictionary.
    """
    def __init__(self, f_update, instance=None):
        self.f_update = f_update
        self.instance = instance
        self._dict = {}

    def __getitem__(self, key):
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            try:
                if self.instance is None:
                    f_key = self.f_update(key)
                else:
                    f_key = self.f_update(self.instance, key)
                self[key] = f_key
                dict.__setitem__(self, key, f_key)
                return dict.__getitem__(self, key)
            except ValueError:
                raise ValueError('"f_update(' + str(key) + ')" not defined '
                                 + 'for key')


def update_and_substitute(expr1, expr2, func, mul_dict):
    """
    Linear expand expr1 and expr2 to get (summation convention)
        expr1 = coefs1[i]*bases1[i]
        expr2 = coefs2[j]*bases2[j]
    where coefs1 and coefs2 are lists of are commutative expressions and
    bases1 and bases2 are lists of bases for the geometric algebra.
    Then evaluate
        expr = coefs1[i]*coefs2[j]*F(bases1[i],bases2[j])
    where F(bases1[i],bases2[j]) is a function that returns the appropriate
    product of bases1[i]*bases2[j] as a linear combination of scalars and
    bases of the geometric algebra.
    """
    if (isinstance(expr1, numbers.Number) or expr1.is_commutative) \
        or (isinstance(expr2, numbers.Number) or expr2.is_commutative):
        return expr1 * expr2
    (coefs1, bases1) = metric.linear_expand(expr1)
    (coefs2, bases2) = metric.linear_expand(expr2)
    expr = S(0)
    for (coef1, base1) in zip(coefs1, bases1):
        for (coef2, base2) in zip(coefs2, bases2):
            #Special cases where base1 and/or base2 is scalar
            if base1 == 1 and base2 == 1:
                expr += coef1 * coef2
            elif base1 == 1:
                expr += coef1 * coef2 * base2
            elif base2 == 1:
                expr += coef1 * coef2 * base1
            else:
                key = (base1, base2)
                #Update mul dictionary for future
                if key not in mul_dict:
                    mul_dict[key] = func(key)
                expr += coef1 * coef2 * mul_dict[key]
    return expr


def nc_subs(expr, base_keys, base_values=None):
    """
    See if expr contains nc keys in base_keys and substitute corresponding
    value in base_values for nc key.  This was written since standard
    sympy subs was very slow in performing this operation for non-commutative
    keys for long lists of keys.
    """
    if base_values is None:
        [base_keys, base_values] = zip(*base_keys)

    if expr.is_commutative:
        return expr
    if isinstance(expr, Add):
        args = expr.args
    else:
        args = [expr]
    s = zero
    for term in args:
        if term.is_commutative:
            s += term
        else:
            c, nc = term.args_cnc(split_1=False)
            key = Mul._from_args(nc)
            coef = Mul._from_args(c)
            if key in base_keys:
                base = base_values[base_keys.index(key)]
                s += coef * base
            else:
                s += term
    return s


class Ga(metric.Metric):
    """
    The vector space (basis, metric, derivatives of basis vectors) is
    defined by the base class 'Metric'.

    The instantiating the class 'Ga' constructs the geometric algebra of
    the vector space defined by the metric.

    The construction includes the multivector bases, multiplication
    tables or functions for the geometric (*), inner (|), outer (^)
    products, plus the left (<) and right (>) contractions.  The
    geometric derivative operator and any required connections for the
    derivative are also calculated.

    Except for the geometric product in the case of a non-orthogonal
    set of basis vectors all products and connections (if needed) are
    calculated when needed and place in dictionaries (lists of tuples)
    to be used when needed.  This greatly speeds up evaluations of
    multivector expressions over previous versions of this code since
    the products of multivector bases and connection are not calculated
    unless they are actually needed in the current calculation.

    Only instantiate the Ga class via the Mv class or any use of enhanced
    printing (text or latex) will cause the bases and multiplication
    table entries to be incorrectly labeled .
    """

    restore = False

    a = []

    presets = {'o3d': 'x,y,z:[1,1,1]:[1,1,0]',
               'cyl3d': 'r,theta,z:[1,r**2,1]:[1,1,0]:norm=True',
               'sph3d': 'r,theta,phi:[1,X[0]**2,X[0]**2*cos(X[1])**2]:[1,1,0]:norm=True',
               'para3d': 'u,v,z:[u**2+v**2,u**2+v**2,1]:[1,1,0]:norm=True'}

    @staticmethod
    def build(*kargs, **kwargs):
        """
        Static method to instanciate geometric algebra and return geometric
        algebra, basis vectors, and grad operator as a tuple.
        """
        GA = Ga(*kargs, **kwargs)
        basis = list(GA.mv())
        return tuple([GA] + basis)

    @staticmethod
    def preset(setting, root='e', debug=False):

        if setting not in Ga.presets:
            raise ValueError(str(setting) + 'not in Ga.presets.')
        set_lst = Ga.presets[setting].split(':')
        X = symbols(set_lst[0], real=True)
        g = eval(set_lst[1])
        simps = eval(set_lst[2])
        kargs = [root]
        kwargs = {'g': g, 'coords': X, 'debug': debug, 'I': True, 'gsym': False}

        if len(set_lst) > 3:
            args_lst = set_lst[-1].split(';')
            for arg in args_lst:
                [name, value] = arg.split('=')
                kwargs[name] = eval(value)

        Ga.set_simp(*simps)
        return Ga(*kargs, **kwargs)

    def __eq__(self, ga):
        if self.name == self.ga:
            return True
        return False

    def __init__(self, bases, **kwargs):

        # Each time a geometric algebra is intialized in setup of append
        # the printer must be restored to the simple text mode (not
        # enhanced text of latex printing) so that when 'str' is used to
        # create symbol names the names are not mangled.

        kwargs = metric.test_init_slots(metric.Metric.init_slots, **kwargs)

        if printer.GaLatexPrinter.latex_flg:
            printer.GaLatexPrinter.restore()
            Ga.restore = True

        metric.Metric.__init__(self, bases, **kwargs)

        self.par_coords = None
        self.build_bases()
        self.dot_mode = '|'
        self.basis_product_tables()

        self.inorm = None
        self.Inorm = None

        if self.coords is not None:
            self.coords = list(self.coords)

        if self.coords is not None:
            self.coord_vec = sum([coord * base for (coord, base) in zip(self.coords, self.basis)])
            self.build_reciprocal_basis(kwargs['gsym'])
            self.Pdop_identity = mv.Pdop({},ga=self)  # Identity Pdop = 1
            self.Pdiffs = {}
            self.sPds = {}
            for x in self.coords:  # Partial derivative operator for each coordinate
                self.Pdiffs[x] = mv.Pdop({x:1}, ga=self)
                self.sPds[x] = mv.Sdop([(S(1), self.Pdiffs[x])], ga=self)
            self.grad, self.rgrad = self.grads()
            self.TSgrads =[]

        if self.connect_flg:
            self.build_connection()

        self.lt_flg = False
        self.i = mv.Mv(self.iobj, ga=self)

        if Ga.restore:  # restore printer to appropriate enhanced mode after ga is instanciated
            printer.GaLatexPrinter.redirect()

        if self.coords is not None:
            self.grads()

        if self.debug:
            print 'Exit Ga.__init__()'

        self.a = []  # List of dummy vectors for Mlt calculations
        self.acoefs = []  # List of dummy vectors coefficients
        self.pdiffs = []  # List of lists dummy vector coefficients
        self.dslot = -1  # kargs slot for dervative, -1 for coordinates

    def make_grad(self, a, cmpflg=False):  # make gradient operator with respect to vector a
        if isinstance(a, mv.Mv):
            ai = a.get_coefs(1)
        else:
            ai = a
        coefs = []
        pdiffs = []
        for (base, coord) in zip(self.r_basis_mv, ai):
            coefs.append(base)
            pdiffs.append(mv.Pdop({coord: 1}, ga=self))
        return mv.Dop(coefs, pdiffs, ga=self, cmpflg=cmpflg)

    def __str__(self):
        return self.name

    def I(self):
        if self.Inorm is None:
            if self.inorm is None:
                self.inorm = (self.i * self.i).scalar()
            self.Inorm = self.i / sqrt(Abs(self.inorm))
        return self.Inorm

    def X(self):
        return self.mv(sum([coord*base for (coord, base) in zip(self.coords, self.basis)]))

    def sdop(self, coefs, pdiffs=None):
        return metric.Sdop(coefs, pdiffs, ga=self)

    def mv(self, root=None, *kargs, **kwargs):
        """
        Instanciate and return a multivector for this, 'self',
        geometric algebra.
        """
        (self.mv_I, self.mv_basis, self.mv_x) = mv.Mv.setup(ga=self)

        if root is None:  # Return ga basis and compute grad and rgrad
            if self.coords is not None:
                self.grads()
            return self.mv_basis

        kwargs['ga'] = self

        if not isinstance(root,str):
            return mv.Mv(root, *kargs, **kwargs)

        if ' ' in root and ' ' not in kargs[0]:
            root_lst = root.split(' ')
            mv_lst = []
            for root in root_lst:
                mv_lst.append(mv.Mv(root, *kargs, **kwargs))
            return tuple(mv_lst)

        if ' ' in root and ' ' in kargs[0]:
            root_lst = root.split(' ')
            mvtype_lst = kargs[0].split(' ')
            if len(root_lst) != len(mvtype_lst):
                raise ValueError('In Ga.mv() for multiple multivectors and ' +
                                  'multivector types incompatible args ' +
                                  str(root_lst) + ' and ' + str(mvtype_lst))

            mv_lst = []
            for (root, mv_type) in zip(root_lst, mvtype_lst):
                kargs = list(kargs)
                kargs[0] = mv_type
                kargs = tuple(kargs)
                mv_lst.append(mv.Mv(root, *kargs, **kwargs))
            return tuple(mv_lst)

        return mv.Mv(root, *kargs, **kwargs)

    def grads(self):
        if not self.is_ortho:
            r_basis = [x / self.inorm for x in self.r_basis_mv]
        else:
            r_basis = self.r_basis_mv
        if self.norm:
            r_basis = [x / e_norm for (x, e_norm) in zip(self.r_basis_mv, self.e_norm)]

        pdx = [self.Pdiffs[x] for x in self.coords]

        self.grad = mv.Dop(r_basis, pdx, ga=self)
        self.rgrad = mv.Dop(r_basis, pdx, ga=self, cmpflg=True)
        return self.grad, self.rgrad

    def dop(self, *kargs, **kwargs):
        """
        Instanciate and return a multivector differential operator for
        this, 'self', geometric algebra.
        """
        kwargs['ga'] = self
        return mv.Dop(*kargs, **kwargs)

    def lt(self, *kargs, **kwargs):
        """
        Instanciate and return a linear transformation for this, 'self',
        geometric algebra.
        """
        if not self.lt_flg:
            self.lt_flg = True
            (self.lt_coords, self.lt_x) = lt.Lt.setup(ga=self)

        kwargs['ga'] = self
        return lt.Lt(*kargs, **kwargs)

    def sm(self, *kargs, **kwargs):
        """
        Instanciate and return a submanifold for this, 'self',
        geometric algebra.  See 'Sm' class for instanciation inputs.
        """
        kwargs['ga'] = self
        SM = Sm(*kargs, **kwargs)
        return SM

    def parametric(self, coords):
        if not isinstance(coords, list):
            raise TypeError('In Ga.parametric coords = ' + str(coords) +
                             ' is not a list.')
        if len(coords) != self.n:
            raise ValueError('In Ga.parametric number of parametric functions' +
                              ' not equal to number of coordinates.')

        self.par_coords = {}

        for (coord, par_coord) in zip(self.coords, coords):
            self.par_coords[coord] = par_coord
        return

    def basis_vectors(self):
        return tuple(self.basis)

    def build_bases(self):
        """
        The bases for the multivector (geometric) algebra are formed from
        all combinations of the bases of the vector space and the scalars.

        Each base is represented as a non-commutative symbol of the form -

                    e_{i_{1}}e_{i_{2}}...e_{i_{r}}

        where 0 < i_{1} < i_{2} < ... < i_{r} and 0 < r <= n the
        dimension of the vector space and 0 < i_{j} <= n. The total
        number of all symbols of this form plus the scalars is 2^{n}.
        Any multivector can be represented as a linear combination
        of these bases and the scalars.

        If the basis vectors are not orthogonal a second set of symbols
        is required given by -

                    e_{i_{1}}^e_{i_{2}}^...^e_{i_{r}}.

        These are called the blade basis for the geometric algebra and
        and multivector can also be represented by a linears combination
        of these blades and the scalars.  The number of basis vectors
        that are in the symbol for the blade is call the grade of the
        blade.

        Representing the multivector as a linear combination of blades
        gives a blade decomposition of the multivector.

        There is a linear mapping from bases to blades and blades to
        bases so that one can easily convert from one representation to
        another.  For the case of an orthogonal set of basis vectors the
        bases and blades are identical.
        """

        # index list for multivector bases and blades by grade
        basis_indexes = tuple(self.n_range)
        self.indexes = [()]
        self.indexes_lst = []
        for i in basis_indexes:
            base_tuple = tuple(combinations(basis_indexes, i + 1))
            self.indexes.append(base_tuple)
            self.indexes_lst += list(base_tuple)
        self.indexes = tuple(self.indexes)

        # list of non-commutative symbols for multivector bases and blades
        # by grade and as a flattened list

        self.blades = []
        self.blades_lst = []
        for grade_index in self.indexes:
            blades = []
            super_scripts = []
            for base_index in grade_index:
                symbol_str = (''.join([str(self.basis[i]) + '^' for i in base_index]))[:-1]
                blade_symbol = Symbol(symbol_str, commutative=False)
                blades.append(blade_symbol)
                self.blades_lst.append(blade_symbol)
            self.blades.append(blades)

        self.blades_lst0 = [S(1)] + self.blades_lst

        self.iobj = self.blades_lst[-1]

        self.blades_to_indexes = []
        self.indexes_to_blades = []
        for (index, blade) in zip(self.indexes_lst, self.blades_lst):
            self.blades_to_indexes.append((blade, index))
            self.indexes_to_blades.append((index, blade))
        self.blades_to_indexes_dict = OrderedDict(self.blades_to_indexes)
        self.indexes_to_blades_dict = OrderedDict(self.indexes_to_blades)

        self.blades_to_grades_dict = {}
        igrade = 0
        for grade in self.blades:
            for blade in grade:
                self.blades_to_grades_dict[blade] = igrade
            igrade += 1

        if not self.is_ortho:

            self.bases = []
            self.bases_lst = []
            for grade_index in self.indexes:
                bases = []
                for base_index in grade_index:
                    symbol_str = (''.join([str(self.basis[i]) + '*' for i in base_index]))[:-1]
                    base_symbol = Symbol(symbol_str, commutative=False)
                    bases.append(base_symbol)
                    self.bases_lst.append(base_symbol)
                self.bases.append(bases)

            self.pseudoI = self.bases_lst[-1]

            self.bases_to_indexes = []
            self.indexes_to_bases = []
            for (index, base) in zip(self.indexes_lst, self.bases_lst):
                self.bases_to_indexes.append((base, index))
                self.indexes_to_bases.append((index, base))
            self.bases_to_indexes_dict = OrderedDict(self.bases_to_indexes)
            self.indexes_to_bases_dict = OrderedDict(self.indexes_to_bases)

            self.bases_to_grades_dict = {}
            igrade = 0
            for grade in self.bases:
                for base in grade:
                    self.bases_to_grades_dict[base] = igrade
                igrade += 1

        if self.coords is None:
            base0 = str(self.basis[0])
            if '_' in base0:
                sub_index = base0.index('_')
                self.basis_super_scripts = [str(base)[sub_index + 1:] for base in self.basis]
            else:
                self.basis_super_scripts = [str(i + 1) for i in self.n_range]
        else:
            self.basis_super_scripts = [str(coord) for coord in self.coords]

        self.blade_super_scripts = []

        for grade_index in self.indexes:
            super_scripts = []
            for base_index in grade_index:
                super_scripts.append(''.join([self.basis_super_scripts[i]
                                     for i in base_index]))
            self.blade_super_scripts.append(super_scripts)

        if self.debug:
            printer.oprint('indexes', self.indexes, 'list(indexes)', self.indexes_lst,
                            'blades', self.blades, 'list(blades)', self.blades_lst,
                            'blades_to_indexes_dict', self.blades_to_indexes_dict,
                            'indexes_to_blades_dict', self.indexes_to_blades_dict,
                            'blades_to_grades_dict', self.blades_to_grades_dict,
                            'blade_super_scripts', self.blade_super_scripts)
            if not self.is_ortho:
                printer.oprint('bases', self.bases, 'list(bases)', self.bases_lst,
                                'bases_to_indexes_dict', self.bases_to_indexes_dict,
                                'indexes_to_bases_dict', self.indexes_to_bases_dict,
                                'bases_to_grades_dict', self.bases_to_grades_dict)
        return

    def basis_product_tables(self):
        """
        For the different products of geometric algebra bases/blade
        initialize auto-updating of bases/blades product lists.  For
        orthogonal bases all basis product lists are generated on the
        fly using functions and the base and blade representations
        are identical.  For a non-orthogonal basis the multiplication
        table for the geometric product is pre-calcuated for base pairs.
        The tables for all other products (including the geometric
        product) are calulated on the fly and updated and are for blade
        pairs.

        All tables are of the form

            [ (blade1*blade2,f(blade1,blade1)),... ]
        """

        self.mul_table = []  # Geometric product (*) of blades
        self.mul_table_dict = {}

        if not self.is_ortho:
            self.non_orthogonal_mul_table()  # Fully populated geometric product (*) multiplication table
            self.base_blade_conversions()  # Generates conversion dictionaries between bases and blades

        self.wedge_table = []  # Outer product (^)
        self.wedge_table_dict = {}

        # All three (|,<,>) types of contractions use the same generation function
        # self.dot_product_basis_blades.  The type of dictionary entry generated depend
        # on self.dot_mode = '|', '<', or '>' as set in self.dot.

        self.dot_table = []  # Inner product (|)
        self.dot_table_dict = {}
        self.left_contract_table = []  # Left contraction (<)
        self.left_contract_table_dict = {}
        self.right_contract_table = []  # Right contraction (>)
        self.right_contract_table_dict = {}

        self.dot_mode = '|'
        if self.debug:
            print 'Exit basis_product_tables.\n'
        return

    def build_connection(self):
        # Partial derivatives of multivector bases multiplied (*,^,|,<,>)
        # on left and right (True and False) by reciprocal basis vectors.
        self.connect = {('*', True): [], ('^', True): [], ('|', True): [],
                        ('<', True): [], ('>', True): [], ('*', False): [],
                        ('^', False): [], ('|', False): [], ('<', False): [],
                        ('>', False): []}
        # Partial derivatives of multivector bases
        self.dbases = {}

        return

    ######## Functions for Calculation products of blades/bases ########

    #******************** Geometric Product (*) ***********************#

    def geometric_product_basis_blades(self, blade12):
        # geometric (*) product for orthogonal basis
        if self.is_ortho:
            (blade1, blade2) = blade12
            index1 = self.blades_to_indexes_dict[blade1]
            index2 = self.blades_to_indexes_dict[blade2]
            blade_index = list(index1 + index2)
            repeats = []
            sgn = 1
            for i in range(1, len(blade_index)):
                save = blade_index[i]
                j = i
                while j > 0 and blade_index[j - 1] > save:
                    sgn = -sgn
                    blade_index[j] = blade_index[j - 1]
                    j -= 1
                blade_index[j] = save
                if blade_index[j] == blade_index[j - 1]:
                    repeats.append(save)
            result = S(sgn)
            for i in repeats:
                blade_index.remove(i)
                blade_index.remove(i)
                result *= self.g[i, i]
            if len(blade_index) > 0:
                result *= self.indexes_to_blades_dict[tuple(blade_index)]
            return result
        else:
            (blade1, blade2) = blade12
            base1 = self.blade_to_base_rep(blade1)
            base2 = self.blade_to_base_rep(blade2)
            base12 = expand(base1 * base2)
            base12 = nc_subs(base12, self.basic_mul_keys, self.basic_mul_values)
            return self.base_to_blade_rep(base12)

    def reduce_basis(self, blst):
        """
        Repetitively applies reduce_basis_loop to blst
        product representation until normal form is
        realized for non-orthogonal basis
        """
        blst = list(blst)
        if blst == []:  # blst represents scalar
            blst_coef = [1]
            blst_expand = [[]]
            return blst_coef, blst_expand
        blst_expand = [blst]
        blst_coef = [1]
        blst_flg = [False]
        # reduce untill all blst revise flgs are True
        while not reduce(operator.and_, blst_flg):
            for i in range(len(blst_flg)):
                if not blst_flg[i]:  # keep revising if revise flg is False
                    tmp = Ga.reduce_basis_loop(self.g, blst_expand[i])
                    if isinstance(tmp, bool):
                        blst_flg[i] = tmp  # revision of blst_expand[i] complete
                    elif len(tmp) == 3:  # blst_expand[i] contracted
                        blst_coef[i] = tmp[0] * blst_coef[i]
                        blst_expand[i] = tmp[1]
                        blst_flg[i] = tmp[2]
                    else:  # blst_expand[i] revised
                        blst_coef[i] = -blst_coef[i]
                        #if revision force one more pass in case revision
                        #causes repeated index previous to revised pair of
                        #indexes
                        blst_flg[i] = False
                        blst_expand[i] = tmp[3]
                        blst_coef.append(-blst_coef[i] * tmp[0])
                        blst_expand.append(tmp[1])
                        blst_flg.append(tmp[2])
        new_blst_coef = []
        new_blst_expand = []
        for (coef, xpand) in zip(blst_coef, blst_expand):
            if xpand in new_blst_expand:
                i = new_blst_expand.index(xpand)
                new_blst_coef[i] += coef
            else:
                new_blst_expand.append(xpand)
                new_blst_coef.append(coef)
        return new_blst_coef, new_blst_expand

    @staticmethod
    def reduce_basis_loop(g, blst):
        """
        blst is a list of integers [i_{1},...,i_{r}] representing the geometric
        product of r basis vectors a_{{i_1}}*...*a_{{i_r}}.  reduce_basis_loop
        searches along the list [i_{1},...,i_{r}] untill it finds i_{j} == i_{j+1}
        and in this case contracts the list, or if i_{j} > i_{j+1} it revises
        the list (~i_{j} means remove i_{j} from the list)

        Case 1: If i_{j} == i_{j+1}, return a_{i_{j}}**2 and
                [i_{1},..,~i_{j},~i_{j+1},...,i_{r}]

        Case 2: If i_{j} > i_{j+1}, return a_{i_{j}}.a_{i_{j+1}},
                [i_{1},..,~i_{j},~i_{j+1},...,i_{r}], and
                [i_{1},..,i_{j+1},i_{j},...,i_{r}]
        """
        nblst = len(blst)  # number of basis vectors
        if nblst <= 1:
            return True  # a scalar or vector is already reduced
        jstep = 1
        while jstep < nblst:
            istep = jstep - 1
            if blst[istep] == blst[jstep]:  # basis vectorindex is repeated
                i = blst[istep]  # save basis vector index
                if len(blst) > 2:
                    blst = blst[:istep] + blst[jstep + 1:]  # contract blst
                else:
                    blst = []
                if len(blst) <= 1 or jstep == nblst - 1:
                    blst_flg = True  # revision of blst is complete
                else:
                    blst_flg = False  # more revision needed
                return g[i, i], blst, blst_flg
            if blst[istep] > blst[jstep]:  # blst not in normal order
                blst1 = blst[:istep] + blst[jstep + 1:]  # contract blst
                a1 = 2 * g[blst[jstep], blst[istep]]  # coef of contraction
                blst = blst[:istep] + [blst[jstep]] + [blst[istep]] + blst[jstep + 1:]  # revise blst
                if len(blst1) <= 1:
                    blst1_flg = True  # revision of blst is complete
                else:
                    blst1_flg = False  # more revision needed
                return a1, blst1, blst1_flg, blst
            jstep += 1
        return True  # revision complete, blst in normal order

    #******************* Outer/wedge (^) product **********************#

    @staticmethod
    def blade_reduce(lst):
        sgn = 1
        for i in range(1, len(lst)):
            save = lst[i]
            j = i
            while j > 0 and lst[j - 1] > save:
                sgn = -sgn
                lst[j] = lst[j - 1]
                j -= 1
            lst[j] = save
            if lst[j] == lst[j - 1]:
                return 0, None
        return sgn, lst

    def wedge_product_basis_blades(self, blade12):  # blade12 = blade1*blade2
        # outer (^) product of basis blades
        # this method works for both orthogonal and non-orthogonal basis
        (blade1, blade2) = blade12
        index1 = self.blades_to_indexes_dict[blade1]
        index2 = self.blades_to_indexes_dict[blade2]
        index12 = list(index1 + index2)

        if len(index12) > self.n:
            return 0
        (sgn, wedge12) = Ga.blade_reduce(index12)
        if sgn != 0:
            return(sgn * self.indexes_to_blades_dict[tuple(wedge12)])
        else:
            return 0

    #****** Dot (|) product, reft (<) and right (>) contractions ******#

    def dot_product_basis_blades(self, blade12):
        # dot (|), left (<), and right (>) products
        # dot product for orthogonal basis
        (blade1, blade2) = blade12
        index1 = self.blades_to_indexes_dict[blade1]
        index2 = self.blades_to_indexes_dict[blade2]
        index = list(index1 + index2)
        grade1 = len(index1)
        grade2 = len(index2)

        if self.dot_mode == '|':
            grade = abs(grade1 - grade2)
        elif self.dot_mode == '<':
            grade = grade2 - grade1
            if grade < 0:
                return 0
        elif self.dot_mode == '>':
            grade = grade1 - grade2
            if grade < 0:
                return 0
        n = len(index)
        sgn = 1
        result = 1
        ordered = False
        while n > grade:
            ordered = True
            i2 = 1
            while i2 < n:
                i1 = i2 - 1
                index1 = index[i1]
                index2 = index[i2]
                if index1 == index2:
                    n -= 2
                    if n < grade:
                        return 0
                    result *= self.g[index1, index1]
                    index = index[:i1] + index[i2 + 1:]
                elif index1 > index2:
                    ordered = False
                    index[i1] = index2
                    index[i2] = index1
                    sgn = -sgn
                    i2 += 1
                else:
                    i2 += 1
            if ordered:
                break
        if n > grade:
            return 0
        else:
            if index == []:
                return sgn * result
            else:
                return sgn * result * self.indexes_to_blades_dict[tuple(index)]

    def non_orthogonal_dot_product_basis_blades(self, blade12):  # blade12 = (blade1,blade2)
        # dot product of basis blades if basis vectors are non-orthogonal
        # inner (|), left (<), and right (>) products of basis blades
        # blade12 is the sympy product of two basis blades
        (blade1, blade2) = blade12
        # Need base rep for blades since that is all we can multiply
        base1 = self.blade_expansion_dict[blade1]
        base2 = self.blade_expansion_dict[blade2]
        # geometric product of basis blades
        base12 = self.mul(base1, base2)
        # blade rep of geometric product
        blade12 = self.base_to_blade_rep(base12)
        # decompose geometric product by grades
        grade_dict = self.grade_decomposition(blade12)
        # grades of input blades
        grade1 = self.blades_to_grades_dict[blade1]
        grade2 = self.blades_to_grades_dict[blade2]
        if self.dot_mode == '|':
            grade_dot = abs(grade2 - grade1)
            if grade_dot in grade_dict:
                return grade_dict[grade_dot]
            else:
                return zero
        elif self.dot_mode == '<':
            grade_contract = grade2 - grade1
            if grade_contract in grade_dict:
                return grade_dict[grade_contract]
            else:
                return zero
        elif self.dot_mode == '>':
            grade_contract = grade1 - grade2
            if grade_contract in grade_dict:
                return grade_dict[grade_contract]
            else:
                return zero
        else:
            raise ValueError('"' + str(self.dot_mode) + '" not allowed '
                             'dot mode in non_orthogonal_dot_basis')

    ############# Non-Orthogonal Tables and Dictionaries ###############

    def non_orthogonal_mul_table(self):
        mul_table = []
        self.basic_mul_keys = []
        self.basic_mul_values = []
        for base1 in self.bases_lst:
            for base2 in self.bases_lst:
                key = base1 * base2
                value = self.non_orthogonal_bases_products((base1, base2))
                mul_table.append((key, value))
                self.basic_mul_keys.append(key)
                self.basic_mul_values.append(value)
        self.basic_mul_table = mul_table
        self.basic_mul_table_dict = OrderedDict(mul_table)

        if self.debug:
            print 'basic_mul_table =\n', self.basic_mul_table
        return

    def non_orthogonal_bases_products(self, base12):  # base12 = (base1,base2)
        # geometric product of bases for non-orthogonal basis vectors
        (base1, base2) = base12
        index = self.bases_to_indexes_dict[base1] + self.bases_to_indexes_dict[base2]

        (coefs, indexes) = self.reduce_basis(index)

        s = 0

        if [] in indexes:  # extract scalar part from multivector expansion
            iscalar = indexes.index([])
            s += coefs[iscalar]
            del indexes[iscalar]
            del coefs[iscalar]

        for (coef, index) in zip(coefs, indexes):
            s += coef * self.indexes_to_bases_dict[tuple(index)]

        return s

    def base_blade_conversions(self):

        blade_expansion = []
        blade_index = []

        # expand blade basis in terms of base basis
        for blade in self.blades_lst:
            index = self.blades_to_indexes_dict[blade]
            grade = len(index)
            if grade == 1:
                blade_expansion.append(blade)
                blade_index.append(index)
            else:
                a = self.indexes_to_blades_dict[(index[0],)]
                Aexpand = blade_expansion[blade_index.index(index[1:])]
                # Formula for outer (^) product of a vector and grade-r multivector
                # a^A_{r} = (a*A + (-1)^{r}*A*a)/2
                # The folowing evaluation takes the most time for setup it is the due to
                # the substitution required for the multiplications
                a_W_A = half * (self.basic_mul(a, Aexpand) - ((-1) ** grade) * self.basic_mul(Aexpand, a))
                blade_index.append(index)
                blade_expansion.append(expand(a_W_A))

        self.blade_expansion = blade_expansion
        self.blade_expansion_dict = OrderedDict(zip(self.blades_lst, blade_expansion))

        if self.debug:
            print 'blade_expansion_dict =', self.blade_expansion_dict

        # expand base basis in terms of blade basis

        base_expand = []

        for (base, blade, index) in zip(self.bases_lst, self.blades_lst, self.indexes_lst):
            grade = len(index)
            if grade == 1:
                base_expand.append((base, base))
            else:  # back substitution of tridiagonal system
                tmp = self.blade_expansion_dict[blade]
                tmp = tmp.subs(base, -blade)
                tmp = -tmp.subs(base_expand)
                base_expand.append((base, expand(tmp)))

        self.base_expand = base_expand
        self.base_expansion_dict = OrderedDict(base_expand)

        if self.debug:
            print 'base_expansion_dict =', self.base_expansion_dict

        return

    def base_to_blade_rep(self, A):

        if self.is_ortho:
            return A
        else:
            #return(expand(A).subs(self.base_expansion_dict))
            return nc_subs(expand(A), self.base_expand)

    def blade_to_base_rep(self, A):

        if self.is_ortho:
            return A
        else:
            #return(expand(A).subs(self.blade_expansion_dict))
            return nc_subs(expand(A), self.blades_lst, self.blade_expansion)

    ###### Products (*,^,|,<,>) for multivector representations ########

    def basic_mul(self, A, B):  # geometric product (*) of base representations
        # only multiplicative operation to assume A and B are in base representation
        AxB = expand(A * B)
        AxB = nc_subs(AxB, self.basic_mul_keys, self.basic_mul_values)
        return expand(AxB)

    def Mul(self, A, B, mode='*'):  # Unifies all products into one function
        if mode == '*':
            return self.mul(A, B)
        elif mode == '^':
            return self.wedge(A, B)
        else:
            self.dot_mode = mode
            return self.dot(A, B)

    def mul(self, A, B):  # geometric (*) product of blade representations
        if A == 0 or B == 0:
            return 0
        return update_and_substitute(A, B, self.geometric_product_basis_blades, self.mul_table_dict)

    def wedge(self, A, B):
        # wedge assumes A and B are in blade rep
        # wedge product is same for both orthogonal and non-orthogonal for A and B in blade rep
        if A == 0 or B == 0:
            return 0
        return update_and_substitute(A, B, self.wedge_product_basis_blades, self.wedge_table_dict)

    def dot(self, A, B):  # inner products |, <, and >
        if A == 0 or B == 0:
            return 0
        if self.is_ortho:
            dot_product_basis_blades = self.dot_product_basis_blades
        else:
            dot_product_basis_blades = self.non_orthogonal_dot_product_basis_blades

        if self.dot_mode == '|':  # Hestenes dot product
            A = self.remove_scalar_part(A)
            B = self.remove_scalar_part(B)
            return update_and_substitute(A, B, dot_product_basis_blades, self.dot_table_dict)
        else:
            if self.dot_mode == '<':  # Left contraction
                return update_and_substitute(A, B, dot_product_basis_blades, self.left_contract_table_dict)
            elif self.dot_mode == '>':  # Right contraction
                return update_and_substitute(A, B, dot_product_basis_blades, self.right_contract_table_dict)
            else:
                raise ValueError('"' + str(self.dot_mode) + '" not a legal mode in dot')

    ######################## Helper Functions ##########################

    def grade_decomposition(self, A):
        """
        Returns dictionary with grades as keys of grades of A.  For example
        if A is a rotor the dictionary keys would be 0 and 2. For a vector
        the single key would be 1
        """
        coefs = collect(expand(A), self.blades_lst, evaluate=False)
        grade_dict = {0: zero}
        for blade in coefs:
            if blade.is_commutative:
                grade_dict[0] += coefs[one]
            else:
                grade = self.blades_to_grades_dict[blade]
                if grade in grade_dict:
                    grade_dict[grade] += coefs[blade] * blade
                else:
                    grade_dict[grade] = coefs[blade] * blade
        return grade_dict

    def remove_scalar_part(self, A):
        """
        Return non-commutative part (sympy object) of A.obj.
        """

        if isinstance(A, mv.Mv):
            return self.remove_scalar_part(A.obj)
        else:
            if isinstance(A, Add):
                A = expand(A)
                return(sum([x for x in A.args if not x.is_commutative]))
            elif isinstance(A, Symbol):
                if A.is_commutative:
                    return 0
                else:
                    return A
            else:
                if A.is_commutative:
                    return 0
                else:
                    return A

    def scalar_part(self, A):
        A = expand(A)
        if isinstance(A, Add):
            return(sum([x for x in A.args if x.is_commutative]))
        else:
            if A.is_commutative:
                return A
            else:
                return zero

    def grades(self, A):  # Return list of grades present in A
        A = self.blade_to_base_rep(A)
        A = expand(A)
        blades = []
        if isinstance(A, Add):
            args = A.args
        else:
            args = [A]
        for term in args:
            blade = term.args_cnc()[1]
            l_blade = len(blade)
            if l_blade > 0:
                if blade[0] not in blades:
                    blades.append(blade[0])
            else:
                if one not in blades:
                    blades.append(one)
        grade_lst = []
        if one in blades:
            grade_lst.append(0)
        for blade in blades:
            if blade != one:
                grade = self.blades_to_grades_dict[blade]
                if grade not in grade_lst:
                    grade_lst.append(grade)
        grade_lst.sort()
        return(grade_lst)

    def reverse(self, A):  # Calculates reverse of A (see documentation)
        A = expand(A)
        blades = {}
        if isinstance(A, Add):
            args = A.args
        else:
            if A.is_commutative:
                return A
            else:
                args = [A]
        for term in args:
            if term.is_commutative:
                if 0 in blades:
                    blades[0] += term
                else:
                    blades[0] = term
            else:
                c, nc = term.args_cnc()
                blade = nc[0]
                grade = self.blades_to_grades_dict[blade]
                if grade in blades:
                    blades[grade] += term
                else:
                    blades[grade] = term
        s = zero
        for grade in blades:
            if (grade * (grade - 1)) / 2 % 2 == 0:
                s += blades[grade]
            else:
                s -= blades[grade]
        return s

    def get_grade(self, A, r):  # Return grade r of A, <A>_{r}
        if r == 0:
            return self.scalar_part(A)
        coefs, bases = metric.linear_expand(A)
        s = zero
        for (coef, base) in zip(coefs, bases):
            if base != one and self.blades_to_grades_dict[base] == r:
                s += coef * base
        return s

    def even_odd(self, A, even=True):  # Return even or odd part of A
        A = expand(A)
        if A.is_commutative and even:
            return A
        if isinstance(A, Add):
            args = A.args
        else:
            args = [A]
        s = zero
        for term in args:
            if term.is_commutative:
                if even:
                    s += term
            else:
                c, nc = term.args_cnc(split_1=False)
                blade = nc[0]
                grade = self.blades_to_grades_dict[blade]
                if even and grade % 2 == 0:
                    s += Mul._from_args(c) * blade
                elif not even and grade % 2 == 1:
                    s += Mul._from_args(c) * blade
        return s

    ##################### Multivector derivatives ######################

    def build_reciprocal_basis(self,gsym):
        """
        Calculate reciprocal basis vectors e^{j} where
                e^{j}\cdot e_{k} = \delta_{k}^{j}
        and \delta_{k}^{j} is the kronecker delta.  We use the formula
        from Doran and Lasenby 4.94 -
            e^{j} = (-1)**{j-1}e_{1}^...e_{j-1}^e_{j+1}^...^e_{n}*E_{n}**{-1}
        where E_{n} = e_{1}^...^e_{n}.

        For non-orthogonal basis e^{j} is not normalized and must be
        divided by E_{n}**2 (self.inorm) in any relevant calculations.

        If gsym = True then (E_{n})**2 is not evaluated, but is represented
        as (E_{n})**2 = (-1)**(n*(n-1)/2)*det(g) where det(g) the determinant
        of the metric tensor is a general scalar function of the coordinates.
        """

        if self.debug:
            print 'Enter build_reciprocal_basis.\n'

        if self.is_ortho:
            self.r_basis = [self.basis[i] / self.g[i, i] for i in self.n_range]
        else:
            self.iobj = self.blade_to_base_rep(self.blades_lst[-1])
            if gsym:
                n = len(self.coords)
                self.inorm_sq = (-1) ** (n*(n - 1)/2) * Function('det(g)')(*self.coords)
            else:
                self.inorm = expand(self.basic_mul(self.iobj, self.iobj))
                self.inorm_sq = metric.Simp.apply(self.inorm * self.inorm)
            if self.debug:
                print 'I**2 =', self.inorm_sq
            duals = list(self.blades_lst[-(self.n + 1):-1])
            duals.reverse()

            sgn = 1
            self.r_basis = []
            for dual in duals:
                dual_base_rep = self.blade_to_base_rep(dual)
                tmp = collect(expand(self.base_to_blade_rep(self.mul(sgn * dual_base_rep, self.iobj))), self.blades_lst)
                self.r_basis.append(tmp)
                sgn = -sgn

            if self.debug:
                printer.oprint('I', self.iobj, 'Inorm', self.inorm, 'unnormalized reciprocal basis =\n', self.r_basis)
                self.dot_mode = '|'
                print 'reciprocal basis test ='
                for ei in self.basis:
                    for ej in self.r_basis:
                        ei_dot_ej = self.dot(ei, ej)
                        if ei_dot_ej == zero:
                            print 'e_{i}|e_{j} = ' + str(ei_dot_ej)
                        else:
                            print 'e_{i}|e_{j} = ' + str(expand(ei_dot_ej / self.inorm))

        self.iobj = self.blades_lst[-1]

        # Dictionary to represent reciprocal basis vectors as expansions
        # in terms of basis vectors.

        self.r_basis_dict = {}
        self.r_basis_mv = []
        for (r_symbol, r_base) in zip(self.r_symbols, self.r_basis):
            self.r_basis_dict[r_symbol] = r_base
            self.r_basis_mv.append(mv.Mv(r_base, ga=self))

        # Replace reciprocal basis vectors with expansion in terms of
        # basis vectors in derivatives of basis vectors

        if self.connect_flg:
            for x_i in self.n_range:
                for jb in self.n_range:
                    if not self.is_ortho:
                        self.de[x_i][jb] = metric.Simp.apply(self.de[x_i][jb].subs(self.r_basis_dict) / self.inorm)
                    else:
                        self.de[x_i][jb] = metric.Simp.apply(self.de[x_i][jb].subs(self.r_basis_dict))

        g_inv = eye(self.n)
        self.dot_mode = '|'

        # Calculate inverse of metric tensor, g^{ij}

        for i in self.n_range:
            rx_i = self.r_symbols[i]
            for j in self.n_range:
                rx_j = self.r_symbols[j]
                if j >= i:
                    g_inv[i, j] = self.dot(self.r_basis_dict[rx_i], self.r_basis_dict[rx_j])
                    if not self.is_ortho:
                        g_inv[i, j] /= self.inorm
                else:
                    g_inv[i, j] = g_inv[j, i]

        self.g_inv = g_inv

        if self.debug:
            print 'reciprocal basis dictionary =\n', self.r_basis_dict

        # True is for left derivative and False is for right derivative
        self.deriv = {('*', True): [], ('^', True): [], ('|', True): [],
                      ('<', True): [], ('>', True): [], ('*', False): [],
                      ('^', False): [], ('|', False): [], ('<', False): [],
                      ('>', False): []}
        return

    def er_blade(self, er, blade, mode='*', left=True):
        """
        Product (*,^,|,<,>) of reciprocal basis vector 'er' and basis
        blade 'blade' needed for application of derivatives to
        multivectors.  left is 'True' means 'er' is multiplying 'blade'
        on the left, 'False' is for 'er' multiplying 'blade' on the
        right.  Symbolically for left geometric product -
            e^{j}*(e_{i_{1}}^...^e_{i_{r}})
        """
        if mode == '*':
            base = self.blade_to_base_rep(blade)
            if left:
                return self.base_to_blade_rep(self.mul(er, base))
            else:
                return self.base_to_blade_rep(self.mul(base, er))
        elif mode == '^':
            if left:
                return self.wedge(er, blade)
            else:
                return self.wedge(blade, er)
        else:
            self.dot_mode = mode
            if left:
                return self.dot(er, blade)
            else:
                return self.dot(blade, er)

    def blade_derivation(self, blade, ib):
        """
        Calculate derivatives of basis blade 'blade' using derivative of
        basis vectors calculated by metric. 'ib' is the index of the
        coordinate the derivation is with respect to or the coordinate
        symbol.  These are requried for the calculation of the geometric
        derivatives in curvilinear coordinates or for more general
        manifolds.

        'blade_derivation' saves the results in a dictionary, 'self.dbases',
        so that the derivation for a given blade and coordinate is never
        calculated more that once.
        """

        if isinstance(ib, int):
            coord = self.coords[ib]
        else:
            coord = ib
            ib = self.coords.index(coord)

        key = (coord, blade)
        if key in self.dbases:
            return self.dbases[key]

        index = self.blades_to_indexes_dict[blade]
        grade = len(index)
        if grade == 1:
            db = self.de[ib][index[0]]
        elif grade == 2:
            db = self.wedge(self.de[ib][index[0]], self.basis[index[1]]) + \
                self.wedge(self.basis[index[0]], self.de[ib][index[1]])
        else:
            db = self.wedge(self.de[ib][index[0]], self.indexes_to_blades[index[1:]]) + \
                self.wedge(self.indexes_to_blades[index[:-1]], self.de[ib][index[-1]])
            for i in range(1, grade - 1):
                db += self.wedge(self.wedge(self.indexes_to_blades[index[:i]], self.de[ib][index[i]]),
                                 self.indexes_to_blades[index[i + 1:]])
        self.dbases[key] = db
        return db

    def pDiff(self, A, coord):
        """
        Compute partial derivative of multivector function 'A' with
        respect to coordinate 'coord'.
        """

        if isinstance(coord, list):
            # Perform multiple partial differentiation where coord is
            # a list of differentiation orders for each coordinate and
            # the coordinate is determinded by the list index.  If the
            # element in the list is zero no differentiation is to be
            # performed for that coordinate index.

            dA = copy.copy(A)  # Make copy of A

            for i in self.n_range:
                x = self.coords[i]
                xn = coord[i]
                if xn > 0:  # Differentiate with respect to coordinate x
                    for j in range(xn):  # xn > 1 multiple differentiation
                        dA = self.pDiff(dA, x)

            return dA

        # Simple partial differentiation, once with respect to a single
        # variable, but including case of non-constant basis vectors

        dA = self.mv(expand(diff(A.obj, coord)))

        if self.connect_flg and self.dslot == -1 and not A.is_scalar():  # Basis blades are function of coordinates
            B = self.remove_scalar_part(A)
            if B != zero:
                if isinstance(B, Add):
                    args = B.args
                else:
                    args = [B]
                for term in args:
                    if not term.is_commutative:
                        c, nc = term.args_cnc(split_1=False)
                        x = self.blade_derivation(nc[0], coord)
                        if x != zero:
                            if len(c) == 1:
                                dA += c[0] * x
                            elif len(c) == 0:
                                dA += x
                            else:
                                dA += reduce(operator.mul, c, one) * x

        return dA

    def grad_sqr(self, A, grad_sqr_mode, mode, left):
        """
        Caclulate '(grad *_{1} grad) *_{2} A' or 'A *_{2} (grad *_{1} grad)'
        where grad_sqr_mode = *_{1} = *, ^, or | and
        mode = *_{2} = *, ^, or |.
        """
        (Sop, Bop) = Ga.DopFop[(grad_sqr_mode, mode)]
        print '(Sop, Bop) =', Sop, Bop

        print 'grad_sqr:A =', A

        self.dot_mode == '|'
        s = zero

        if Sop is False and Bop is False:
            return s

        dA_i = []
        for coord_i in self.coords:
            dA_i.append(self.pDiff(A, coord_i))

        print 'dA_i =', dA_i

        if Sop:
            for i in self.n_range:
                coord_i = self.coords[i]
                if self.connect_flg:
                    s += self.grad_sq_scl_connect[coord_i] * dA_i[i]

                for j in self.n_range:
                    d2A_j = self.pDiff(dA_i[i], self.coords[j])
                    s += self.g_inv[i, j] * d2A_j

        if Bop and self.connect_flg:
            for i in self.n_range:
                coord_i = self.coords[i]
                print 'mode =', mode
                print 'dA_i[i] =', dA_i[i]
                if left:
                    if mode == '|':
                        s += self.dot(self.grad_sq_mv_connect[coord_i], dA_i[i])
                    if mode == '^':
                        s += self.wedge(self.grad_sq_mv_connect[coord_i], dA_i[i])
                    if mode == '*':
                        s += self.mul(self.grad_sq_mv_connect[coord_i], dA_i[i])
                else:
                    if mode == '|':
                        s += self.dot(dA_i[i], self.grad_sq_mv_connect[coord_i])
                    if mode == '^':
                        s += self.wedge(dA_i[i], self.grad_sq_mv_connect[coord_i])
                    if mode == '*':
                        s += self.mul(dA_i[i], self.grad_sq_mv_connect[coord_i])
        return s

    def connection(self, rbase, key_base, mode, left):
        """
        Compute required multivector connections of the form
        (Einstein summation convention) e^{j}*(D_{j}e_{i_{1}...i_{r}})
        and (D_{j}e_{i_{1}...i_{r}})*e^{j} where * could be *, ^, |,
        <, or > depending upon the mode and e^{j} are reciprocal
        basis vectors
        """
        mode_key = (mode, left)
        keys = [i for i, j in self.connect[mode_key]]
        if left:
            key = rbase * key_base
        else:
            key = key_base * rbase
        if key not in keys:
            keys.append(key)
            C = zero
            for ib in self.n_range:
                x = self.blade_derivation(key_base, ib)
                if self.norm:
                    x /= self.e_norm[ib]
                C += self.er_blade(self.r_basis[ib], x, mode, left)
            # Update connection dictionaries
            self.connect[mode_key].append((key, C))
        return C

    def ReciprocalFrame(self, basis, mode='norm'):
        dim = len(basis)

        indexes = tuple(range(dim))
        index = [()]

        for i in indexes[-2:]:
            index.append(tuple(combinations(indexes, i + 1)))

        MFbasis = []

        for igrade in index[-2:]:
            grade = []
            for iblade in igrade:
                blade = self.mv('1', 'scalar')
                for ibasis in iblade:
                    blade ^= basis[ibasis]
                blade = blade.trigsimp()
                grade.append(blade)
            MFbasis.append(grade)
        E = MFbasis[-1][0]
        E_sq = trigsimp((E * E).scalar())

        duals = copy.copy(MFbasis[-2])

        duals.reverse()
        sgn = 1
        rbasis = []
        for dual in duals:
            recpv = (sgn * dual * E).trigsimp()
            rbasis.append(recpv)
            sgn = -sgn

        if mode != 'norm':
            rbasis.append(E_sq)
        else:
            for i in range(dim):
                rbasis[i] = rbasis[i] / E_sq

        return tuple(rbasis)

    def Mlt(self,*kargs,**kwargs):
        return lt.Mlt(kargs[0], self, *kargs[1:], **kwargs)

class Sm(Ga):
    """
    Submanifold is a geometric algebra defined on a submanifold of a
    base geometric algebra defined on a manifold.  The submanifold is
    defined by a mapping from the coordinates of the base manifold to
    the coordinates of the submanifold. The inputs required to define
    the submanifold are:

        u      {kargs[0]} The coordinate map defining the submanifold
               which is a list of functions of coordinates of the base
               manifold in terms of the coordinates of the submanifold.
               for example if the manifold is a unit sphere then -
               'u = [sin(u)*cos(v),sin(u)*sin(v),cos(u)]'.

        coords {kargs[1]} The coordinate list for the submanifold, for
               example '[u,v]'.

    See 'init_slots' for possible other inputs.  The 'Ga' member function
    'sm' can be used to instanciate the submanifold via (o3d is the base
    manifold)

        coords = (u,v) = symbols(',v',real=True)
        sm_example = o3d.sm([sin(u)*cos(v),sin(u)*sin(v),cos(u)],coords)

        (eu,ev) = sm_example.mv()
        sm_grad = sm_example.grad
    """
    init_slots = {'debug': (False, 'True for debug output'),
                  'root': ('e', 'Root symbol for basis vectors'),
                  'name': (None, 'Name of submanifold'),
                  'norm': (False, 'Normalize basis if True'),
                  'ga': (None, 'Base Geometric Algebra')}

    def __init__(self, *kargs, **kwargs):

        if printer.GaLatexPrinter.latex_flg:
            printer.GaLatexPrinter.restore()
            Ga.restore = True

        kwargs = metric.test_init_slots(Sm.init_slots, **kwargs)
        u = kargs[0]  # Coordinate map to define submanifold
        coords = kargs[1]  # List of cordinates
        ga = kwargs['ga']  # base geometric algebra
        if ga is None:
            raise ValueError('Base geometric algebra must be specified for submanifold.')

        g_base = ga.g_raw
        n_base = ga.n
        n_sub = len(coords)

        if len(u) != n_base:
            raise ValueError('In submanifold dimension of base manifold' +
                              ' not equal to dimension of mapping.')

        dxdu = []
        for x_i in u:
            tmp = []
            for u_j in coords:
                tmp.append(diff(x_i, u_j))
            dxdu.append(tmp)

        sub_pairs = zip(ga.coords, u)

        g = eye(n_sub)
        n_range = range(n_sub)
        for i in n_range:
            for j in n_range:
                s = zero
                for k in ga.n_range:
                    for l in ga.n_range:
                        s += dxdu[k][i] * dxdu[l][j] * g_base[k, l].subs(sub_pairs)
                #g[i, j] = trigsimp(s, method='combined')
                g[i, j] = trigsimp(s)

        root = kwargs['root']

        basis_str = ''
        for x in coords:
            basis_str += root + '_' + str(x) + ' '
        basis_str = basis_str[:-1]

        norm = kwargs['norm']
        debug = kwargs['debug']

        if Ga.restore:  # restore printer to appropriate enhanced mode after sm is instantiated
            printer.GaLatexPrinter.redirect()

        Ga.__init__(self, root, g=g, coords=coords, norm=norm, debug=debug)

        self.ga = ga
        self.u = u

        if debug:
            print 'Exit Sm.__init__()'

if __name__ == "__main__":
    pass
