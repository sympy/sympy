#!/usr/bin/python

"""
The module symbolicGA implements symbolic Geometric Algebra in python.
The relevant references for this module are:

    1. "Geometric Algebra for Physicists" by C. Doran and A. Lazenby,
       Cambridge University Press, 2003.

    2. "Geometric Algebra for Computer Science" by Leo Dorst, Daniel Fontijne,
       and Stephen Mann, Morgan Kaufmann Publishers, 2007

    3. Sympy Tutorial, http://docs.sympy.org/
"""
import sys
import os, string, types, copy
import numpy, sympy
import re as regrep
import sympy.galgebra.latex_ex

NUMPAT = regrep.compile( '([\-0-9])|([\-0-9]/[0-9])')
"""Re pattern for rational number"""

ZERO = sympy.Rational(0)
ONE  = sympy.Rational(1)
TWO  = sympy.Rational(2)
HALF = sympy.Rational(1,2)

sym_type = sympy.core.symbol.Symbol
pow_type = sympy.core.power.Pow
abs_type = sympy.abs
mul_type = sympy.core.mul.Mul
add_type = sympy.core.add.Add

global MAIN_PROGRAM

MAIN_PROGRAM = ''

@sympy.vectorize(0)
def substitute_array(array,*args):
    return(array.subs(*args))

def is_quasi_unit_numpy_array(array):
    """
    Determine if a array is square and diagonal with
    entries of +1 or -1.
    """
    shape = numpy.shape(array)
    if len(shape) == 2 and (shape[0] == shape[1]):
        n = shape[0]
        ix = 0
        while ix < n:
            iy = 0
            while iy < ix:
                if array[ix][iy] != ZERO:
                    return(False)
                iy += 1
            if sympy.abs(array[ix][ix]) != ONE:
                return(False)
            ix += 1
        return(True)
    else:
        return(False)

def set_main(main_program):
    global MAIN_PROGRAM
    MAIN_PROGRAM = main_program
    return

def plist(lst):
    if type(lst) == types.ListType:
        for x in lst:
            plist(x)
    else:
        sys.stderr.write(lst+'\n')
    return

def numeric(num_str):
    """
    Returns rational numbers compatible with symbols.
    Input is a string representing a fraction or integer
    or a simple integer.
    """
    if type(num_str) == types.IntType:
        a = num_str
        b = 1
    else:
        tmp = num_str.split('/')
        if len(tmp) == 1:
            a = int(tmp[0])
            b = 1
        else:
            a = int(tmp[0])
            b = int(tmp[1])
    return(sympy.Rational(a,b))

def collect(expr,lst):
    """
    Wrapper for sympy.collect.
    """
    lst = MV.scalar_to_symbol(lst)
    return(sympy.collect(expr,lst))

def sqrt(expr):
    return(sympy.sqrt(expr))

def isint(a):
    """
    Test for integer.
    """
    return(type(a) == types.IntType)

def make_null_array(n):
    """
    Return list of n empty lists.
    """
    a = []
    for i in range(n):
        a.append([])
    return(a)

def test_int_flgs(lst):
    """
    Test if all elements in list are 0.
    """
    for i in lst:
        if i:
            return(1)
    return(0)

def comb(N,P):
    """
    Calculates the combinations of the integers [0,N-1] taken P at a time.
    The output is a list of lists of integers where the inner lists are
    the different combinations. Each combination is sorted in ascending
    order.
    """
    def rloop(n,p,combs,comb):
        if p:
            for i in range(n):
                newcomb = comb+[i]
                np = p-1
                rloop(i,np,combs,newcomb)
        else:
            combs.append(comb)
    combs = []
    rloop(N,P,combs,[])
    for comb in combs:
        comb.sort()
    return(combs)

def diagpq(p,q=0):
    """
    Return string equivalent metric tensor for signature (p,q).
    """
    n = p+q
    D = ''
    rn = range(n)
    for i in rn:
        for j in rn:
            if i ==j:
                if i < p:
                    D += '1 '
                else:
                    D += '-1 '
            else:
                D += '0 '
        D = D[:-1]+','
    return(D)

def make_scalars(symnamelst):
    """
    make_symbols takes a string of symbol names separated by
    blanks and converts them to MV scalars separately
    accessible by the main program and addition returns a list
    of the symbols.
    """
    global MAIN_PROGRAM
    if type(symnamelst) == types.StringType:
        symnamelst = symnamelst.split()
    scalar_lst = []
    isym = 0
    for name in symnamelst:
        tmp = sympy.Symbol(name)
        tmp = MV(tmp,'scalar')
        scalar_lst.append(tmp)
        setattr(MAIN_PROGRAM,name,tmp)
        isym += 1
    return(scalar_lst)

def make_symbols(symnamelst):
    """
    make_symbols takes a string of symbol names separated by
    blanks and converts them to MV scalars separately
    accessible by the main program and addition returns a list
    of the symbols.
    """
    global MAIN_PROGRAM
    if type(symnamelst) == types.StringType:
        symnamelst = symnamelst.split()
    sym_lst = []
    isym = 0
    for name in symnamelst:
        tmp = sympy.Symbol(name)
        sym_lst.append(tmp)
        setattr(MAIN_PROGRAM,name,tmp)
        isym += 1
    return(sym_lst)

def israt(numstr):
    """
    Test if string represents a rational number.
    """
    global NUMPAT
    if NUMPAT.match(numstr):
        return(1)
    return(0)

def dualsort(lst1, lst2):
    """
    Inplace dual sort of lst1 and lst2 keyed on sorted lst1.
    """
    _indices = range(len(lst1))
    _indices.sort(key=lst2.__getitem__)
    lst1[:] = map(lst1.__getitem__, _indices)
    lst2[:] = map(lst2.__getitem__, _indices)
    return

def cp(A,B):
    """
    Calculates the comutator product (A*B-B*A)/2 for
    the objects A and B.
    """
    return(HALF*(A*B-B*A))

def reduce_base(k,base):
    """
    If base is a list of sorted integers [i_1,...,i_R] then reduce_base
    sorts the list [k,i_1,...,i_R] and calculates whether an odd or even
    number of permutations is required to sort the list. The sorted list
    is returned and +1 for even permutations or -1 for odd permutations.
    """
    if k in base:
        return(0,base)
    grade = len(base)
    if grade == 1:
        if k < base[0]:
            return(1,[k,base[0]])
        else:
            return(-1,[base[0],k])
    ilo = 0
    ihi = grade-1
    if k < base[0]:
        return(1,[k]+base)
    if k > base[ihi]:
        if grade%2 == 0:
            return(1,base+[k])
        else:
            return(-1,base+[k])
    imid = ihi+ilo
    if grade == 2:
        return(-1,[base[0],k,base[1]])
    while True:
        if ihi-ilo == 1:
            break
        if base[imid] > k:
            ihi = imid
        else:
            ilo = imid
        imid = (ilo+ihi)/2
    if ilo%2 == 1:
        return(1,base[:ihi]+[k]+base[ihi:])
    else:
        return(-1,base[:ihi]+[k]+base[ihi:])

def sub_base(k,base):
    """
    If base is a list of sorted integers [i_1,...,i_R] then sub_base returns
    a list with the k^th element removed. Note that k=0 removes the first
    element.  The is no test to see if k is in the range of the list.
    """
    n = len(base)
    if n == 1:
        return([])
    if n == 2:
        if k == base[0]:
            return([base[1]])
        else:
            return([base[0]])
    return(base[:k]+base[k+1:])

def magnitude(vector):
    """
    Calculate magnitude of vector containing trig expressions
    and simplify.  This is a hack because of way sign of
    magsq is determined and because of the way absoluted
    values are removed.
    """
    magsq = sympy.expand((vector|vector)())
    magsq = sympy.trigsimp(magsq,deep=True,recursive=True)
    #print magsq
    magsq_str = sympy.galgebra.latex_ex.LatexPrinter()._print(magsq)
    if magsq_str[0] == '-':
        magsq = -magsq
    mag = unabs(sqrt(magsq))
    #print mag
    return(mag)

def LaTeX_lst(lst,title=''):
    """
    Ouput a list in LaTeX format.
    """
    if title != '':
        LaTeX(title)
    for x in lst:
        LaTeX(x)
    return

def unabs(x):
    """
    Remove absolute values from expressions so a = sqrt(a**2).
    This is a hack.
    """
    if type(x) == mul_type:
        y = unabs(x.args[0])
        for yi in x.args[1:]:
            y *= unabs(yi)
        return(y)
    if type(x) == pow_type:
        if x.args[1] == HALF and type(x.args[0]) == add_type:
            return(x)
        y = 1/unabs(x.args[0])
        return(y)
    if len(x.args) == 0:
        return(x)
    if type(x) == abs_type:
        return(x.args[0])
    return(x)

def function_lst(fstr,xtuple):
    sys.stderr.write(fstr+'\n')
    fct_lst = []
    for xstr in fstr.split():
        f = sympy.Function(xstr)(*xtuple)
        fct_lst.append(f)
    return(fct_lst)

def vector_fct(Fstr,x):
    """
    Create a list of functions of arguments x.  One function is
    created for each variable in x.  Fstr is a string that is
    the base name of each function while each fuction in the
    list is given the name Fstr+'__'+str(x[ix]) so that if
    Fstr = 'f' and str(x[1]) = 'theta' then the LaTeX output
    of the second element in the output list would be 'f^{\theta}'.
    """
    nx = len(x)
    Fvec = []
    for ix in range(nx):
        ftmp = sympy.Function(Fstr+'__'+sympy.galgebra.latex_ex.LatexPrinter.str_basic(x[ix]))(*tuple(x))
        Fvec.append(ftmp)
    return(Fvec)

def print_lst(lst):
    for x in lst:
        print x
    return

def normalize(elst,nname_lst):
    """
    Normalize a list of vectors and rename the normalized vectors. 'elist' is the list
    (or array) of vectors to be normalized and nname_lst is a list of the names for the
    normalized vectors.  The function returns the numpy arrays enlst and mags containing
    the normalized vectors (enlst) and the magnitudes of the original vectors (mags).
    """
    i = 0
    mags = numpy.array(MV.n*[ZERO],dtype=numpy.object)
    enlst= numpy.array(MV.n*[ZERO],dtype=numpy.object)
    for (e,nname) in zip(elst,nname_lst):
        emag = magnitude(e)
        emaginv = 1/emag
        mags[i] = emag
        enorm = emaginv*e
        enorm.name = nname
        enlst[i] = enorm
        i += 1
    return(enlst,mags)

def build_base(base_index,base_vectors,reverse=False):
    base = base_vectors[base_index[0]]
    if len(base_index) > 1:
        for i in base_index[1:]:
            base = base^base_vectors[i]
    if reverse:
        base = base.rev()
    return(base)

class MV(object):

    is_setup = False
    basislabel_lst = 0
    curvilinear_flg = False
    coords = None

    @staticmethod
    def pad_zeros(value,n):
        """
        Pad list with zeros to lenght n. If length is > n
        truncate list.  Return padded list.
        """
        nvalue = len(value)
        if nvalue < n:
            value = value+(n-nvalue)*[ZERO]
        if nvalue > n:
            value = value[:n]
        return(value)

    @staticmethod
    def define_basis(basis):
        """
        Calculates all the MV static variables needed for
        basis operations.  See reference 5 section 2.
        """
        MV.vbasis     = basis
        MV.vsyms      = make_symbols(MV.vbasis)
        MV.n          = len(MV.vbasis)
        MV.nrg        = range(MV.n)
        MV.n1         = MV.n+1
        MV.n1rg       = range(MV.n1)
        MV.npow       = 2**MV.n
        MV.index      = range(MV.n)
        MV.gabasis    = [[]]
        MV.basis      = (MV.n+1)*[0]
        MV.basislabel = (MV.n+1)*[0]
        MV.basis[0]   = []
        MV.basislabel[0] = '1'
        MV.basislabel_lst = [['1']]
        MV.nbasis = numpy.array((MV.n+1)*[1],dtype=numpy.object)
        for igrade in range(1,MV.n+1):
            tmp = comb(MV.n,igrade)
            MV.gabasis += [tmp]
            ntmp = len(tmp)
            MV.nbasis[igrade] = ntmp
            MV.basis[igrade] = tmp
            gradelabels = []
            gradelabel_lst = []
            for i in range(ntmp):
                tmp_lst = []
                bstr = ''
                for j in tmp[i]:
                    bstr += MV.vbasis[j]
                    tmp_lst.append(MV.vbasis[j])
                gradelabel_lst.append(tmp_lst)
                gradelabels.append(bstr)
            MV.basislabel_lst.append(gradelabel_lst)
            MV.basislabel[igrade] = gradelabels
        MV.basis_map = [{'':0}]
        igrade = 1
        while igrade <= MV.n:
            tmpdict = {}
            bases = MV.gabasis[igrade]
            nbases = len(bases)
            ibases = 0
            for base in bases:
                tmpdict[str(base)] = ibases
                ibases += 1
            MV.basis_map.append(tmpdict)
            igrade += 1

        if MV.debug:
            print 'basis strings =',MV.vbasis
            print 'basis symbols =',MV.vsyms
            print 'basis labels  =',MV.basislabel
            print 'basis         =',MV.basis
            print 'grades        =',MV.nbasis
            print 'index         =',MV.index
        return

    @staticmethod
    def define_metric(metric):
        """
        Calculates all the MV static variables needed for
        metric operations.  See reference 5 section 2.
        """
        if MV.metric_str:
            name_flg = False
            MV.g = []
            MV.metric = numpy.array(MV.n*[MV.n*[ZERO]],dtype=numpy.object)
            if metric == '':
                metric = numpy.array(MV.n*[MV.n*['#']],dtype=numpy.object)
            for i in MV.index:
                for j in MV.index:
                    gij = metric[i][j]
                    if israt(gij):
                        MV.metric[i][j] = numeric(gij)
                    else:
                        if gij == '#':
                            name_flg = True
                            if i == j:
                                gij = '('+MV.vbasis[j]+'**2)'
                                name = MV.vbasis[j]+'sq'
                            else:
                                gij = '('+MV.vbasis[min(i,j)]+'.'+MV.vbasis[max(i,j)]+')'
                                name = MV.vbasis[min(i,j)]+'dot'+MV.vbasis[max(i,j)]
                        tmp = sympy.Symbol(gij)
                        MV.metric[i][j] = tmp
                        if i <= j:
                            MV.g.append(tmp)
                            if name_flg:
                                setattr(MAIN_PROGRAM,name,tmp)
                                name_flg = False
        else:
            MV.metric = metric
            MV.g = []
            for row in metric:
                g_row = []
                for col in metric:
                    g_row.append(col)
                    MV.g.append(g_row)
        if MV.debug:
            print 'metric =',MV.metric
        return

    @staticmethod
    def define_reciprocal_frame():
        """
        Calculates unscaled reciprocal vectors (MV.brecp) and scale
        factor (MV.Esq). The ith scaled reciprocal vector is
        (1/MV.Esq)*MV.brecp[i].  The pseudoscalar for the set of
        basis vectors is MV.E.
        """
        if MV.tables_flg:
            MV.E = MV.bvec[0]
            MV.brecp = []
            for i in range(1,MV.n):
                MV.E = MV.E^MV.bvec[i]
            for i in range(MV.n):
                tmp = ONE
                if i%2 != 0:
                    tmp = -ONE
                for j in range(MV.n):
                    if i != j:
                        tmp = tmp^MV.bvec[j]
                tmp = tmp*MV.E
                MV.brecp.append(tmp)
            MV.Esq = (MV.E*MV.E)()
            MV.Esq_inv = ONE/MV.Esq
            for i in range(MV.n):
                MV.brecp[i] = MV.brecp[i]*MV.Esq_inv
        return

    @staticmethod
    def reduce_basis_loop(blst):
        """
        Makes one pass through basis product representation for
        reduction of representation to normal form.
        See reference 5 section 3.
        """
        nblst = len(blst)
        if nblst <= 1:
            return(1)
        jstep = 1
        while jstep <nblst:
            istep = jstep-1
            if blst[istep] == blst[jstep]:
                i = blst[istep]
                if len(blst) >2:
                    blst = blst[:istep]+blst[jstep+1:]
                else:
                    blst = []
                if len(blst) <= 1 or jstep == nblst-1:
                    blst_flg = 0
                else:
                    blst_flg = 1
                return(MV.metric[i][i],blst,blst_flg)
            if blst[istep] > blst[jstep]:
                blst1 = blst[:istep]+blst[jstep+1:]
                a1 = TWO*MV.metric[blst[jstep]][blst[istep]]
                blst = blst[:istep]+[blst[jstep]]+[blst[istep]]+blst[jstep+1:]
                if len(blst1) <= 1:
                    blst1_flg = 0
                else:
                    blst1_flg = 1
                return(a1,blst1,blst1_flg,blst)
            jstep +=1
        return(1)

    @staticmethod
    def reduce_basis(blst):
        """
        Repetively applies reduce_basis_loop to basis
        product representation untill normal form is
        realized.  See reference 5 section 3.
        """
        if blst == []:
            blst_coef   = []
            blst_expand = []
            for i in MV.n1rg:
                blst_coef.append([])
                blst_expand.append([])
            blst_expand[0].append([])
            blst_coef[0].append(ONE)
            return(blst_coef,blst_expand)
        blst_expand = [blst]
        blst_coef      = [ONE]
        blst_flg         = [1]
        while test_int_flgs(blst_flg):
            for i in range(len(blst_flg)):
                if blst_flg[i]:
                    tmp = MV.reduce_basis_loop(blst_expand[i])
                    if tmp ==1:
                        blst_flg[i] = 0
                    else:
                        if len(tmp) == 3:
                            blst_coef[i] = tmp[0]*blst_coef[i]
                            blst_expand[i] = tmp[1]
                            blst_flg[i] = tmp[2]
                        else:
                            blst_coef[i]      = -blst_coef[i]
                            blst_flg[i]        = 1
                            blst_expand[i] = tmp[3]
                            blst_coef.append(-blst_coef[i]*tmp[0])
                            blst_expand.append(tmp[1])
                            blst_flg.append(tmp[2])
        (blst_coef,blst_expand) = MV.combine_common_factors(blst_coef,blst_expand)
        return(blst_coef,blst_expand)

    @staticmethod
    def combine_common_factors(blst_coef,blst_expand):
        new_blst_coef = []
        new_blst_expand = []
        for i in range(MV.n1):
            new_blst_coef.append([])
            new_blst_expand.append([])
        nfac = len(blst_coef)
        for ifac in range(nfac):
            blen = len(blst_expand[ifac])
            new_blst_coef[blen].append(blst_coef[ifac])
            new_blst_expand[blen].append(blst_expand[ifac])
        for i in range(MV.n1):
            if len(new_blst_coef[i]) > 1:
                MV.contract(new_blst_coef[i],new_blst_expand[i])
        return(new_blst_coef,new_blst_expand)

    @staticmethod
    def contract(coefs,bases):
        dualsort(coefs,bases)
        n = len(bases)-1
        i = 0
        while i < n:
            j = i+1
            if bases[i] == bases[j]:
                coefs[i] += coefs[j]
                bases.pop(j)
                coefs.pop(j)
                n -= 1
            else:
                i += 1
        n = len(coefs)
        i = 0
        while i < n:
            if coefs[i] == ZERO:
                coefs.pop(i)
                bases.pop(i)
                n -= 1
            else:
                i +=1
        return

    @staticmethod
    def convert(coefs,bases):
        mv = MV()
        mv.bladeflg = 0
        for igrade in MV.n1rg:
            coef = coefs[igrade]
            base = bases[igrade]
            if len(coef) > 0:
                nbases = MV.nbasis[igrade]
                mv.mv[igrade] = numpy.array(nbases*[ZERO],dtype=numpy.object)
                nbaserg = range(len(base))
                for ibase in nbaserg:
                    if igrade > 0:
                        k = MV.basis[igrade].index(base[ibase])
                        mv.mv[igrade][k] = coef[ibase]
                    else:
                        mv.mv[igrade] = numpy.array([coef[0]],dtype=numpy.object)
        return(mv)

    @staticmethod
    def set_str_format(str_mode=0):
        MV.str_mode = str_mode
        return

    @staticmethod
    def str_rep(mv):
        """
         Converts internal representation of a multivector to a string
         for outputing.  If lst_mode = 1, str_rep outputs a list of
         strings where each string contains one multivector coefficient
         concatenated with the corresponding base or blade symbol.

            MV.str_mode     Effect
                0           Print entire multivector on one line (default)
                1           Print each grade on a single line
                2           Print each base on a single line
         """

        if MV.bladeprint:
            mv.convert_to_blades()
            labels = MV.bladelabel
        else:
            if not mv.bladeflg:
                labels = MV.basislabel
            else:
                labels = MV.bladelabel
        mv.compact()
        if isinstance(mv.mv[0],types.IntType):
            value = ''
        else:
            value = (mv.mv[0][0]).__str__()
            value = value.replace(' ','')
        dummy = sympy.Symbol('dummy')
        for igrade in MV.n1rg[1:]:
            if isinstance(mv.mv[igrade],numpy.ndarray):
                j = 0
                for x in mv.mv[igrade]:
                    if x != ZERO:
                        xstr = (x*dummy).__str__()
                        xstr = xstr.replace(' ','')
                        if xstr[0] != '-' and len(value) > 0:
                            xstr = '+'+xstr
                        if xstr.find('dummy') < 2 and xstr[-5:] != 'dummy':
                            xstr = xstr.replace('dummy*','')+'*'+labels[igrade][j]
                        else:
                            xstr = xstr.replace('dummy',labels[igrade][j])
                        if MV.str_mode == 2:
                            xstr += '\n'
                        value += xstr
                    j += 1
                if MV.str_mode == 1:
                    value += '\n'
        if value == '':
            value = '0'
        #value = value.replace(' ','')
        value = value.replace('dummy','1')
        return(value)

    @staticmethod
    def xstr_rep(mv,lst_mode=0):
        """
         Converts internal representation of a multivector to a string
         for outputing.  If lst_mode = 1, str_rep outputs a list of
         strings where each string contains one multivector coefficient
         concatenated with the corressponding base or blade symbol.
         """
        if lst_mode:
            outlst = []
        if MV.bladeprint:
            mv.convert_to_blades()
            labels = MV.bladelabel
        else:
            if not mv.bladeflg:
                labels = MV.basislabel
            else:
                labels = MV.bladelabel
        value = ''
        for igrade in MV.n1rg:
            tmp = []
            if isinstance(mv.mv[igrade],numpy.ndarray):
                j = 0
                xsum = 0
                for x in mv.mv[igrade]:
                    if x != ZERO:
                        xstr = x.__str__()
                        if xstr == '+1' or xstr == '1' or xstr == '-1':
                           if xstr == '+1' or xstr == '1':
                               xstr = '+'
                           else:
                               xstr = '-'
                        else:
                           if xstr[0] != '+':
                               xstr = '+('+xstr+')'
                           else:
                               xstr = '+('+xstr[1:]+')'
                        value += xstr+labels[igrade][j]
                        if MV.str_mode and not lst_mode:
                            value += value+'\n'
                        if lst_mode:
                            tmp.append(value)
                    j += 1
            if lst_mode:
                if len(tmp) > 0:
                    outlst.append(tmp)
                value = ''
        if not lst_mode:
            if len(value) > 1 and value[0] == '+':
                value = value[1:]
            if len(value) == 0:
                value = '0'
        else:
            value = outlst
        return(value)

    @staticmethod
    def setup(basis,metric='',rframe=False,coords=None,debug=False,offset=0):
        """
        MV.setup initializes the MV class by calculating the static
        multivector tables required for geometric algebra operations
        on multivectors.  See reference 5 section 2 for details on
        basis and metric arguments.
        """
        global MAIN_PROGRAM
        MV.is_setup = True
        MV.metric_str = False
        MV.debug = debug
        MV.bladeprint = 0
        MV.tables_flg = 0
        MV.str_mode  = 0
        MV.basisroot = ''
        MV.index_offset = offset
        if coords == None:
            MV.coords = None
        else:
            MV.coords = tuple(coords)
            rframe= True
        if type(basis) == types.StringType:
            basislst = basis.split()
            if len(basislst) == 1:
                MV.basisroot = basislst[0]
                basislst = []
                for coord in coords:
                    basislst.append(MV.basisroot+'_'+str(coord))
            MV.define_basis(basislst)
        if type(metric) == types.StringType:
            MV.metric_str = True
            if len(metric) > 0:
                if metric[0] == '[' and metric[-1] == ']':
                    tmps = metric[1:-1].split(',')
                    N = len(tmps)
                    metric = []
                    itmp = 0
                    for tmp in tmps:
                        xtmp = N*['0']
                        xtmp[itmp] = tmp
                        itmp += 1
                        metric.append(xtmp)
                else:
                    tmps = metric.split(',')
                    metric = []
                    for tmp in tmps:
                        xlst = tmp.split()
                        xtmp = []
                        for x in xlst:
                            xtmp.append(x)
                        metric.append(xtmp)

        MV.define_metric(metric)
        MV.multiplication_table()
        MV.blade_table()
        MV.inverse_blade_table()
        MV.tables_flg = 1
        isym = 0
        MV.bvec = []
        for name in MV.vbasis:
            bvar = MV(value=isym,mvtype='basisvector',mvname=name)
            bvar.bladeflg = 1
            MV.bvec.append(bvar)
            setattr(MAIN_PROGRAM,name,bvar)
            isym += 1
        if rframe:
            MV.define_reciprocal_frame()
        MV.I = MV(ONE,'pseudo','I')
        MV.ZERO = MV()
        Isq = (MV.I*MV.I)()
        MV.Iinv = (1/Isq)*MV.I
        return('Setup of '+basis+' complete!')

    @staticmethod
    def set_coords(coords):
        MV.coords = coords
        return

    @staticmethod
    def scalar_fct(fct_name):
        """
        Create multivector scalar function with name fct_name (string) and
        independent varibles coords (list of variable).  Default variables are
        those associated with each dimension of vector space.
        """
        phi = sympy.Function(fct_name)(*MV.coords)
        Phi = MV(phi,'scalar')
        Phi.name = fct_name
        return(Phi)

    @staticmethod
    def vector_fct(fct_name,vars=''):
        """
        Create multivector vector function with name fct_name (string) and
        independent varibles coords (list of variable).  Default variables are
        those associated with each dimension of vector space.
        """
        if isinstance(vars,types.StringType):
            Acoefs = vector_fct(fct_name,MV.coords)
        else:
            Acoefs =numpy.array(MV.n*[ZERO],dtype=numpy.object)
            x = MV.coords
            if isinstance(vars,sympy.core.symbol.Symbol):
                for icoef in MV.nrg:
                    Acoefs[icoef] = sympy.Function(fct_name+'__'+\
                                   sympy.galgebra.latex_ex.LatexPrinter.str_basic(x[icoef]))(vars)
            else:
                for icoef in MV.nrg:
                    Acoefs[icoef] = sympy.Function(fct_name+'__'+\
                                   sympy.galgebra.latex_ex.LatexPrinter.str_basic(x[icoef]))(*tuple(vars))
        A = MV(Acoefs,'vector',fct_name)
        return(A)

    @staticmethod
    def rebase(x,coords,base_name='',debug=False,debug_level=0):
        """
        Define curvilinear coordinates for previously defined vector (multivector) space (MV.setup has been run)
        with position vector, x, that is a vector function of the independendent coordinates, coords (list of
        sympy variables equal in length to dimension of vector space), and calculate:

            1. Frame (basis) vectors
            2. Normalized frame (basis) vectors.
            3. Metric tensor
            4. Reciprocal frame vectors
            5. Reciprocal metric tensor
            6. Connection multivectors

        The basis vectors are named with the base_name (string) and a subscript derived from the name of each
        coordinate.  So that if the base name is 'e' and the coordinated are [r,theta,z] the variable names
        of the frame vectors would be e_r, e_theta, and e_z.  For LaTeX output the names of the frame vectors
        would be e_{r}, e_{\theta}, and e_{z}.  Everthing needed to compute the geometric, outer, and inner
        derivatives of multivector functions in curvilinear coordinates is calculated.

        If debug is True all the quatities in the above list are output in LaTeX format.

        Currently rebase works with cylindrical and spherical coordinates in any dimension.  The limitation is the
        ability to automatically simplify complex sympy expressions generated while calculating the quantities in
        the above list.  This is why the debug option is included.  The debug_level can equal 0,1,2, or 3 and
        determines how far in the list to calculate (imput 0 to do the entire list) while debugging.
        """
        global MAIN_PROGRAM

        #Form root names for basis, reciprocal basis, normalized basis, and normalized reciprocal basis

        if base_name == '':
            base_name = MV.basisroot+'prm'

        LaTeX_base = sympy.galgebra.latex_ex.LatexPrinter.extended_symbol(base_name)
        bm = '\\bm{'+LaTeX_base+'}'
        bmhat = '\\hat{'+bm+'}'
        bstr = bmhat+'_{[i_{1},\dots, i_{R}]}'
        base_name += 'bm'
        base_name_hat = base_name+'hat'

        base_name_lst   = []
        nbase_name_lst  = []
        rbase_name_lst  = []
        rnbase_name_lst = []
        coords_lst = []

        for coord in coords:
            coord_str = sympy.galgebra.latex_ex.LatexPrinter.str_basic(coord)
            coords_lst.append(coord_str)
            base_name_lst.append(base_name+'_'+coord_str)
            rbase_name_lst.append(base_name+'__'+coord_str)
            nbase_name_lst.append(base_name_hat+'_'+coord_str)
            rnbase_name_lst.append(base_name_hat+'__'+coord_str)

        if not (MV.n == len(coords) == len(base_name_lst)):
            print 'rebaseMV inputs not congruent:'
            print 'MV.n =',MV.n
            print 'coords =',coords
            print 'bases =',base_name
            sys.exit(1)

        if isinstance(x,MV):

            #Calculate basis vectors from derivatives of position vector x

            bases = numpy.array(MV.n*[ZERO],dtype=numpy.object)
            i = 0
            for coord in coords:
                ei = x.diff(coords[i])
                ei.set_name(base_name_lst[i])
                bases[i] = ei
                i += 1

            #Calculate normalizee basis vectors and basis vector magnitudes

            if debug:
                print 'Coordinate Generating Vector'
                print x
                print 'Basis Vectors'
                for base in bases:
                    print base

        else:

            #Input basis vectors as N vector fields

            bases = x

            for (base,name) in zip(bases,base_name_lst):
                base.set_name(name)

            if debug:
                print 'Basis Vectors'
                for base in bases:
                    print base
                if debug_level == 1:
                    return

        if debug_level == 1:
            return

        #Calculate normalized basis vectors and magnitudes of
        #unormalized basis vectors

        (nbases,mags) = normalize(bases,nbase_name_lst)

        if debug:
            print 'Magnitudes'
            print '\\abs{'+LaTeX_base+'_{i}} = ',mags
            print 'Normalized Basis Vectors'
            for nbase in nbases:
                print nbase

        g =  numpy.array(MV.n*[MV.n*[ZERO]],dtype=numpy.object)

        for irow in MV.nrg:
            for icol in MV.nrg:
                magsq = sympy.expand((nbases[irow]|nbases[icol])())
                g[irow][icol]  = sympy.simplify(sympy.trigsimp(magsq,deep=True,recursive=True))

        if debug:
            print 'Metric $\\hat{g}_{ij} = \\hat{'+LaTeX_base+'}_{i}\\cdot \\hat{'+\
                  LaTeX_base+'}_{j}$'
            print r'\hat{g}_{ij} =',sympy.galgebra.latex_ex.LaTeX(g)

        if debug_level == 2:
            return

        #Calculate reciprocal normalized basis vectors

        rnbases = []

        if is_quasi_unit_numpy_array(g):
            ibasis = 0
            while ibasis < MV.n:
                base = g[ibasis][ibasis]*nbases[ibasis]
                base.set_name(rnbase_name_lst[ibasis])
                base.simplify()
                base.trigsimp()
                rnbases.append(base)
                ibasis += 1
        else:
            rnbases = reciprocal_frame(nbases,rnbase_name_lst)
            ibase = 0
            for base in rnbases:
                base.simplify()
                base.trigsimp()
                rnbases[ibase] = base
                ibase += 1

        if debug:
            if debug_level != 0:
                sympy.galgebra.latex_ex.MV_format(1)
            print 'Reciprocal Normalized Basis Vectors'
            for rnbase in rnbases:
                print rnbase

        if debug_level == 3:
            return

        #Calculate components of inverse vectors

        Acoef = []

        for ibasis in MV.nrg:
            evec = numpy.array(MV.n*[ZERO],dtype=numpy.object)
            for jbasis in MV.nrg:
                evec[jbasis] = (MV.bvec[ibasis]|rnbases[jbasis])()
            Acoef.append(evec)

        #Calculat metric tensors

        gr = numpy.array(MV.n*[MV.n*[ZERO]],dtype=numpy.object)

        for irow in MV.nrg:
            for icol in MV.nrg:
                magsq = sympy.expand((rnbases[irow]|rnbases[icol])())
                gr[irow][icol] = sympy.simplify(sympy.trigsimp(magsq,deep=True,recursive=True))

        if debug:
            print 'Metric $\\hat{g}^{ij} = \\hat{'+LaTeX_base+'}^{i}\\cdot \\hat{'+\
                  LaTeX_base+'}^{j}$'
            print r'\hat{g}^{ij} =',sympy.galgebra.latex_ex.LaTeX(gr)

        if debug_level == 4:
            return

        #Calculate bases and reciprocal bases for curvilinear mulitvectors

        MV_bases = [[ONE]]
        MV_rbases = [[ONE]]
        igrade = 1
        while igrade <= MV.n:
            base_index = MV.gabasis[igrade]
            nbase_index = len(base_index)
            grade_bases = []
            rgrade_bases = []
            for index in base_index:
                base = build_base(index,nbases)
                base.simplify()
                base.trigsimp()
                rbase = build_base(index,rnbases,True)
                rbase.simplify()
                rbase.trigsimp()
                grade_bases.append(base)
                rgrade_bases.append(rbase)
            igrade += 1
            MV_bases.append(grade_bases)
            MV_rbases.append(rgrade_bases)

        #Calculate connection multivectors for geometric derivative

        MV_connect = [[ZERO]]
        igrade = 1
        while igrade <= MV.n:
            grade_connect = []
            ibase = 0
            for base in MV_bases[igrade]:
                sum = MV()
                itheta = 0
                for (theta,etheta) in zip(coords,rnbases):
                    psum = (1/mags[itheta])*etheta*base.diff(theta)
                    psum.trigsimp()
                    sum += psum
                    itheta += 1
                sum.simplify()
                sum.trigsimp()
                grade_connect.append(sum)
                ibase += 1
            MV_connect.append(grade_connect)
            igrade += 1

        if debug:
            print 'Curvilinear Bases: $'+bstr+' = '+bmhat+'_{i_{1}}\\W\\dots\\W'+bmhat+'_{i_{R}}$'
            igrade = 1
            for grade in MV_bases[1:]:
                ibase = 0
                for base in grade:
                    index = MV.gabasis[igrade][ibase]
                    sub_str = ''
                    for i in index:
                        sub_str += sympy.galgebra.latex_ex.LatexPrinter.extended_symbol(coords_lst[i])
                    base_str = bmhat+'_{['+sub_str+']} = '
                    print base_str,base
                    ibase += 1
                igrade += 1

        if debug_level == 5:
            return

        #Calculate representation of connection multivectors in curvilinear system

        MV_Connect = [[ZERO]]
        igrade = 1
        while igrade <= MV.n:
            grade_connect = []
            ibase = 0
            nbase = len(MV_bases[igrade])

            if igrade < MV.n:
                ibase = 0
                p1base = len(MV_bases[igrade+1])
                m1base = len(MV_bases[igrade-1])
                while ibase < nbase:
                    Cm1 =  numpy.array(m1base*[ZERO],dtype=numpy.object)
                    Cp1 =  numpy.array(p1base*[ZERO],dtype=numpy.object)
                    C = MV_connect[igrade][ibase]
                    if igrade == 1:
                        X = C(0)
                    else:
                        X = MV()
                    jbase = 0
                    while jbase < m1base:
                        Cm1[jbase] = sympy.trigsimp((MV.inner_product(MV_rbases[igrade-1][jbase],C))(),deep=True,recursive=True)
                        jbase += 1
                    jbase = 0
                    while jbase < p1base:
                        Cp1[jbase] = sympy.trigsimp((MV.inner_product(MV_rbases[igrade+1][jbase],C))(),deep=True,recursive=True)
                        jbase += 1
                    X += MV((igrade-1,Cm1),'grade')+MV((igrade+1,Cp1),'grade')
                    X.simplify()
                    X.trigsimp()
                    grade_connect.append(X)
                    ibase += 1
            else:
                ibase = 0
                m1base = len(MV_bases[igrade-1])
                while ibase < nbase:
                    Cm1 =  numpy.array(m1base*[ZERO],dtype=numpy.object)
                    C = MV_connect[igrade][ibase]
                    jbase = 0
                    while jbase < m1base:
                        Cm1[jbase] = sympy.trigsimp((MV.inner_product(MV_rbases[igrade-1][jbase],C))(),deep=True,recursive=True)
                        jbase += 1
                    X = MV()
                    X.mv[MV.n-1] = Cm1
                    X.simplify()
                    X.trigsimp()
                    grade_connect.append(X)
                    ibase += 1

            MV_Connect.append(grade_connect)
            igrade += 1

        base_str = ''
        for coord in coords:
            base_str += base_name+'_'+sympy.galgebra.latex_ex.LatexPrinter.str_basic(coord)+' '
        base_str = base_str[:-1]

        old_names = MV.vbasis

        MV.setup(base_str,g,True,coords)

        MV.curvilinear_flg = True
        MV.Connect = MV_Connect
        sympy.galgebra.latex_ex.LatexPrinter.latex_bases()
        MV.Rframe = numpy.array(MV.n*[ZERO],dtype=numpy.object)
        ibasis = 0
        while ibasis < MV.n:
            base = MV()
            jbasis = 0
            while jbasis < MV.n:
                base.add_in_place(gr[ibasis][jbasis]*MV.bvec[jbasis])
                jbasis += 1
            base.scalar_mul_inplace(1/mags[ibasis])
            MV.Rframe[ibasis] = base
            ibasis += 1

        MV.org_basis = []
        for ibasis in MV.nrg:
            evec = MV(Acoef[ibasis],'vector',old_names[ibasis])
            setattr(MAIN_PROGRAM,evec.name,evec)
            MV.org_basis.append(evec)

        if MV.coords[0] == sympy.Symbol('t'):
            MV.dedt = []
            for coef in dedt_coef:
                MV.dedt.append(MV(coef,'vector'))
        else:
            MV.dedt = None

        if debug:
            print 'Representation of Original Basis Vectors'
            for evec in MV.org_basis:
                 print evec

            print 'Renormalized Reciprocal Vectors '+\
                  '$\\bfrac{'+bmhat+'^{k}}{\\abs{\\bm{'+LaTeX_base+'}_{k}}}$'

            ibasis = 0
            while ibasis < MV.n:
                c_str = sympy.galgebra.latex_ex.LatexPrinter.extended_symbol(coords_lst[ibasis])
                print '\\bfrac{\\bm{\\hat{'+LaTeX_base+\
                       '}}^{'+c_str+'}}{\\abs{\\bm{'+LaTeX_base+\
                       '}_{'+c_str+'}}} =',MV.Rframe[ibasis]
                ibasis += 1

            title_str = 'Connection Multivectors: $C\\lbrc'+bstr+'\\rbrc = '+\
                        '\\bfrac{'+bmhat+'^{k}}{\\abs{'+bmhat+\
                        '_{k}}}\\pdiff{'+bstr+'}{\\theta^{k}}$'

            print title_str
            igrade = 1
            for grade in MV.Connect[1:]:
                ibase = 0
                for base in grade:
                    index = MV.gabasis[igrade][ibase]
                    sub_str = ''
                    for i in index:
                        sub_str += sympy.galgebra.latex_ex.LatexPrinter.extended_symbol(coords_lst[i])

                    base_str = 'C\\lbrc\\hat{'+LaTeX_base+'}_{['+sub_str+']}\\rbrc = '
                    print base_str,base
                    ibase += 1
                igrade += 1
        return

    @staticmethod
    def print_blades():
        """
        Set multivector output to blade representation.
        """
        MV.bladeprint = 1
        return

    @staticmethod
    def print_bases():
        """
        Set multivector output to base representation.
        """
        MV.bladeprint = 0
        return

    @staticmethod
    def multiplication_table():
        """
        Calculate geometric product base multiplication table.
        See reference 5 section 3 for details.
        """
        MV.mtable = []
        for igrade in MV.n1rg:
            MV.mtable.append([])
            for ibase in range(MV.nbasis[igrade]):
                MV.mtable[igrade].append([])
                if igrade == 0:
                    base1 = []
                else:
                    base1 = MV.basis[igrade][ibase]
                for jgrade in MV.n1rg:
                    MV.mtable[igrade][ibase].append([])
                    for jbase in range(MV.nbasis[jgrade]):
                        if jgrade == 0:
                            base2 = []
                        else:
                            base2 = MV.basis[jgrade][jbase]
                        base = base1+base2
                        (coefs,bases) = MV.reduce_basis(base)
                        product = MV.convert(coefs,bases)
                        product.name = '('+MV.basislabel[igrade][ibase]+')('+\
                                       MV.basislabel[jgrade][jbase]+')'
                        MV.mtable[igrade][ibase][jgrade].append(product)
        if MV.debug:
            print 'Multiplication Table:'
            for level1 in MV.mtable:
                for level2 in level1:
                    for level3 in level2:
                        for mv in level3:
                            mv.printmv()
        return

    @staticmethod
    def geometric_product(mv1,mv2):
        """
        MV.geometric_product(mv1,mv2) calculates the geometric
        product the multivectors mv1 and mv2 (mv1*mv2). See
        reference 5 section 3.
        """
        product = MV()
        if isinstance(mv1,MV) and isinstance(mv2,MV):
            bladeflg1 = mv1.bladeflg
            bladeflg2 = mv2.bladeflg
            if bladeflg1:
                mv1.convert_from_blades()
            if bladeflg2:
                mv2.convert_from_blades()
            mul = MV()
            for igrade in MV.n1rg:
                gradei = mv1.mv[igrade]
                if isinstance(gradei,numpy.ndarray):
                    for ibase in range(MV.nbasis[igrade]):
                        xi = gradei[ibase]
                        if xi != ZERO:
                            for jgrade in MV.n1rg:
                                gradej = mv2.mv[jgrade]
                                if isinstance(gradej,numpy.ndarray):
                                    for jbase in range(MV.nbasis[jgrade]):
                                        xj = gradej[jbase]
                                        if xj != ZERO:
                                            xixj = MV.mtable[igrade][ibase][jgrade][jbase].scalar_mul(xi*xj)
                                            product.add_in_place(xixj)
            product.bladeflg = 0
            if bladeflg1:
                mv1.convert_to_blades()
            if bladeflg2:
                mv1.convert_to_blades()
            if bladeflg1 and bladeflg2:
                product.convert_to_blades()
        else:
            if isinstance(mv1,MV):
                product = mv1.scalar_mul(mv2)
            else:
                product = mv2.scalar_mul(mv1)
        return(product)

    @staticmethod
    def wedge(igrade1,blade1,vector2,name=''):
        """
        Calculate the outer product of a multivector blade
        and a vector.  See reference 5 section 5 for details.
        """
        w12 = blade1*vector2
        w21 = vector2*blade1
        if igrade1%2 != 0:
            w = w12-w21
        else:
            w = w12+w21
        w.name = name
        return(w*HALF)

    @staticmethod
    def blade_table():
        """
        Calculate basis blades in terms of bases. See reference 5
        section 5 for details. Used to convert from blade to base
        representation.
        """
        MV.btable = []
        MV.bladelabel = []
        basis_str = MV.basislabel[1]
        for igrade in MV.n1rg:
            MV.bladelabel.append([])
            if igrade == 0:
                MV.bladelabel[0].append('1')
                tmp = [MV(value=ONE,mvtype='scalar',mvname='1')]
            if igrade == 1:
                tmp = []
                for ibase in range(MV.nbasis[1]):
                    MV.bladelabel[1].append(basis_str[ibase])
                    tmp.append(MV(value=ibase,mvtype='basisvector',mvname=basis_str[ibase]))
            if igrade >= 2:
                tmp = []
                basis = MV.basis[igrade]
                nblades = MV.nbasis[igrade]
                iblade = 0
                for blade in basis:
                    name = ''
                    for i in blade:
                        name += basis_str[i]+'^'
                    name = name[:-1]
                    MV.bladelabel[igrade].append(name)
                    lblade = MV.basis[igrade-1].index(blade[:-1])
                    rblade = blade[-1]
                    igrade1 = igrade-1
                    blade1  = MV.btable[igrade1][lblade]
                    vector2 = MV.btable[1][rblade]
                    b1Wv2 = MV.wedge(igrade1,blade1,vector2,name)
                    tmp.append(b1Wv2)
            MV.btable.append(tmp)
        if MV.debug:
            print 'Blade Tabel:'
            for grade in MV.btable:
                for mv in grade:
                    print mv
            print 'Blade Labels:'
            print MV.bladelabel
        return

    @staticmethod
    def inverse_blade_table():
        """
        Calculate bases in terms of basis blades. See reference 5
        section 5 for details. Used to convert from base to blade
        representation.
        """
        MV.ibtable = []
        for igrade in MV.n1rg:
            if igrade == 0:
                tmp = [MV(value=ONE,mvtype='scalar',mvname='1')]
            if igrade == 1:
                tmp = []
                for ibase in range(MV.nbasis[1]):
                    tmp.append(MV(value=ibase,mvtype='basisvector'))
            if igrade >= 2:
                tmp = []
                iblade = 0
                for blade in MV.btable[igrade]:
                    invblade = MV()
                    for i in range(igrade-1):
                        invblade.mv[i] = -blade.mv[i]
                        invblade.mv[igrade] = +blade.mv[igrade]
                        invblade.bladeflg = 1
                    if igrade >= 4:
                        jgrade = igrade-2
                        while jgrade > 1:
                            for ibase in range(MV.nbasis[jgrade]):
                                invblade.substitute_base(jgrade,ibase,MV.ibtable[jgrade][ibase])
                            jgrade -= 2
                    invblade.name = MV.basislabel[igrade][iblade]
                    iblade += 1
                    tmp.append(invblade)
            MV.ibtable.append(tmp)
        if MV.debug:
            print 'Inverse Blade Tabel:'
            for grade in MV.ibtable:
                for mv in grade:
                    mv.printmv()
        return

    @staticmethod
    def outer_product(mv1,mv2):
        """
        MV.outer_product(mv1,mv2) calculates the outer (exterior,wedge)
        product of the multivectors mv1 and mv2 (mv1^mv2). See
        reference 5 section 6.
        """
        if type(mv1) == type(MV) and type(mv2) == type(MV):
            if mv1.is_scalar() and mv2.is_scalar():
                return(mv1*mv2)

        if isinstance(mv1,MV) and isinstance(mv2,MV):
            product = MV()
            product.bladeflg = 1
            mv1.convert_to_blades()
            mv2.convert_to_blades()
            for igrade1 in MV.n1rg:
                if isinstance(mv1.mv[igrade1],numpy.ndarray):
                    pg1 = mv1.project(igrade1)
                    for igrade2 in MV.n1rg:
                        igrade = igrade1+igrade2
                        if igrade <= MV.n:
                            if isinstance(mv2.mv[igrade2],numpy.ndarray):
                                pg2 = mv2.project(igrade2)
                                pg1pg2 = pg1*pg2
                                product.add_in_place(pg1pg2.project(igrade))
        else:
            if isinstance(mv1,MV):
                product = mv1.scalar_mul(mv2)
            if isinstance(mv2,MV):
                product = mv2.scalar_mul(mv1)
        return(product)

    @staticmethod
    def inner_product(mv1,mv2,mode='s'):
        """
        MV.inner_product(mv1,mv2) calculates the inner

        mode = 's' - symmetic (Doran & Lasenby)
        mode = 'l' - left contraction (Dorst)
        mode = 'r' - right contraction (Dorst)
        """
        if type(mv1) == type(MV) and type(mv2) == type(MV):
            if mv1.is_scalar() and mv2.is_scalar():
                return(mv1*mv2)

        if isinstance(mv1,MV) and isinstance(mv2,MV):
            product = MV()
            product.bladeflg = 1
            mv1.convert_to_blades()
            mv2.convert_to_blades()
            for igrade1 in range(MV.n1):
                if isinstance(mv1.mv[igrade1],numpy.ndarray):
                    pg1 = mv1.project(igrade1)
                    for igrade2 in range(MV.n1):
                        igrade = igrade1-igrade2
                        if mode == 's':
                            igrade = igrade.__abs__()
                        else:
                            if mode == 'l':
                                igrade = -igrade
                        if igrade >= 0:
                            if isinstance(mv2.mv[igrade2],numpy.ndarray):
                                pg2 = mv2.project(igrade2)
                                pg1pg2 = pg1*pg2
                                product.add_in_place(pg1pg2.project(igrade))
            return(product)
        else:
            if mode == 's':
                if isinstance(mv1,MV):
                    product = mv1.scalar_mul(mv2)
                if isinstance(mv2,MV):
                    product = mv2.scalar_mul(mv1)
            else:
                product = None
        return(product)


    @staticmethod
    def addition(mv1,mv2):
        """
        MV.addition(mv1,mv2) calculates the sum
        of the multivectors mv1 and mv2 (mv1+mv2).
        """
        sum = MV()
        if isinstance(mv1,MV) and isinstance(mv2,MV):
            if mv1.bladeflg or mv2.bladeflg:
                mv1.convert_to_blades()
                mv2.convert_to_blades()
                sum.bladeflg = 1
            for i in MV.n1rg:
                if isinstance(mv1.mv[i],numpy.ndarray) and isinstance(mv2.mv[i],numpy.ndarray):
                    sum.mv[i] = mv1.mv[i]+mv2.mv[i]
                else:
                    if isinstance(mv1.mv[i],numpy.ndarray) and not isinstance(mv2.mv[i],numpy.ndarray):
                        sum.mv[i] = +mv1.mv[i]
                    else:
                        if isinstance(mv2.mv[i],numpy.ndarray) and not isinstance(mv1.mv[i],numpy.ndarray):
                            sum.mv[i] = +mv2.mv[i]
            return(sum)
        else:
            if isinstance(mv1,MV):
                return(mv1+MV(mv2,'scalar'))
            else:
                return(MV(mv1,'scalar')+mv2)

    @staticmethod
    def subtraction(mv1,mv2):
        """
        MV.subtraction(mv1,mv2) calculates the difference
        of the multivectors mv1 and mv2 (mv1-mv2).
        """
        diff = MV()
        if isinstance(mv1,MV) and isinstance(mv2,MV):
            if mv1.bladeflg or mv2.bladeflg:
                mv1.convert_to_blades()
                mv2.convert_to_blades()
                diff.bladeflg = 1
            for i in MV.n1rg:
                if isinstance(mv1.mv[i],numpy.ndarray) and isinstance(mv2.mv[i],numpy.ndarray):
                    diff.mv[i] = mv1.mv[i]-mv2.mv[i]
                else:
                    if isinstance(mv1.mv[i],numpy.ndarray) and not isinstance(mv2.mv[i],numpy.ndarray):
                        diff.mv[i] = +mv1.mv[i]
                    else:
                        if not isinstance(mv1.mv[i],numpy.ndarray) and isinstance(mv2.mv[i],numpy.ndarray):
                            diff.mv[i] = -mv2.mv[i]
            return(diff)
        else:
            if isinstance(mv1,MV):
                return(mv1-MV(mv2,'scalar'))
            else:
                return(MV(mv1,'scalar')-mv2)

    @staticmethod
    def vdiff(vec,x):
        dvec = numpy.array(len(vec)*[ZERO])
        ivec = 0
        for veci in vec:
            dvec[ivec] = sympy.diff(veci,x)
            ivec += 1
        return(dvec)

    @staticmethod
    def scalar_to_symbol(scalar):
        if isinstance(scalar,MV):
            return(scalar.mv[0][0])
        if type(scalar) == types.ListType:
            sym = []
            for x in scalar:
                sym.append(MV.scalar_to_symbol(x))
            return(sym)
        return(scalar)

    def __init__(self,value='',mvtype='',mvname='',fct=False,vars=None):
        """
        Initialization of multivector X. Inputs are as follows

        mvtype           value                       result

        default         default                  Zero multivector
        'basisvector'   int i                    ith basis vector
        'basisbivector' int i                    ith basis bivector
        'scalar'        symbol x                 scalar of value x
                        string s
        'grade'         [int i, symbol array A]  X.grade(i) = A
                        [int i, string s]
        'vector'        symbol array A           X.grade(1) = A
                        string s
        'grade2'        symbol array A           X.grade(2) = A
                        string s
        'pseudo'        symbol x                 X.grade(n) = x
                        string s
        'spinor'        string s                 spinor with coefficients
                                                 s__indices and name sbm

        mvname is name of multivector.
        If fct is 'True' and MV.coords is defined in MV.setup then a
        multivector field of MV.coords is instanciated.
        """

        self.name      = mvname
        self.mv        = MV.n1*[0]
        self.bladeflg  = 0  #1 for blade expansion
        self.puregrade = 1
        if mvtype == 'basisvector':
            self.mv[1] = numpy.array(MV.nbasis[1]*[ZERO],dtype=numpy.object)
            self.mv[1][value] = ONE
        if mvtype == 'basisbivector':
            self.mv[2] = numpy.array(MV.nbasis[2]*[ZERO],dtype=numpy.object)
            self.mv[2][value] = ONE
        if mvtype == 'scalar':
            if isinstance(value,types.StringType):
                value = sympy.Symbol(value)
            if isinstance(value,types.IntType):
                value = sympy.Rational(value)
            self.mv[0] = numpy.array([value],dtype=numpy.object)
        if mvtype == 'pseudo':
            if isinstance(value,types.StringType):
                value = sympy.Symbol(value)
            self.mv[MV.n] = numpy.array([value],dtype=numpy.object)
        if mvtype == 'vector':
            if isinstance(value,types.StringType): #Most general vector
                symbol_str = ''
                for ibase in MV.nrg:
                    if MV.coords == None:
                        symbol = value+'__'+str(ibase+MV.index_offset)
                        symbol_str += symbol+' '
                    else:
                        symbol = value+'__'+(MV.coords[ibase]).name
                        symbol_str += symbol+' '
                symbol_lst = make_symbols(symbol_str)
                self.mv[1] = numpy.array(symbol_lst,dtype=numpy.object)
                self.name = value
            else:
                value = MV.pad_zeros(value,MV.nbasis[1])
                self.mv[1] = numpy.array(value,dtype=numpy.object)
        if mvtype == 'grade2':
            if isinstance(value,types.StringType): #Most general grade-2 multivector
                if value != '':
                    symbol_str = ''
                    for base in MV.basis[2]:
                        symbol = value+MV.construct_index(base)
                        symbol_str += symbol+' '
                    symbol_lst = make_symbols(symbol_str)
                    self.mv[2] = numpy.array(symbol_lst,dtype=numpy.object)
            else:
                value = MV.pad_zeros(value,MV.nbasis[2])
                self.mv[2] = numpy.array(value,dtype=numpy.object)
        if mvtype == 'grade':
            igrade = value[0]
            coefs  = value[1]
            if isinstance(coefs,types.StringType): #Most general pure grade multivector
                base_symbol = coefs
                coefs = []
                bases = MV.basis[igrade]
                if igrade == 0:
                    self.mv[0] = numpy.array([sympy.Symbol(base_symbol)],dtype=numpy.object)
                else:
                    for base in bases:
                        coef = base_symbol+MV.construct_index(base)
                        coef = sympy.Symbol(coef)
                        coefs.append(coef)
                    self.mv[igrade] = numpy.array(coefs,dtype=numpy.object)
            else:
                self.mv[igrade] = coefs
        if mvtype == 'base':
            self.mv[value[0]] = numpy.array(MV.nbasis[value[0]]*[ZERO],dtype=numpy.object)
            self.mv[value[0]][value[1]] = ONE
        if mvtype == 'spinor':
            if isinstance(value,types.StringType): #Most general spinor
                for grade in MV.n1rg:
                    if grade%2 == 0:
                        symbol_str = ''
                        if grade != 0:
                            for base in MV.basis[grade]:
                                symbol = value+MV.construct_index(base)
                                symbol_str += symbol+' '
                            symbol_lst = make_symbols(symbol_str)
                            self.mv[grade] = numpy.array(symbol_lst,dtype=numpy.object)
                        else:
                            self.mv[0] = numpy.array([sympy.Symbol(value)],dtype=numpy.object)
                self.name = value+'bm'
        if isinstance(value,types.StringType) and mvtype == '': #Most general multivector
            if value != '':
                for grade in MV.n1rg:
                    symbol_str = ''
                    if grade != 0:
                        for base in MV.basis[grade]:
                            symbol = value+MV.construct_index(base)
                            symbol_str += symbol+' '
                        symbol_lst = make_symbols(symbol_str)
                        self.mv[grade] = numpy.array(symbol_lst,dtype=numpy.object)
                    else:
                        self.mv[0] = numpy.array([sympy.Symbol(value)],dtype=numpy.object)
                self.name = value+'bm'
        if fct:
            if vars != None:
                vars = tuple(vars)
            for grade in MV.n1rg:
                if not isinstance(self.mv[grade],types.IntType):
                    if grade == 0:
                        coef = sympy.galgebra.latex_ex.LatexPrinter.str_basic(self.mv[0][0])
                        if vars == None and MV.coords != None:
                            self.mv[0]= numpy.array([sympy.Function(coef)(*MV.coords)],dtype=numpy.object)
                        else:
                            self.mv[0]= numpy.array([sympy.Function(coef)(*vars)],dtype=numpy.object)
                    else:
                        for base in range(MV.nbasis[grade]):
                            coef = sympy.galgebra.latex_ex.LatexPrinter.str_basic(self.mv[grade][base])
                            if vars == None and MV.coords != None:
                                self.mv[grade][base] = sympy.Function(coef)(*MV.coords)
                            else:
                                self.mv[grade][base] = sympy.Function(coef)(*vars)

    @staticmethod
    def construct_index(base):
        index_str = ''
        if len(base) == 0:
            return('')
        if MV.coords == None:
            for ix in base:
                index_str += str(ix+MV.index_offset)
        else:
            for ix in base:
                index_str += (MV.coords[ix]).name
        return('__'+index_str)

    def set_name(self,namestr):
        self.name = namestr
        return

    def max_grade(self):
        """
        X.max_grade() is maximum grade of non-zero grades of X.
        """
        for i in range(MV.n,-1,-1):
            if isinstance(self.mv[i],numpy.ndarray):
                return(i)
        return(-1)

    @staticmethod
    def coord(xname,offset=0):
        xi_str = ''
        for i in MV.nrg:
            xi_str += xname+str(i+offset)+' '
        xi = make_symbols(xi_str)
        x = MV(xi,'vector')
        return(x)

    def x(self,i):
        if isint(self.mv[1]):
            return(ZERO)
        return(self.mv[1][i])

    def set_coef(self,grade,base,value):
        if isinstance(self.mv[grade],types.IntType):
            self.mv[grade] = numpy.array(MV.nbasis[grade]*[ZERO],dtype=numpy.object)
        self.mv[grade][base] = value
        return

    @staticmethod
    def named(mvname,value='',mvtype=''):
        name = mvname
        tmp = MV(value=value,mvtype=mvtype,mvname=name)
        setattr(sys.modules[__name__],name,tmp)
        return

    @staticmethod
    def printnm(tpl):
        for a in tpl:
            print a.name,' =',a.mv
        return

    def __str__(self):
        return(MV.str_rep(self))

    def printmv(self,name=''):
        title = ''
        if name:
            title += name+' = '
        else:
            if self.name:
                title += self.name+' = '
        print title+MV.str_rep(self)
        return

    def set_value(self,igrade,ibase,value):
        if isinstance(self.mv[igrade],numpy.ndarray):
            self.mv[igrade][ibase] = value
        else:
            self.mv[igrade] = numpy.array(MV.nbasis[igrade]*[ZERO],dtype=numpy.object)
            self.mv[igrade][ibase] = value
        return

    def add_in_place(self,mv):
        """
        X.add_in_place(mv) increments multivector X by multivector
        mv.
        """
        if self.bladeflg or mv.bladeflg:
            self.convert_to_blades()
            mv.convert_to_blades()
        for i in MV.n1rg:
            if isinstance(self.mv[i],numpy.ndarray) and isinstance(mv.mv[i],numpy.ndarray):
                self.mv[i] += mv.mv[i]
            else:
                if not isinstance(self.mv[i],numpy.ndarray) and isinstance(mv.mv[i],numpy.ndarray):
                    self.mv[i] = +mv.mv[i]
        return

    def sub_in_place(self,mv):
        """
        X.sub_in_place(mv) decrements multivector X by multivector
        mv.
        """
        if self.bladeflg or mv.bladeflg:
            self.convert_to_blades()
            mv.convert_to_blades()
        for i in MV.n1rg:
            if isinstance(self.mv[i],numpy.ndarray) and isinstance(mv.mv[i],numpy.ndarray):
                self.mv[i] -= mv.mv[i]
            else:
                if not isinstance(self.mv[i],numpy.ndarray) and isinstance(mv.mv[i],numpy.ndarray):
                    self.mv[i] = +mv.mv[i]
        return

    def __pos__(self):
        p = MV()
        p.bladeflg = self.bladeflg
        for i in MV.n1rg:
            if isinstance(self.mv[i],numpy.ndarray):
                p.mv[i] = +self.mv[i]
        return(p)

    def __neg__(self):
        n = MV()
        n.bladeflg = self.bladeflg
        for i in MV.n1rg:
            if isinstance(self.mv[i],numpy.ndarray):
                n.mv[i] = -self.mv[i]
        return(n)

    def __add_ab__(self,mv):
        self.add_in_place(mv)
        return

    def __sub_ab__(self,mv):
        self.sub_in_place(mv)
        return

    def __add__(self,mv):
        """See MV.addition(self,mv)"""
        return(MV.addition(self,mv))

    def __radd__(self,mv):
        """See MV.addition(mv,self)"""
        return(MV.addition(mv,self))

    def __sub__(self,mv):
        """See MV.subtraction(self,mv)"""
        return(MV.subtraction(self,mv))

    def __rsub__(self,mv):
        """See MV.subtraction(mv,self)"""
        return(MV.subtraction(mv,self))

    def __xor__(self,mv):
        """See MV.outer_product(self,mv)"""
        return(MV.outer_product(self,mv))

    def __pow__(self,mv):
        """See MV.outer_product(self,mv)"""
        return(MV.outer_product(self,mv))

    def __rxor__(self,mv):
        """See MV.outer_product(mv,self)"""
        return(MV.outer_product(mv,self))

    def __or__(self,mv):
        """See MV.inner_product(self,mv)"""
        return(MV.inner_product(self,mv))

    def __ror__(self,mv):
        """See MV.inner_product(mv,self)"""
        return(MV.inner_product(mv,self))

    def __lt__(self,mv):
        """See MV.inner_product(self,mv)"""
        return(MV.inner_product(self,mv,'l'))

    def __lshift__(self,mv):
        """See MV.inner_product(self,mv)"""
        return(MV.inner_product(self,mv,'l'))

    def __rlshift__(self,mv):
        """See MV.inner_product(self,mv)"""
        return(MV.inner_product(mv,self,'l'))

    def lc(self,mv):
        return(MV.inner_product(self,mv,'l'))

    def __gt__(self,mv):
        """See MV.inner_product(self,mv)"""
        return(MV.inner_product(self,mv,'r'))

    def __rshift__(self,mv):
        """See MV.inner_product(self,mv)"""
        return(MV.inner_product(self,mv,'r'))

    def __rrshift__(self,mv):
        """See MV.inner_product(self,mv)"""
        return(MV.inner_product(mv,self,'r'))

    def rc(self,mv):
        return(MV.inner_product(self,mv,'r'))

    def scalar_mul(self,c):
        """
        Y = X.scalar_mul(c), multiply multivector X by scalar c and return
        result.
        """
        mv = MV()
        mv.bladeflg = self.bladeflg
        mv.puregrade = self.puregrade
        for i in MV.n1rg:
            if isinstance(self.mv[i],numpy.ndarray):
                #print self.mv[i]
                #print c,type(c)
                mv.mv[i] = self.mv[i]*c
        return(mv)

    def scalar_mul_inplace(self,c):
        """
        X.scalar_mul_inplace(c), multiply multivector X by scalar c and save
        result in X.
        """
        for i in MV.n1rg:
            if isinstance(self.mv[i],numpy.ndarray):
                self.mv[i] = self.mv[i]*c
        return

    def __mul__(self,mv):
        """See MV.geometric_product(self,mv)"""
        return(MV.geometric_product(self,mv))

    def __rmul__(self,mv):
        """See MV.geometric_product(mv,self)"""
        return(MV.geometric_product(mv,self))

    def __div__(self,scalar):
        div = MV()
        div.bladeflg = self.bladeflg
        for i in MV.n1rg:
            if isinstance(self.mv[i],numpy.ndarray):
                div.mv[i] = self.mv[i]/scalar
        return(div)

    def __div_ab__(self,scalar):
        for i in MV.n1rg:
            if isinstance(self.mv[i],numpy.ndarray):
                self.mv[i] /= scalar
        return

    def __call__(self,igrade=0,ibase=0):
        """
        X(i,j) returns symbol in ith grade and jth base or blade of
        multivector X.
        """
        if not isinstance(self.mv[igrade],numpy.ndarray):
            return(ZERO)
        return(self.mv[igrade][ibase])

    @staticmethod
    def equal(mv1,mv2):
        mv1.compact()
        if isinstance(mv2,MV):
            mv2.compact()
        pure_grade = mv1.is_pure()
        if not isinstance(mv2,MV) and pure_grade != 0:
            return(False)
        if not isinstance(mv2,MV) and pure_grade == 0:
            if isinstance(mv1.mv[0],types.IntType):
                return(mv2 == 0)
            else:
                return(mv1.mv[0][0] == mv2)
        for (mvi,mvj) in zip(mv1.mv,mv2.mv):
            if isint(mvi) ^ isint(mvj):
                return(False)
            if isinstance(mvi,numpy.ndarray) and isinstance(mvj,numpy.ndarray):
                for (x,y) in zip(mvi,mvj):
                    if x != y:
                        return(False)
        return(True)

    def __eq__(self,mv):
        return(MV.equal(self,mv))

    def copy(self,sub=0):
        """
        Y = X.copy(), make a deep copy of multivector X in multivector
        Y so that Y can be modified without affecting X.
        """
        cpy = MV()
        cpy.name      = self.name
        cpy.bladeflg  = self.bladeflg
        cpy.puregrade = self.puregrade
        for i in MV.n1rg:
            if sub:
                if isinstance(self.mv[i],numpy.ndarray):
                    cpy.mv[i] = -self.mv[i]
            else:
                if isinstance(self.mv[i],numpy.ndarray):
                    cpy.mv[i] = +self.mv[i]
        return(cpy)

    def substitute_base(self,igrade,base,mv):
        if not isinstance(self.mv[igrade],numpy.ndarray):
            return
        if isinstance(base,numpy.ndarray):
            ibase = MV.basis[igrade].index(base)
        else:
            ibase = base
        coef = self.mv[igrade][ibase]
        if coef == ZERO:
            return
        self.mv[igrade][ibase] = ZERO
        self.add_in_place(mv*coef)
        return

    def convert_to_blades(self):
        """
        X.convert_to_blades(), inplace convert base representation
        to blade representation.  See reference 5 section 5.
        """
        if self.bladeflg:
            return
        self.bladeflg = 1
        for igrade in range(2,MV.n1):
            if isinstance(self.mv[igrade],numpy.ndarray):
                for ibase in range(MV.nbasis[igrade]):
                    coef = self.mv[igrade][ibase]
                    if (not coef == ZERO):
                        self.mv[igrade][ibase] = ZERO
                        self.add_in_place(MV.ibtable[igrade][ibase]*coef)
        return

    def convert_from_blades(self):
        """
        X.convert_from_blades(), inplace convert blade representation
        to base representation.  See reference 5 section 5.
        """
        if not self.bladeflg:
            return
        self.bladeflg = 0
        for igrade in range(2,MV.n1):
            if isinstance(self.mv[igrade],numpy.ndarray):
                for ibase in range(MV.nbasis[igrade]):
                    coef = self.mv[igrade][ibase]
                    if (not coef == ZERO):
                        self.mv[igrade][ibase] = ZERO
                        self.add_in_place(MV.btable[igrade][ibase]*coef)
        return

    def project(self,r):
        """
        Grade projection operator. For multivector X, X.project(r)
        returns multivector of grade r components of X if r is an
        integer. If r is a multivector X.project(r) returns a
        mutivector consisting of the grade of X for which r has non-
        zero grades. For example if X is a general multvector and
        r is a general spinor then X.project(r) will return the even
        grades of X.
        """
        if isinstance(r,types.IntType):
            grade_r = MV()
            if r > MV.n:
                return(grade_r)
            self.convert_to_blades()
            if not isinstance(self.mv[r],numpy.ndarray):
                return(grade_r)
            grade_r.bladeflg = 1
            grade_r.puregrade = 1
            grade_r.mv[r] = +self.mv[r]
            return(grade_r)
        if isinstance(r,MV):
            self.convert_to_blades()
            r.convert_to_blades()
            proj = MV()
            for i in MV.n1rg:
                if not isinstance(r.mv[i],types.IntType):
                    proj.mv[i] = self.mv[i]
            proj.bladeflg = self.bladeflg
            return(proj)
        return(None)

    def even(self):
        """
        Even grade projection operator. For multivector X, X.even()
        returns multivector of even grade components of X.
        """
        egrades = MV()
        self.convert_to_blades()
        egrades.bladeflg  = self.bladeflg
        egrades.puregrade = self.puregrade
        for igrade in range(0,MV.n1,2):
            egrades.mv[igrade] = +self.mv[igrade]
        return(egrades)

    def odd(self):
        """
        Odd grade projection operator. For multivector X, X.odd()
        returns multivector of odd grade components of X.
        """
        ogrades = MV()
        self.convert_to_blades()
        ogrades.bladeflg  = self.bladeflg
        ogrades.puregrade = self.puregrade
        for igrade in range(1,MV.n1,2):
            ogrades.mv[igrade] = +self.mv[igrade]
        return(ogrades)

    def rev(self):
        """
        Revisioin operator. For multivector X, X.rev()
        returns reversed multivector of X.
        """
        revmv = MV()
        self.convert_to_blades()
        revmv.bladeflg  = self.bladeflg
        revmv.puregrade = self.puregrade
        for igrade in MV.n1rg:
            if isinstance(self.mv[igrade],numpy.ndarray):
                if igrade < 2 or (not (((igrade*(igrade-1))/2)%2)):
                    revmv.mv[igrade] = +self.mv[igrade]
                else:
                    revmv.mv[igrade] = -self.mv[igrade]
        return(revmv)

    def cse(self,grade):
        cse_lst = []
        if isinstance(self.mv[grade],numpy.ndarray):
            for ibase in range(MV.nbasis[grade]):
                if self.mv[grade][ibase] != ZERO:
                    cse_lst.append(sympy.cse(self.mv[grade][ibase]))
        return(cse_lst)

    def div(self,grade,divisor):
        div_lst = []
        if isinstance(self.mv[grade],numpy.ndarray):
            for ibase in range(MV.nbasis[grade]):
                if self.mv[grade][ibase] != ZERO:
                    div_lst.append(self.mv[grade][ibase].as_coefficient(divisor))
        return(div_lst)

    def collect(self,lst):
        """
        Applies sympy collect function to each component
        of multivector.
        """
        for igrade in MV.n1rg:
            if isinstance(self.mv[igrade],numpy.ndarray):
                for ibase in range(MV.nbasis[igrade]):
                    if self.mv[igrade][ibase] != ZERO:
                        self.mv[igrade][ibase] = \
                        sympy.collect(self.mv[igrade][ibase],lst)
        return

    def sqrfree(self,lst):
        """
        Applies sympy sqrfree function to each component
        of multivector.
        """
        for igrade in MV.n1rg:
            if isinstance(self.mv[igrade],numpy.ndarray):
                for ibase in range(MV.nbasis[igrade]):
                    if self.mv[igrade][ibase] != ZERO:
                        self.mv[igrade][ibase] = \
                        sympy.sqrfree(self.mv[igrade][ibase],lst)
        return

    def flatten(self):
        flst = []
        for igrade in MV.n1rg:
            if isinstance(self.mv[igrade],types.IntType):
                flst += MV.nbasis[igrade]*[ZERO]
            else:
                for coef in self.mv[igrade]:
                    flst.append(coef)
        return(flst)

    def subs(self,*args):
        X = MV()
        X.bladeflg  = self.bladeflg
        X.puregrade = self.puregrade
        for igrade in MV.n1rg:
            if isinstance(self.mv[igrade],numpy.ndarray):
                X.mv[igrade] = numpy.array(substitute_array(self.mv[igrade],*args))
        return(X)

    def sub_mv(self,mv1,mv2):
        mv1_flat = mv1.flatten()
        mv2_flat = mv2.flatten()
        self.sub_scalar(mv1_flat,mv2_flat)
        return

    def sub_scalar(self,expr1,expr2):
        if (isinstance(expr1,types.ListType) and isinstance(expr2,types.ListType)) or \
           (isinstance(expr1,types.TupleType) and isinstance(expr2,types.TupleType)):
            for (var1,var2) in zip(expr1,expr2):
                self.sub_scalar(var1,var2)
        else:
            for igrade in MV.n1rg:
                if not isinstance(self.mv[igrade],types.IntType):
                    for ibase in range(MV.nbasis[igrade]):
                        if expr1 != ZERO:
                            self.mv[igrade][ibase] = self.mv[igrade][ibase].subs(expr1,expr2)
        return

    def simplify(self):
        """
        Applies sympy simplify function
        to each component of multivector.
        """
        for igrade in MV.n1rg:
            if isinstance(self.mv[igrade],numpy.ndarray):
                for ibase in range(MV.nbasis[igrade]):
                    if self.mv[igrade][ibase] != ZERO:
                        self.mv[igrade][ibase] = \
                        sympy.simplify(self.mv[igrade][ibase])
        return

    def trigsimp(self):
        """
        Applies sympy trigsimp function
        to each component of multivector.
        """
        for igrade in MV.n1rg:
            if isinstance(self.mv[igrade],numpy.ndarray):
                for ibase in range(MV.nbasis[igrade]):
                    if self.mv[igrade][ibase] != ZERO:
                        self.mv[igrade][ibase] = \
                        sympy.trigsimp(self.mv[igrade][ibase],deep=True,recursive=True)
        return

    def cancel(self):
        """
        Applies sympy cancel function
        to each component of multivector.
        """
        for igrade in MV.n1rg:
            if isinstance(self.mv[igrade],numpy.ndarray):
                for ibase in range(MV.nbasis[igrade]):
                    if self.mv[igrade][ibase] != ZERO:
                        self.mv[igrade][ibase] = \
                        sympy.cancel(self.mv[igrade][ibase])
        return

    def trim(self):
        """
        Applies sympy trim function
        to each component of multivector.
        """
        for igrade in MV.n1rg:
            if isinstance(self.mv[igrade],numpy.ndarray):
                for ibase in range(MV.nbasis[igrade]):
                    if self.mv[igrade][ibase] != ZERO:
                        self.mv[igrade][ibase] = \
                        sympy.trim(self.mv[igrade][ibase])
        return

    def expand(self):
        """
        Applies sympy expand function to each component
        of multivector.
        """
        for igrade in MV.n1rg:
            if isinstance(self.mv[igrade],numpy.ndarray):
                for ibase in range(MV.nbasis[igrade]):
                    if self.mv[igrade][ibase] != ZERO:
                        self.mv[igrade][ibase] = \
                        self.mv[igrade][ibase].expand()
        return

    def is_pure(self):
        igrade = -1
        ngrade = 0
        self.compact()
        for i in MV.n1rg:
            if isinstance(self.mv[i],numpy.ndarray):
                for base in self.mv[i]:
                    if base != 0:
                        igrade = i
                        ngrade += 1
                        break
                if ngrade > 1:
                    return(-1)
        if igrade == -1:
            return(0)
        return(igrade)

    def is_scalar(self):
        if self.is_pure() == 0:
            return(True)
        return(False)

    def compact(self):
        """
        Convert zero numpy arrays to single interge zero place holder
        in grade list for instanciated multivector. For example if
        numpy array of grade one components is a zero array then replace
        with single integer equal to zero.
        """
        for i in MV.n1rg:
            if not isint(self.mv[i]):
                zero_flg = True
                for x in self.mv[i]:
                    if x != 0:
                        zero_flg = False
                        break
                if zero_flg:
                    self.mv[i] = 0
        return

    def diff(self,x):
        """
        Calculate partial derivative of multivector with respect to
        argument x.
        """
        D = MV()
        igrade = 0
        for grade in self.mv:
            if not isinstance(grade,types.IntType):
                D.mv[igrade] = MV.vdiff(grade,x)
            igrade += 1
        return(D)

    def ddt(self):
        if MV.coords[0] != sympy.Symbol('t'):
            return(MV())
        dxdt = self.diff(MV.coords[0])
        for ibase in MV.nrg:
            dxdt += self.mv[1][ibase]*MV.dedt[ibase]
        #dxdt.simplify()
        #dxdt.trigsimp()
        return(dxdt)

    def grad(self):
        """
        Calculate geometric (grad) derivative of multivector function
        """
        D = []
        dD = MV()
        for theta in MV.coords:
            D.append(self.diff(theta))
        if MV.curvilinear_flg:
            recp = MV.Rframe
        else:
            recp = MV.brecp
        for (rbase,iD) in zip(recp,D):
            dD.add_in_place(rbase*iD)
        if MV.curvilinear_flg: #Add Connection
            igrade = 1
            while igrade <= MV.n:
                coefs = self.mv[igrade]
                if type(coefs) != types.IntType:
                    for (coef,connect) in zip(coefs,MV.Connect[igrade]):
                        dD.add_in_place(coef*connect)
                igrade += 1
        return(dD)

    def grad_ext(self):
        """
        Calculate outer (exterior,curl) derivative of multivector function.
        """
        D = []
        dD = MV()
        for ix in MV.coords:
            D.append(self.diff(ix))
        if MV.curvilinear_flg:
            recp = MV.Rframe
        else:
            recp = MV.brecp
        for (irbase,iD) in zip(recp,D):
            dD.add_in_place(irbase^iD)
        if MV.curvilinear_flg: #Add Connection
            igrade = 1
            while igrade <= MV.n:
                coefs = self.mv[igrade]
                if type(coefs) != types.IntType:
                    for (coef,connect) in zip(coefs,MV.Connect[igrade]):
                        if igrade < MV.n:
                            dD.add_in_place(coef*connect.project(igrade+1))
                igrade += 1
        return(dD)

    def curl(self):
        return(self.grad_ext())

    def grad_int(self):
        """
        Calculate inner (interior,div) derivative of multivector function.
        """
        D = []
        dD = MV()
        for ix in MV.coords:
            D.append(self.diff(ix))
        if MV.curvilinear_flg:
            recp = MV.Rframe
        else:
            recp = MV.brecp
        for (irbase,iD) in zip(recp,D):
            dD.add_in_place(irbase|iD)
        if MV.curvilinear_flg: #Add Connection
            igrade = 1
            while igrade <= MV.n:
                coefs = self.mv[igrade]
                if type(coefs) != types.IntType:
                    for (coef,connect) in zip(coefs,MV.Connect[igrade]):
                        dD.add_in_place(coef*connect.project(igrade-1))
                igrade += 1
        return(dD)

    def div(self):
        return(self.grad_int())

    def mag2(self):
        """
        Calculate scalar component of square of multivector.
        """
        return((self|self)())

    def Name(self):
        """
        Get LaTeX name of multivector.
        """
        return(sympy.galgebra.latex_ex.LatexPrinter.extended_symbol(self.name))

def set_names(var_lst,var_str):
    """
    Set the names of a list of multivectors (var_lst) for a space delimited
    string (var_str) containing the names.
    """
    var_str_lst = var_str.split()
    if len(var_lst) == len(var_str_lst):
        for (var,var_name) in zip(var_lst,var_str_lst):
            var.set_name(var_name)
        return
    sys.stderr.write('Error in set_names. Lists incongruent!\n')
    sys.exit()
    return

def reciprocal_frame(vlst,names=''):
    """
    Calculate reciprocal frame of list (vlst) of vectors.  If desired name each
    vector in list of reciprocal vectors with names in space delimited string
    (names).
    """
    E = vlst[0]
    recp = []
    if type(names) != types.StringType:
        name_lst = names
    else:
        if names != '':
            name_lst = names.split()
    for i in range(1,MV.n):
        E = E^vlst[i]
    for i in range(MV.n):
        tmp = ONE
        if i%2 != 0:
            tmp = -ONE
        for j in range(MV.n):
            if i != j:
                tmp = tmp^vlst[j]
        tmp = tmp*E
        recp.append(tmp)
    Esq = sympy.trigsimp(E.mag2(),deep=True,recursive=True)
    print Esq
    print sympy.simplify(Esq)
    Esq_inv = ONE/Esq
    i = 0
    for i in range(MV.n):
        recp[i].trigsimp()
        recp[i] = recp[i]*Esq_inv
        if names != '':
            recp[i].set_name(name_lst[i])
            i += 1
    return(recp)

def S(value):
    return(MV(value,'scalar'))



