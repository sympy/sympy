#!/usr/bin/python

"""
The module symbolicGA implements symbolic Geometric Algebra in python.
The relevant references for this module are:

    1. "Geometric Algebra for Physicists" by C. Doran and A. Lazenby,
       Cambridge University Press, 2003.

    2. GiNaC 1.4.1, An open framwork for symbolic computation with the
       C++ programming language, 9/5/2007, http://www.ginac.de

    3. Symbolic Computation with GiNaC and Python 1.3.5, 2007-02-01,
       http://swiginac.berlios.de/ginac-tutorial.py.html

    4. Sympy Tutorial, http://code.google.com/p/sympy/wiki/Tutorial

    5. "Design of a Python Module for Symbolic Geometric Algebra
       Calculations" by Alan Bromborsky, included as symbolGA.pdf
"""

import numpy
import os,sys,string,types,re,copy
import sympy
#from latex_out import *

NUMPAT = re.compile( '([\-0-9])|([\-0-9]/[0-9])')
"""Re pattern for rational number"""

ZERO = sympy.Rational(0)
ONE  = sympy.Rational(1)
TWO  = sympy.Rational(2)
HALF = sympy.Rational(1,2)

global MAIN_PROGRAM

MAIN_PROGRAM = ''

def set_main(main_program):
    global MAIN_PROGRAM
    MAIN_PROGRAM = main_program
    return

def plist(lst):
    if type(lst) == types.ListType:
        for x in lst:
            plist(x)
    else:
        print lst
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


def symbol(sym_str):
    """
    Symbol converts a string to a sympy/sympy symbol.
    """
    sym = sympy.Symbol(sym_str)
    return(sym)

def expand(expr):
    return(sympy.expand(expr))


def collect(expr,lst):
    """
    Wrapper for sympy.collect and sympy.collect.
    See references 2, 3, and 4.
    """
    lst = MV.scalar_to_symbol(lst)
    return(sympy.collect(expr,lst))

def sqrfree(expr,lst):
    """
    Wrapper for sympy.sqrfree and sympy.factor.
    See references 2, 3, and 4.
    """
    return(sympy.sqrfree(expr,lst))

def collect_common_factors(expr):
    """
    Wrapper for sympy.collect_common_factors and sympy.factor.
    See references 2, 3, and 4.
    """
    return(sympy.collect_common_factors(expr))

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
        tmp = symbol(name)
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
        tmp = symbol(name)
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

class MV:
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
    pad_zeros = staticmethod(pad_zeros)

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
        MV.nbasis = numpy.array((MV.n+1)*[1],dtype=numpy.object)
        for igrade in range(1,MV.n+1):
            tmp = comb(MV.n,igrade)
            MV.gabasis += [tmp]
            ntmp = len(tmp)
            MV.nbasis[igrade] = ntmp
            MV.basis[igrade] = tmp
            gradelabels = []
            for i in range(ntmp):
                bstr = ''
                for j in tmp[i]:
                    bstr += MV.vbasis[j]
                gradelabels.append(bstr)
            MV.basislabel[igrade] = gradelabels
        if MV.debug:
            print 'basis strings =',MV.vbasis
            print 'basis symbols =',MV.vsyms
            print 'basis labels  =',MV.basislabel
            print 'basis         =',MV.basis
            print 'grades        =',MV.nbasis
            print 'index         =',MV.index
        return
    define_basis = staticmethod(define_basis)

    def define_metric(metric):
        """
        Calculates all the MV static variables needed for
        metric operations.  See reference 5 section 2.
        """
        name_flg = False
        MV.g = []
        MV.metric = numpy.array(MV.n*[MV.n*[ZERO]],dtype=numpy.object)
        if not metric:
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
                    tmp = symbol(gij)
                    MV.metric[i][j] = tmp
                    if i <= j:
                        MV.g.append(tmp)
                        if name_flg:
                            setattr(MAIN_PROGRAM,name,tmp)
                            name_flg = False
        if MV.debug:
            print 'metric =',MV.metric
        return
    define_metric = staticmethod(define_metric)

    def define_reciprocal_frame():
        """
        Calculates unscaled reciprocal vectors (MV.brecp) and scale
        factor (MV.Esq). The ith scaled reciprocal vector is
        (1/MV.Esq)*MV.brecp[i].  The pseudoscalar for the set of
        basis vectors is MV.E.
        """
        if MV.tables_flg:
            if MV.rframe_flg:
                return
            MV.rframe_flg = True
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
    define_reciprocal_frame = staticmethod(define_reciprocal_frame)

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
    reduce_basis_loop = staticmethod(reduce_basis_loop)

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
    reduce_basis = staticmethod(reduce_basis)

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
    combine_common_factors = staticmethod(combine_common_factors)

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
    contract = staticmethod(contract)

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
    convert = staticmethod(convert)

    def set_str_format(str_mode=0):
        MV.str_mode = str_mode
        return
    set_str_format = staticmethod(set_str_format)

    def str_rep(mv,lst_mode=0):
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
                                xstr = '+{'+xstr+'}'
                            else:
                                xstr = '+{'+xstr[1:]+'}'
                        value += xstr+labels[igrade][j]
                        if MV.str_mode and not lst_mode:
                            value += '\n'
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
        if MV.latexflg:
            value = LaTeXstr(value)
        return(value)
    str_rep = staticmethod(str_rep)

    def LaTeX():
        MV.latexflg = not MV.latexflg
        return
    LaTeX = staticmethod(LaTeX)

    def setup(basis,metric='',debug=0):
        """
        MV.setup initializes the MV class by calculating the static
        multivector tables required for geometric algebra operations
        on multivectors.  See reference 5 section 2 for details on
        basis and metric arguments.
        """
        global MAIN_PROGRAM
        MV.latexflg = False
        MV.debug = debug
        MV.bladeprint = 0
        MV.tables_flg = 0
        MV.str_mode  = 0
        MV.rframe_flg = False
        if type(basis) == types.StringType:
            basislst = basis.split()
            MV.define_basis(basislst)
        if len(metric) > 0 and type(metric) == types.StringType:
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
        return('Setup of '+basis+' complete!')
    setup = staticmethod(setup)

    def print_blades():
        """
        Set multivector output to blade representation.
        """
        MV.bladeprint = 1
        return
    print_blades=staticmethod(print_blades)

    def print_bases():
        """
        Set multivector output to base representation.
        """
        MV.bladeprint = 0
        return
    print_bases=staticmethod(print_bases)

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
                        product.name = '('+MV.basislabel[igrade][ibase]+')('+MV.basislabel[jgrade][jbase]+')'
                        MV.mtable[igrade][ibase][jgrade].append(product)
        if MV.debug:
            print 'Multiplication Table:'
            for level1 in MV.mtable:
                for level2 in level1:
                    for level3 in level2:
                        for mv in level3:
                            mv.printmv()
        return
    multiplication_table = staticmethod(multiplication_table)

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
    geometric_product = staticmethod(geometric_product)

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
    wedge = staticmethod(wedge)

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
    blade_table = staticmethod(blade_table)

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
    inverse_blade_table = staticmethod(inverse_blade_table)

    def outer_product(mv1,mv2):
        """
        MV.outer_product(mv1,mv2) calculates the outer (exterior,wedge)
        product of the multivectors mv1 and mv2 (mv1^mv2). See
        reference 5 section 6.
        """
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
            else:
                product = mv2.scalar_mul(mv1)
        return(product)
    outer_product = staticmethod(outer_product)

    def inner_product(mv1,mv2):
        """
        MV.inner_product(mv1,mv2) calculates the inner (scalar,dot)
        product of the multivectors mv1 and mv2 (mv1|mv2).  See
        reference 5 section 6.
        """
        if isinstance(mv1,MV) and isinstance(mv2,MV):
            product = MV()
            product.bladeflg = 1
            mv1.convert_to_blades()
            mv2.convert_to_blades()
            for igrade1 in range(1,MV.n1):
                if isinstance(mv1.mv[igrade1],numpy.ndarray):
                    pg1 = mv1.project(igrade1)
                    for igrade2 in range(1,MV.n1):
                        igrade = (igrade1-igrade2).__abs__()
                        if isinstance(mv2.mv[igrade2],numpy.ndarray):
                            pg2 = mv2.project(igrade2)
                            pg1pg2 = pg1*pg2
                            product.add_in_place(pg1pg2.project(igrade))
            return(product)
        return(MV())
    inner_product = staticmethod(inner_product)

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
        else:
            if isinstance(mv1,MV):
                sum = mv1.copy()
                if isinstance(sum.mv[0],numpy,ndarray):
                    sum.mv[0] += numpy.array([mv2],dtype=numpy.object)
                else:
                    sum.mv[0] = numpy.array([mv2],dtype=numpy.object)
            else:
                sum = mv2.copy()
                if isinstance(sum.mv[0],numpy.ndarray):
                    sum.mv[0] += numpy.array([mv1],dtype=numpy.object)
                else:
                    sum.mv[0] = numpy.array([mv1],dtype=numpy.object)
        return(sum)
    addition = staticmethod(addition)

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
        else:
            if isinstance(mv1,MV):
                diff = mv1.copy()
                if isinstandce(diff.mv[0],numpy.ndarray):
                    diff.mv[0] -= numpy.array([mv2],dtype=numpy.object)
                else:
                    diff.mv[0] = -numpy.array([mv2],dtype=numpy.object)
            else:
                diff = mv2.copy(1)
                if isinstance(diff.mv[0],numpy.ndarray):
                    diff.mv[0] += numpy.array([mv1],dtype=numpy.object)
                else:
                    diff.mv[0] = +numpy.array([mv1],dtype=numpy.object)
        return(diff)
    subtraction = staticmethod(subtraction)

    def scalar_to_symbol(scalar):
        if isinstance(scalar,MV):
            return(scalar.mv[0][0])
        if type(scalar) == types.ListType:
            sym = []
            for x in scalar:
                sym.append(MV.scalar_to_symbol(x))
            return(sym)
        return(scalar)
    scalar_to_symbol = staticmethod(scalar_to_symbol)

    def __init__(self,value='',mvtype='',mvname=''):
        """
        Initialization of multivector X. Inputs are as follows

        mvtype           value                       result

        default         default                  Zero multivector
        'basisvector'   int i                    ith basis vector
        'basisbivector' int i                    ith basis bivector
        'scalar'        symbol x                 scalar of value x
        'grade'         [int i, symbol array A]  X.grade(i) = A
        'vector'        symbol array A           X.grade(1) = A
        'grade2'        symbol array A           X.grade(2) = A

        mvname is name of multivector.
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
            self.mv[0] = numpy.array([value],dtype=numpy.object)
        if mvtype == 'grade':
            igrade          = value[0]
            coefs           = value[1]
            coefs = MV.pad_zeros(coefs,MV.nbasis[igrade])
            self.mv[igrade] = numpy.array(coefs,dtype=numpy.object)
        if mvtype == 'vector':
            value = MV.pad_zeros(value,MV.nbasis[1])
            self.mv[1]    = numpy.array(value,dtype=numpy.object)
        if mvtype == 'grade2':
            value = MV.pad_zeros(value,MV.nbasis[2])
            self.mv[2]    = numpy.array(value,dtype=numpy.object)

    def max_grade(self):
        """
        X.max_grade() is maximum grade of non-zero grades of X.
        """
        for i in range(MV.n,-1,-1):
            if isinstance(self.mv[i],numpy.ndarray):
                return(i)
        return(-1)

    def coord(xname,offset=0):
        xi_str = ''
        for i in MV.nrg:
            xi_str += xname+str(i+offset)+' '
        xi = make_symbols(xi_str)
        x = MV(xi,'vector')
        return(x)
    coord = staticmethod(coord)

    def x(self,i):
        if isint(self.mv[1]):
            return(ZERO)
        return(self.mv[1][i])

    def named(mvname,value='',mvtype=''):
        name = mvname
        tmp = MV(value=value,mvtype=mvtype,mvname=name)
        setattr(sys.modules[__name__],name,tmp)
        return
    named=staticmethod(named)

    def printnm(tpl):
        for a in tpl:
            print a.name,' =',a.mv
        return
    printnm = staticmethod(printnm)

    def __str__(self):
        """See MV.str_rep(self)"""
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
        return(mv)

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

    def __eq__(self,mv):
        if not isinstance(mv,MV):
            return(False)
        for (mvi,mvj) in zip(self.mv,mv.mv):
            if isint(mvi) ^ isint(mvj):
                return(False)
            if isinstance(mvi,numpy.ndarray) and isinstance(mvj,numpy.ndarray):
                for (x,y) in zip(mvi,mvj):
                    if x != y:
                        return(False)
        return(True)

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
        returns multivector of grade r components of X.
        """
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

    def collect(self,lst):
        """
        Applies sympy/sympy collect function to each component
        of multivector.
        """
        for igrade in MV.n1rg:
            if isinstance(self.mv[igrade],numpy.ndarray):
                for ibase in range(MV.nbasis[igrade]):
                    if self.mv[igrade][ibase] != ZERO:
                        self.mv[igrade][ibase] = \
                        collect(self.mv[igrade][ibase],lst)
                        self.mv[igrade][ibase] = \
                        collect(self.mv[igrade][ibase],lst)
        return

    def sqrfree(self,lst):
        """
        Applies sympy/sympy sqrfree function to each component
        of multivector.
        """
        for igrade in MV.n1rg:
            if isinstance(self.mv[igrade],numpy.ndarray):
                for ibase in range(MV.nbasis[igrade]):
                    if self.mv[igrade][ibase] != ZERO:
                        self.mv[igrade][ibase] = \
                        sqrfree(self.mv[igrade][ibase],lst)
        return

    def collect(self,faclst):
        """
        Applies sympy collect_common_factors function
        to each component of multivector.
        """
        for igrade in MV.n1rg:
            if isinstance(self.mv[igrade],numpy.ndarray):
                for ibase in range(MV.nbasis[igrade]):
                    if self.mv[igrade][ibase] != ZERO:
                        self.mv[igrade][ibase] = \
                        sympy.collect(self.mv[igrade][ibase],faclst)
        return

    def subs(self,var,substitute):
        """
        Applies sympy subs function
        to each component of multivector.
        """
        for igrade in MV.n1rg:
            if isinstance(self.mv[igrade],numpy.ndarray):
                for ibase in range(MV.nbasis[igrade]):
                    if self.mv[igrade][ibase] != ZERO:
                        self.mv[igrade][ibase] = \
                        self.mv[igrade][ibase].subs(var,substitute)
        return

    def simplify(self):
        """
        Applies sympy subs function
        to each component of multivector.
        """
        for igrade in MV.n1rg:
            if isinstance(self.mv[igrade],numpy.ndarray):
                for ibase in range(MV.nbasis[igrade]):
                    if self.mv[igrade][ibase] != ZERO:
                        self.mv[igrade][ibase] = \
                        sympy.simplify(self.mv[igrade][ibase])
        return

    def cancel(self):
        """
        Applies sympy cancle function
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
        Applies sympy cancle function
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
        Applies sympy/sympy expand function to each component
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
        for i in MV.n1rg:
            if isinstance(self.mv[i],numpy.ndarray):
                igrade = i
                ngrade += 1
                if ngrade > 1:
                    return(-1)
        return(igrade)

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

##class LT:
##    def __init__(self,fct_lst):
##        """
##        Initialize linear transformation.  fct_lst is a list of vectors A_{i}
##        such that A_{i} = LT(a_{i}) where the a_{i} are the basis vectos of
##        the MV class.
##        """
##        self.fct_lst = fct_lst
##        self.norm = ONE
##        self.ltblades = [[],self.fct_lst]
##        for igrade in range(2,MV.n1):
##            tmp = []
##            for ibasis in range(MV.nbasis[igrade]):
##                ilst = MV.basis[igrade][ibasis][:-1]
##                jlst = MV.basis[igrade-1]
##                iblade = jlst.index(ilst)
##                jblade = MV.basis[igrade][ibasis][-1]
##                tmp.append(self.ltblades[igrade-1][iblade]^self.fct_lst[jblade])
##            self.ltblades.append(tmp)
##        return
##
##    def Fvector(self,mv):
##        """
##        Return the linear transformation of an arbitrary vector.
##        """
##        mv.convert_to_blades()
##        if mv.is_pure() == 1:
##            fofa = MV()
##            fofa.bladeflg = 1
##            for i in range(MV.nbasis[1]):
##                fofa += mv.mv[1][i]*self.fct_lst[i]
##            fofa.collect_common_factors()
##            return(fofa/self.norm)
##        print 'Error in LT.F, mv =',mv,' not a vector!'
##        return
##
##    def __call__(self,mv):
##        """
##        Return the linear transformation of an arbitrary multivector.
##        """
##        if mv.is_pure() == 1:
##            return(self.Fvector(mv))
##        fofx = MV()
##        mv.convert_to_blades()
##        fofx.bladeflg = 1
##        for igrade in range(1,MV.n1):
##            if isinstance(mv.mv[igrade],numpy.ndarray):
##                for ibase in range(MV.nbasis[igrade]):
##                    if mv.mv[igrade][ibase] != ZERO:
##                        fofx.add_in_place(mv.mv[igrade][ibase]*self.ltblades[igrade][ibase])
##        fofx.collect_common_factors()
##        return(fofx/self.norm)
##
##    def __str__(self):
##        LTstr = ''
##        for ibasis in range(MV.n):
##            if self.norm != ONE:
##                LTstr += 'F('+MV.basislabel[1][ibasis]+') = ('
##            else:
##                LTstr += 'F('+MV.basislabel[1][ibasis]+') = '
##            LTstr += MV.str_rep(self.fct_lst[ibasis])
##            if self.norm != ONE:
##                LTstr += ')/('+str(self.norm)+')\n'
##            else:
##                LTstr += '\n'
##        return(LTstr[:-1])
##
##    def adj(self):
##        MV.define_reciprocal()
##        adj = []
##        norm = MV.Esq
##        norm = norm*self.norm
##        for ibasis in range(MV.n):
##            tmp = MV()
##            tmp.bladeflg = 1
##            for jbasis in range(MV.n):
##                tmp.add_in_place(MV.brecp[jbasis]*(MV.bvec[ibasis]|self.fct_lst[jbasis]))
##            tmp.expand()
##            tmp.sqrfree([])
##            adj.append(tmp)
##        adj = LT(adj)
##        adj.norm = norm
##        return(adj)
