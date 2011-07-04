#TODO:
# -Implement Clebsch-Gordan symmetries
# -Improve simplification method
# -Implement new simpifications
"""Clebsch-Gordon Coefficients."""

from sympy import Expr, Add, Function, Mul, Pow, sqrt, Sum, symbols, sympify, Wild
from sympy.printing.pretty.stringpict import prettyForm, stringPict

from sympy.physics.quantum.kronecker import KroneckerDelta
from sympy.physics.wigner import wigner_3j, clebsch_gordan


__all__ = [
    'Wigner3j',
    'CG',
    'cg_simp'
]

#-----------------------------------------------------------------------------
# CG Coefficients
#-----------------------------------------------------------------------------

class Wigner3j(Expr):
    """Class for the Wigner-3j symbols

    Wigner 3j-symbols are coefficients determined by the coupling of
    two angular momenta. When created, they are expressed as symbolic
    quantities that can be evaluated using the doit() method.

    Parameters
    ==========

    j1, m1, j2, m2, j3, m3 : Number, Symbol
        Terms determining the angular momentum of coupled angular momentum
        systems.

    Examples
    ========

    Declare a Wigner-3j coefficient and calcualte its value

        >>> from sympy.physics.quantum.cg import Wigner3j
        >>> w3j = Wigner3j(6,0,4,0,2,0)
        >>> w3j
        (6, 4, 2)
        (0, 0, 0)
        >>> w3j.doit()
        715**(1/2)/143

    References
    ==========

    [1] Varshalovich, D A, Quantum Theory of Angular Momentum. 1988.
    """
    def __new__(cls, j1, m1, j2, m2, j3, m3):
        j1,m1,j2,m2,j3,m3 = map(sympify, (j1,m1,j2,m2,j3,m3))
        return Expr.__new__(cls, j1, m1, j2, m2, j3, m3)

    @property
    def j1(self):
        return self.args[0]

    @property
    def m1(self):
        return self.args[1]

    @property
    def j2(self):
        return self.args[2]

    @property
    def m2(self):
        return self.args[3]

    @property
    def j3(self):
        return self.args[4]

    @property
    def m3(self):
        return self.args[5]

    @property
    def is_symbolic(self):
        return not (self.j1.is_number and self.j2.is_number and self.j3.is_number and
            self.m1.is_number and self.m2.is_number and self.m3.is_number)

    # This is modified from the _print_Matrix method
    def _sympystr(self, printer, *args):
        res = [[printer._print(self.j1), printer._print(self.j2), printer._print(self.j3)], \
            [printer._print(self.m1), printer._print(self.m2), printer._print(self.m3)]]
        maxw = [-1] * 3
        for j in range(3):
            maxw[j] = max([ len(res[i][j]) for i in range(2) ])
        for i, row in enumerate(res):
            for j, elem in enumerate(row):
                row[j] = elem.rjust(maxw[j])
            res[i] = "(" + ", ".join(row) + ")"
        return '\n'.join(res)

    # This is modified from the _print_Matrix method
    def _pretty(self, printer, *args):
        m = ((printer._print(self.j1), printer._print(self.m1)), \
            (printer._print(self.j2), printer._print(self.m2)), \
            (printer._print(self.j3), printer._print(self.m3)))
        hsep = 2
        vsep = 1
        maxw = [-1] * 3
        for j in range(3):
            maxw[j] = max([ m[j][i].width() for i in range(2) ])
        D = None
        for i in range(2):
            D_row = None
            for j in range(3):
                s = m[j][i]
                wdelta = maxw[j] - s.width()
                wleft  = wdelta //2
                wright = wdelta - wleft

                s = prettyForm(*s.right(' '*wright))
                s = prettyForm(*s.left(' '*wleft))

                if D_row is None:
                    D_row = s
                    continue
                D_row = prettyForm(*D_row.right(' '*hsep))
                D_row = prettyForm(*D_row.right(s))
            if D is None:
                D = D_row
                continue
            for _ in range(vsep):
                D = prettyForm(*D.below(' '))
            D = prettyForm(*D.below(D_row))
        D = prettyForm(*D.parens())
        return D

    def _latex(self, printer, *args):
        return r'\left(\begin{array}{ccc} %s & %s & %s \\ %s & %s & %s \end{array}\right)' % \
            (printer._print(self.j1), printer._print(self.j2), printer._print(self.j3), \
            printer._print(self.m1), printer._print(self.m2), printer._print(self.m3))

    def doit(self, **hints):
        if self.is_symbolic:
            raise ValueError("Coefficients must be numerical")
        return wigner_3j(self.j1, self.j2, self.j3, self.m1, self.m2, self.m3)


class CG(Wigner3j):
    """Class for Clebsch-Gordan coefficient

    Clebsch-Gordan coefficients describe the angular momentum coupling between
    two systems. The coefficients give the expansion of a coupled total angular
    momentum state and an uncoupled tensor product state. The Clebsch-Gordan
    coefficients are defined as:
    CG(j1,m1,j2,m2,j3,m3) = <j1,m1; j2,m2 | j3,m3>

    Parameters
    ==========

    j1, m1, j2, m2, j3, m3 : Number, Symbol
        Terms determining the angular momentum of coupled angular momentum
        systems.

    Examples
    ========

    Define a Clebsch-Gordan coefficient and evaluate its value

        >>> from sympy.physics.quantum.cg import CG
        >>> from sympy import S
        >>> cg = CG(S(3)/2, S(3)/2, S(1)/2, -S(1)/2, 1, 1)
        >>> cg
        CG(3/2, 3/2, 1/2, -1/2, 1, 1)
        >>> cg.doit()
        3**(1/2)/2

    References
    ==========

    [1] Varshalovich, D A, Quantum Theory of Angular Momentum. 1988.
    """

    def doit(self, **hints):
        if self.is_symbolic:
            raise ValueError("Coefficients must be numerical")
        return clebsch_gordan(self.j1,self.j2, self.j3, self.m1, self.m2, self.m3)

    def _sympystr(self, printer, *args):
        return 'CG(%s, %s, %s, %s, %s, %s)' % \
            (printer._print(self.j1), printer._print(self.m1), printer._print(self.j2), \
            printer._print(self.m2), printer._print(self.j3), printer._print(self.m3))

    def _pretty(self, printer, *args):
        bot = printer._print(self.j1)
        bot = prettyForm(*bot.right(','))
        bot = prettyForm(*bot.right(printer._print(self.m1)))
        bot = prettyForm(*bot.right(','))
        bot = prettyForm(*bot.right(printer._print(self.j2)))
        bot = prettyForm(*bot.right(','))
        bot = prettyForm(*bot.right(printer._print(self.m2)))
        top = printer._print(self.j3)
        top = prettyForm(*top.right(','))
        top = prettyForm(*top.right(printer._print(self.m3)))

        pad = max(top.width(), bot.width())

        bot = prettyForm(*bot.left(' '))
        top = prettyForm(*top.left(' '))
        if not pad == bot.width():
            bot = prettyForm(*bot.right(' ' * (pad-bot.width())))
        if not pad == top.width():
            top = prettyForm(*top.right(' ' * (pad-top.width())))
        s = stringPict('C' + ' '*pad)
        s = prettyForm(*s.below(bot))
        s = prettyForm(*s.above(top))
        return s

    def _latex(self, printer, *args):
        return r'C^{%s,%s}_{%s,%s,%s,%s}' % \
            (printer._print(self.j3), printer._print(self.m3),
            printer._print(self.j1), printer._print(self.m1),
            printer._print(self.j2), printer._print(self.m2))


def cg_simp(e):
    """Simplify and combine CG coefficients

    This function uses various symmetry and properties of sums and
    products of Clebsch-Gordan coefficients to simplify statements
    involving these terms

    Examples
    ========

    Simplify the sum over CG(a,alpha,0,0,a,alpha) for all alpha to
    2*a+1

        >>> from sympy.physics.quantum.cg import CG, cg_simp
        >>> a = CG(1,1,0,0,1,1)
        >>> b = CG(1,0,0,0,1,0)
        >>> c = CG(1,-1,0,0,1,-1)
        >>> cg_simp(a+b+c)
        3

    References
    ==========

    [1] Varshalovich, D A, Quantum Theory of Angular Momentum. 1988.
    """
    if isinstance(e, Add):
        return _cg_simp_add(e)
    elif isinstance(e, Sum):
        return _cg_simp_sum(e)
    elif isinstance(e, Mul):
        return Mul(*[cg_simp(arg) for arg in e.args])
    elif isinstance(e, Pow):
        return Pow(cg_simp(e.base), e.exp)
    else:
        return e


def _cg_simp_add(e):
    #TODO: Improve simplification method
    """Takes a sum of terms involving Clebsch-Gordan coefficients and
    simplifies the terms.

    This method creates lists of the terms in each of the arguments of the
    summation. Each term is then parsed for the Clebsch-Gordan coefficients in
    it and, where applicable, apply symmetries to the simplify the terms in
    the sum.

    References
    ==========

    [1] Varshalovich, D A, Quantum Theory of Angular Momentum. 1988.
    """
    cg_part = []
    other_part = []
    for arg in e.args:
        if arg.has(CG):
            if isinstance(arg, CG):
                terms = [arg]
            elif isinstance(arg, Sum):
                other_part.append(_cg_simp_sum(arg))
                continue
            elif isinstance(arg, Mul):
                terms = []
                for term in arg.args:
                    if isinstance(term, Pow) and sympify(term.exp).is_number:
                        [ terms.append(term.base) for _ in range(term.exp) ]
                    elif isinstance(term, Sum):
                        terms.append(_cg_simp_sum(term))
                    else:
                        terms.append(term)
            elif isinstance(arg, Pow):
                terms = []
                if sympify(arg.exp).is_number:
                    [ terms.append(arg.base) for i in range(arg.exp) ]
                else:
                    other_part.append(arg)
                    continue
            else:
                other_part.append(arg)
                continue
            cg_part.append(terms)
        else:
            other_part.append(arg)
    i = 0
    while i < len(cg_part):
        cg, coeff, sign = _cg_list(cg_part[i])
        cg_count = len(cg)
        if cg_count == 1:
            if not cg[0].j1.is_number:
                i += 1
                continue
            # Simplifies Varshalovich 8.7.1 Eq 1
            # Sum( CG(a,alpha,b,0,a,alpha), (alpha,-a,a) ) = (2*a+1) * delta_b,0
            cg_index = [None]*(2*cg[0].j1+1)
            for term in cg_part:
                cg2, coeff2, sign2 = _cg_list(term)
                if not len(cg2) == 1:
                    continue
                cg2, coeff2, sign2 = _has_cg(_check_871_1, cg, sign, cg2, coeff2, sign2)
                if not len(cg2) == 1:
                    continue
                if not cg2[0].m1.is_number or not cg2[0].j1.is_number:
                    continue
                cg_index[cg2[0].m1+cg2[0].j1] = term,cg2[0],coeff2
            if cg_index.count(None) == 0:
                min_coeff = min([abs(term[2]) for term in cg_index])
                for term in cg_index:
                    cg_part.pop(cg_part.index(term[0]))
                    if not term[2] == min_coeff*sign:
                        cg_part.append((term[1],term[2]-min_coeff*sign))
                other_part.append((2*cg[0].j1+1)*min_coeff*sign*KroneckerDelta(cg[0].j2,0))
                continue
            # Simplifies Varshalovich 8.7.1 Eq 2
            # Sum( (-1)**(a-alpha) * CG(a,alpha,a,-alpha,c,0), (alpha,-a,a) ) = sqrt(2*a+1) * delta_c,0
            cg_index = [None]*(2*cg[0].j1+1)
            for term in cg_part:
                cg2, coeff2, sign2 = _cg_list(term)
                if not len(cg2) == 1:
                    continue
                cg2, coeff2, sign2 = _has_cg(_check_871_2, cg, sign, cg2, coeff2, sign2)
                if not len(cg2) == 1:
                    continue
                if not cg2[0].m1.is_number or not cg2[0].j1.is_number:
                    continue
                cg_index[cg2[0].m1+cg2[0].j1] = term,cg2[0],coeff2
            if cg_index.count(None) == 0:
                min_coeff = min([abs(term[2]) for term in cg_index])
                for term in cg_index:
                    sign2 = (-1)**((term[1].j1-term[1].m1)+(cg[0].j1-cg[0].m1))
                    cg_part.pop(cg_part.index(term[0]))
                    if not term[2] == min_coeff*sign*sign2:
                        cg_part.append((term[1],term[2]-min_coeff*sign*sign2))
                other_part.append(sqrt(2*cg[0].j1+1)*min_coeff*sign*sign2*KroneckerDelta(cg[0].j3,0))
                continue
        if cg_count == 2:
            if not (cg[0].j1.is_number and cg[0].j2.is_number):
                i += 1
                continue
            # Simplifies Varshalovich 8.7.2 Eq 9
            # Sum( CG(a,alpha,b,beta,c,gamma) * CG(a,alpha',b,beta',c,gamma), (c, abs(a-b), a+b), (gamma, -c, c) ) =
            #    delta_alpha,alpha' * delta_beta,beta'
            if cg[0].m1.is_number and cg[0].m2.is_number:
                cg_index = [None]*(cg[0].j1+cg[0].j2-max(abs(cg[0].j1-cg[0].j2),abs(cg[0].m1+cg[0].m2))+1)
            else:
                # TODO: Symbolic simplification for this case
                #cg_index = [None]*(cg[0].j1+cg[0].j2-abs(cg[0].j1-cg[0].j2)+1)
                continue
            for term in cg_part:
                cg2, coeff2, sign2 = _cg_list(term)
                if not len(cg2) == 2:
                    continue
                cg2, coeff2, sign2 = _has_cg(_check_872_9, cg, sign, cg2, coeff2, sign2)
                if not len(cg2) == 2:
                    continue
                if cg[0].m1.is_number and cg[0].m2.is_number:
                    cg_index[cg2[0].j3-max(abs(cg[0].j1-cg[0].j2),abs(cg[0].m1+cg[0].m2))] = term,cg2[0]*cg2[1],coeff2
                else:
                    #cg_index[cg2[0].j3-abs(cg[0].j1-cg[0].j2)] = term,cg2[0]*cg2[1],coeff2
                    continue
            if cg_index.count(None) == 0:
                min_coeff = min([abs(term[2]) for term in cg_index])
                for term in cg_index:
                    cg_part.pop(cg_part.index(term[0]))
                    if not term[2] == min_coeff*sign:
                        cg_part.append((term[1],term[2]-min_coeff*sign))
                other_part.append(min_coeff*sign*KroneckerDelta(cg[0].m1,cg[1].m1)*KroneckerDelta(cg[0].m2,cg[1].m2))
                continue
        i += 1
    return Add(*[Mul(*i) for i in cg_part])+Add(*other_part)

#TODO: Implement symbolic simplification of Sum objects
def _cg_simp_sum(e):
    e = _check_varsh_sum_871_1(e)
    e = _check_varsh_sum_871_2(e)
    e = _check_varsh_sum_872_4(e)
    return e

def _check_varsh_sum_871_1(e):
    a = Wild('a')
    alpha = symbols('alpha')
    b = Wild('b')
    match = e.match(Sum(CG(a,alpha,b,0,a,alpha),(alpha,-a,a)))
    if not match is None and len(match) == 2:
        return (2*match.get(a)+1)*KroneckerDelta(match.get(b),0)
    return e

def _check_varsh_sum_871_2(e):
    a = Wild('a')
    alpha = symbols('alpha')
    c = Wild('c')
    match = e.match(Sum((-1)**(a-alpha)*CG(a,alpha,a,-alpha,c,0),(alpha,-a,a)))
    if not match is None and len(match) == 2:
        return sqrt(2*match.get(a)+1)*KroneckerDelta(match.get(c),0)
    return e

def _check_varsh_sum_872_4(e):
    a = Wild('a')
    alpha = Wild('alpha')
    b = Wild('b')
    beta = Wild('beta')
    c = Wild('c')
    cp = Wild('cp')
    gamma = Wild('gamma')
    gammap = Wild('gammap')
    match1 = e.match(Sum(CG(a,alpha,b,beta,c,gamma)*CG(a,alpha,b,beta,cp,gammap),(alpha,-a,a),(beta,-b,b)))
    if not match1 is None and len(match1) == 8:
        return KroneckerDelta(match1.get(c),match1.get(cp))*KroneckerDelta(match1.get(gamma),match1.get(gammap))
    match2 = e.match(Sum(CG(a,alpha,b,beta,c,gamma)**2,(alpha,-a,a),(beta,-b,b)))
    if not match2 is None and len(match2) == 6:
        return 1
    return e

def _check_871_1(cg1, cg2, sign1, sign2):
    cg1 = cg1[0]
    cg2 = cg2[0]
    if cg1.j1 == cg2.j1 and cg1.j2 == cg2.j2 and cg1.j3 == cg2.j3 and \
        cg2.m2 == 0 and cg2.m1 == cg2.m3 and cg2.j1 == cg2.j3 and \
        sign1 == sign2:
            return True
    return False

def _check_871_2(cg1, cg2, sign1, sign2):
    cg1 = cg1[0]
    cg2 = cg2[0]
    if cg1.j1 == cg2.j1 and cg1.j2 == cg2.j2 and cg1.j3 == cg2.j3 and \
        cg2.j1 == cg2.j2 and cg2.m3 == 0 and cg2.m1 == -cg2.m2 and \
        sign1*(-1)**(cg1.j1-cg1.m1) == sign2*(-1)**(cg2.j1-cg2.m1):
            return True
    return False

def _check_872_9(cg1, cg2, sign1, sign2):
    cg1_1 = cg1[0]
    cg1_2 = cg1[1]
    cg2_1 = cg2[0]
    cg2_2 = cg2[1]
    if cg1_1.j1 == cg2_1.j1 and cg1_1.j2 == cg2_1.j2 and cg1_1.m1 == cg2_1.m1 and cg1_1.m2 == cg2_1.m2 and \
            cg1_2.j1 == cg2_2.j1 and cg1_2.j2 == cg2_2.j2 and cg1_2.m1 == cg2_2.m1 and cg1_2.m2 == cg2_2.m2 and \
            cg2_1.j1 == cg2_2.j1 and cg2_1.j2 == cg2_2.j2 and cg2_1.j3 == cg2_2.j3 and cg2_1.m3 == cg2_2.m3 and \
            sign1 == sign2:
        return True
    return False

def _cg_list(arg_list):
    coeff = 1
    cg = []
    for term in arg_list:
        if isinstance(term, CG):
            cg.append(term)
        elif isinstance(term, Pow):
            terms = []
            if isinstance(term.base, CG) and sympify(term.exp).is_number:
                [ cg.append(term.base) for i in range(term.exp) ]
        else:
            coeff *= term
    return cg, coeff, coeff/abs(coeff)


def _has_cg(f, cg1, sign1, cg2, coeff2, sign2):
    # This should be extended to check if symmetries can be applied to give the necessary cg terms
    if f(cg1, cg2, sign1, sign2):
        return cg2, coeff2, sign2
    return [], [], []
