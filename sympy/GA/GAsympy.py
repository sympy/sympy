#GAsympy.py

import sys
from sympy import expand,Mul,Add,Symbol,S,Expr,Wild,Pow,diff,trigsimp,\
                  simplify,Matrix

try:
    from numpy import matrix
    #from numpy.linalg import det,eig
    numpy_loaded = True
except ImportError:
    numpy_loaded = False

ONE  = S(1)
ZERO = S(0)

def get_commutative_coef(expr):
    if isinstance(expr,Mul):
        (coefs,bases) = expr.args_cnc()
        return(Mul(*coefs))

def half_angle_reduce(expr,theta):
    (s,c) = symbols('s,c')
    sub_dict = {sin(theta/2):s,cos(theta/2):c}
    new_expr = expr.subs(sub_dict)
    sub_dict = {s*c:sin(theta)/2,s**2:(1-cos(theta))/2,c**2:(1+cos(theta))/2}
    #print new_expr
    new_expr = trigsimp(simplify(new_expr.subs(sub_dict)),recursive=True)
    #print expand(new_expr)
    return(new_expr)

def linear_expand(expr):
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
        bases = bases[0]
    elif isinstance(expr,Symbol): #term is Symbol
        coefs = ONE
        bases = expr
    elif isinstance(expr,Add): #expr has multiple terms
        coefs = []
        bases = []
        for arg in expr.args:
            term = arg.args_cnc()
            coef = Mul(*term[0])
            base = term[1][0]
            if base in bases: #increment coefficient of base
                ibase = bases.index(base)
                coefs[ibase] += coef
            else: #add base to list
                coefs.append(coef)
                bases.append(base)

    if not isinstance(coefs,list): #convert single coef to list
        coefs = [coefs]
    if not isinstance(bases,list): #convert single base to list
        bases = [bases]
    coefs = tuple(coefs)
    bases = tuple(bases)
    return(coefs,bases)

def linear_projection(expr,plist=None):
    """
    If a sympy 'Expr' is of the form:

    expr = expr_0+expr_1*a_1+...+expr_n*a_n

    where all the a_j are noncommuting symbols in basis then

    proj(expr) returns the sum of those terms where a_j is in plist
    """
    if expr.is_commutative and plist == None: #return scalar projection
        return(expr)
    expr = expand(expr)
    if isinstance(expr,Mul): #expr has single term
        x = (coefs,bases) = expr.args_cnc()
        if bases[0] in plist: #vector term to be projected
            return(Mul(*coefs)*bases[0])
        else:
            return(ZERO)
    elif isinstance(expr,Symbol): #base vector to be projected
        if expr in plist:
            return(expr)
        else:
            return(ZERO)
    elif isinstance(expr,Add): #expr has multiple terms
        result = ZERO
        for arg in expr.args:
            term = arg.args_cnc()
            if term[1] == [] and plist == None: #scalar term to be projected
                result += Mul(*term[0])
            elif term[1] != [] and plist != None and term[1][0] in plist: #vector term to be projected
                    result += Mul(*term[0])*term[1][0]
        return(result)

def non_scalar_projection(expr):
    """
    If a sympy 'Expr' is of the form:

    expr = expr_0*ONE+expr_1*a_1+...+expr_n*a_n

    where all the a_j are noncommuting symbols in basis then

    proj(expr) returns the sum of those terms where a_j is in plist
    """
    if expr.is_commutative: #return scalar projection
        return(ZERO)
    expr = expand(expr)
    if isinstance(expr,Mul): #expr has single term
        x = (coefs,bases) = expr.args_cnc()
        if bases[0] != MV.ONE: #vector term to be projected
            return(Mul(*coefs)*bases[0])
        else:
            return(ZERO)
    elif isinstance(expr,Symbol): #base vector to be projected
        if expr != MV.ONE:
            return(expr)
        else:
            return(ZERO)
    elif isinstance(expr,Add): #expr has multiple terms
        result = ZERO
        for arg in expr.args:
            term = arg.args_cnc()
            if term[1] != MV.ONE: #vector term to be projected
                    result += Mul(*term[0])*term[1][0]
        return(result)

def nc_substitue(expr,sub_dict):
    (coefs,bases) = linear_expand(expr)
    result = ZERO
    for (coef,base) in zip(coefs,bases):
        if base != 1:
            result += coef*sub_dict[base]
    return(result)

def linear_function(expr,fct):
    """
    If a sympy 'Expr' is of the form:

    expr = expr_0+expr_1*a_1+...+expr_n*a_n

    where all the a_j are noncommuting symbols in basis then

    f(expr) = expr_0+expr_1*f(a_1)+...+expr_n*f(a_n)

    is returned
    """
    if expr.is_commutative:
        return(expr)
    expr = expand(expr)
    if isinstance(expr,Mul):
        x = (coefs,bases) = expr.args_cnc()
        return(Mul(*coefs)*fct(bases[0]))
    elif isinstance(expr,Symbol):
        return(fct(expr))
    elif isinstance(expr,Add):
        result = ZERO
        for arg in expr.args:
            term = arg.args_cnc()
            if term[1] == []:
                result += Mul(*term[0])
            else:
                result += Mul(*term[0])*fct(term[1][0])
    return(result)

def coef_function(expr,fct):
    """
    If a sympy 'Expr' is of the form:

    expr = expr_0+expr_1*a_1+...+expr_n*a_n

    where all the a_j are noncommuting symbols in basis then

    f(expr) = fct(expr_0)+fct(expr_1)*a_1+...+fct(expr_n)*a_n

    is returned
    """
    expr = expand(expr)
    if isinstance(expr,Mul):
        x = (coefs,bases) = expr.args_cnc()
        return(fct(Mul(*coefs))*bases[0])
    elif isinstance(expr,Symbol):
        if expr.is_commutative:
            return(fct(expr))
        else:
            return(expr)
    elif isinstance(expr,Add):
        result = ZERO
        for arg in expr.args:
            term = arg.args_cnc()
            if term[1] == []:
                result += fct(Mul(*term[0]))
            else:
                result += fct(Mul(*term[0]))*fct(term[1][0])
    return(result)

def bilinear_product(expr,fct):
    """
    If a sympy 'Expr' is of the form:

    expr = expr_ij*a_i*a_j or expr_0 or expr_i*a_i

    where all the a_i are noncommuting symbols in basis and the expr's

    are commuting expressions then

    bilinear_product(expr) = expr_ij*fct(a_i,a_j)

    bilinear_product(expr_0) = expr_0

    bilinear_product(expr_i*a_i) = expr_i*a_i
    """

    def bilinear_term(expr,fct):
        if expr == ZERO:
            return(expr)
        if isinstance(expr,Mul): #bases in expr
            tmp = (coefs,bases) = expr.args_cnc()
            coef = Mul(*tuple(coefs))
            if isinstance(bases[0],Pow): #base is a_i**2
                args = bases[0].args
                return(coef*fct(args[0],args[0]))
            elif len(bases) == 1: #base is a_i
                return(expr)
            else: #base is a_i*a_j
                return(coef*fct(bases[0],bases[1]))
        elif isinstance(expr,Pow): # expr is a_i*a_i
            args = expr.args
            return(fct(args[0],args[0]))
        elif isinstance(expr,Symbol):
            return(expr)
        else:
            sys.stderr.write('!!!!Cannot compute bilinear_product for '+str(expr)+'!!!!\n')
            sys.exit(1)

    expr = expand(expand(expr))

    if not isinstance(expr,Add):
        return(bilinear_term(expr,fct))
    else:
        result = ZERO
        for term in expr.args:
            tmp = bilinear_term(term,fct)
            result += tmp
        return(result)

def multilinear_product(expr,fct):
    """
    If a sympy 'Expr' is of the form:

    expr = expr_i1i2...irj*a_i1*a_i2*...*a_ir or expr_0

    where all the a_i are noncommuting symbols in basis and the expr's

    are commuting expressions then

    multilinear_product(expr) = expr_i1i2...ir*fct(a_i1,a_i2,...,a_ir)

    bilinear_product(expr_0) = expr_0

    where fct() is defined for r <= n the total number of bases
    """
    if expr.is_commutative: #no bases in expr
        return(expr)
    if isinstance(expr,Mul): #bases in expr
        (coefs,bases) = expr.args_cnc()
        if len(coefs) == 0: #expr_ij = 1
            coefs = [ONE]
        coef = Mul(*tuple(coefs))
        new_bases = []
        for base in bases:
            if isinstance(base,Pow):
                args = base.args
                new_bases += args[1]*[args[0]]
            else:
                new_bases.append(base)
        return(coef*fct(new_bases))

def bilinear_function(expr,fct):
    """
    If a sympy 'Expr' is of the form:

    expr = expr_0+expr_1*a_1+...+expr_n*a_n+expr_11*a_1*a_1
           +...+expr_ij*a_i*a_j+...+expr_nn*a_n*a_n

    where all the a_j are noncommuting symbols in basis then

    bilinear_function(expr) = bilinear_product(expr_0)+bilinear_product(expr_1*a_1)+...+bilinear_product(expr_n*a_n)
                              +bilinear+product(expr_11*a_1*a_1)+...+bilinear_product(expr_nn*a_n*a_n)
    """
    if expr.is_commutative:
        return(expr)
    expr = expand(expr)
    if isinstance(expr,(Mul,Pow,Symbol)): #only one additive term
        return(bilinear_product(expr,fct))
    elif isinstance(expr,Add): #multiple additive terms
        result = ZERO
        for arg in expr.args:
            result += bilinear_product(arg,fct)
    return(result)

def multilinear_function(expr,fct):
    """
    If a sympy 'Expr' is of the form summation convention):

    expr = expr_0+Sum{0<r<=n}{expr_i1i2...ir*a_i1*a_i2*...*a_ir}

    where all the a_j are noncommuting symbols in basis then and the
    dimension of the basis in n then

    bilinear_function(expr) = multilinear_product(expr_0)
        +Sum{0<r<=n}multilinear_product(expr_i1i2...ir*a_i1*a_i2*...*a_ir)
    """
    if expr.is_commutative:
        return(expr)
    expr = expand(expr)
    if isinstance(expr,(Mul,Pow,Symbol)): #only one additive term
        return(bilinear_product(expr,fct))
    elif isinstance(expr,Add): #multiple additive terms
        result = ZERO
        for arg in expr.args:
            result += bilinear_product(arg,fct)
    return(result)

def linear_derivation(expr,fct,x):
    """
    If a sympy 'Expr' is of the form:

    expr = expr_0+expr_1*a_1+...+expr_n*a_n

    where all the a_j are noncommuting symbols in basis then

    linear_drivation(expr) = diff(expr_0,x)+diff(expr_1,x)*a_1+...
                            +diff(expr_n,x)*a_n+expr_1*fct(a_1,x)+...
                            +expr_n*fct(a_n,x)
    """
    if expr.is_commutative:
        return(diff(expr,x))
    expr = expand(expr)
    if isinstance(expr,Mul):
        x = (coefs,bases) = expr.args_cnc()
        coef = Mul(*coefs)
        return(diff(coef,x)*bases[0]+coef*fct(bases[0],x))
    elif isinstance(expr,Symbol):
        return(fct(expr,x))
    elif isinstance(expr,Add):
        result = ZERO
        for arg in expr.args:
            term = arg.args_cnc()
            coef = Mul(*term[0])
            if term[1] == []:
                result += diff(coef,x)
            else:
                result += diff(coef,x)*term[1][0]+coef*fct(term[1][0],x)
    return(result)

def product_derivation(F,fct,x):
    """
    If a sympy 'Expr' is of the form:

    expr = expr_0*a_1*...*a_n

    where all the a_j are noncommuting symbols in basis then

    product_derivation(expr) = diff(expr_0,x)*a_1*...*a_n
                              +expr_0*(fct(a_1,x)*a_2*...*a_n+...
                              +a_1*...*a_(i-1)*fct(a_i,x)*a_(i+1)*...*a_n+...
                              +a_1*...*a_(n-1)*fct(a_n,x))
    """
    if F.is_commutative:
        return(diff(F,x))
    elif isinstance(F,Mul):
        (coefs,bases) = F.args_cnc()
        coef = Mul(*coefs)
        dcoef = diff(coef,x)
        if len(bases) == 1:
            return(dcoef*bases[0]+coef*fct(bases[0],x))
        else:
            result = dcoef*Mul(*bases)
            for ib in range(len(bases)):
                result += coef*Mul(*bases[:ib])*fct(bases[ib],x)*Mul(*bases[ib+1:])
            return(result)
    elif isinstance(F,Symbol):
        return(fct(F,x))

def multilinear_derivation(F,fct,x):
    """
    If a sympy 'Expr' is of the form (summation convention):

    expr = expr_0 + expr_i1i2...ir*a_i1*...*a_ir

    where all the a_j are noncommuting symbols in basis then

    dexpr = diff(expr_0,x) + d(expr_i1i2...ir*a_i1*...*a_ir)

    is returned where d() is the product derivation
    """
    if F.is_commutative:
        return(diff(F,x))
    elif isinstance(F,Mul) or isinstance(F,Symbol):
        return(product_derivation(F,fct,x))
    elif isinstance(F,Add):
        result = ZERO
        for term in F.args:
            result += product_derivation(term,fct,x)
        return(result)

def numpy_matrix(M):
    if not numpy_loaded:
        print '!!!Cannot use "numpy_matrix" since "numpy" is not loaded!!!'
        sys.exit(1)
    Mlst = M.tolist()
    nrows = len(Mlst)
    ncols = len(Mlst[0])
    for irow in range(nrows):
        for icol in range(ncols):
            try:
                Mlst[irow][icol] = float(Mlst[irow][icol])
            except ValueError:
                print 'In Matrix:'
                print M
                print 'Cannot convert '+str(Mlst[irow][icol])+' to python float.'
                sys.exit(1)
    return(matrix(Mlst))
"""
M = Matrix([[0,1,2,3],\
            [1,0,4,5],\
            [2,4,0,6],\
            [3,5,6,0]])
print M
Mnp = numpy_matrix(M)
print Mnp
(e,v) = eig(Mnp)
print Mnp*v[:,3]-e[3]*v[:,3]
"""
