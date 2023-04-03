"""
#------------------------------------------------------------------------#
#                                                                        #
#                         z_transform of f[n]                            #
#                         z_transform_inverse of F[z]                    #
#                                                                        #
#------------------------------------------------------------------------#
# http://blog.espol.edu.ec/telg1001/transformada-z-pares-fn-y-fz-con-sympy-python/
# Edison Del Rosario edelros@espol.edu.ec
# 2023/03/31
"""
### import for z_transform
from sympy.core import diff, S
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.numbers import Rational
from sympy.core.power import Pow
from sympy.core.symbol import Wild, Symbol
from sympy.core.function import ( expand, WildFunction)
from sympy.simplify import powsimp
from sympy.simplify.simplify import simplify
from sympy.polys.polytools import pdiv,real_roots,degree,LC,factor,Poly
from sympy.polys.partfrac import apart
from sympy.functions.elementary.trigonometric import cos, sin, acos, atan
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.exponential import log
from sympy.functions.combinatorial.factorials import factorial
from sympy.concrete.summations import summation, Sum
from sympy.functions.special.delta_functions import DiracDelta, Heaviside
from sympy import Float, preorder_traversal
from sympy.utilities.misc import debug

##########################################################################
# Helpers / Utilities
##########################################################################

class ZTransformError(NotImplementedError):
    """
    Exception raised in relation to problems computing transforms.

    Explanation
    ===========

    This class is mostly used internally; if integrals cannot be computed
    objects representing unevaluated transforms are usually returned.

    The hint ``needeval=True`` can be used to disable returning transform
    objects, and instead raise this exception if an z_transform cannot be
    computed.
    """
    def __init__(self, transform, function, msg):
        super().__init__(
            "%s could not be computed: %s." % (transform, msg))
        self.function = function

#------------------------------------------------------------------------#
#                         z_transform of f[n]                            #
#------------------------------------------------------------------------#

def _z_pairs_table(n, z):
    """
    z_transform basic pairs table of discrete time 'n_domain' and
    'z_domain', complemented with convergence conditions.

    Explanation
    ===========
    Table as a List with simplified arguments formatted as:
    z_transform_table = [(time_domain,
                          z_domain,condition,
                          convergence plane,
                          preparation function),]
    Only basic table is used in order to test and verify
    function z_properties to obtain other z_transforms like:
    n*Heaviside(n), a**n*Heaviside(n). Heaviside(n-1), etc.

    Examples
    ========
    >>> from sympy.discrete.z_transforms import _z_pairs_table
    >>> from sympy.abc import n,z
    >>> table = _z_pairs_table(n, z)
    >>> len(table)
    5
    >>> print(table[0][0:4])
    (DiracDelta(n), 1, True, 0)
    >>> print(table[1][0:4])
    (Heaviside(n), z/(z - 1), True, Abs(z) > 1)
    >>> print(table[3][0:4])
    (sin(n*a_), z*sin(a_)/(z**2 - 2*z*cos(a_) + 1), Abs(a_) > 0, Abs(z) > 1)
    >>>

    basic pairs table:
    (DiracDelta(n), 1, True, 0)
    (DiracDelta(n*a_), 1, Abs(a_) > 0, 0)
    (Heaviside(n), z/(z - 1), True, Abs(z) > 1)
    (cos(n*a_), z*(z - cos(a_))/(z**2 - 2*z*cos(a_) + 1), Abs(a_) > 0, Abs(z) > 1)
    (sin(n*a_), z*sin(a_)/(z**2 - 2*z*cos(a_) + 1), Abs(a_) > 0, Abs(z) > 1)

    Parameters
    ==========
    n : variable on discrete n_domain
    z : variable on z_domain

    Reference
    =========
    http://blog.espol.edu.ec/telg1001/transformada-z-tabla/
    """
    a = Wild('a', exclude=[n,z])
    b = Wild('b', exclude=[n,z])
    c = Wild('c', exclude=[n,z])
    dco = lambda f: _z_args_collect(f,n)
    # ( time domain,
    #   z domain,
    #   condition, convergence plane, preparation function )
    z_transform_table = [
    # DiracDelta
    (DiracDelta(n),
     S.One,
     S.true, S.Zero, dco),
    # Heaviside
    (Heaviside(n),
     z/(z-1),
     S.true, Abs(z) > 1, dco),
    # cos, sin
    (cos(a*n),
     (z*(z-cos(a)))/(z**2-(2*cos(a))*z+1),
     Abs(a)>0, Abs(z) > 1, dco),
    (sin(a*n),
     (sin(a)*z+0)/(z**2-(2*cos(a))*z+1),
     Abs(a)>0, Abs(z) > 1, dco),
    (cos(b*n+c),
     z*(z*cos(c)-cos(b-c))/(z**2-(2*cos(b))*z+1),
     Abs(b)>0, Abs(z) > 1, dco),
    ]

    def _z_args_collect(f, n):
        """
        Collects 'f(n)' arguments from tree to match 'f(a*n+b)'.
        'f(w*n-1*n-c)' is simplified as 'f((w-1)*n-c)'
        """
        f_dco = f
        if len(f.args) > 0:
            fn = f.func
            arg_n = list(f.args)
            arg_n = [_z_args_collect(arg_k, n) for arg_k in arg_n]
            f_dco = fn(*arg_n)
            if fn.is_Add:
                f_dco = f_dco.collect(n)
        return f_dco
    return z_transform_table

def z_pairs_properties(f,n,z, apply_properties=True):
    """
    z_transform based on pairs table and properties,
    otherwise F[z] pass to try with summation.

    Parameters
    ==========
    f : f[n] expression at discrete n_domain
    n : variable on discrete n_domain
    z : variable on z_domain
    apply_properties:
       True: check with properties if not found at table
       False: check only with pairs table
    debug_on : show operations on f[n]

    Examples
    ========
    >>> from sympy.discrete.z_transforms import z_pairs_properties
    >>> from sympy.abc import n,z
    >>> import sympy as sym
    >>> f = sym.Heaviside(n)
    >>> z_pairs_properties(f,n,z)
    (z/(z - 1), Abs(z) > 1, True)
    >>> f = n*sym.Heaviside(n)
    >>> z_pairs_properties(f,n,z)
    (z/(z - 1)**2, Abs(z) > 1, True)
    >>>
    >>> f = sym.DiracDelta(n)
    >>> z_pairs_properties(f,n,z)
    (1, 0, True)
    >>>
    >>> f = sym.cos(3*n)*sym.Heaviside(n)
    >>> z_pairs_properties(f,n,z)
    (z*(z - cos(3))/(z**2 - 2*z*cos(3) + 1), Abs(z) > 1, True)
    >>>
    >>> f = sym.sin(2*n)*sym.Heaviside(n)
    >>> z_pairs_properties(f,n,z)
    (z*sin(2)/(z**2 - 2*z*cos(2) + 1), Abs(z) > 1, True)
    """

    # Try _z_pairs_table
    Fz = None

    # extract constant factors from Fz
    k, fn = f.as_independent(n, as_Add=False)
    z_pairs = _z_pairs_table(n, z) # get table
    z_pairs_len = len(z_pairs)
    i = 0

    # lookup on table
    while (i<z_pairs_len) and (Fz is None):
        n_dom, z_dom, check, plane, prep = z_pairs[i]
        ma = prep(fn).match(n_dom)
        if ma or ma=={}: # expression equivalent found
            debug('\n_z_pairs_table match:')
            debug('  k:',k,' ; fn:',fn)
            debug('  z_pair f[n]:',n_dom)
            debug('  z_pair F[z]:',z_dom)
            try:
                c = check.xreplace(ma)
                debug('  try,check: %s -> %s'%(check, c))
                if c:
                    Fz = k*z_dom.xreplace(ma)
                    plane_z = plane.xreplace(ma)
                    cond_z  = S.true
                    Fz = (Fz, plane_z, cond_z)
            except Exception:
                debug(' _z_pairs_table did not match.')
        i = i+1 # next pair

    # f[n] did not match z_pairs, try z_properties()
    if Fz is None and apply_properties is True and f.has(n):
        Fz = z_properties(f, n, z)

    return Fz

def z_properties(f,n,z):
    """
    z_transform properties check:
      Multiplication by n    , n*f[n] <--> z*diff(Fz)
      Multiplication by a**n , a**n*f[n] <--> F[z/a]
      Time shifting, f[n-b] <--> z**(-b)*F +Fsum
      Time reversal, f[-n] <--> F[1/z]
    returns None if not found by pairs and properties
    Parameters
    ==========
    f : f[n] expression at discrete n_domain
    n : variable on discrete n_domain
    z : variable on z_domain
    debug_on : show operations on f[n]

    Examples
    ========
    >>> from sympy.discrete.z_transforms import z_properties
    >>> from sympy.abc import n,z
    >>> import sympy as sym
    >>> f = n*sym.Heaviside(n)
    >>> z_properties(f,n,z)
    (z/(z - 1)**2, Abs(z) > 1, True)
    >>>
    >>> f = n*n*sym.Heaviside(n)
    >>> z_properties(f,n,z)
    (z*(z + 1)/(z - 1)**3, Abs(z) > 1, True)
    >>>
    >>> f = sym.Heaviside(n-2)
    >>> z_properties(f,n,z)
    (1/(z*(z - 1)), Abs(z) > 1, True)
    >>>
    >>> f = ((sym.S.One/2)**2)*sym.Heaviside(n-2)
    >>> z_properties(f,n,z)
    (1/(4*z*(z - 1)), Abs(z) > 1, True)
    >>>
    >>> f = ((1/2)**2)*sym.Heaviside(n-2)
    >>> z_properties(f,n,z)
    (0.25/(z*(z - 1)), Abs(z) > 1, True)
    >>>
    >>> f = sym.Heaviside(-n)
    >>> z_properties(f,n,z)
    (-1/(z - 1), Abs(z) < 1, True)
    >>>
    >>> f = -sym.Heaviside(-n-1)
    >>> z_properties(f,n,z)
    (z/(z - 1), Abs(z) < 1, True)
    >>>
    """
    if not n in f.free_symbols:
        msg = 'n variable is not in f[n] expression'
        raise ZTransformError('z_properties',f,msg)
    Fz = None
    f = expand(powsimp(f))
    a = Wild('a', exclude=[n,z])
    b = Wild('b', exclude=[n,z])
    y = Wild('y')
    g = WildFunction('g', nargs=1)
    k, fn = f.as_independent(n, as_Add=False)
    debug('\nz_transform properties_check -------------')
    debug(' k:',k,' ; fn:',fn)
    if fn.is_Add:
        msg = 'f[n] is an Additive expression, suggested to use z_transform(fn,n,z) instead'
        raise ZTransformError('z_properties',fn,msg)
    # check fn por pow(n,a) or pow(a,n)
    f_powna = S.One
    f_powan = S.One
    factor_Mul = Mul.make_args(fn)
    for factor_k in factor_Mul:
        ma_powna = factor_k.match(Pow(n,a))
        ma_powan = factor_k.match(Pow(a,n))
        ma_powan_ = factor_k.match(Pow(a,-n))
        # debug(factor_k,ma_powna,ma_powan_)
        if (ma_powna and not(ma_powna[a]==S.Zero) and
            ma_powna[a].is_integer):
            f_powna = f_powna*factor_k
            debug('  f_powna ; n**a:',factor_k,' ; ',ma_powna)
        if ma_powan or ma_powan_:
            if ma_powan:
                f_powan = f_powan*factor_k
            if ma_powan_:
                f_powan = (1/f_powan)*factor_k
            debug('  f_powan ; a**n:',factor_k,' ; ',ma_powan,ma_powan_)
    ma_un = fn.match(Heaviside(y))
    ma_gn = fn.match(g)
    ma_gu = fn.match(g*Heaviside(y))

    # z_property a**n*f[n] <--> F[z/a]
    if not f_powan==S.One : # a**n factor match
        ma_powan  = f_powan.match(Pow(a,n))
        ma_powan_ = f_powan.match(Pow(a,-n))
        fn_na = fn/f_powan
        debug(' _z_property a**n*f[n] ----------')
        debug('   a**n ; fn/f_pow :',f_powan,' ; ',fn_na)
        Fz = z_pairs_properties(fn_na,n,z)
        if not Fz is None:
            an_value =1
            if ma_powan:
                an_value = ma_powan[a]
            if ma_powan_:
                an_value = Rational(S.One,ma_powan_[a])
            #debug('Fz',Fz,ma_powan,ma_powan_,an_value)
            Fz_k = k*factor(Fz[0].subs(z,z/an_value))
            cond_z = list(Fz[1].args)
            num_position = None
            fun_position = None
            for i in range(len(cond_z)):
                cond_z[i] = cond_z[i].subs(z,z*an_value)
                if not cond_z[i].has(z):
                    num_position = i
                if cond_z[i].has(z):
                    fun_position = i
                    k_Rel, f_Rel = cond_z[i].as_independent(z, as_Add=False)
            if num_position is not None:
                cond_z[fun_position] = cond_z[fun_position]/k_Rel
                cond_z[num_position] = cond_z[num_position]/k_Rel
            cond_z = Fz[1].func(cond_z[0],cond_z[1])
            Fz = (Fz_k,cond_z,Fz[2])

            debug('\n  _z_property a**n*f[n] <--> F[z/a]:\n  ',Fz)

    # z_property n*f[n] <--> -z*dF(z)/dz --------
    elif not f_powna==S.One : # n**a factor match
        ma_powna = f_powna.match(Pow(n,a))
        fn_n = fn/n
        debug('\n _z_property n*f[n]')
        debug('  n**a:',n,' ; fn/n, :',fn_n)
        Fz = z_pairs_properties(fn_n,n,z)
        if not Fz is None:
            Fz = (k*(-z)*factor(diff(Fz[0],z,1)),
                  Fz[1], Fz[2])
            debug('\n  _z_property n*f[n] <--> -z*diff(Fz):\n  ',Fz)
    # match Heaviside(n)
    elif ma_un:
        arg_n = ma_un[y].collect(n)
        ma_u = arg_n.match(a*n-b)
        debug(' ma_un ; ma_u: ',ma_un,' ; ',ma_u)
        # time shift only
        if ma_u[a]>0 and Abs(ma_u[b])>=0:
            fn_g = Heaviside(n)
            debug(' _z_property time_shift u[n-b]--------')
            debug('  fn_gn:', fn_g)
            b_shift = ma_u[b]
            Fz = z_pairs_properties(fn_g, n, z)
            if not Fz is None:
                Fz = (k*factor(z**(-b_shift)*Fz[0]),
                      Fz[1], Fz[2])
            msg = ' _z_property time_shift'
            debug( msg,'f[n-b] <--> z**(-b)*F +Fsum:\n ', Fz)
        # u time reversal and shift
        if ma_u[a]<0 and Abs(ma_u[b])>=0:
            fn_g = Heaviside(n)
            b_shift = -ma_u[b]
            msg ='  _z_property time_reversal and shift'
            debug(msg,' u[-(n+b)]---')
            debug(' fn_g ; b_shift:',fn_g,' ; ',b_shift)
            Fz = z_pairs_properties(fn_g,n, z)
            if not Fz is None:
                plane_reversal = simplify(Fz[1].func(Fz[1].args[1],Fz[1].args[0]))
                Fz = (k*factor((Fz[0].subs(z,1/z))*z**(-b_shift)),
                      plane_reversal, Fz[2])
            debug(msg,' u[-(n+b)]:\n ',Fz)

    # match g function
    elif ma_gn:
        arg_n = ma_gn[g].args[0].collect(n)
        ma_g = arg_n.match(a*n-b)
        debug(' ma_gn ; ma_g : ',ma_gn,ma_g)
        # g time shift only
        if not ma_g[a]==S.Zero and not Abs(ma_g[b])==S.Zero:
            b_shift = ma_g[b]
            fn_g = ma_gn[g].subs(arg_n,ma_g[a]*n)
            debug(' _z_property time_shift --------')
            debug('  fn_gn ; b_shift:', fn_g,' ; ',b_shift)
            Fz = z_pairs_properties(fn_g, n, z)
            Fsum = 0
            m = b_shift # Sumation term x[n-m]u[n]
            # DiracDelta does not apply
            # right shift n
            if not ma_g[b]==S.Zero and not fn_g.has(DiracDelta):
                fn_gz = fn_g.subs(arg_n,-n)*z**n
                Fsum = summation(fn_gz,(n,1,m))
            # left shift n
            if ma_g[b]<S.Zero and not fn_g.has(DiracDelta):
                fn_gz = fn_g.subs(arg_n,n)*z**(-n)
                Fsum = -summation(fn_gz,(n,0,-m-1))
            Fsum = expand(z**(-m)*Fsum)
            # apply shift
            if not Fz is None:
                Fz = (Fsum + k*factor((z**(-b_shift))*Fz[0]),
                      Fz[1], Fz[2])
                msg = ' _z_property time_shift'
                debug(msg,' z**(-b)*F[z]+Fsum:\n ', Fz)
        # g time reversal and shift
        if ma_g[a]<0 and Abs(ma_g[b])>0:
            b_shift = -ma_g[b]
            fn_g = ma_gn[g].subs(arg_n,ma_g[a]*n)
            debug(' _z_property time reversal and shift ---')
            debug('  fn_gn ; b_shift:', fn_g,' ; ',b_shift)
            Fz = z_pairs_properties(fn_g, n, z)
            Fsum = 0
            m = b_shift # Sumation term x[n-m]u[n]
            # DiracDelta does not apply
            # right shift n
            if ma_g[b]>0 and not fn_g.has(DiracDelta):
                fn_gz = fn_g.subs(arg_n,-n)*z**n
                Fsum = summation(fn_gz,(n,1,m))
            # left shift n
            if ma_g[b]<0 and not fn_g.has(DiracDelta):
                fn_gz = fn_g.subs(arg_n,n)*z**(-n)
                Fsum = -summation(fn_gz,(n,0,-m-1))
            Fsum = expand(z**(-m)*Fsum)
            if not Fz is None: # apply shift
                Fz = (Fsum+k*factor((z**(-b_shift))*(Fz[0].subs(z,1/z))),
                      Fz[1], Fz[2])
                msg = ' _z_property time_shift'
                debug(msg,' z**(-b)*F[z]+Fsum:\n ', Fz)

    # match g*u
    if ma_gu: # match g*u
        debug(' _z_property g[n]*u[n] ---------')
        debug('  ma_gu: ',ma_gu)
        # Heaviside factor
        u_arg = ma_gu[y].collect(n)
        u_k, u_arg = (factor(u_arg)).as_independent(n, as_Add=False)
        p_invert = False  # plain invertion
        if u_arg == (-n-1): # u = -u(-n-1)
            k = - k
            u_arg = n
            p_invert = True
        ma_u = u_arg.match(a*n-b)
        # g function factor
        g_arg = ma_gu[g].args[0].collect(n)
        g_k, g_arg = (factor(g_arg)).as_independent(n, as_Add=False)
        ma_g  = g_arg.match(a*n-b)
        # check args g(a1*n-b1)*u(a2*n-b2)
        PQ = pdiv(g_arg,u_arg)
        debug('  g_k*(g_arg): ',g_k,'*(',g_arg,')')
        debug('  u_k*(u_arg): ',u_k,'*(',u_arg,')')
        debug('  g_arg/u_arg poly division:',
              g_k,'*(',PQ[0],'*(',u_arg,')+',PQ[1],')')

        if ma_u[a]==S.One and ma_u[b]==S.Zero: # g*u[n]
            fn_g = ma_gu[g]
            debug('  fn_g:',fn_g,)
            Fz = z_pairs_properties(fn_g,n,z)
            if not Fz is None:
                plane_z =  Fz[1]
                if p_invert:
                    plane_z = simplify(Fz[1].func(Fz[1].args[1],
                                                  Fz[1].args[0]))
                Fz = (k*Fz[0], plane_z, Fz[2])
                debug('\n  _z_property g*u[n] Fz:\n ',Fz)

        else: # g[q(n-k)+p]*u[n-k]
            g_arg_m = PQ[0]*n+PQ[1]
            fn_g  = ma_gu[g].subs(ma_gu[g].args[0],g_k*g_arg_m)
            debug('  fn_g:',fn_g,)
            Fz = z_pairs_properties(fn_g,n,z)
            if not Fz is None:
                Fz = (k*Fz[0], Fz[1], Fz[2])
                debug(' _z_property g[q(n-k)+p]*u[n-k] Fz:\n ',Fz)

    # clean floats as integers, 4z/(4z-1) to z/(z-1/4)
    if not Fz is None:
        Fz0 = _simplify_z(Fz[0],z)
        Fz0 = _round_float_is_int(Fz0)
        Fz1 = _round_float_is_int(Fz[1])
        Fz2 = Fz[2]
        Fz = (Fz0,Fz1,Fz2)
    debug(' F[z]:',Fz)
    return Fz

def _z_transform_summation(f, n, z):
    """
    z_transform unilateral summation,for a simple term (auxiliar function)
    """
    f = f+0*n
    Fz = None
    if f.has(n):
        fzn = f*(z**(-n))
        # unilateral summation
        Fz_sum = summation(fzn,(n,0,S.Infinity))
        debug(' _z_transform_summation: ',Fz)
        if Fz_sum.is_Piecewise:
            Fz_eq = Fz_sum.args[0]  # first interval
            Fz = (Fz_eq[0].simplify(),Fz_eq[1]) # expression found
        else:
            Fz = Fz_sum
        debug(' _z_transform_summation: ',Fz)
    return Fz

def z_transform(f,n,z):
    """
    z_transform of f[n] to F[z], using _z_pairs_table() and
    z_properties(), otherwise still not found, then apply _z_summation.

    Explanation
    ===========
    f[n] is expanded as Add terms to simplify the process.
    Returns F[z] expression is all done, otherwise return
    a tuple(F[z], f[n] not processed)

    Examples
    ========
    >>> from sympy.discrete.z_transforms import z_transform
    >>> from sympy.abc import n,z
    >>> import sympy as sym
    >>> f = sym.DiracDelta(n)
    >>> z_transform(f,n,z)
    (1, 0, True)
    >>>
    >>> f = sym.Heaviside(n-1)
    >>> z_transform(f,n,z)
    (1/(z - 1), Abs(z) > 1, True)
    >>>
    >>> f = sym.cos(2*n)
    >>> z_transform(f,n,z)
    (z*(z - cos(2))/(z**2 - 2*z*cos(2) + 1), Abs(z) > 1, True)
    >>>
    >>> f = (n**2)*sym.Heaviside(n)
    >>> z_transform(f,n,z)
    (z*(z + 1)/(z - 1)**3, Abs(z) > 1, True)
    >>>
    >>> f = ((S.One/2)**n)*sym.Heaviside(n)
    >>> z_transform(f,n,z)
    (z/(z - 1/2), Abs(z) > 2, True)
    >>>
    >>> f = ((1/2)**n)*sym.Heaviside(n)
    >>> z_transform(f,n,z)
    (z/(z - 0.5), Abs(z) > 2, True)
    >>> f = ((S.One/2)**n)*n*sym.Heaviside(n)
    >>> z_transform(f,n,z)
    (z/(2*(z - 1/2)**2), Abs(z) > 2, True)
    >>>
    >>> f = DiracDelta(n) + cos(4*n)
    >>> z_transform(f,n,z)
    z*(z - cos(4))/(z**2 - 2*z*cos(4) + 1) + 1
    >>>

    Parameters
    ==========
    f : f[n] expression at discrete n_domain
    n : variable on discrete n_domain
    z : variable on z_domain
    """
    f = expand(powsimp(f,n))
    FT = 0*n
    F_noeval = []
    term_sum = Add.make_args(f)
    debug('\nterm_sum:',len(term_sum))
    # to_do express properties with a single Term, no Add
    if len(term_sum)==1:
        FT = z_pairs_properties(f, n, z) #[0]
        if FT is None: # not done by _z_pairs_table or z_properties
            FT = _z_transform_summation(f, n, z)
            debug('  FT_summation:',FT)
    else:
        for term_k in term_sum:
            k, fn = term_k.as_independent(n, as_Add=False)
            debug('  k:',k,' ; term_k:',fn)
            if not fn.has(n):
                msg = 'is a constant term, consider to make it causal by *u[n]'
                raise ZTransformError('z_transform',term_k,msg)
            term_z = z_pairs_properties(fn, n, z)
            if term_z is None: # not done by _z_pairs_table or z_properties
                term_z = _z_transform_summation(fn, n, z)
                debug('  k',k,' ; term_z_summation:',term_z)
                if term_z is None or isinstance(term_z,Sum):
                    F_noeval.append(term_k)
            if not term_z is None :
                debug('\n  k:',k,' ; term_z:',term_z,'\n')
                if not isinstance(term_z,Sum):
                    term_z = _round_float_is_int(k*term_z[0])
                FT = FT + term_z
##        if len(F_noeval)>0:
##            FT = (FT,F_noeval)
    debug('\n F[z]:',FT)

    return FT

#------------------------------------------------------------------------#
#                         z_transform_inverse of F[z]                    #
#------------------------------------------------------------------------#
def z_pairs_prop_inverse(F,z,n,apply_properties=True):
    """
    z_transform based on pairs table and properties,
    otherwise returns None .

    Explanation
    ===========
    z_transform pairs table
    http://blog.espol.edu.ec/telg1001/transformada-z-tabla/

    Examples
    ========
    >>> from sympy.discrete.z_transforms import z_pairs_prop_inverse
    >>> from sympy.abc import z,n
    >>> F = z/(z - 1)
    >>> z_pairs_prop_inverse(F,z,n)
    (Heaviside(n), Abs(z) > 1, True)
    >>> F = z/(z - 1)**2
    >>> z_pairs_prop_inverse(F,z,n)
    (n*Heaviside(n), Abs(z) > 1, True)
    >>>

    Parameters
    ==========
    F : F[z] expression at z_domain
    z : variable on z_domain
    n : variable on discrete n_domain
    apply_properties:
       True: check with properties if not found at table
       False: check only with pairs table
    """
    fn = None
    # check expression F[z] = P/Q
    F = simplify(F+0*z) # F+0*z if F=constant
    k, Fz = F.as_independent(z, as_Add=False)
    a = Wild('a', exclude=[z])
    b = Wild('b', exclude=[z])
    c = Wild('c', exclude=[z])
    r = Wild('r', exclude=[z])

    # Fz poles and zeros
    PQ = _Fz_as_PQ_roots(Fz,z)
    P = PQ['P']
    Q = PQ['Q']
    P_degree = PQ['P_degree']
    Q_poles  = PQ['Q_poles']
    Q_degree = PQ['Q_degree']

    # Q denominator check
    ma_Qza = Q.match(z**a)
    ma_Q2  = None # Q degree 2
    if (fn is None and not ma_Qza) and len(Q_poles)<Q_degree:
        #debug('--ma_Q2 block')
        ma_P1 = None
        ma_P2 = None
        if P_degree==1:
            ma_P1 = P.match(a*z+ b)
        if P_degree==2:
            ma_P2 = simplify(P).match(r*z*(a*z+b))
        ma_Q2 = Q.match(a*z**2+ b*z + r**2)
        k_r = None
        if ma_Q2 and not(ma_Q2[r]==S.Zero or ma_Q2[a]==S.Zero):
            k_r = ma_Q2[r]
            k_r = Rational(k_r).limit_denominator(100)
            k_w = acos(ma_Q2[b]/(-2*k_r)).evalf()
            k_w = Rational(k_w).limit_denominator(100)
            if ma_P2: # cos[a*n] ; P.match(r*z*(a*z+b))
                debug('  -- ma_P2 block at z_pairs_prop_inverse')
                k_Q2 = ma_Q2[b]/ma_P2[b]
                k_Q2 = Rational(k_Q2.evalf()).limit_denominator(100)
                k_P2 = ma_P2[b]/cos(k_w)#/ma_P2[a]/cos(k_w)
                k_P2 = Abs(Rational(k_P2.evalf()).limit_denominator(100))
                k = (k*ma_P2[a]*ma_P2[r]).evalf()
                k = Abs(Rational(k).limit_denominator(100))
                a_Q = Rational(ma_P2[a]).limit_denominator(100)
                if ma_P2[r]==1 and ma_P2[a]==1 and k_Q2==2:
                    F_numer = z*(a_Q*z-k_P2*cos(k_w))
                    F_denom = z**2-k_Q2*k_r*cos(k_w)*z+k_r**2
                    Fz = F_numer/F_denom
                    debug('  -- Fz: ',Fz)
            if ma_P1 and ma_P1[b]==S.Zero: # sin[a*n] ; P.match(a*z+ b)
                #debug(' -- ma_P1 block')
                k_Q2 = ma_Q2[b]/Abs(cos(k_w).evalf())
                k_Q2 = Rational(k_Q2.evalf()).limit_denominator(100)
                k_P1 = (ma_P1[a]/sin(k_w)).evalf()
                k_P1 = Abs(Rational(k_P1).limit_denominator(100))
                if not(ma_P1[a]==S.Zero) and ma_P1[b]==S.Zero:
                    k = k/sin(k_w)
                    F_numer = z*sin(k_w)
                    F_denom = z**2-k_Q2*k_r*cos(k_w)*z+k_r**2
                    Fz = F_numer/F_denom
    msg = 'pair'
    # Try _z_pairs_table
    z_pairs = _z_pairs_table(n, z)
    z_pairs_len = len(z_pairs)
    i=0
    while i<z_pairs_len and fn is None:
        n_dom, z_dom, check, plane, prep = z_pairs[i]
        ma_z = Fz.match(z_dom)
        #debug('z_pair F[z]:',z_dom, ma_z,Fz)
        if ma_z or ma_z=={}:
            debug('\n_z_pairs_table match:')
            debug('  k:',k,' ; Fz:',Fz)
            debug('  z_pair F[z]:',z_dom,' ; ma_z:',ma_z)
            debug('  z_pair f[n]:',n_dom)
            try:
                c = check.xreplace(ma_z)
                debug('  try,check  : %s -> %s'%(check, c))
                if c or F==1:
                    fn = n_dom.xreplace(ma_z)
                    if not fn.has(Heaviside) and not fn.has(DiracDelta):
                        fn = fn*Heaviside(n)
                    plane_z = plane.xreplace(ma_z)
                    cond_z  = S.true
                    fn = (k*fn, plane_z, cond_z)
            except Exception:
                debug(' _z_pairs_table did not match.')
        i = i+1 # next pair

    # f[n] did not match z_pairs, try z_properties_inverse()
    if (fn is None and apply_properties is True and F.has(z)):
        fn = z_properties_inverse(F, z, n)
        msg = 'prop'

    # clean floats as integers
    if not fn is None:
        fn0 = _round_float_is_int(fn[0])
        fn1 = _round_float_is_int(fn[1])
        fn2 = _round_float_is_int(fn[2])
        fn = (fn0,fn1,fn2)
        debug('   '+msg+' f[n]:',fn)
    return fn

def _Fz_as_PQ_roots(Fz,z):
    """
    auxiliary function to analyze F[z] as P/Q with numerical
    and symbolic parts
    """
    a = Wild('a', exclude=[z])
    # Fz poles and zeros
    PQ = _simplify_z(Fz,z,PQ_args=True)
    P = PQ['P']
    Q = PQ['Q']
    k = PQ['k']
    Fz = PQ['F']
    P_sign = 1
    if k<0:
        P_sign = -1
        k = Abs(k)
    # P numerator review
    P_zeros = []
    P_degree = None # P is not z**a symbol
    ma_za  = P.match(z**a)
    if P.has(z): # numerator with pole
        if ma_za: # P_zeros with 0 only or z**(symbol)
            if len(P.free_symbols)==1:
                P_zeros = [0]*ma_za[a]
            if len(P.free_symbols)==2: # z**(symbol)
                P_zeros = [0] # at least a Zero
                P_degree = ma_za[a]
        elif len(P.free_symbols)==1 and not ma_za:
            P_zeros  = real_roots(P.evalf())
        else:
            raise ZTransformError('_Fz_as_PQ_roots',Fz,
                                  'unable to determine P_zeros from expression F[z] = P/Q')
    if P_degree is None:
        P_degree = degree(P,z)

    # Q denominator review
    Q_poles = []
    Q_degree = None # Q is not z**(symbol)
    ma_za  = Q.match(z**a)
    if Q.has(z): # denominator with pole
        if ma_za: # 0 poles only or z**(symbol)
            if len(Q.free_symbols)==1:
                Q_poles = [0]*ma_za[a]
            elif len(Q.free_symbols)==2: #z**(symbol)
                Q_poles  = [0] # at least a Zero
                Q_degree = ma_za[a]
        elif len(Q.free_symbols)==1:
            Q_poles  = real_roots(Q.evalf())
        elif len(Q.free_symbols)>1 and not ma_za:
            Q_poles  = real_roots(Poly(Q,z))
        else:
            raise ZTransformError('_Fz_as_PQ_roots',Fz,
                                  'unable to determine Q_poles from expression F[z] = P/Q')
    # unique poles
    Q_unique = []
    for pole in Q_poles:
        if pole not in Q_unique:
            Q_unique.append(pole)
    if Q_degree is None:
        Q_degree = degree(Q,z)
    Fz = P/Q
    PQ = {'P':P,
          'P_degree': P_degree, 'P_zeros':P_zeros,
          'Q':Q,
          'Q_degree': Q_degree, 'Q_poles':Q_poles,
          'Q_unique': Q_unique,
          'k':k, 'P_sign':P_sign, 'F':Fz
          }
    return PQ

def _split_f_kgu(fn,n):
    """
    auxiliary function split f[n]:
    constant k , funcion g[n], Heaviside u[n]
    """
    a = Wild('a', exclude=[n])
    y = Wild('y')
    g = WildFunction('g', nargs=1)
    ma_un = fn.match(a*Heaviside(y))
    ma_gn = fn.match(a*g)
    ma_gu = fn.match(a*g*Heaviside(y))
    f_kgu = (1,1,1)
    if ma_un:
        f_kgu = (ma_un[a],1,ma_un[y])
    if ma_gn:
        f_kgu = (ma_gn[a],ma_gn[g],1)
    if ma_gu:
        f_kgu = (ma_gu[a],ma_gu[g],ma_gu[y])
    return f_kgu

def z_properties_inverse(F,z,n):
    """
    z_transform inverse properties check:

    Explanation
    ===========
    z_transform table properties
    http://blog.espol.edu.ec/telg1001/transformada-z-tabla-de-propiedades/

    Multiplication by n    , n*f[n] <--> z*diff(Fz)
    Multiplication by a**n , a**n*f[n] <--> F[z/a]
    Time shifting , f[n-b] <--> z**(-b)*F +Fsum
    Time reversal , f[-n] <--> F[1/z]
    case, r*(|gamma|**n)*cos(b*n+c)
         <--> z*(A*z+B)/(z**2-2*a*z+gamma**2)
    case, n*(n-1)*(n-2)...(n-m+1)*(gamma**n)*u[n]/(m!*gamma**m)')
         <--> z/(z-gamma)**m+1')

    Examples
    ========
    >>> from sympy.discrete.z_transforms import z_properties_inverse
    >>> from sympy.abc import z,n
    >>> F = z/(z - 1)**2
    >>> z_properties_inverse(F,z,n)
    (n*Heaviside(n), Abs(z) > 1, True)
    >>>
    >>> F = z/(z - S.Half)
    >>> z_properties_inverse(F,z,n)
    (Heaviside(n)/2**n, Abs(z) > 1/2, True)
    >>>
    >>> z_properties_inverse(1/z**2,z,n)
    (DiracDelta(n - 2), 0, True)
    >>> z_properties_inverse(z**2,z,n)
    (DiracDelta(n + 2), 0, True)
    >>>
    >>> F = z*2/(z-S.One/2)**2
    >>> z_properties_inverse(F,z,n)
    (4*n*Heaviside(n)/2**n, Abs(z) > 1/2, True)
    >>>
    >>> from sympy.functions.elementary.trigonometric import cos
    >>> F = z*(0.968912421710645*z + 0.178246055649492)/(z**2 - 2*z*cos(2) + 1)
    >>> z_properties_inverse(F,z,n)
    (cos(2*n + 0.25)*Heaviside(n), Abs(z) > 1, True)
    >>>

    Parameters
    ==========
    F : F[z] expression at z_domain
    z : variable on z_domain
    n : variable on discrete n_domain
    """
    if not z in F.free_symbols:
        msg = 'z variable is not in F[z] expression'
        raise ZTransformError('z_properties_inverse',F,msg)
    fn = None
    F = simplify(F+0*z) # F+0*z if F=constant
    debug('\nz_properties_inverse_z_transform')
    a = Wild('a', exclude=[n,z])
    b = Wild('b', exclude=[n,z])
    c = Wild('c', exclude=[n,z])
    r = Wild('r', exclude=[n,z])

    # Fz poles and zeros
    PQ = _Fz_as_PQ_roots(F,z)
    P = PQ['P']
    Q = PQ['Q']
    P_zeros  = PQ['P_zeros']
    P_degree = PQ['P_degree']
    Q_poles  = PQ['Q_poles']
    Q_degree = PQ['Q_degree']
    #Q_unique = PQ['Q_unique']
    P_sign = PQ['P_sign']
    k = PQ['k']
    Fz = PQ['F']

    # P numerator check
    ma_Pza = P.match(z**a)
    ma_P1  = P.match(a*z+ b)
    # Q denominator check
    ma_Qza = Q.match(z**a)
    ma_Q1  = Q.match((c*z-a)**b)
    if ma_Q1 and ma_Q1[b]==0:
        ma_Q1 = None # no z in Q denominator
    ma_Q2 = None # Q degree 2
    # Imag poles
    if (fn is None and not ma_Qza) and len(Q_poles)<Q_degree:
        ma_P1 = None
        ma_P2 = None
        if P_degree==1: # sin[a*n]
            ma_P1 = P.match(a*z+ b)
        if P_degree==2: # cos[a*n]
            ma_P2 = simplify(P).match(r*z*(a*z+b))
        ma_Q2 = Q.match(a*z**2+ b*z + r**2)
        #debug('    ma_Q2:',ma_Q2)
        k_a = None
        k_r = None
        k_Q2 = None
        # Q.match(a*z**2+ b*z + r**2)
        if (ma_Q2 and not(ma_Q2[r]==S.Zero or ma_Q2[a]==S.Zero)):
            k_r = ma_Q2[r]/ma_Q2[a]
            k_r = Rational(k_r).limit_denominator(100)
            k_w = acos(ma_Q2[b]/ma_Q2[a]/(-2*k_r))
            k_w = Rational(k_w.evalf()).limit_denominator(100)
            if ma_P2: # cos[a*n]
                #debug('  -- ma_P2 block cos[a*n]--')
                k_Q2 = ma_Q2[b]/ma_Q2[a]/k_r
                if ma_P2[r]==1 and ma_P2[a]==1 and k_Q2==2:
                    F_numer = z*(z-k_r*cos(k_w))
                    F_denom = z**2-k_Q2*k_r*cos(k_w)*z+k_r**2
                    Fz = F_numer/F_denom
                    Fz = k*Fz
            if ma_P1 and ma_P1[b]==S.Zero: # sin[a*n] ; P.match(a*z+ b)
                #debug(' -- ma_P1 block sin[a*n]--')
                # ma_Q2 = Q.match(a*z**2+ b*z + r**2)
                k_Q2 = ma_Q2[b]/ma_Q2[a]/Abs(cos(k_w))/k_r
                if not(ma_P1[a]==S.Zero) and ma_P1[b]==0:
                    #k = (k/sin(k_w)).evalf()
                    F_numer = z*sin(k_w)
                    F_denom = z**2-k_Q2*k_r*cos(k_w)*z+k_r**2
                    Fz = F_numer/F_denom
                    Fz = k*Fz

    #  ma_Q1 = Q.match(c*(z-a)**b)
    if ma_Q2 and (ma_P2 or ma_P1):
        if ma_P2:
            ma_P2 = P.match(r*z*(a*z+b))
        if ma_P1:
            ma_P1 = P.match(a*z+ b)
    # debug partial results
    debug(' P:',P,'\n Q:',Q)
    debug(' k:',k,' ; P_signo:',P_sign,
          ' ; P_degree:',P_degree,' ; P_zeros:',P_zeros)
    debug(' Q_degree:',Q_degree,' ; Q_poles (real):',Q_poles)
    debug(' ma_Q1 (c*z-a)**b :',ma_Q1)
    if ma_Q2:
        debug(' ma_P1 (a*z+ b)    :',ma_P1)
        debug(' ma_P2 r*z*(b*z+a) :',ma_P2)
        debug(' ma_Q2 a*z**2+ b*z + r**2:\n ',ma_Q2)
        debug('  k_r',k_r,' ; k_w:',k_w,' ; k_Q2 :',k_Q2)
    debug(' F[z]: ',Fz)

    cond1 = not(ma_Q1 and ma_Q1[b]==S.Zero)
    # avoid z/(z-a)**3
    cond2 = not(ma_Q1 and ma_Q1[a]>0 and ma_Q1[b]>2 and
            ma_P1 and ma_P1[a]>S.Zero and ma_P1[b]==S.Zero)
    # avoid z/(z-1)**3
    cond3 = not(ma_Q1 and ma_Q1[a]>0 and ma_Q1[b]>2 and
            ma_P1 and ma_P1[a]>S.Zero and ma_P1[b]==S.Zero)

    # _z_property nf[n] <--> -z*diff(F[z])
    if (fn is None and ma_Q1 and ma_Q1[a]==1
        and P_degree>=1 and P_degree<=Q_degree and
        cond1 and cond2):
        FunI = factor(Fz/(-z),z)
        from sympy.integrals.integrals import integrate
        Fz = factor(integrate(FunI,z))
        F0 = 0 # Integral constant
        C = -Fz.subs(z,0)+F0
        FunC = factor(Fz+C)
        debug('\n _z_property multiply nf[n] <--> -z*diff(F[z])')
        debug(' Fz = integrate(factor(-Fz/z),z):\n\t=',Fz,C)
        debug(' Fz = integrate(factor(Fz)/z,z):\n\t=',FunC)
        if not FunC.has(log):
            fn = z_pairs_prop_inverse(FunC, z, n)
        else:
            fn = None
        if not fn is None:
            fn = (P_sign*k*n*fn[0],fn[1],fn[2])
            debug('\n _z_property multiply nf[n]:\n  ',fn)

    #  _z_property (a**n)*f[n] <--> F(z/a)
    # ma_Q1 = Q.match((z-a)**b) # ma_Q2 = Q.match(a*z**2+ b*z + r**2)
    elif ((ma_Q1 and not ma_Q1[a]==S.One  and not ma_Q1[a]==S.Zero) or
          (ma_Q2 and not(ma_Q2[r]==S.One) and
            ((ma_P2 and ma_P2[a]==S.One) or ma_P1)) and cond2):
        debug('\n _z_property multiply (a**n)*f[n] <--> F(z/a) ')
        k_a = None
        # Q1 is not a/z
        if (ma_Q1 and not ma_Q1[a]==S.One and not ma_Q1[a]==S.Zero and
            not ma_Q1[c]==S.Zero):
            k_a = ma_Q1[a]
            Fz = factor(Fz.subs(z,k_a*z))
        # quadratic denominator
        if ma_Q2 and not ma_Q2[r]==S.One:
            k_a = ma_Q2[r]/ma_Q2[a]
            Fz = factor(Fz.subs(z,k_a*z))
            [P2,Q2] = Fz.as_numer_denom()
            if ma_P2 and not ma_P2[a]==S.Zero:
                P2 = P2/k_a
            if ma_P1 and ma_P1[b]==0:
                P2 = P2/k_a
                k  = 1
            Q2 = Q2/k_a
            Fz = P2/Q2
        debug(' P_sign:',P_sign,'; k:',k,' ; k_a:',k_a,type(k_a))
        debug(' Fz:',Fz)
        fn = z_pairs_prop_inverse(Fz, z, n)
        if not fn is None:
            fnk = P_sign*k*(k_a**n)*fn[0]
            cond_z = list(fn[1].args)
            num_position = None
            fun_position = None
            for i in range(len(cond_z)):
                cond_z[i] = cond_z[i].subs(z,z/k_a)
                if not cond_z[i].has(z):
                    num_position = i
                if cond_z[i].has(z):
                    fun_position = i
                    k_Rel,f_rel = cond_z[i].as_independent(z, as_Add=False)
            if num_position is not None:
                cond_z[fun_position] = cond_z[fun_position]/k_Rel
                cond_z[num_position] = cond_z[num_position]/k_Rel
            cond_z = fn[1].func(cond_z[0],cond_z[1])
            fn = (fnk,cond_z,fn[2])
            debug('\n  _z_property multiply (a**n)*f[n]:\n ',fn)

    # _z_property time_shift f[n-b] <--> z**(-b)*F +Fsum
    elif ((0 in Q_poles) or
          (P_sign>0 and not(0 in Q_poles) and not P.has(z)) or
          (P_degree>Q_degree) or
          (ma_Pza and not ma_Pza[a]==S.Zero) and cond3):
        debug('\n _z_property time_shift f[n-b] <--> z**(-b)*F +Fsum ')
        debug(' Q_roots:',Q_poles)
        b_shift = 0
        if (not P.has(z) and ma_Q1 and not(ma_Q1[a]==S.Zero) and
            not ma_Q1[b]==S.Zero):
            b_shift = ma_Q1[a]
        if (0 in Q_poles) and not ma_Qza:
            b_shift = Q_poles[0]+1
        # only z**-a with a positive
        if (0 in Q_poles) and ma_Qza:
            debug(' only z**-a with a positive')
            b_shift = Q_degree
        # only z**a with a positive
        if (ma_Qza and ma_Qza[a]==S.Zero and ma_Pza):
            debug(' only z**a with a positive')
            b_shift = -1*P_degree
        # z**(1-a)
        if (0 in P_zeros) and not P_degree.is_Number:
            debug(' P_degree has symbol')
            P_symbol = list(P_degree.free_symbols)
            P_numero = P_degree.subs(P_symbol[0],0)
            b_shift = -(P_degree-P_numero)
        if not ma_Qza and P_degree.is_Number and P_degree>Q_degree:
            debug(' P_degree>Q_degree')
            b_shift = -1

        Fz = simplify(Fz*z**(b_shift))
        debug(' b_shift, Q_roots:',b_shift,' ; ', Q_poles)
        debug(' Fz:',Fz)
        debug(' k:',k)

        fn = z_pairs_prop_inverse(Fz, z, n)
        if not fn is None:
            f_kgu = _split_f_kgu(fn[0],n)
            f_g=1
            f_u=1
            f_k = P_sign*k*f_kgu[0]
            if f_kgu[1]!=1:
                f_g = f_kgu[1].subs(f_kgu[1].args[0],
                                    f_kgu[1].args[0]-b_shift)
            if f_kgu[2]!=1:
                f_u = Heaviside(f_kgu[2]-b_shift)
            fn = (f_k*f_g*f_u,fn[1],fn[2])
        debug('\n _z_property time_shift f[n-b] :\n ',fn)

    # _z_property time_reversal f[-n] <--> F[1/z]
    elif P_sign<0 and not P.has(z):
        Fz = simplify(Fz.subs(z,1/z))
        debug('\n _z_property time_reversal F[1/z]-- ')
        debug(' Fz : ', Fz)
        fn = z_pairs_prop_inverse(Fz, z, n)
        debug('  -----fn : ', fn,P_sign,k,fn[0])
        if not fn is None:
            f_kgu = _split_f_kgu(fn[0],n)
            f_g=1
            f_u=1
            f_k = f_kgu[0]*k*P_sign
            if f_kgu[1]!=1:
                f_g = f_kgu[1].subs(f_kgu[1].args[0],-f_kgu[1].args[0])
            if f_kgu[2]!=1:
                f_u = Heaviside(-f_kgu[2])
            plane_reversal = fn[1].func(fn[1].args[1],fn[1].args[0])
            plane_reversal = simplify(plane_reversal)
            fn = (f_k*f_g*f_u,plane_reversal,fn[2])
            debug('  _z_property time_reversal f[n]:\n ',fn)

    # case: r*(|gamma|**n)*cos(beta*n+theta)
    #          <--> z*(A*z+B)/(z**2-2*a*z+gamma**2)
    elif fn is None and ma_Q2 and ma_P2: # denominador cuadratico
        # ma_P2 : r*z*(a*z+b)  ; ma_Q2 : a*z**2+ b*z + r**2
        A = ma_P2[a]*ma_P2[r]
        B = ma_P2[b]*ma_P2[r]
        aq= ma_Q2[b]/ma_Q2[a]/2
        R = sqrt((ma_Q2[r]**2)/ma_Q2[a])
        R = Rational(R).limit_denominator(1000)
        r_P = ((A**2)*(R**2)+(B**2)-2*A*aq*B).evalf()
        r_Q = ((R**2)-(aq**2)).evalf()
        k_r = sqrt(r_P/r_Q)
        beta = acos(-aq/Abs(R)).evalf()
        theta_Q = A*sqrt(r_Q)
        theta = atan((A*aq-B)/theta_Q).evalf()
        fn = (P_sign*k*k_r*(Abs(R)**n*cos(beta*n+theta)*Heaviside(n)),
              Abs(z)>R,True)
        debug('\n  case : r*(|gamma|**n)*cos(beta*n+theta)')
        debug('          <--> z*(A*z+B)/(z**2-2*a*z+gamma**2)')
        debug('   F[z] : z(Az+B)/(z**2+2az+|r|**2)')
        debug('   k    :',k,' ; ','k_r:',k_r)
        debug('   gamma:',R)
        debug('   beta :',beta)
        debug('   theta:',theta)
        debug('   case f[n]:',fn)

    #case: n*(n-1)*(n-2)...(n-m+1)*(gamma**n)*u[n]/(m!*gamma**m)')
    #      <--> z/(z-gamma)**m+1')
    # z/(z-a)**m+1 # ma_Q1 = Q.match(c*(z-a)**b)
    elif (fn is None and ma_Q1 and not(ma_Q1[a]==S.Zero)
          and ma_Q1[b]>2 and P_degree==S.One):
        m = ma_Q1[b]-1
        ma_P1=None
        debug('\n case: n*(n-1)*(n-2)...(n-m+1)*(gamma**n)*u[n]/(m!*gamma**m)')
        debug('       <--> z/(z-gamma)**m+1')
        debug('   P:',P,' ; Q:',Q,
              ' ; m:',m,' ; gamma:',ma_Q1[a])
        if P_degree==1:
            ma_P1 = P.match(a*z+ b)
        if ma_P1 and ma_P1[b]==0:
            k_a = ma_P1[a]
            term = 1+0*z
            for i in range(0,m,1):
                term = term*(n-i)
            debug('   factors:',term)
            fn = term*(ma_Q1[a]**n)/((ma_Q1[a]**m)*factorial(m))
            fn = (P_sign*k*fn*Heaviside(n),Abs(z) > 1,0)
            debug('   case f[n]:',fn)
    if not fn is None:
        fn0 = _round_float_is_int(fn[0])
        fn1 = _round_float_is_int(fn[1])
        fn2 = fn[2]
        fn = (fn0,fn1,fn2)
        debug('   f[n]:',fn)
    return fn

def inverse_z_transform(F,z,n):
    """
    inverse_z_transform of F[z] to f[n], using _z_pairs_table()
    and z_properties_inverse().

    Explanation
    ===========
    F[z] is expanded with modified partial fractions as Add terms
    to simplify the process.
    Returns f[n] expression with for each Add term
    otherwise returns tuple(f[n], F[z] terms not processed)
    "Many of the transforms X[z] of practical interest
    are rational functions (ratio of polynomials in z), which
    can be expressed as a sum of partial fractions."[1]
    Inverse transforms can be found in a table of transform
    and applying some properties as done by:
    z_pairs_prop_inverse(F, z, n, apply_properties=True )
    "The partial fraction method works because for every
    transformable x[n] defined for n >= 0,
    there is a corresponding unique X[z] defined for |z| > r0
    (where r0 is some constant),and vice versa."[1]

    Parameters
    ==========
    F : F[z] expression at z_domain
    z : variable on z_domain
    n : variable on discrete n_domain
    debug_on : show operations on F[z]

    Examples
    ========
    >>> from sympy.discrete.z_transforms import inverse_z_transform
    >>> from sympy.abc import z,n
    >>> F = 5*z/(z-1)
    >>> inverse_z_transform(F,z,n)
    5*Heaviside(n)
    >>>
    >>> F = (8*z-19)/((z-2)*(z-3))
    >>> inverse_z_transform(F,z,n)
    (9*2**n + 10*3**n)*Heaviside(n)/6 - 19*DiracDelta(n)/6
    >>>

    See Also: apart_z(F)

    References
    ==========
    .. [1] [Lathi]5.1-1 p495 Lathi B.P and Green R.A.(2018).
           Linear Systems and Signals Third Edition.
           Oxford University Press.
           ISBN-13: 978-0190200176, ISBN-10: 0190200170
    """
    debug('\ninverse_z_transform(F,z,n):',F)
    F = factor(F,z)
    F = apart_z(F)

    fT = 0*n
    f_noeval = []
    term_sum = Add.make_args(F)
    if len(term_sum)>1:
        debug(' apart_z(F) as Add:',term_sum)
    for term_k in term_sum:
        if len(term_sum)>1:
            debug('  term_k[z]:',term_k)
        term_n = None # not paired yet
        if not term_k==S.Zero: # avoid term_k=0
            term_n = z_pairs_prop_inverse(term_k, z, n)
        if term_n is None and not term_k==S.Zero:
            f_noeval.append(term_k)
            term_n = None
        if not term_n is None:
            fT = fT + _round_float_is_int(term_n[0])
        if len(term_sum)>1:
            debug('  term_k[n]:',term_n,)
    debug('\n f[n]:',fT)
    if len(term_sum)>1:
        fT = fT.collect(Heaviside(n))
        fT = fT.collect(DiracDelta(n))
        term_sum = Add.make_args(fT)
        fT = 0
        for term_k in term_sum:
            term_k = simplify(term_k)
            fT = fT+term_k
        debug(' f[n]:',fT)
    if len(f_noeval)>0:
        fT = (fT,f_noeval)
    return fT

def apart_z(F):
    """
    modified partial fractions expansion at z_domain.
    Apply partial fractions to F[z]/z and later restore F[z]*z
    to get expressions with u[n] rather than u[n-1]

    Parameters
    ==========
    F : F[z] expression at z_domain
    debug_on : show operations on F[z]

    Explanation
    ===========
    F[z] expanded into partial fractions directly, always obtain
    an answer is multiplied by u[n-1]. It is prefered to have
    expressions with u[n] rather than u[n-1].
    This goal can be accomplished by expanding x[z]/z into
    partial fractions and then multiplying both sides by z. [1]

    Examples
    ========
    >>> from sympy.discrete.z_transforms import apart_z
    >>> from sympy.abc import z
    >>> F = (8*z-19)/((z-2)*(z-3))
    >>> apart_z(F)
    3*z/(2*(z - 2)) + 5*z/(3*(z - 3)) - 19/6
    >>>
    >>> from sympy.polys.partfrac import apart
    >>> apart(F)
    3/(z - 2) + 5/(z - 3)
    >>>
    >>> F = -z/(z-1)
    >>> apart(F)
    -1 - 1/(z - 1)
    >>> apart_z(F)
    -z/(z - 1)
    >>>

    References
    ==========
    .. [1] [Lathi]5.3a p495 Lathi B.P and Green R.A.(2018).
           Linear Systems and Signals Third Edition.
           Oxford University Press.
           ISBN-13: 978-0190200176, ISBN-10: 0190200170
    """
    z = Symbol('z')
    F = F + 0*z # check if only a constant
    if not z in F.free_symbols:
        return F

    debug('\napart_z(F)     :',F)
    a = Wild('a', exclude=[z])
    b = Wild('b', exclude=[z])
    c = Wild('c', exclude=[z])

    # separate symbols z**a , a
    F = factor(F)
    k_a = 1 # factor constant and symbol
    z_a = 1 # z**a , a is symbol
    z_b = 1 # z**a , a is integer
    factor_Mul = Mul.make_args(F)
    for factor_k in factor_Mul:
        ma_za = factor_k.match(z**a)
        ma_ka = factor_k.match(a)
        debug('  z**a, ma_za:',ma_za,' ; a ,ma_ka:',ma_ka)
        if ma_za and not ma_za[a]==S.Zero and not ma_za[a].is_integer is True:
            za_symbol = list(ma_za[a].free_symbols)
            debug('  za_symbol:',za_symbol)
            for a_symbol in za_symbol: # z**(a+b+c-1)
                za_integer = ma_za[a].subs(a_symbol,0)
            za_symbol = ma_za[a] - za_integer
            z_a = z**za_symbol #z**(a+b+c)
            F = F/z_a
            debug('  z_a:',z_a) # z**(-1) only integer pow
            if Abs(za_integer)>0:
                z_b = z_b**(za_integer)
                F = F/z_b
        # symbolic factors, as constant
        if ma_ka and not ma_ka[a].is_integer:
            k_a = k_a*ma_ka[a]
            F = F/ma_ka[a]
            debug('  k_a:',k_a)

    # apply partial fractions to F[z]/z
    Fzz = factor(F/z,z)

    # Fz poles and zeros
    PQ = _Fz_as_PQ_roots(Fzz,z)
    P = PQ['P']
    Q = PQ['Q']
    P_zeros  = PQ['P_zeros']
    P_degree = PQ['P_degree']
    Q_degree = PQ['Q_degree']
    Q_unique = PQ['Q_unique']

    debug('  simplify(F/z):',Fzz)
    debug('  P:',P,'\n  P_degree:',P_degree,' ; P_zeros :',P_zeros)
    debug('  Q:',Q,'\n  Q_degree:',Q_degree,' ; Q_unique:',Q_unique)

    # avoid: z**3/(z-1) ; z**2/(z-a) ; 1/(z-1)**2 ; z**2/(z-a)**2
    if (P.has(z) and (len(Q_unique)>1 or P_degree==(Q_degree-1))):
        Fzz = apart(Fzz,z)
    debug('  apart(F/z):',Fzz)

    #  restore F[z]*z
    term_suma = Add.make_args(Fzz)
    Fzp = 0*z
    for term_k in term_suma:
        if term_k.has(z): # not a constant
            # make Q lead coefficient=1
            Pk,Qk = term_k.as_numer_denom()
            ma_Qk = Qk.match((a*z+b)**c)
            ma_Qk2 = Qk.match(a*z**2+b*z+c)
            debug('  ma_Qk:',ma_Qk,' ; ma_Qk2:',ma_Qk2)
            # Q_degree = 1, not constant
            if (ma_Qk and not ma_Qk[c]==S.Zero and
                not ma_Qk[a]==S.One and
                not ma_Qk[a]==S.Zero):
                Pk = Pk/ma_Qk[a].doit()
                Qk = Qk/ma_Qk[a].doit()
                term_k = Pk/Qk
            # Q_degree = 2
            elif (ma_Qk2 and not(ma_Qk2[a]==S.One) and
                  not(ma_Qk2[a]==S.Zero) and
                  not(ma_Qk and ma_Qk[c]==S.Zero)):
                Pk = Pk/ma_Qk2[a].doit()
                Qk = Qk/ma_Qk2[a].doit()
                term_k = Pk/Qk
        #  restore F[z]*z and symbolic factor z**a or a
        term_k = _simplify_z(term_k,z)
        Fzp = Fzp + term_k*z*z_a*k_a*z_b

    Fzp = _round_float_is_int(Fzp)
    debug('\n apart z*(Fz/z)  :',Fzp)
    return Fzp

def _round_float_is_int(F, nearzero=1e-10):
    """
    auxiliary function to round floats in F expresion with no decimal part,
    ending with .0 to simplify reading expression
    """
    for a in preorder_traversal(F):
        if isinstance(a, Float):
            r = a%1
            if r<nearzero:
                F = F.subs(a, int(a))
    return F

def _round_float(F, precision=15):
    """
    auxiliary function to round floats F[Z] in expression to simplify
    reading with a given precision or to match "assert" at test files.
    """
    for a in preorder_traversal(F):
        if isinstance(a, Float):
            F = F.subs(a, round(a,precision))
    return F

def _simplify_z(F,z,PQ_args = False):
    """
    auxiliary function to simplify z expressions
    4z/4z+3 is simplified as z/(z - 3/4)
    PQ_args = True, returns{P,Q,k,F}
    """
    a = Wild('a', exclude=[z])
    b = Wild('b', exclude=[z])
    c = Wild('c', exclude=[z])
    # review leadcoef 4z/4z+3 is simplified as z/(z - 3/4)
    P,Q = F.as_numer_denom()
    debug('\n_simplify_z:',F,'\n P:',P,'\n Q:',Q)

    # P numerator review
    count_vars = P.free_symbols
    ma = P.match((a*z-b)**c)
    if P.has(z) and len(count_vars)==1:
        k = LC(P)
    elif P.has(z) and len(count_vars)>1:
        k = 1
        for term_k in Mul.make_args(P):
            if not term_k.has(z):
                k = k*term_k
    else: # numerator is a constant
        k = P

    # repeated zeros
    ma = P.match((a*z-b)**c)
    cond = ma and len(ma[c].free_symbols)==0 \
        and ma[c]>S.One and not ma[c]==S.Zero \
        and not(ma[a]==S.One) and Abs(ma[a])>S.Zero
    if cond:
        P = (z-ma[b]/ma[a])**ma[c]
    else:
        P = (P/k).doit()
    P = _round_float_is_int(P)
    k = _round_float_is_int(k) #constant factor
    debug('  P match (a*z-b)**c:',ma,
          '\n  P_leadcoef:',k,' : P:',P)

    # Q denominator review
    count_vars = Q.free_symbols
    ma = Q.match((a*z-b)**c)
    if Q.has(z) and len(count_vars)==1:
        Q_leadcoef = LC(Q)
    elif Q.has(z) and len(count_vars)>1:
        Q_leadcoef = 1
        for term_k in Mul.make_args(Q):
            if not term_k.has(z):
                Q_leadcoef = Q_leadcoef*term_k
    else: # denominator is a constant
        Q_leadcoef = Q

    # repeated poles
    ma = Q.match((a*z-b)**c)
    cond = ma and len(ma[c].free_symbols)==0 \
        and ma[c]>S.One and not ma[c]==S.Zero and \
        not(ma[a]==S.One) and Abs(ma[a])>S.Zero
    if cond:
        Q = (z-ma[b]/ma[a])**ma[c]
    else:
        Q = (Q/Q_leadcoef).doit()
    Q = _round_float_is_int(Q)
    Q_leadcoef =_round_float_is_int(Q_leadcoef)
    k = _round_float_is_int(k/Q_leadcoef)   # constant factor
    debug('  Q match (a*z-b)**c:',ma,
          '\n  Q_leadcoef:',Q_leadcoef,'  Q:',Q)

    F = k*P/Q
    debug('  k:',k,' ; Fz:',P/Q,'\n  simplify_z F:',F)
    PQ = {'P':P,'Q':Q,
          'k':k,'F':P/Q}
    if PQ_args is True:
        F = PQ
    return F
