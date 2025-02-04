"""
#------------------------------------------------------------------------#
#                         z_transform of f[n]                            #
#                         z_transform_inverse of F[z]                    #
#------------------------------------------------------------------------#
# http://blog.espol.edu.ec/telg1001/transformada-z-pares-fn-y-fz-con-sympy-python/
# Edison Del Rosario edelros@espol.edu.ec
# 2023/04/27
"""
from sympy.core import diff, S
from sympy.core.numbers import Integer
from sympy.core.add import Add
from sympy.core.mul import Mul
#from sympy.core.numbers import Rational
from sympy.core.power import Pow
from sympy.core.symbol import Wild, Symbol
from sympy.core.function import expand, WildFunction
from sympy.simplify import powsimp
from sympy.simplify.simplify import simplify
from sympy.polys.polytools import pdiv,degree,LC,factor,Poly
from sympy.polys.partfrac import apart
from sympy.polys.polyroots import roots
from sympy.functions.elementary.trigonometric import cos, sin, acos, atan
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.exponential import log
from sympy.functions.combinatorial.factorials import factorial
from sympy.concrete.summations import summation, Sum
from sympy.functions.special.delta_functions import DiracDelta, Heaviside
from sympy.logic.boolalg import And
from sympy import Float, preorder_traversal, sympify,monic
from sympy.utilities.misc import debug

#------------------------------------------------------------------------#
# Helpers / Utilities
#------------------------------------------------------------------------#

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
        super().__init__("%s could not be computed: %s." % (transform, msg))
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
    >>> from sympy.integrals.ztransforms import _z_pairs_table
    >>> from sympy.abc import n,z
    >>> table = _z_pairs_table(n, z)
    >>> len(table)
    9
    >>> print(table[0][0:4])
    (DiracDelta(n), 1, True, True)
    >>> print(table[3][0:4])
    (Heaviside(n), z/(z - 1), True, Abs(z) > 1)
    >>> print(table[5][0:4])
    (sin(n*a_), z*sin(a_)/(z**2 - 2*z*cos(a_) + 1), Abs(a_) > 0, Abs(z) > 1)
    >>>

    basic pairs table:
    (DiracDelta(n), 1, True, True)
    (DiracDelta(n-a),z**(-a),a>0, Abs(z) > 0)
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
    d = Wild('d', exclude=[n,z])
    r = Wild('r', exclude=[n,z])
    dco = lambda f: _z_args_collect(f,n)
    #nearzero = 1e-10
    # ( time domain,
    #   z domain,
    #   condition, convergence plane, preparation function )
    z_transform_table = [
    # DiracDelta
    (DiracDelta(n),
     S.One,
     S.true, S.true, dco),
    (DiracDelta(n-a),
     z**(-a),
     a>0, Abs(z) > 0, dco),
    (DiracDelta(n+a),
     z**a,
     a>0, Abs(z) < S.Infinity, dco),
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
    (cos(b*n+a),
     z*(z*cos(a)-cos(b - a))/(z**2-(2*cos(b))*z+1),
     And(Abs(b)>0,Abs(a)>0), Abs(z) > 1, dco),
    # [1] fn to Fz
    (c*(r**n)*cos(b*n+a),
     c*z*(z*cos(a)-cos(b-a))/(z**2-2*r*cos(b)*z+r**2),
     And(Abs(b)>0,Abs(a)>0,Abs(r)>0,Abs(c)>0),
     Abs(z) > 1, dco),
    # [1] Fz to Fn
    (c*(r**n)*cos(b*n+a),
     c*z*(z*cos(a)-d)/(z**2-2*r*cos(b)*z+r**2),
     And(Abs(b)>0,Abs(a)>0,Abs(c)>0,Abs(d-cos(b-a))>0),
     Abs(z) > 1, dco)
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

def _z_pairs_properties(f,n,z, apply_properties=True,itera = 0):
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
    itera : pairs and properties table reference operation

    Examples
    ========
    >>> from sympy.integrals.ztransforms import _z_pairs_properties
    >>> from sympy.abc import n,z
    >>> import sympy as sym
    >>> f = sym.Heaviside(n)
    >>> _z_pairs_properties(f,n,z)
    (z/(z - 1), Abs(z) > 1, True)
    >>> f = n*sym.Heaviside(n)
    >>> _z_pairs_properties(f,n,z)
    (z/(z - 1)**2, Abs(z) > 1, True)
    >>>
    >>> f = sym.DiracDelta(n)
    >>> _z_pairs_properties(f,n,z)
    (1, True, True)
    >>>
    >>> f = sym.cos(3*n)*sym.Heaviside(n)
    >>> _z_pairs_properties(f,n,z)
    (z*(z - cos(3))/(z**2 - 2*z*cos(3) + 1), Abs(z) > 1, True)
    >>>
    >>> f = sym.sin(2*n)*sym.Heaviside(n)
    >>> _z_pairs_properties(f,n,z)
    (z*sin(2)/(z**2 - 2*z*cos(2) + 1), Abs(z) > 1, True)
    """

    # Try _z_pairs_table
    Fz = None

    # extract constant factors from Fz
    k, f = f.as_independent(n, as_Add=False)
    z_pairs = _z_pairs_table(n, z) # get table
    debug('['+str(itera)+'] _z_pairs_properties -----')
    debug('  k:',k,' ; f:',f)
    # lookup on table z_pairs
    i = 0
    while (i<len(z_pairs)) and (Fz is None):
        n_dom, z_dom, check, plane, prep = z_pairs[i]
        ma = prep(f).match(n_dom)
        if ma or ma=={}: # pair found
            debug('\n['+str(itera)+']_z_pairs_table match:')
            debug('  k:',k,' ; f:',f,
                  '\n  z_pair f[n]:',n_dom,
                  '\n  z_pair F[z]:',z_dom)
            #try:
            if check.xreplace(ma):
                Fz = (k*z_dom.xreplace(ma),
                      plane.xreplace(ma),
                      S.true)
            #except Exception:
            #    debug(' _z_pairs_table did not match.')
        i = i+1 # next pair

    # f[n] did not match z_pairs, try z_properties()
    if Fz is None and apply_properties is True and f.has(n):
        Fz = _z_properties(f, n, z,itera=itera+1)
        if Fz is not None:
            Fz = (k*Fz[0],Fz[1],Fz[2])
    debug('\n['+str(itera)+'] _z_pairs_properties F[z]:',Fz)
    return Fz

def _z_properties(f,n,z,itera=0):
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
    itera : pairs and properties table reference operation

    Examples
    ========
    >>> from sympy.integrals.ztransforms import _z_properties
    >>> from sympy.abc import n,z
    >>> import sympy as sym
    >>> f = n*sym.Heaviside(n)
    >>> _z_properties(f,n,z)
    (z/(z - 1)**2, Abs(z) > 1, True)
    >>>
    >>> f = n*n*sym.Heaviside(n)
    >>> _z_properties(f,n,z)
    (z*(z + 1)/(z - 1)**3, Abs(z) > 1, True)
    >>>
    >>> f = sym.Heaviside(n-2)
    >>> _z_properties(f,n,z)
    (1/(z*(z - 1)), Abs(z) > 1, True)
    >>>
    >>> f = ((sym.S.One/2)**2)*sym.Heaviside(n-2)
    >>> _z_properties(f,n,z)
    (1/(4*z*(z - 1)), Abs(z) > 1, True)
    >>>
    >>> f = ((1/2)**2)*sym.Heaviside(n-2)
    >>> _z_properties(f,n,z)
    (0.25/(z*(z - 1)), Abs(z) > 1, True)
    >>>
    >>> f = sym.Heaviside(-n)
    >>> _z_properties(f,n,z)
    (-1/(z - 1), Abs(z) < 1, True)
    >>>
    >>> f = -sym.Heaviside(-n-1)
    >>> _z_properties(f,n,z)
    (z/(z - 1), Abs(z) < 1, True)
    >>>
    >>> f = sym.cos(2*n + 0.25)*sym.Heaviside(n)
    >>> _z_properties(f,n,z)
    (z*(z*cos(0.25) - cos(1.75))/(z**2 - 2*z*cos(2) + 1), Abs(z) > 1, True)
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

    debug('\n['+str(itera)+'] z_transform properties check -------------')
    if f.is_Add:
        msg = 'f[n] is an Additive expression, suggested to use z_transform(f,n,z) instead'
        raise ZTransformError('z_properties',f,msg)

    # split f[n]: (k,n**a,a**n,g[n],Heaviside[n])
    # Heaviside do not match WildFunction
    k, f_powna, f_powan, gn, un = _split_f(f,n)
    f = f/k

    # z_property a**n*f[n] <--> F[z/a]
    if not f_powan==S.One : # a**n factor match
        ma_powan  = f_powan.match(Pow(a,n))
        ma_powan_ = f_powan.match(Pow(a,-n))
        f = f/f_powan
        debug(' _z_property a**n*f[n] ----------')
        debug('   k:',k,' ; a**n:',f_powan,' ; f:',f)
        Fz = _z_pairs_properties(f,n,z,itera=itera+1)
        if Fz is not None:
            a_k = 1 # from a**n
            if ma_powan:
                a_k = ma_powan[a]
            if ma_powan_:
                a_k = S.One/ma_powan_[a]
            #debug('Fz',Fz,ma_powan,ma_powan_,a_k)
            Fz_k = k*factor(Fz[0].subs(z,z/a_k))

            cond_z = list(Fz[1].args)
            num_i = None
            fun_i = None
            k_cond = None
            for i in range(len(cond_z)):
                cond_z[i] = cond_z[i].subs(z,z/a_k)
                if not cond_z[i].has(z):
                    num_i = i
                else: # cond_z[i].has(z)
                    fun_i = i
                    k_cond = (cond_z[i].as_independent(z, as_Add=False))[0]
            if k_cond is not None:
                cond_z[fun_i] = cond_z[fun_i]/k_cond
                cond_z[num_i] = cond_z[num_i]/k_cond
            cond_z = Fz[1].func(cond_z[0],cond_z[1])
            cond_z = simplify(cond_z)
            Fz = (Fz_k,cond_z,Fz[2])
            debug('\n['+str(itera)+']  _z_property a**n*f[n] <--> F[z/a]:\n  ',Fz)

    # z_property n*f[n] <--> -z*dF(z)/dz --------
    elif Fz is None and not f_powna==S.One : # n**a factor match
        #ma_powna = f_powna.match(Pow(n,a))
        f = f/n
        debug('\n _z_property n*f[n]',
              '\n  n**a:',n,' ; f:',f)
        Fz = _z_pairs_properties(f,n,z,itera=itera+1)
        if not Fz is None:
            Fz = (k*(-z)*factor(diff(Fz[0],z,1)),
                  Fz[1], Fz[2])
            debug('\n['+str(itera)+']  _z_property n*f[n] <--> -z*diff(Fz):\n  ',Fz)

    # match Heaviside(n) only
    elif Fz is None and not(un==S.One) and gn==S.One:
        ma_un = f.match(a*Heaviside(y))
        args_n = ma_un[y].collect(n)
        ma_args = args_n.match(a*n-b)
        debug(' ma_un:',ma_un,' ; ma_args:',ma_args)
        f = Heaviside(n)
        Fz = _z_pairs_properties(f,n,z,itera=itera+1)
        # time shift only
        if  ma_args[a]>0 and Abs(ma_args[b])>=0:
            b_shift = ma_args[b]
            debug(' _z_property time_shift u[n-b]--------')
            debug('  b_shift:',b_shift,' ; f:', f,)
            if not Fz is None:
                Fz = (k*factor(z**(-b_shift)*Fz[0]),
                      Fz[1], Fz[2])
            msg = '['+str(itera)+']  _z_property time_shift'
            debug(msg,'f[n-b] <--> z**(-b)*F +Fsum:\n ', Fz)
        # u time reversal and shift
        if ma_args and ma_args[a]<0 and Abs(ma_args[b])>=0:
            b_shift = -ma_args[b]
            debug(' _z_property time_reversal and shift u[-(n+b)]---')
            debug('  b_shift:',b_shift,' ; f:',f,)
            if Fz is not None:
                Fz_k = k*factor((Fz[0].subs(z,1/z))*z**(-b_shift))
                plane_reversal = Fz[1].func(Fz[1].args[1],
                                            Fz[1].args[0])
                plane_reversal = simplify(plane_reversal)
                Fz = (Fz_k,plane_reversal, Fz[2])
            msg = '['+str(itera)+']  _z_property time_reversal and shift'
            debug(msg,' u[-(n+b)]:\n ',Fz)

    # match g function only, no Heaviside
    elif Fz is None and not(gn==S.One) and un==S.One:
        ma_gn = f.match(g)
        args_n = ma_gn[g].args[0].collect(n)
        ma_args = args_n.match(a*n-b)
        debug(' ma_gn ; ma_args : ',ma_gn,ma_args)
        # g time shift only
        if not ma_args[a]==S.Zero and not Abs(ma_args[b])==S.Zero:
            b_shift = ma_args[b]
            f = ma_gn[g].subs(args_n,ma_args[a]*n)
            debug(' _z_property time_shift --------')
            debug('  b_shift:',b_shift,' ; f:',f,'\n')
            Fz = _z_pairs_properties(f,n,z,itera=itera+1)
            Fsum = 0
            m = b_shift # Sumation term x[n-m]u[n]
            # DiracDelta does not apply
            # right shift n
            if not ma_args[b]==S.Zero and not f.has(DiracDelta):
                fn_gz = f.subs(args_n,-n)*z**n
                Fsum = summation(fn_gz,(n,1,m))
            # left shift n
            if ma_args[b]<S.Zero and not f.has(DiracDelta):
                fn_gz = f.subs(args_n,n)*z**(-n)
                Fsum = -summation(fn_gz,(n,0,-m-1))
            Fsum = expand(z**(-m)*Fsum)
            # apply shift
            if not Fz is None:
                Fz = (Fsum + k*factor((z**(-b_shift))*Fz[0]),
                      Fz[1], Fz[2])
                msg = '\n['+str(itera)+'] _z_property time_shift'
                debug(msg,' z**(-b)*F[z]+Fsum:\n ', Fz)
        # g time reversal and shift
        if ma_args[a]<0 and Abs(ma_args[b])>0:
            b_shift = -ma_args[b]
            f = ma_gn[g].subs(args_n,ma_args[a]*n)
            debug(' _z_property time reversal and shift ---')
            debug('  b_shift:',b_shift,' ; f:',f)
            Fz = _z_pairs_properties(f,n,z,itera=itera+1)
            Fsum = 0
            m = b_shift # Sumation term x[n-m]u[n]
            # DiracDelta does not apply
            # right shift n
            if ma_args[b]>0 and not f.has(DiracDelta):
                fn_gz = f.subs(args_n,-n)*z**n
                Fsum = summation(fn_gz,(n,1,m))
            # left shift n
            if ma_args[b]<0 and not f.has(DiracDelta):
                fn_gz = f.subs(args_n,n)*z**(-n)
                Fsum = -summation(fn_gz,(n,0,-m-1))
            Fsum = expand(z**(-m)*Fsum)
            if not Fz is None: # apply shift
                Fz = (Fsum+k*factor((z**(-b_shift))*(Fz[0].subs(z,1/z))),
                      Fz[1], Fz[2])
                msg = '\n['+str(itera)+'] _z_property time reversal and shift'
                debug(msg,' z**(-b)*F[z]+Fsum:\n ', Fz)

    # match g*Heaviside
    if Fz is None and not gn==S.One and not un==S.One:
        ma_gu = f.match(g*Heaviside(y))
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
        # ma_g  = g_arg.match(a*n-b)
        # check args g(a1*n-b1)*u(a2*n-b2)
        PQ = pdiv(g_arg,u_arg)
        debug('  g_k*(g_arg): ',g_k,'*(',g_arg,')')
        debug('  u_k*(u_arg): ',u_k,'*(',u_arg,')')
        debug('  g_arg/u_arg poly division:',
              g_k,'*(',PQ[0],'*(',u_arg,')+',PQ[1],')')

        if ma_u[a]==S.One and ma_u[b]==S.Zero: # g*u[n]
            f = ma_gu[g]
            debug('  f:',f,)
            Fz = _z_pairs_properties(f,n,z,itera=itera+1)
            if not Fz is None:
                plane_z =  Fz[1]
                if p_invert:
                    plane_z = simplify(Fz[1].func(Fz[1].args[1],
                                                  Fz[1].args[0]))
                Fz = (k*Fz[0], plane_z, Fz[2])
                debug('\n['+str(itera)+']  _z_property g*u[n] Fz:\n ',Fz)

        else: # g[q(n-k)+p]*u[n-k]
            g_arg_m = PQ[0]*n+PQ[1]
            f = ma_gu[g].subs(ma_gu[g].args[0],g_k*g_arg_m)
            debug('  f:',f,)
            Fz = _z_pairs_properties(f,n,z,itera=itera+1)
            if not Fz is None:
                Fz = (k*Fz[0], Fz[1], Fz[2])
                debug('\n['+str(itera)+'] _z_property g[q(n-k)+p]*u[n-k] Fz:\n ',Fz)

    # clean floats as integers, 4z/(4z-1) to z/(z-1/4)
    if not Fz is None:
        Fz0 = _simplify_z(Fz[0],z)
        Fz0 = _round_float_is_int(Fz0)
        Fz1 = _round_float_is_int(Fz[1])
        Fz2 = Fz[2]
        Fz = (Fz0,Fz1,Fz2)
    debug('\n['+str(itera)+'] _z_properties F[z]:',Fz)

    return Fz

def _split_f(f,n):
    """
    auxiliary function to split f[n] to:
    constant k , a**n, n**a,function g[n], Heaviside u[n]
    """
    a = Wild('a', exclude=[n])
    y = Wild('y')
    g = WildFunction('g', nargs=1)
    k, f = f.as_independent(n, as_Add=False)
    f_un = S.One
    f_gn = S.One
    f_powna = S.One
    f_powan = S.One
    for factor_k in Mul.make_args(f):
        fk_ma = False
        ma = factor_k.match(Pow(n,a))
        if ma and not ma[a]==S.Zero and isinstance(ma[a],Integer):
            f_powna = f_powna*factor_k
            fk_ma = True
        ma = factor_k.match(Pow(a,n)) # n positive
        ma_ = factor_k.match(Pow(a,-n)) # n negative
        if ma or ma_:
            if ma:
                f_powan = f_powan*factor_k
                fk_ma = True
            if ma_:
                f_powan = (1/f_powan)*factor_k
                fk_ma = True
        ma = factor_k.match(g)
        if not fk_ma and ma:
            f_gn = f_gn*factor_k
            fk_ma = True
        ma = factor_k.match(Heaviside(y))
        if not fk_ma and ma:
            f_un = f_un*factor_k
            fk_ma = True
        if fk_ma is False:
            f_gn = f_gn*factor_k

    f_parts = (k,f_powna,f_powan,f_gn,f_un)
    debug('  f[n]: (k,n**a,a**n,g[n],Heaviside[n])')
    debug('  f[n]:',f_parts)
    return f_parts

def _z_transform_summation(f, n, z):
    """
    z_transform unilateral summation,for a simple term (auxiliar function)
    """
    f = f+0*n
    Fz = None
    if f.has(n):
        fzn = f*(z**(-n))
        # unilateral summation
        Fz_sum = summation(fzn,(n,S.Zero,S.Infinity))
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
    >>> from sympy.integrals.ztransforms import z_transform
    >>> from sympy.abc import n,z
    >>> from sympy.core import S
    >>> import sympy as sym
    >>> f = sym.DiracDelta(n)
    >>> z_transform(f,n,z)
    (1, True, True)
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
    (z/(z - 1/2), Abs(z) > 1/2, True)
    >>>
    >>> f = ((1/2)**n)*sym.Heaviside(n)
    >>> z_transform(f,n,z)
    (z/(z - 0.5), Abs(z) > 0.5, True)
    >>> f = ((S.One/2)**n)*n*sym.Heaviside(n)
    >>> z_transform(f,n,z)
    (z/(2*(z - 1/2)**2), Abs(z) > 1/2, True)
    >>>
    >>> f = sym.DiracDelta(n) + sym.cos(4*n)
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
        FT = _z_pairs_properties(f, n, z) #[0]
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
            term_z = _z_pairs_properties(fn, n, z)
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

def _z_pairs_prop_inverse(F,z,n,apply_properties=True,itera=0):
    """
    z_transform based on pairs table and properties,
    otherwise returns None .

    Explanation
    ===========
    z_transform pairs table
    http://blog.espol.edu.ec/telg1001/transformada-z-tabla/

    Examples
    ========
    >>> from sympy.integrals.ztransforms import _z_pairs_prop_inverse
    >>> from sympy.abc import z,n
    >>> F = z/(z - 1)
    >>> _z_pairs_prop_inverse(F,z,n)
    (Heaviside(n), Abs(z) > 1, True)
    >>> F = z/(z - 1)**2
    >>> _z_pairs_prop_inverse(F,z,n)
    (n*Heaviside(n), Abs(z) > 1, True)
    >>>
    >>> from sympy.functions.elementary.trigonometric import cos
    >>> F = z*(z*cos(0.25) - cos(1.75))/(z**2 - 2*z*cos(2) + 1)
    >>> _z_pairs_prop_inverse(F,z,n)
    (cos(2*n + 0.25)*Heaviside(n), Abs(z) > 1, True)

    Parameters
    ==========
    F : F[z] expression at z_domain
    z : variable on z_domain
    n : variable on discrete n_domain
    apply_properties:
       True: check with properties if not found at table
       False: check only with pairs table
    itera : pairs and properties table reference operation
    """
    fn = None
    F = sympify(F) # if F is constant
    # check expression F[z] = P/Q
    debug('['+str(itera)+'] _z_pairs_prop_inverse ----')
    PQ = _simplify_z(F,z,PQ_args=True)
    k_sign = PQ['k_sign']
    k = sympify(PQ['k'])
    Fz = PQ['F']
    P_zc = PQ['P_zc'] # z**(symbol)
    Q_zc = PQ['Q_zc'] # z**(symbol)
    zc = not(P_zc==S.One and Q_zc==S.One)

    # Try _z_pairs_table
    msg = 'pair'
    z_pairs = _z_pairs_table(n, z)
    z_pairs_len = len(z_pairs)
    i = 0
    while i<z_pairs_len and fn is None and not zc:
        # [n_dom, z_dom, check, plane, prep] from zpairs
        zp = z_pairs[i]
        n_dom = zp[0]
        z_dom = zp[1]
        check = zp[2]
        plane = zp[3]
        ma_z = Fz.match(z_dom)
        #debug('z_pair F[z]:',z_dom, ma_z,Fz)
        if ma_z or ma_z=={}:
            debug(' _z_pairs_table match:','\n  F:',Fz)
            debug('  z_pair F[z]:',z_dom,' ; ma_z:',ma_z)
            debug('  z_pair f[n]:',n_dom)
            #try:
            c = check.xreplace(ma_z)
            debug('  try,check  : %s -> %s'%(check, c))
            if c or F==1:
                fn = n_dom.xreplace(ma_z)
                if not fn.has(Heaviside) and not fn.has(DiracDelta):
                    fn = fn*Heaviside(n)
                plane_z = plane.xreplace(ma_z)
                cond_z  = S.true
                fn = (fn, plane_z, cond_z)
            #except Exception:
            #    debug(' _z_pairs_table did not match.')
        i = i+1 # next pair

    # f[n] did not match z_pairs, try z_properties_inverse()
    if fn is None and apply_properties is True:
        fn = _z_properties_inverse(F, z, n,PQ,itera=itera+1)
        msg = 'prop'

    # clean floats as integers as n
    if not fn is None:
        if msg=='pair':
            fn0 = _round_float_is_int(k_sign*k*fn[0])
        if msg=='prop':
            fn0 = _round_float_is_int(fn[0])
        fn1 = _round_float_is_int(fn[1])
        fn2 = _round_float_is_int(fn[2])
        fn = (fn0,fn1,fn2)
        debug('['+str(itera)+']   '+msg+' f[n]:',fn)
    return fn

def _z_properties_inverse(F,z,n,PQ=None, itera=0):
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
    case, n*(n-1)*(n-2)...(n-m+1)*(gamma**n)*u[n]/(m!*gamma**m)')
         <--> z/(z-gamma)**m+1')

    Examples
    ========
    >>> from sympy.integrals.ztransforms import _z_properties_inverse
    >>> from sympy.abc import z,n
    >>> from sympy.core import S
    >>> F = z/(z - 1)**2
    >>> _z_properties_inverse(F,z,n)
    (n*Heaviside(n), Abs(z) > 1, True)
    >>>
    >>> F = z/(z - S.Half)
    >>> _z_properties_inverse(F,z,n)
    (Heaviside(n)/2**n, Abs(z) > 1/2, True)
    >>>
    >>> F = 1/z**2
    >>> _z_properties_inverse(F,z,n)
    (DiracDelta(n - 2), True, True)
    >>> F = z**2
    >>> _z_properties_inverse(F,z,n)
    (DiracDelta(n + 2), True, True)
    >>>
    >>> F = z*2/(z-S.One/2)**2
    >>> _z_properties_inverse(F,z,n)
    (4*n*Heaviside(n)/2**n, Abs(z) > 1/2, True)
    >>>

    Parameters
    ==========
    F : F[z] expression at z_domain
    z : variable on z_domain
    n : variable on discrete n_domain
    PQ : F parts as numerator P and Q denominator, degree, zeros and poles
    itera : pairs and properties table reference operation
    """

    fn = None
    F = sympify(F) # if F is constant
    debug('\n['+str(itera)+'] z_properties_inverse_z_transform')
    a = Wild('a', exclude=[n,z])
    b = Wild('b', exclude=[n,z])
    c = Wild('c', exclude=[n,z])
    r = Wild('r', exclude=[n,z])

    # F poles and zeros
    if PQ is None:
        PQ = _roots_z(F,z)
    else: # PQ pre-processed by _z_pairs_prop_inverse
        PQ = _roots_z(F,z,PQ)
    P = PQ['P']
    Q = PQ['Q']
    P_zeros  = PQ['P_zeros']
    P_degree = PQ['P_degree']
    Q_poles  = PQ['Q_poles']
    Q_degree = PQ['Q_degree']
    k_sign = PQ['k_sign']
    k = PQ['k']
    F = PQ['F']
    P_zc = PQ['P_zc']
    Q_zc = PQ['Q_zc']

    # z**(c) factor at F[z]
    ma_Pzc = None
    ma_Qzc = None
    if not Q_zc==S.One: # no constant
        ma_Qzc = Q_zc.match(z**c)
    if not P_zc==S.One : # no constant
        ma_Pzc = P_zc.match(z**c)

    # denominator
    ma_Q1 = None # real Q_poles
    ma_Q2 = None # complex Q_poles, Q degree 2
    if Q.has(z):
        ma_Q1 = Q.match((a*z-b)**c)
        ma_Q2 = Q.match(a*z**2+ b*z + r**2)
        if ma_Q2 and ma_Q2[a]==S.Zero:
            ma_Q2 = None

    # numerator
    ma_P1 = None
    ma_P2 = None
    if P.has(z):
        ma_P1 = P.match(a*z+b)
        ma_P2 = P.match(r*z*(a*z+b))
        if ma_P2 and ma_P2[a]==S.Zero:
            ma_P2 = None

    # debug match fromf P an Q
    if ma_P1:
        debug('  ma_P1 (a*z+ b):   ',ma_P1)
    if ma_P2:
        debug('  ma_P2 r*z*(b*z+a):',ma_P2)
    if ma_Q1:
        debug('  ma_Q1 (a*z-b)**c :',ma_Q1)
    if ma_Q2:
        debug('  ma_Q2 a*z**2+ b*z + r**2:\n ',ma_Q2)

    # no Factor P_zc or Q_zc, both are S.One
    c_not_zc = P_zc==S.One and Q_zc==S.One
    # Q1 has no symbols and Q1_degree<=2, P1 is not only z
    # avoid z/(z-a)**3 is symbolic and Q_degree>2
    c_nf1 = not(ma_Q1 and len(ma_Q1[b].free_symbols)==1
                and ma_Q1[c]>2
                and ma_P1 and ma_P1[b]==S.Zero)
    # ma_Q1 and z/(z-a)**c
    c_nf2 = ma_Q1 and (1 <= P_degree <Q_degree)
    # _z_property nf[n] <--> -z*diff(F[z])
    #  ma_Q1 = Q.match(a*(z-b)**c)
    if fn is None and c_not_zc and c_nf1 and c_nf2:
        FunI = factor(F/(-z),z)
        from sympy.integrals.integrals import integrate
        F = factor(integrate(FunI,z))
        F0 = 0 # Integral constant
        C = -F.subs(z,0)+F0
        FunC = factor(F+C)
        debug('\n _z_property multiply nf[n] <--> -z*diff(F[z])')
        debug(' F = integrate(factor(-F/z,z):\n\t=',F,'+',C)
        debug(' \t=',FunC)
        if not FunC.has(log):
            fn = _z_pairs_prop_inverse(FunC, z, n,itera=itera+1)
        else:
            fn = None
        if not fn is None:
            fn = (k_sign*k*n*fn[0],fn[1],fn[2])
            debug('\n _z_property multiply nf[n]:\n  ',fn)

    #  _z_property (a**n)*f[n] <--> F(z/a)
    # ma_Q1 = Q.match((a*z-b)**c)
    c_an1 = ma_Q1 and not ma_Q1[b]==S.One # Q = (z-b)**c
    c_an2 = ma_Q2 and not ma_Q2[r]==S.One # to be used for cos[an] or sin[an]
    if fn is None and c_not_zc and c_an1 or c_an2:
        debug('\n _z_property multiply (a**n)*f[n] <--> F(z/a) ')
        k_b = None
        # Q1 is (z-b)
        if c_an1:
            k_b = ma_Q1[b]
            F = factor(F.subs(z,k_b*z))

        # Q2 quadratic a*z**2+ b*z + r**2
        # P2 r*z*(b*z+a)
        if ma_Q2 and not ma_Q2[r]==S.One and ma_Q2[a]==S.One:
            k_b = ma_Q2[r]
            F = factor(F.subs(z,k_b*z))

        debug(' k_sign:',k_sign,'; k:',k,' ; k_b:',k_b)
        debug(' F[z]:',F)

        fn = _z_pairs_prop_inverse(F, z, n, itera=itera+1)

        if fn is not None:
            fnk = k_sign*k*(k_b**n)*fn[0]
            # conditions for z
            cond_z = list(fn[1].args)
            num_position = None
            fun_position = None
            for i in range(len(cond_z)):
                cond_z[i] = cond_z[i].subs(z,z/k_b)
                if not cond_z[i].has(z):
                    num_position = i
                if cond_z[i].has(z):
                    fun_position = i
                    k_Rel = (cond_z[i].as_independent(z, as_Add=False))[0]
            if num_position is not None:
                cond_z[fun_position] = cond_z[fun_position]/k_Rel
                cond_z[num_position] = cond_z[num_position]/k_Rel
            cond_z = fn[1].func(cond_z[0],cond_z[1])

            fn = (fnk,cond_z,fn[2])
            debug('\n  _z_property multiply (a**n)*f[n]:\n ',fn)


    # _z_property time_shift f[n-b] <--> z**(-b)*F +Fsum
    c_shift0 = 0 in Q_poles # at least one pole
    c_shift1 = not P_zc==S.One or not Q_zc==S.One # z**c symbolic
    c_shift2 = (P_degree>Q_degree) # check P_zeros example: z*(z/(z-1))
    if fn is None and (c_shift0 or c_shift1 or c_shift2):
        debug(' _z_property time_shift f[n-b] <--> z**(-b)*F +Fsum ')
        b_shift = S.Zero
        msg = 'None'
        # only z**a with a positive
        if (0 in Q_poles) and ma_Qzc:
            b_shift = Q_poles[0]
            msg =' Q_zc = z**c with c positive'
        # only z**a with a positive
        if (0 in P_zeros) and ma_Pzc: #and P_zeros[0]>S.One:
            if F.has(z):
                b_shift = -(P_zeros[0]-1) # one z at numerator
            else:
                b_shift = -(P_zeros[0]) # F is a constant
            msg =' P_zc = z**-c with c negative'

        b_shift = sympify(b_shift)
        debug('   b_shift:',b_shift,' ;',msg,'\n')

        fn = _z_pairs_prop_inverse(F, z, n,itera=itera+1)
        if not fn is None:
            f_k = k_sign*k*(fn[0].subs(n,n-b_shift))
            fn = (f_k,fn[1],fn[2])
        debug('\n _z_property time_shift f[n-b] :\n ',fn)

    # _z_property time_reversal f[-n] <--> F[1/z]
    if fn is None  and k_sign<0 and not P.has(z):
        F = simplify(F.subs(z,1/z))
        debug('\n _z_property time_reversal F[1/z]-- ')
        debug(' F[z] : ', F)
        fn = _z_pairs_prop_inverse(F, z, n, itera=itera+1)
        debug('  -----fn : ', fn,k_sign,k,fn[0])
        if not fn is None:
            f_kgu = _split_f(fn[0],n)
            f_g=1
            f_u=1
            f_k = f_kgu[0]*k*k_sign
            if f_kgu[3]!=S.One:
                f_g = f_kgu[1].subs(f_kgu[3].args[0],-f_kgu[3].args[0])
            if f_kgu[4]!=S.One:
                f_u = Heaviside(-f_kgu[4].args[0])
            plane_reversal = fn[1].func(fn[1].args[1],fn[1].args[0])
            plane_reversal = simplify(plane_reversal)
            fn = (f_k*f_g*f_u,plane_reversal,fn[2])
            debug('  _z_property time_reversal f[n]:\n ',fn)

    #case: n*(n-1)*(n-2)...(n-m+1)*(gamma**n)*u[n]/(m!*gamma**m)')
    #      <--> z/(z-gamma)**m+1')
    # ma = Q.match(a*(z-b)**c) # z/(z-a)**m+1
    c_case1 = ma_Q1 and not(ma_Q1[b]==S.Zero) and ma_Q1[c]>2 \
              and P_degree==S.One
    if fn is None and c_not_zc and c_case1:
        m = ma_Q1[c]-1
        ma_P1 = None
        debug('\n case: n*(n-1)*(n-2)...(n-m+1)*(gamma**n)*u[n]/(m!*gamma**m)')
        debug('       <--> z/(z-gamma)**m+1')
        debug('   P:',P,' ; Q:',Q,
              ' ; m:',m,' ; gamma:',ma_Q1[b])
        ma = P.match(a*z+ b)
        if ma and ma[b]==S.Zero and not ma[a]==S.Zero:
            term = S.One
            for i in range(0,m,1):
                term = term*(n-i)
            debug('   factors:',term)
            fn = term*(ma_Q1[b]**n)/((ma_Q1[b]**m)*factorial(m))
            fn = (k_sign*k*fn*Heaviside(n),Abs(z) > 1,0)
            debug('   case f[n]:',fn)

    if not fn is None:
        fn0 = _round_float_is_int(fn[0])
        fn1 = _round_float_is_int(fn[1])
        fn2 = fn[2]
        fn = (fn0,fn1,fn2)
        debug('['+str(itera)+']  _z_properties_inverse f[n]:',fn)
    return fn

def inverse_z_transform(F,z,n,H0=S.Half):
    """
    inverse_z_transform of F[z] to f[n], pairs and properties tables.

    Explanation
    ===========
    F[z] is expanded with modified partial fractions as Add terms
    to simplify the process.
    Returns f[n] expression with each Add term transformed to z domain
    otherwise returns tuple(f[n], F[z] terms not processed)
    "Many of the transforms F[z] of practical interest
    are rational functions (ratio of polynomials in z), which
    can be expressed as a sum of partial fractions."[1]
    Inverse transforms can be found in a table of transform
    and applying some properties as done by:
    _z_pairs_prop_inverse(F, z, n, apply_properties=True )
    "The partial fraction method works because for every
    transformable f[n] defined for n >= 0,
    there is a corresponding unique F[z] defined for |z| > r0
    (where r0 is some constant),and vice versa."[1]

    Parameters
    ==========
    F : F[z] expression at z_domain
    z : variable on z_domain
    n : variable on discrete n_domain
    H0 : change Heaviside[0] to H0=S.One, default H0=S.Half

    Examples
    ========
    >>> from sympy.integrals.ztransforms import inverse_z_transform
    >>> from sympy.abc import z,n
    >>> F = 5*z/(z-1)
    >>> inverse_z_transform(F,z,n)
    (5*Heaviside(n), Abs(z) > 1, True)
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
        debug(' apart_z(F) as partial fractions:',term_sum)
    for term_k in term_sum:
        if len(term_sum)>1:
            debug('  term_k[z]:',term_k)
        term_n = None # not paired yet
        if not term_k==S.Zero: # avoid term_k=0
            term_n = _z_pairs_prop_inverse(term_k, z, n)
        if term_n is None and not term_k==S.Zero:
            f_noeval.append(term_k)
            term_n = None
        if term_n is not None and len(term_sum)>1:
            fT = fT + _round_float_is_int(term_n[0])
        else:
            fT = _round_float_is_int(term_n)
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
        if H0==S.One:
            fT = Heaviside_edge_toOne(fT,n)
    if len(f_noeval)>0:
        fT = [fT,f_noeval]
    return fT

def Heaviside_edge_toOne(f,n):
    fT = 0*n
    term_sum = Add.make_args(f)
    print(term_sum)
    for term_k in term_sum:
        if term_k.has(Heaviside):
            term_kH = 1
            factor_Mul = Mul.make_args(term_k)
            for factor_k in factor_Mul:
                if factor_k.has(Heaviside):
                    H_args = factor_k.args
                    H = Heaviside(H_args[0],1)
                    term_kH = term_kH*H
                else:
                    term_kH = term_kH*factor_k
            fT = fT + term_kH
        else:
            fT = fT + term_k
    return(fT)

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
    >>> from sympy.integrals.ztransforms import apart_z
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
    F = sympify(F) # check if only a constant
    if not z in F.free_symbols:
        return F

    debug('\napart_z(F) :',F)
    a = Wild('a', exclude=[z])
    b = Wild('b', exclude=[z])
    c = Wild('c', exclude=[z])

    # Fz poles and zeros
    F = factor(F)
    PQ = _simplify_z(F,z,PQ_args=True)
    # P = PQ['P']
    # Q = PQ['Q']
    k = PQ['k']
    k_sign = PQ['k_sign']
    F = PQ['F']
    P_zc = PQ['P_zc'] # z**(const+symbol)
    Q_zc = PQ['Q_zc'] # z**(const+symbol)

    if P_zc is not S.One:
        P_zc = P_zc/z
        F = F*z
    if Q_zc is not S.One:
        Q_zc = Q_zc/z
        F = F/z

    # apply partial fractions to F[z]/z
    F_z = factor(F/z,z)

    # Fz poles and zeros
    PQ_z = _roots_z(F_z,z)
    P_z = PQ_z['P']
    #Q_z = PQ_z['Q']
    #P_z_zeros  = PQ_z['P_zeros']
    P_z_degree = PQ_z['P_degree']
    Q_z_degree = PQ_z['Q_degree']
    Q_z_unique = PQ_z['Q_unique']

    # avoid: z**3/(z-1) ; z**2/(z-a) ; z**2/(z-a)**2
    if (P_z.has(z) and (len(Q_z_unique)>1 or P_z_degree==(Q_z_degree-1))):
        F_z = apart(F_z,z)
    debug('  apart(F/z):',F_z)

    #  restore F[z]*z
    term_suma = Add.make_args(F_z)
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
        #  restore F[z]*z and constant k and z**c factors
        term_k = _simplify_z(term_k,z)
        Fzp = Fzp + (term_k*z)*(k*k_sign*P_zc/Q_zc)

    Fzp = _round_float_is_int(Fzp)
    debug('\n apart z*(Fz/z)  :',Fzp,'\n')
    return Fzp

def _round_float_is_int(F, nearzero=1e-10):
    """
    auxiliary function to round floats in F expresion with no decimal part,
    ending with .0 to simplify reading expression
    """
    for a in preorder_traversal(F):
        if isinstance(a, Float):
            r = a%1
            if Abs(r)<nearzero:
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
    F = powsimp(F,z)
    P,Q = F.as_numer_denom()
    debug('_simplify_z:',F)

    # P numerator, set apart z**c symbolic for leadcoef
    P_zc = S.One
    Q_zc = S.One
    for factor_k in Mul.make_args(P):
        ma = factor_k.match(z**c)
        if ma and not ma[c]==S.Zero: # z**c found
            # check z**(1-c)=(z)*(z**-c)
            factor_zc = expand(factor_k)
            P_factor, Q_factor = expand(factor_zc).as_numer_denom()
            Q_zc = Q_zc*Q_factor
            if Q.has(z): # example: (z**c)*z/(z-1)
                P_zc = P_zc*(P_factor)/z
            else: # only Q == 1 denominator
                P_zc = P_zc*factor_k
    P = powsimp(P/P_zc*Q_zc)
    # leadcoef
    k = LC(P,z)
    k = _round_float_is_int(k) #constant factor
    P = powsimp(P/k)
    # P as monic, use rationals at expressions
    if P.has(z):
        P_rational = 1
        for factor_k in Mul.make_args(P):
            if factor_k.is_Pow and factor_k.has(z):
                factor_args = factor_k.args
                P_rational = P_rational*Pow(P_rational*monic(factor_args[0]),factor_args[1])
            elif not(factor_k.is_Pow) and factor_k.has(z): # not a constant coefficient
                P_rational = P_rational*monic(factor_k)
            elif not factor_k.has(z):
                P_rational = P_rational*1
        P = P_rational.as_expr()

    if not k==S.One or not P_zc==S.One:
        debug('  P_leadcoef:',k,' ; P_zc:',P_zc,' ; P:',P)
    if not Q_zc==S.One:
        debug('      Q_zc:',Q_zc,' ; from P')

    # Q denominator review, # set apart z**c symbolic for leadcoef
    Q_zc0 = S.One
    for factor_k in Mul.make_args(Q):
        ma = factor_k.match(z**c)
        if ma and not ma[c]==S.Zero:
            Q_zc0 = Q_zc0*factor_k

    Q = powsimp(Q/Q_zc0)
    Q_zc = Q_zc*Q_zc0
    # leadcoef
    Q_leadcoef = LC(Q,z)
    Q_leadcoef = _round_float_is_int(Q_leadcoef) #constant factor
    Q = Q/Q_leadcoef
    # Q as monic, Divide all coefficients of Q by LC(Q)
    if Q.has(z):
        Q_rational = 1
        for factor_k in Mul.make_args(Q):
            if factor_k.is_Pow:
                factor_args = factor_k.args
                Q_rational = Q_rational*Pow(Q_rational*monic(factor_args[0]),
                                            factor_args[1])
            elif factor_k.has(z): # not a constant coefficient
                Q_rational = Q_rational*monic(factor_k)
        Q = Q_rational.as_expr()
    k = k/Q_leadcoef
    if not Q_leadcoef==S.One or not Q_zc==S.One:
        debug('  Q_leadcoef:',Q_leadcoef,' ; Q_zc:',Q_zc,' ; Q:',Q)

    # denominator Q = z**2 + 2*r*cos(a)z + r**2
    # pair 12c at http://blog.espol.edu.ec/telg1001/transformada-z-tabla/
    ma = Q.match(z**2+a*z + b)
    if ma and not ma[a]==S.Zero and not ma[b]==S.Zero:
        P,Q,k = _simplify_Q2(P,Q,k,z,nearzero=1e-10)

    # separate sign from k,
    k = sympify(k) # k is not symbolic
    k_sign = 1
    if len(k.free_symbols)==0 and k<0:
        k_sign = -1
        k = Abs(k)

    F = _round_float_is_int(P/Q)
    k = _round_float_is_int(k)
    PQ = {'P':P,'Q':Q,
          'k_sign':k_sign,'k':k,'F':F,
          'P_zc': P_zc,'Q_zc': Q_zc}
    if not k==S.One or not P_zc==S.One or not Q_zc==S.One:
        debug('  _simplify_z k_sign:',k_sign,' ; k:',k,' ; P_zc:', P_zc,' ; Q_zc:', Q_zc)
    debug('    F[z]:',F)
    # returns
    if PQ_args:
        F = PQ
    else:
        F = _round_float_is_int(k_sign*k*(P_zc/Q_zc)*F)
    return F

def _simplify_Q2(P,Q,k,z,nearzero=1e-10):
    '''
    denominator Q = z**2 + 2*r*cos(a)z + r**2
    pair 12c at http://blog.espol.edu.ec/telg1001/transformada-z-tabla/
    '''
    a = Wild('a', exclude=[z])
    b = Wild('b', exclude=[z])
    c = Wild('c', exclude=[z])
    ma = Q.match(z**2+a*z + b)
    check_Q2 = False
    if ma and not ma[a]==S.Zero and not ma[b]<S.Zero:
        k_r = sqrt(ma[b])
        if len(ma[a].free_symbols)==0:
            k_b = acos(ma[a]/(-2*k_r)).evalf()
            Q = z**2 - 2.0*k_r*cos(k_b, evaluate=False)*z+k_r**2
            check_Q2 = True
        else: # a is symbolic
            # case: b*z*sin(a)/(z**2 - 2*z*cos(a) + 1)
            ma = Q.match(z**2 - c*z*cos(a) + b)
            if ma:
                k_r = sqrt(ma[b])
                k_b = ma[a]
                if len(k_r.free_symbols)==0 and Abs(ma[c]/2-k_r)<nearzero:
                    Q = z**2 - (2.0)*k_r*cos(k_b)*z+k_r**2
                    check_Q2 = True
        if check_Q2:
            debug('  ma_Q2 (z**2+a*z + b):',ma,
                  '\n  k_b:',k_b,' ; k_r:',k_r)
        debug('  Q:',Q,' ; check_Q2:',check_Q2)

        check_P = False
        #check P = a*z+b for sin[beta*n]
        ma = P.match(a*z+b)
        if check_Q2 and ma and not ma[a]==S.Zero and ma[b]==S.Zero:
            debug('  ma_P1 (a*z+b):',ma)
            if len(k.free_symbols)==0:
                k = (k/sin(k_b)/k_r).evalf()
            else:
                # case: b*z*sin(a)/(z**2 - 2*z*cos(2) + 1)
                k = (k/sin(k_b,evaluate=False)/k_r)
            P = z*k_r*sin(k_b, evaluate=False)
            check_P = True

        #check P = c*z*(a*z+b) for cos[beta*n]
        ma = P.match(c*z*(a*z+b))
        if check_Q2 and not check_P and ma and ma[a]==S.One \
           and not ma[b]==S.One:
            k_P2 = (ma[b]/cos(k_b)/k_r).evalf()
            debug('  P:',P,' ' ,ma,'\n  k_P2:',k_P2)
            k = (k*ma[a]*ma[c]).evalf()
            a_Q = ma[a]
            # P match(z*(z-cos(a))
            if Abs(Abs(k_P2)-S.One)<nearzero: # P match(z*(z-cos(a))
                P = z*(a_Q*z-k_r*cos(k_b,evaluate=False))
                debug('  P match(z*(z-cos(a)):',P)
            else: # case: cos[beta*n+theta]
                P = P*k
                k = 1
                # denominator
                ma = Q.match(z**2+a*z+b)
                aq = ma[a]/2
                R = sqrt(ma[b])
                # numerator
                ma = P.match(c*z*(a*z+b))
                A = ma[c]*ma[a]
                B = ma[c]*ma[b]
                # calc theta
                rP = ((A**2)*(R**2)+(B**2)-2*A*aq*B).evalf()
                rQ = ((R**2)-(aq**2)).evalf()
                kr = sqrt(rP/rQ)
                beta = acos(-aq/Abs(R)).evalf()
                theta = atan((A*aq-B)/(A*sqrt(rQ))).evalf()
                debug('  cos[beta*n+theta]',
                      '  beta:',beta,' ; theta:',theta,' ; kr:',kr)
                P = kr*z*(z*cos(theta,evaluate=False)-R*cos(beta-theta,evaluate=False))

        debug('  P:',P)
    return P,Q,k

def _roots_z(F,z,PQ=None):
    """
    auxiliary function to analyze F[z] as P/Q with numerical
    and symbolic parts
    """
    a = Wild('a', exclude=[z])
    b = Wild('b', exclude=[z])
    c = Wild('c', exclude=[z])
    # Fz poles and zeros
    if PQ is None:
        PQ = _simplify_z(F,z,PQ_args=True)
    P = PQ['P']
    Q = PQ['Q']
    #k_sign = PQ['k_sign']
    #k = PQ['k']
    #Fz = PQ['F']
    #P_zc = PQ['P_zc'] # z**(symbol)
    #Q_zc = PQ['Q_zc'] # z**(symbol)

    # P numerator review
    P_zeros = {}
    # P_zeros by symbols,
    #check roots((z-0.5)**3) numeric vs match
    ma = P.match((a*z+b)**c)
    cond = ma and not ma[c]==S.Zero # not a constant
    if cond and not ma[a]==S.Zero and not ma[b]==S.Zero:
        P_zeros = {-ma[b]:ma[c]}
    else:
        P_zeros = roots(Poly(P,z)) # consider symbols z/(z-a)

    # P_degree
    P = factor(P)
    P_degree = degree(Poly(P,z)) # consider symbols z/(z-a)
    if not PQ['P_zc']==S.One: # P has z**(symbol+number)
        ma = PQ['P_zc'].match(z**a)
        if S.Zero in P_zeros:
            P_zeros[0] = P_zeros[0] + ma[a]
        else:
            P_zeros[0] = ma[a]
            P_degree = P_degree + ma[a]
    if P_degree is None:
        P_degree = degree(P,z)

    # Q denominator review
    Q_poles = {}
    # Q_poles by symbols,
    #check roots((z-0.5)**3) numeric vs match
    ma = Q.match((a*z+b)**c)
    cond = ma and not ma[c]==S.Zero # not a constant
    if cond and not ma[a]==S.Zero and not ma[b]==S.Zero:
        Q_poles = {-ma[b]:ma[c]}
    else:
        Q_poles = roots(Poly(Q,z)) # consider symbols z/(z-a)

    # Q_degree
    Q = factor(Q)
    Q_degree = degree(Poly(Q,z))
    if not PQ['Q_zc']==S.One: # Q has z**(symbol+number)
        ma  = PQ['Q_zc'].match(z**a)
        Q_degree = Q_degree + ma[a]
        if 0 in Q_poles:
            Q_poles[0] = Q_poles[0] + ma[a]
        else:
            Q_poles[0] = ma[a]

    debug('  P:',P,'\n  Q:',Q,
          '\n  P_degree:',P_degree,' ; P_zeros:',P_zeros,
          '\n  Q_degree:',Q_degree,' ; Q_poles:',Q_poles)

    PQ['P_degree'] = P_degree
    PQ['P_zeros'] = P_zeros
    PQ['Q_degree'] = Q_degree
    PQ['Q_poles'] = Q_poles
    PQ['Q_unique'] = list(Q_poles.keys())

    return PQ
