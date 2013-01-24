"""
Implementation of the trigsimp algorithm by Fu et al

work started by:
Dimitar Vlahovski @ Technological School "Electronic systems"
30.11.2011

References
==========
http://rfdz.ph-noe.ac.at/fileadmin/Mathematik_Uploads/ACDCA/
DESTIME2006/DES_contribs/Fu/simplification.pdf
"""

from sympy.simplify.simplify import simplify, powsimp, ratsimp, combsimp
from sympy.core.sympify import sympify
from sympy.functions.elementary.trigonometric import cos, sin, tan, cot
from sympy.core.core import C
from sympy.core.mul import Mul
from sympy.core.function import expand_mul
from sympy.core.add import Add
from sympy.core.symbol import Wild
from sympy.core.exprtools import factor_terms

#-------Transformation rule 0 - simplification of rational polynomials
#--Trying to simplify the expression, combine things like 3*x + 2*x, etc.

def TR0(rv):
    rv = simplify(rv)
    rv = powsimp(rv)
    rv = ratsimp(rv)
    rv = combsimp(rv)
    rv = expand_mul(rv)
    return rv

#-------Transformation rule 1 - replace sec, csc with 1/cos, 1/sin

def TR1(rv):
    try:
        rv = rv.replace(sec, lambda x: 1/cos(x))
        rv = rv.replace(csc, lambda x: 1/sin(x))
    except:
        pass
    return rv

#-------Transformation rule 2 - replace tan and cot with sin/cos and cos/sin

def TR2(rv):
    """
    Examples
    ========
    >>> TR2(tan(x))
    sin(x)/cos(x)
    >>> TR2(cot(x))
    cos(x)/sin(x)

    """
    rv = rv.replace(tan, lambda x: sin(x)/cos(x))
    rv = rv.replace(cot, lambda x: cos(x)/sin(x))
    return rv

#-------Transformation rule 3 - induced formula : example sin(-a) = -sin(a)

def TR3(rv):
    """
    Examples
    ========
    >>> TR3(cos(y-x*(y-x)))
    cos(x*(x - y) + y)
    >>> cos(pi/2+x)
    -sin(x)
    >>> cos(30*pi/2+x)
    -cos(x)
    """
    from sympy.simplify.simplify import signsimp

    # Negative argument (already automatic for funcs like sin(-x) -> -sin(x)
    # but more complicated expressions can use it, too
    return rv.replace(
        lambda x: isinstance(x, C.TrigonometricFunction),
        lambda x: x.func(signsimp(x.args[0])))
    #The following are automatically handled
    #Argument of type: pi/2 +/- angle
    #Argument of type: pi +/- angle
    #Argument of type : 2k*pi +/- angle

#-------Transformation rule 4 - values of special angles

def TR4(rv):
    """
    Examples
    ========
    >>> for s in (0, pi/6, pi/4, pi/3, pi/2):
    ...    for f in (cos, sin, tan, cot):
    ...      print f(s),
    ...    print
    ...
    1 0 0 zoo
    sqrt(3)/2 1/2 sqrt(3)/3 sqrt(3)
    sqrt(2)/2 sqrt(2)/2 1 1
    1/2 sqrt(3)/2 sqrt(3) sqrt(3)/3
    0 1 zoo 0
    """
    # special values at 0, pi/6, pi/4, pi/3, pi/2 already handled
    return rv

#-------Transformation rule 5 - substitution of sin**2 or sin**4

def TR5(rv):
    """
    Examples
    ========
    >>> TR5(sin(x)**2)
    -cos(x)**2 + 1
    >>> TR5(sin(x)**4)
    (-cos(x)**2 + 1)**2
    >>> TR5(sin(x)**6)
    sin(x)**6
    """
    rv = rv.replace(
        lambda x: x.is_Pow and x.base.func == sin and x.exp in (2, 4),
        lambda x: (1 - cos(x.base.args[0])**2)**(x.exp//2))
    return rv

#-------Transformation rule 6 - substitution of cos**2 or cos**4

def TR6(rv):
    """
    Examples
    ========
    >>> TR6(cos(x)**4)
    (-sin(x)**2 + 1)**2
    >>> TR6(cos(x)**2)
    -sin(x)**2 + 1
    >>> TR6(cos(x)**6)
    cos(x)**6
    """
    rv = rv.replace(
        lambda x: x.is_Pow and x.base.func == cos and x.exp in (2, 4),
        lambda x: (1 - sin(x.base.args[0])**2)**(x.exp//2))
    return rv

#-------Transformation rule 7 - Lowering the degree of cos square:

def TR7(rv):
    """
    Examples
    ========
    >>> TR7(cos(x)**2)
    cos(2*x)/2 + 1/2
    >>> TR7(cos(x)**2 + 1)
    cos(2*x)/2 + 3/2

    """
    rv = rv.replace(
        lambda x: x.is_Pow and x.base.func == cos and x.exp == 2,
        lambda x: (1 + cos(2*x.base.args[0]))/2)
    return rv

#-------Transformation rule 8 - Converting cos*cos, cos*sin or sin*sin to
#-------a sum or difference of cos and or sin terms

def TR8(rv):
    """
    Examples
    ========
    >>> TR8(cos(2)*cos(3))
    cos(5)/2 + cos(1)/2
    >>> TR8(cos(2)*sin(3))
    sin(5)/2 - sin(1)/2
    >>> TR8(sin(2)*sin(3))
    -cos(5)/2 + cos(1)/2
    """
    if not rv.is_Mul:
        return rv.replace(lambda x: x.is_Mul, lambda x: TR8(x))
    args = {cos: [], sin: [], None: []}
    for a in Mul.make_args(rv):
        if a.func in (cos, sin):
            args[a.func].append(a.args[0])
        else:
            args[None].append(a)
    c = args[cos]  # do these need to be ordered?
    s = args[sin]
    if not (c and s or len(c) > 1 or len(s) > 1):
        return rv
    args = args[None]
    n = min(len(c), len(s))
    for i in range(n):
        a1 = c.pop()
        a2 = s.pop()
        args.append((sin(a1 + a2) + sin(a1 - a2))/2)
    while len(c) > 1:
        a1 = c.pop()
        a2 = c.pop()
        args.append((cos(a1 + a2) + cos(a1 - a2))/2)
    if c:
        args.append(cos(c.pop()))
    while len(s) > 1:
        a1 = s.pop()
        a2 = s.pop()
        args.append((-cos(a1 + a2) + cos(a1 - a2))/2)
    if s:
        args.append(sin(s.pop()))
    return Mul(*args)

#-------Transformation rule 9 - sum of cos or sin

def TR9(rv):
    """
    Examples
    ========
    >>> TR9(cos(1)+cos(2))
    2*cos(1/2)*cos(3/2)
    >>> TR9(cos(1)-cos(2))
    2*sin(1/2)*sin(3/2)
    >>> TR9(sin(1)-sin(2))
    -2*sin(1/2)*cos(3/2)
    >>> TR9(sin(1)+sin(2))
    2*sin(3/2)*cos(1/2)
    >>> TR9(cos(1) + 2*sin(1) + 2*sin(2))
    cos(1) + 4*sin(3/2)*cos(1/2)
    """
    if not rv.is_Add:
        return rv.replace(lambda x: x.is_Add, lambda x: TR9(x))
    # combine two sines or two cosines from args that
    # have the same coefficients
    a, b = map(Wild, 'ab')
    pat = {
        (cos, cos, 1): 2*cos((a + b)/2)*cos((a - b)/2),
        (cos, cos, -1): -2*sin((a + b)/2)*sin((a - b)/2),
        (sin, sin, 1): 2*sin((a + b)/2)*cos((a - b)/2),
        (sin, sin, -1): 2*cos((a + b)/2)*sin((a - b)/2)}
    from sympy.utilities.iterables import combinations
    for c in combinations(rv.args, 2):
        m = factor_terms(Add(*c))
        for mi in Mul.make_args(m):
            if not (mi.is_Add and len(mi.args) == 2):
                continue
            r = mi.match(cos(a) + cos(b))
            if r:
                return TR9(rv - Add(*c) + m/mi*pat[(cos, cos, 1)].subs(r))
            r = mi.match(cos(a) - cos(b))
            if r:
                return TR9(rv - Add(*c) + m/mi*pat[(cos, cos, -1)].subs(r))
            r = mi.match(sin(a) + sin(b))
            if r:
                return TR9(rv - Add(*c) + m/mi*pat[(sin, sin, 1)].subs(r))
            r = mi.match(sin(a) - sin(b))
            if r:
                return TR9(rv - Add(*c) + m/mi*pat[(sin, sin, -1)].subs(r))
    return rv

#-------Transformation rule 10 - sin or cos of sum

def TR10(rv):
    """
    Examples
    ========
    >>> TR10(cos(a+b))
    -sin(a)*sin(b) + cos(a)*cos(b)
    >>> TR10(sin(a+b))
    sin(a)*cos(b) + sin(b)*cos(a)
    >>> TR10(sin(a+b+c))
    (-sin(a)*sin(b) + cos(a)*cos(b))*sin(c) + (sin(a)*cos(b) + sin(b)*cos(a))*cos(c)
    """
    if rv.func not in (cos, sin):
        return rv.replace(lambda x: x.func in (cos, sin), lambda x: TR10(x))
    f = rv.func
    arg = rv.args[0]  # should expand_mul be used?
    if arg.is_Add:
        a, b = arg.as_two_terms()
        if b.is_Add:
            if f == sin:
                return sin(a)*TR10(cos(b)) + cos(a)*TR10(sin(b))
            else:
                return cos(a)*TR10(cos(b)) - sin(a)*TR10(sin(b))
        else:
            if f == sin:
                return sin(a)*cos(b) + cos(a)*sin(b)
            else:
                return cos(a)*cos(b) - sin(a)*sin(b)
    return rv

#-------Transformation rule 10^-1 - Inverse of Sum or difference of angles:

def TR10i(rv, first=True):
    """
    Examples
    ========
    >>> TR10i(cos(1)*cos(3) + sin(1)*sin(3))
    cos(2)
    >>> TR10i(cos(1)*cos(3) - sin(1)*sin(3))
    cos(4)
    >>> TR10i(cos(1)*sin(3) - sin(1)*cos(3))
    sin(2)
    >>> TR10i(cos(1)*sin(3) + sin(1)*cos(3))
    sin(4)
    >>> TR10i(cos(1)*sin(3) + sin(1)*cos(3) + 7)
    sin(4) + 7
    >>> TR10i(cos(1)*sin(3) + sin(1)*cos(3) + cos(3))
    cos(3) + sin(4)
    >>> TR10i(2*cos(1)*sin(3) + 2*sin(1)*cos(3)+cos(3))
    2*sin(4) + cos(3)
    """
    if not rv.is_Add:
        return rv.replace(lambda x: x.is_Add, lambda x: TR9(x))

    a1, a2 = map(Wild, 'ab')
    from sympy.utilities.iterables import combinations
    for c in combinations(rv.args, 2):
        m = factor_terms(Add(*c))
        for mi in Mul.make_args(m):
            if not (mi.is_Add and len(mi.args) == 2):
                continue
            hit = False
            for i in range(2):
                r = mi.match(sin(a1)*cos(a2) + cos(a1)*sin(a2))
                if r:
                    a, b = r[a1], r[a2]
                    new = sin(a + b)
                    hit = True
                    break
                r = mi.match(sin(a1)*cos(a2) - cos(a1)*sin(a2))
                if r:
                    a, b = r[a1], r[a2]
                    new = sin(a - b)
                    hit = True
                    break
                r = mi.match(cos(a1)*cos(a2) + sin(a1)*sin(a2))
                if r:
                    a, b = r[a1], r[a2]
                    new = cos(a - b)
                    hit = True
                    break
                r = mi.match(cos(a1)*cos(a2) - sin(a1)*sin(a2))
                if r:
                    a, b = r[a1], r[a2]
                    new = cos(a + b)
                    hit = True
                    break
                assert (-(a1 + a2)).is_Add
                mi = -mi
            if hit:
                if i == 1:
                    mi = -mi
                    new = -new
                rv = rv - Add(*c) + m/mi*new
    return rv

#-------Transformation rule 11 - sin or cos of double angle:

def TR11(rv):
    """
    Examples
    ========
    >>> TR11(sin(2))
    2*sin(1)*cos(1)
    >>> TR11(sin(4))
    4*(-sin(1)**2 + cos(1)**2)*sin(1)*cos(1)

    >>> TR11(cos(2))
    -sin(1)**2 + cos(1)**2
    >>> TR11(cos(4))
    -4*sin(1)**2*cos(1)**2 + (-sin(1)**2 + cos(1)**2)**2
    """
    if rv.func not in (cos, sin):
        return rv.replace(
            lambda x: x.func in (sin, cos),
            lambda x: TR11(x))
    c, m = rv.args[0].as_coeff_Mul()
    if c % 2 == 0:
        arg = c//2*m
        c = TR11(cos(arg))
        s = TR11(sin(arg))
        if rv.func == sin:
            rv = 2*s*c
        else:
            rv = c**2 - s**2
    return rv

#-------Transformation rule 12 ---- tan of sum or difference:

def TR12(rv):
    """
    Examples
    ========
    >>> from sympy.simplify.fu import TR12
    >>> TR12(tan(x + y))
    (tan(x) + tan(y))/(-tan(x)*tan(y) + 1)
    """
    if rv.func != tan:
        return rv.replace(lambda x: x.func == tan, lambda x: TR12(x))
    arg = rv.args[0]  # should expand_mul be used?
    if arg.is_Add:
        a, b = arg.as_two_terms()
        if b.is_Add:
            tb = TR12(tan(b))
        else:
            tb = tan(b)
        return (tan(a) + tb)/(1 - tan(a)*tb)
    return rv

#-------Transformation rule 13 ---- Product of tan or cot:

def TR13(rv):
    """
    Examples
    ========
    >>> TR13(tan(3)*tan(2))
    -(tan(2) + tan(3))*cot(5) + 1
    >>> TR13(cot(3)*cot(2))
    1 + (cot(3) + cot(2))*cot(5)
    """
    if not rv.is_Mul:
        return rv.replace(lambda x: x.is_Mul, lambda x: TR13(x))
    args = {tan: [], cot: [], None: []}
    for a in Mul.make_args(rv):
        if a.func in (tan, cot):
            args[a.func].append(a.args[0])
        else:
            args[None].append(a)
    t = args[tan]  # do these need to be ordered?
    c = args[cot]
    if len(t) < 2 and len(c) < 2:
        return rv
    args = args[None]
    while len(t) > 1:
        t1 = t.pop()
        t2 = t.pop()
        args.append(1 - (tan(t1) + tan(t2))*cot(t1 + t2))
    if t:
        args.append(tan(t.pop()))
    while len(c) > 1:
        t1 = c.pop()
        t2 = c.pop()
        args.append(1 + (cot(t1) + cot(t2))*cot(t1 + t2))
    if c:
        args.append(cot(t.pop()))
    return Mul(*args)

#---------Length of trigonometric expression (count occurances of 'sin' , 'cos', etc.)

def L(rv):
    """
    Examples
    ========
    >>> L(cos(x)+sin(x))
    2
    """
    return rv.count(C.TrigonometricFunction)

#--------Combination transformation rule 1 : ---------------------

def CTR1(rv):
    rv1 = TR5(rv)
    rv1 = TR0(rv1)
    rv2 = TR6(rv)
    rv2 = TR0(rv2)
    if( (L(rv1) < L(rv)) and (L(rv1) <= L(rv2)) ):
       return rv1
    elif( (L(rv2) < L(rv)) and (L(rv2) <= L(rv1)) ):
       return rv2
    else:
        return rv

#--------Combination transformation rule 2: ---------------------

def CTR2(rv):
    rv1 = TR11(rv)
    rv1 = TR5(rv1)
    rv2 = TR11(rv)
    rv2 = TR6(rv2)
    rv3 = TR11(rv)
    rv1 = TR0(rv1)
    rv2 = TR0(rv2)
    rv3 = TR0(rv3)
    if( (L(rv1) < L(rv3)) and (L(rv1) <= L(rv2)) ):
        return rv1
    elif( (L(rv2) < L(rv3)) and (L(rv2) <= L(rv1)) ):
        return rv2
    else:
        return rv3

#--------Combination transformation rule 3: ---------------------

def CTR3(rv):
    rv1 = TR8(rv)
    rv1 = TR0(rv1)
    rv2 = TR8(rv)
    rv2 = TR10i(rv2)
    rv2 = TR0(rv2)
    if( (L(rv2) < L(rv)) ):
        return rv2
    elif( (L(rv1) < L(rv)) ):
        return rv1
    else:
        return rv

#--------Combination transformation rule 4: ---------------------

def CTR4(rv):
    rv1 = TR4(rv)
    rv1 = TR0(rv1)
    if( (L(rv1) < L(rv)) ):
        return rv1
    else:
        return rv

#-------------Rule list 1-------------------------

def RL1(rv):
    rv = TR4(rv)
    rv = TR3(rv)
    rv = TR4(rv)
    rv = TR12(rv)
    rv = TR4(rv)
    rv = TR13(rv)
    rv = TR4(rv)
    rv = TR0(rv)
    return rv

#-------------Rule list 2-------------------------

def RL2(rv):
    rv = TR4(rv)
    rv = TR3(rv)
    rv = TR10(rv)
    rv = TR4(rv)
    rv = TR3(rv)
    rv = TR11(rv)
    rv = TR5(rv)
    rv = TR7(rv)
    rv = TR11(rv)
    rv = TR4(rv)
    rv = CTR3(rv)
    rv = TR0(rv)
    rv = CTR1(rv)
    rv = TR9(rv)
    rv = CTR2(rv)
    rv = TR4(rv)
    rv = TR9(rv)
    rv = TR0(rv)
    rv = TR9(rv)
    rv = CTR4(rv)
    return rv


def fu(rv):
    """reduces expression by using transformation rules given in the Fu et al algorithm

    Examples
    ========
    >>> from sympy.simplify.simplify import fu
    >>> a = sin(50)**2 + cos(50)**2 + sin(pi/6)
    >>> fu(a)
    3/2
    >>> a = sin(100)**4 - cos(50)**2 + sin(50)**2 + 2*cos(100)**2
    >>> fu(a)
    2 - 2*cos(50)**2 + cos(100)**4
    """
    rv = sympify(rv)
    rv = TR0(rv)
    rv = TR1(rv)
    if rv.has(tan, cot):
        rv1 = RL1(rv)
        if(L(rv1) < L(rv)):
            rv = rv1

    rv = TR0(rv)
    rv = TR2(rv)
    rv = TR0(rv)
    if rv.has(sin, cos):
        rv = RL2(rv)
    return rv
