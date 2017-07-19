'''
Utility functions for Rubi integration.

See: http://www.apmaths.uwo.ca/~arich/IntegrationRules/PortableDocumentFiles/Integration%20utility%20functions.pdf
'''

from sympy.functions.elementary.integers import floor, frac
from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf)
from sympy.functions.elementary.hyperbolic import acosh, asinh, atanh, acoth, acsch, acsch, cosh, sinh, tanh, coth, sech, csch
from sympy.functions.elementary.trigonometric import atan, acsc, asin, acot, acos, asec
from sympy.polys.polytools import degree, Poly, quo, rem
from sympy.simplify.simplify import fraction, simplify, cancel, factor, nthroot
from sympy.core.expr import UnevaluatedExpr
from sympy.core.sympify import sympify
from sympy.utilities.iterables import postorder_traversal
from sympy.core.expr import UnevaluatedExpr
from sympy.functions.elementary.complexes import im, re, Abs
from sympy.core.exprtools import factor_terms
from sympy import (exp, polylog, N, Wild, factor, gcd, Sum, S, I, Mul, Add, hyper,
    Symbol, symbols, sqf_list, sqf, Max, gcd, hyperexpand, trigsimp, factorint,
    Min, Max, sign, E, expand_trig, poly, apart, lcm, And)
from mpmath import ellippi, ellipe, ellipf, appellf1
from sympy.polys.polytools import poly_from_expr
from sympy.utilities.iterables import flatten
from sympy.strategies import distribute

def Int(expr, var, showsteps=False):
    from .rubi import rubi_integrate
    if expr == None:
        return None
    return rubi_integrate(expr, var, showsteps)

def Set(expr, value):
    return {expr: value}

def With(subs, expr):
    if isinstance(subs, dict):
        k = list(subs.keys())[0]
        expr = expr.subs(k, subs[k])
    else:
        for i in subs:
            k = list(i.keys())[0]
            expr = expr.subs(k, i[k])
    return expr

def Scan(f, expr):
    # evaluates f applied to each element of expr in turn.
    for i in expr:
        yield f(i)

def MapAnd(f, l, x=None):
    # MapAnd[f,l] applies f to the elements of list l until False is returned; else returns True
    if x:
        for i in l:
            if f(i, x) == False:
                return False
        return True
    else:
        for i in l:
            if f(i) == False:
                return False
        return True

def FalseQ(u):
    return u == False

def ZeroQ(expr):
    if isinstance(expr, list):
        return any(ZeroQ(i) for i in expr)
    else:
        return expr == 0

def NegativeQ(u):
    res = u < 0
    if not res.is_Relational:
        return res
    return False

def NonzeroQ(expr):
    return expr != 0

def FreeQ(nodes, var):
    if isinstance(nodes, list):
        return not any(expr.has(var) for expr in nodes)
    else:
        return not nodes.has(var)

def List(*var):
    return list(var)

def Log(e):
    return log(e)

def PositiveQ(var):
    res = var > 0
    if not res.is_Relational:
        return res
    return False

def PositiveIntegerQ(*args):
    return all(var.is_Integer and PositiveQ(var) for var in args)

def NegativeIntegerQ(*args):
    return all(var.is_Integer and NegativeQ(var) for var in args)

def IntegerQ(var):
    if isinstance(var, int):
        return True
    else:
        return var.is_Integer

def IntegersQ(*var):
    return all(IntegerQ(i) for i in var)

def ComplexNumberQ(*var):
    return all((im(i)!=0) for i in var)

def RealNumericQ(u):
    return u.is_real

def PositiveOrZeroQ(u):
    return u.is_real and u >= 0

def NegativeOrZeroQ(u):
    return u.is_real and u <= 0

def FractionOrNegativeQ(u):
    return FractionQ(u) or NegativeQ(u)

def NegQ(var):
    return NegativeQ(var)

def Equal(a, b):
    return a == b

def Unequal(a, b):
    return a != b

def IntPart(u):
    # IntPart[u] returns the sum of the integer terms of u.
    if ProductQ(u):
        if IntegerQ(First(u)):
            return First(u)*IntPart(Rest(u))
    elif IntegerQ(u):
        return u
    elif FractionQ(u):
        return IntegerPart(u)
    elif SumQ(u):
        res = 0
        for i in u.args:
            res += IntPart(i)
        return res
    else:
        return 0

def FracPart(u):
    # FracPart[u] returns the sum of the non-integer terms of u.
    if ProductQ(u):
        if IntegerQ(First(u)):
            return First(u)*FracPart(Rest(u))

    if IntegerQ(u):
        return 0
    elif FractionQ(u):
        return FractionalPart(u)
    elif SumQ(u):
        res = 0
        for i in u.args:
            res += FracPart(i)
        return res
    else:
        return u

def RationalQ(*nodes):
    return all(var.is_Rational for var in nodes)

def ProductQ(expr):
    return expr.is_Mul

def SumQ(expr):
    return expr.is_Add

def NonsumQ(expr):
    return not SumQ(expr)

def Subst(a, x, y):
    if None in [a, x, y]:
        return None
    return a.subs(x, y)

def First(expr, d=None):
    # gives the first element if it exists, or d otherwise.
    if isinstance(expr, list):
        return expr[0]
    else:
        if SumQ(expr) or ProductQ(expr):
            l = Sort(expr.args)
            return l[0]
        else:
            return expr.args[0]

def Rest(expr):
    if isinstance(expr, list):
        return expr[1:]
    else:
        if SumQ(expr) or ProductQ(expr):
            l = Sort(expr.args)
            return expr.func(*l[1:])
        else:
            return expr.args[1]

def SqrtNumberQ(expr):
    # SqrtNumberQ[u] returns True if u^2 is a rational number; else it returns False.
    if expr.is_Pow:
        m = expr.base
        n = expr.exp
        return (IntegerQ(n) and SqrtNumberQ(m)) or (IntegerQ(n-S(1)/2) and RationalQ(m))
    elif expr.is_Mul:
        return all(SqrtNumberQ(i) for i in expr.args)
    else:
        return RationalQ(expr) or expr == I

def SqrtNumberSumQ(u):
    return SumQ(u) & SqrtNumberQ(First(u)) & SqrtNumberQ(Rest(u)) | ProductQ(u) & SqrtNumberQ(First(u)) & SqrtNumberSumQ(Rest(u))

def LinearQ(expr, x):
    if expr.is_polynomial(x):
        if degree(Poly(expr, x), gen=x) == 1:
            return True
    return False

def Sqrt(a):
    return sqrt(a)

def ArcCosh(a):
    return acosh(a)

def Coefficient(expr, var, n=1):
    a = Poly(expr, var)
    if (degree(a) - n) < 0:
        return 0
    else:
        return a.all_coeffs()[degree(a) - n]

def Denominator(var):
    return fraction(var)[1]

def Hypergeometric2F1(a, b, c, z):
    return hyper([a, b], [c], z)

def ArcTan(a):
    return atan(a)

def Not(var):
    if isinstance(var, bool):
        return not var
    elif var.is_Relational:
        var = False
    return not var

def Simplify(expr):
    return simplify(expr)

def FractionalPart(a):
    return frac(a)

def IntegerPart(a):
    return floor(a)

def AppellF1(a, b1, b2, c, x, y):
    return appellf1(a, b1, b2, c, x, y)

def EllipticPi(*args):
    return ellippi(*args)

def EllipticE(*args):
    return ellipe(*args)

def EllipticF(Phi, m):
    return ellipf(Phi, m)

def ArcTan(a):
    return atan(a)

def ArcTanh(a):
    return atanh(a)

def ArcSin(a):
    return asin(a)

def ArcSinh(a):
    return asinh(a)

def ArcCos(a):
    return acos(a)

def ArcCsc(a):
    return acsc(a)

def ArcCsch(a):
    return acsch(a)

def Sinh(u):
    return sinh(u)

def Tanh(u):
    return tanh(u)

def Cosh(u):
    return cosh(u)

def Sech(u):
    return sech(u)

def Csch(u):
    return csch(u)

def Coth(u):
    return coth(u)

def LessEqual(*args):
    for i in range(0, len(args) - 1):
        try:
            if args[i] > args[i + 1]:
                return False
        except:
            return False
    return True

def Less(*args):
    for i in range(0, len(args) - 1):
        try:
            if args[i] >= args[i + 1]:
                return False
        except:
            return False
    return True

def Greater(*args):
    for i in range(0, len(args) - 1):
        try:
            if args[i] <= args[i + 1]:
                return False
        except:
            return False
    return True

def GreaterEqual(*args):
    for i in range(0, len(args) - 1):
        try:
            if args[i] < args[i + 1]:
                return False
        except:
            return False
    return True

def FractionQ(*args):
    return all(i.is_Rational for i in args)

def IntLinearcQ(a, b, c, d, m, n, x):
    # returns True iff (a+b*x)^m*(c+d*x)^n is integrable wrt x in terms of non-hypergeometric functions.
    return IntegerQ(m) | IntegerQ(n) | IntegersQ(3*m, 3*n) | IntegersQ(4*m, 4*n) | IntegersQ(2*m, 6*n) | IntegersQ(6*m, 2*n) | IntegerQ(m + n)

Defer = UnevaluatedExpr

def Expand(expr):
    return expr.expand()

def IndependentQ(u, x):
    return FreeQ(u, x)

def PowerQ(expr):
    return expr.is_Pow

def IntegerPowerQ(u):
    return PowerQ(u) and IntegerQ(u.args[1])

def PositiveIntegerPowerQ(u):
    return PowerQ(u) and IntegerQ(u.args[1]) and PositiveQ(u.args[1])

def FractionalPowerQ(u):
    return PowerQ(u) & FractionQ(u.args[1])

def AtomQ(expr):
    expr = sympify(expr)
    if isinstance(expr, list):
        return False
    if expr in [None, True, False]: # [None, True, False] are atoms in mathematica
        return True
    elif isinstance(expr, list):
        return all(AtomQ(i) for i in expr)
    else:
        return expr.is_Atom

def ExpQ(u):
    return Head(u) == exp

def LogQ(u):
    return u.func == log

def Head(u):
    return u.func

def MemberQ(l, u):
    if isinstance(l, list):
        return u in l
    else:
        return u in l.args

def TrigQ(u):
    if AtomQ(u):
        x = u
    else:
        x = Head(u)
    return MemberQ([sin, cos, tan, cot, sec, csc], x)

def SinQ(u):
    return Head(u) == sin

def CosQ(u):
    return Head(u) == cos

def TanQ(u):
    return Head(u) == tan

def CotQ(u):
    return Head(u) == cot

def SecQ(u):
    return Head(u) == sec

def CscQ(u):
    return Head(u) == csc

def HyperbolicQ(u):
    if AtomQ(u):
        x = u
    else:
        x = Head(u)
    return MemberQ([sinh, cosh, tanh, coth, sech, csch], x)

def SinhQ(u):
    return Head(u) == sinh

def CoshQ(u):
    return Head(u) == cosh

def TanhQ(u):
    return Head(u) == tanh

def CothQ(u):
    return Head(u) == coth

def SechQ(u):
    return Head(u) == sech

def CschQ(u):
    return Head(u) == csch

def InverseTrigQ(u):
    if AtomQ(u):
        x = u
    else:
        x = Head(u)
    return MemberQ([asin, acos, atan, acot, asec, acsc], x)

def SinCosQ(f):
    return MemberQ([sin, cos, sec, csc], Head(f))

def SinhCoshQ(f):
    return MemberQ([sinh, cosh, sech, csch], Head(f))

def Rt(val, n):
    if ComplexNumberQ(nthroot(val, n)):
        return nthroot(val, n)
    if FractionalPart(float(nthroot(val, n))) == 0:
        return sympify(int(nthroot(val, n)))
    else:
        return sympify(float(nthroot(val, n)))

def LeafCount(expr):
    return len(list(postorder_traversal(expr)))

def Numerator(u):
    return fraction(u)[0]

def NumberQ(u):
    return u.is_number

def NumericQ(u):
    return N(u).is_number

def Length(expr):
    # returns number of elements in the experssion
    if isinstance(expr, list):
        return len(expr)
    return len(expr.args)

def ListQ(u):
    return isinstance(u, list)

def Im(u):
    return im(u)

def Re(u):
    return re(u)

def InverseHyperbolicQ(u):
    if not u.is_Atom:
        u = Head(u)
    return u in [acosh, asinh, atanh, acoth, acsch, acsch]

def InverseFunctionQ(u):
    # returns True if u is a call on an inverse function; else returns False.
    return LogQ(u) or InverseTrigQ(u) and Length(u) <= 1 or InverseHyperbolicQ(u) or u.func == polylog

def TrigHyperbolicFreeQ(u, x):
    # If u is free of trig, hyperbolic and calculus functions involving x, TrigHyperbolicFreeQ[u,x] returns true; else it returns False.
    if AtomQ(u):
        return True
    else:
        if TrigQ(u) | HyperbolicQ(u) | CalculusQ(u):
            return FreeQ(u, x)
        else:
            for i in u.args:
                if not TrigHyperbolicFreeQ(i, x):
                    return False
            return True

def InverseFunctionFreeQ(u, x):
    # If u is free of inverse, calculus and hypergeometric functions involving x, InverseFunctionFreeQ[u,x] returns true; else it returns False.
    if AtomQ(u):
        return True
    else:
        if InverseFunctionQ(u) | CalculusQ(u) | u.func == hyper | u.func == appellf1:
            return FreeQ(u, x)
        else:
            for i in u.args:
                if not ElementaryFunctionQ(i):
                    return False
            return True

def RealQ(u):
    if ListQ(u):
        return MapAnd(RealQ, u)
    elif NumericQ(u):
        return ZeroQ(Im(N(u)))
    elif u.is_Pow:
        u = u.base
        v = u.exp
        return RealQ(u) & RealQ(v) & (IntegerQ(v) | PositiveOrZeroQ(u))
    elif u.is_Mul:
        return all(RealQ(i) for i in u.args)
    elif u.is_Add:
        return all(RealQ(i) for i in u.args)
    elif u.is_Function:
        f = u.func
        u = u.args[0]
        if f in [sin, cos, tan, cot, sec, csc, atan, acot, erf]:
            return RealQ(u)
        else:
            if f in [asin, acos]:
                return LE(-1, u, 1)
            else:
                if f == log:
                    return PositiveOrZeroQ(u)
                else:
                    return False
    else:
        return False

def EqQ(u, v):
    return ZeroQ(u - v)

def FractionalPowerFreeQ(u):
    if AtomQ(u):
        return True
    elif FractionalPowerQ(u):
        return False

def ComplexFreeQ(u):
    if AtomQ(u) and Not(ComplexNumberQ(u)):
        return True
    else:
         return False

def PolynomialQ(u, x):
    return u.is_polynomial(x)

def FactorSquareFree(u):
    return sqf(u)

def PowerOfLinearQ(expr, x):
    u = Wild('u')
    w = Wild('w')
    m = Wild('m', exclude=[x])
    n = Wild('n', exclude=[x])
    Match = expr.match(u**m)
    if PolynomialQ(Match[u], x):
        if IntegerQ(Match[m]):
            e = FactorSquareFree(Match[u]).match(w**n)
            return LinearQ(e[w], x)
        else:
            return LinearQ(Match[u], x)

def Exponent(expr, x, *k):
    expr = Expand(S(expr))
    if not k:
        if S(expr).is_number or (not expr.has(x)):
            return 0
        if expr.is_polynomial(x):
            return degree(expr, gen = x)
        else:
            return 0
    else:
        if S(expr).is_number or (not expr.has(x)):
            return [0]
        if expr.is_Add:
            lst = []
            k = 1
            for t in expr.args:
                if t.has(x):
                    lst += [degree(t, gen = x)]
                else:
                    if k == 1:
                        lst += [0]
                        k += 1
            lst.sort()
            return lst
        else:
            return [degree(expr, gen = x)]

def QuadraticQ(u, x):
    # QuadraticQ(u, x) returns True iff u is a polynomial of degree 2 and not a monomial of the form a x^2
    if ListQ(u):
        for expr in u:
            if Not(PolyQ(expr, x, 2) and Not(Coefficient(expr, x, 0) == 0 and Coefficient(expr, x, 1) == 0)):
                return False
        return True
    else:
        return PolyQ(u, x, 2) and Not(Coefficient(u, x, 0) == 0 and Coefficient(u, x, 1) == 0)

def LinearPairQ(u, v, x):
    # LinearPairQ(u, v, x) returns True iff u and v are linear not equal x but u/v is a constant wrt x
    return LinearQ(u, x) and LinearQ(v, x) and NonzeroQ(u-x) and ZeroQ(Coefficient(u, x, 0)*Coefficient(v, x, 1)-Coefficient(u, x, 1)*Coefficient(v, x, 0))

def BinomialParts(u, x):
    if PolynomialQ(u, x):
        if Exponent(u, x) > 0:
            lst = Exponent(u, x, List)
            if len(lst)==1:
                return [0, Coefficient(u, x, Exponent(u, x)), Exponent(u,x)]
            elif len(lst) == 2 and lst[0] == 0:
                return [Coefficient(u, x, 0), Coefficient(u, x, Exponent(u, x)), Exponent(u, x)]
            else:
                return False
        else:
            return False
    elif PowerQ(u):
        if u.args[0] == x and FreeQ(u.args[1], x):
            return [0, 1, u.args[1]]
        else:
            return False
    elif ProductQ(u):
        if FreeQ(First(u), x):
            lst2 = BinomialParts(Rest(u), x)
            if AtomQ(lst2):
                return False
            else:
                return [First(u)*lst2[0], First(u)*lst2[1], lst2[2]]
        elif FreeQ(Rest(u), x):
            lst1 = BinomialParts(First(u), x)
            if AtomQ(lst1):
                return False
            else:
                return [Rest(u)*lst1[0], Rest(u)*lst1[1], lst1[2]]
        lst1 = BinomialParts(First(u), x)
        if AtomQ(lst1):
            return False
        lst2 = BinomialParts(Rest(u), x)
        if AtomQ(lst2):
            return False
        a = lst1[0]
        b = lst1[1]
        m = lst1[2]
        c = lst2[0]
        d = lst2[1]
        n = lst2[2]
        if ZeroQ(a):
            if ZeroQ(c):
                return [0, b*d, m + n]
            elif ZeroQ(m + n):
                return [b*d, b*c, m]
            else:
                return False
        if ZeroQ(c):
            if ZeroQ(m + n):
                return [b*d, a*d, n]
            else:
                return False
        if EqQ(m, n) and ZeroQ(a*d + b*c):
            return [a*c, b*d, 2*m]
        else:
            return False
    elif SumQ(u):
        if FreeQ(First(u),x):
            lst2 = BinomialParts(Rest(u), x)
            if AtomQ(lst2):
                return False
            else:
                return [First(u) + lst2[0], lst2[1], lst2[2]]
        elif FreeQ(Rest(u), x):
            lst1 = BinomialParts(First(u), x)
            if AtomQ(lst1):
                return False
            else:
                return[Rest(u) + lst1[0], lst1[1], lst1[2]]
        lst1 = BinomialParts(First(u), x)
        if AtomQ(lst1):
            return False
        lst2 = BinomialParts(Rest(u),x)
        if AtomQ(lst2):
            return False
        if EqQ(lst1[2], lst2[2]):
            return [lst1[0] + lst2[0], lst1[1] + lst2[1], lst1[2]]
        else:
            return False
    else:
        return False

def TrinomialParts(u, x):
    # If u is equivalent to a trinomial of the form a + b*x^n + c*x^(2*n) where n!=0, b!=0 and c!=0, TrinomialParts[u,x] returns the list {a,b,c,n}; else it returns False.
    u = sympify(u)
    if PolynomialQ(u, x):
        lst = CoefficientList(u, x)
        if len(lst)<3 or EvenQ(sympify(len(lst))) or ZeroQ((len(lst)+1)/2):
            return False
        #Catch(
         #   Scan(Function(if ZeroQ(lst), Null, Throw(False), Drop(Drop(Drop(lst, [(len(lst)+1)/2]), 1), -1];
          #  [First(lst), lst[(len(lst)+1)/2], Last(lst), (len(lst)-1)/2]):
    if PowerQ(u):
        if EqQ(u.args[1], 2):
            lst = BinomialParts(u.args[0], x)
            if not lst or ZeroQ(lst[0]):
                return False
            else:
                return [lst[0]**2, 2*lst[0]*lst[1], lst[1]**2, lst[2]]
        else:
            return False
    if ProductQ(u):
        if FreeQ(First(u), x):
            lst2 = TrinomialParts(Rest(u), x)
            if not lst2:
                return False
            else:
                return [First(u)*lst2[0], First(u)*lst2[1], First(u)*lst2[2], lst2[3]]
        if FreeQ(Rest(u), x):
            lst1 = TrinomialParts(First(u), x)
            if not lst1:
                return False
            else:
                return [Rest(u)*lst1[0], Rest(u)*lst1[1], Rest(u)*lst1[2], lst1[3]]
        lst1 = BinomialParts(First(u), x)
        if not lst1:
            return False
        lst2 = BinomialParts(Rest(u), x)
        if not lst2:
            return False
        a = lst1[0]
        b = lst1[1]
        m = lst1[2]
        c = lst2[0]
        d = lst2[1]
        n = lst2[2]
        if EqQ(m, n) and NonzeroQ(a*d+b*c):
            return [a*c, a*d + b*c, b*d, m]
        else:
            return False
    if SumQ(u):
        if FreeQ(First(u), x):
            lst2 = TrinomialParts(Rest(u), x)
            if not lst2:
                return False
            else:
                return [First(u)+lst2[0], lst2[1], lst2[2], lst2[3]]
        if FreeQ(Rest(u), x):
            lst1 = TrinomialParts(First(u), x)
            if not lst1:
                return False
            else:
                return [Rest(u)+lst1[0], lst1[1], lst1[2], lst1[3]]
        lst1 = TrinomialParts(First(u), x)
        if not lst1:
            lst3 = BinomialParts(First(u), x)
            if not lst3:
                return False
            lst2 = TrinomialParts(Rest(u), x)
            if not lst2:
                lst4 = BinomialParts(Rest(u), x)
                if not lst4:
                    return False
                if EqQ(lst3[2], 2*lst4[2]):
                    return [lst3[0]+lst4[0], lst4[1], lst3[1], lst4[2]]
                if EqQ(lst4[2], 2*lst3[2]):
                    return [lst3[0]+lst4[0], lst3[1], lst4[1], lst3[2]]
                else:
                    return False
            if EqQ(lst3[2], lst2[3]) and NonzeroQ(lst3[1]+lst2[1]):
                return [lst3[0]+lst2[0], lst3[1]+lst2[1], lst2[2], lst2[3]]
            if EqQ(lst3[2], 2*lst2[3]) and NonzeroQ(lst3[1]+lst2[2]):
                return [lst3[0]+lst2[0], lst2[1], lst3[1]+lst2[2], lst2[3]]
            else:
                return False
        lst2 = TrinomialParts(Rest(u), x)
        if AtomQ(lst2):
            lst4 = BinomialParts(Rest(u), x)
            if not lst4:
                return False
            if EqQ(lst4[2], lst1[3]) and NonzeroQ(lst1[1]+lst4[0]):
                return [lst1[0]+lst4[0], lst1[1]+lst4[1], lst1[2], lst1[3]]
            if EqQ(lst4[2], 2*lst1[3]) and NonzeroQ(lst1[2]+lst4[1]):
                return [lst1[0]+lst4[0], lst1[1], lst1[2]+lst4[1], lst1[3]]
            else:
                return False
        if EqQ(lst1[3], lst2[3]) and NonzeroQ(lst1[1]+lst2[1]) and NonzeroQ(lst1[2]+lst2[2]):
            return [lst1[0]+lst2[0], lst1[1]+lst2[1], lst1[2]+lst2[2], lst1[3]]
        else:
            return False
    else:
        return False

def PolyQ(u, x, n):
    # returns True iff u is a polynomial of degree n.
    if u.is_polynomial(x):
        return degree(u, gen=x) == n
    return False

def EvenQ(u):
    # gives True if expr is an even integer, and False otherwise.
    return u.is_Integer and u%2 == 0

def OddQ(u):
    # gives True if expr is an odd integer, and False otherwise.
    return u.is_Integer and u%2 == 1

def PerfectSquareQ(u):
    # (* If u is a rational number whose squareroot is rational or if u is of the form u1^n1 u2^n2 ...
    # and n1, n2, ... are even, PerfectSquareQ[u] returns True; else it returns False. *)
    if RationalQ(u):
        return u > 0 and RationalQ(Sqrt(u))
    elif PowerQ(u):
        return EvenQ(u.args[1])
    elif ProductQ(u):
        return PerfectSquareQ(First(u)) & PerfectSquareQ(Rest(u))
    elif SumQ(u):
        s = Simplify(u)
        return NonsumQ(s) & PerfectSquareQ(s)
    else:
        return False

def NiceSqrtAuxQ(u):
    if RationalQ(u):
        return u > 0
    elif PowerQ(u):
        return EvenQ(u.args[1])
    elif ProductQ(u):
        return NiceSqrtAuxQ(First(u)) & NiceSqrtAuxQ(Rest(u))
    elif SumQ(u):
        s = Simplify(u)
        return  NonsumQ(s) & NiceSqrtAuxQ(s)
    else:
        return False

def NiceSqrtQ(u):
    return Not(NegativeQ(u)) & NiceSqrtAuxQ(u)

def Together(u):
    return factor(u)

def FixSimplify(u):
    return u

def TogetherSimplify(u):
    return simplify(u)
    #return With(Set(v, Together(Simplify(Together(u)))), FixSimplify(v))

def PosAux(u):
    if RationalQ(u):
        return u>0
    elif NumberQ(u):
        if ZeroQ(Re(u)):
            return Im(u) > 0
        else:
            return Re(u) > 0
    elif NumericQ(u):
        v = N(u)
        if ZeroQ(Re(v)):
            return Im(v) > 0
        else:
            return Re(v) > 0
    elif PowerQ(u) and OddQ(u.args[1]):
        return PosAux(u.args[0])
    elif ProductQ(u):
        if PosAux(First(u)):
            return PosAux(Rest(u))
        else:
            return not PosAux(Rest(u))
    elif SumQ(u):
        return PosAux(First(u))
    else:
        return True

def PosQ(u):
    # If u is not 0 and has a positive form, PosQ[u] returns True, else it returns False.
    return PosAux(TogetherSimplify(u))

def CoefficientList(u, x):
    if PolynomialQ(u, x):
        return list(reversed(Poly(u, x).all_coeffs()))
    else:
        return []

def ReplaceAll(expr, args):
    return expr.subs(args)

def ExpandLinearProduct(v, u, a, b, x):
    # If u is a polynomial in x, ExpandLinearProduct[v,u,a,b,x] expands v*u into a sum of terms of the form c*v*(a+b*x)^n.
    if FreeQ([a, b], x) and PolynomialQ(u, x):
        lst = CoefficientList(ReplaceAll(u, {x: (x - a)/b}), x)
        lst = [SimplifyTerm(i, x) for i in lst]
        res = 0
        for k in range(1, len(lst)+1):
            res = res + v*lst[k-1]*(a + b*x)**(k - 1)
        return res
    return u*v

def GCD(*args):
    return gcd(*args)

def ContentFactor(expn):
    return factor_terms(expn)

def NumericFactor(u):
    # returns the real numeric factor of u.
    if NumberQ(u):
        if ZeroQ(Im(u)):
            return u
        elif ZeroQ(Re(u)):
            return Im(u)
        else:
            return S(1)
    elif PowerQ(u):
        if RationalQ(u.base) and RationalQ(u.exp):
            if u.exp > 0:
                return 1/Denominator(u.base)
            else:
                return 1/(1/Denominator(u.base))
        else:
            return S(1)
    elif ProductQ(u):
        return Mul(*[NumericFactor(i) for i in u.args])
    elif SumQ(u):
        if LeafCount(u) < 50:
            c = ContentFactor(u)
            if SumQ(c):
                return S(1)
            else:
                return NumericFactor(c)
        else:
            m = NumericFactor(First(u))
            n = NumericFactor(Rest(u))
            if m < 0 and n < 0:
                return -GCD(-m, -n)
            else:
                return GCD(m, n)
    return S(1)

def NonnumericFactors(u):
    if NumberQ(u):
        if ZeroQ(Im(u)):
            return S(1)
        elif ZeroQ(Re(u)):
            return I
        return u
    elif PowerQ(u):
        if RationalQ(u.base) and FractionQ(u.exp):
            return u/NumericFactor(u)
        return u
    elif ProductQ(u):
        result = 1
        for i in u.args:
            result *= NonnumericFactors(i)
        return result
    elif SumQ(u):
        if LeafCount(u) < 50:
            i = ContentFactor(u)
            if SumQ(i):
                return u
            else:
                return NonnumericFactors(i)
        n = NumericFactor(u)
        result = 0
        for i in u.args:
            result += i/n
        return result
    return u

def MakeAssocList(u, x, alst=[]):
    # (* MakeAssocList[u,x,alst] returns an association list of gensymed symbols with the nonatomic
    # parameters of a u that are not integer powers, products or sums. *)
    if AtomQ(u):
        return alst
    elif IntegerPowerQ(u):
        return MakeAssocList(u.args[0], x, alst)
    elif ProductQ(u) or SumQ(u):
        return MakeAssocList(Rest(u), x, MakeAssocList(First(u), x, alst))
    elif FreeQ(u, x):
        tmp = []
        for i in alst:
            if i.args[1] == u:
                tmp.append(i)
                break
        if tmp == []:
            # Append[alst,{Unique["Rubi"],u}],
            alst.append(u)
        return alst
    return alst

def GensymSubst(u, x, alst):
    # (* GensymSubst[u,x,alst] returns u with the kernels in alst free of x replaced by gensymed names. *)
    if AtomQ(u):
        return u
    elif IntegerPowerQ(u):
        return GensymSubst(u.args[0], x, alst)**u.exp
    elif ProductQ(u) or SumQ(u):
        return u.func(*[GensymSubst(i, x, alst) for i in u.args])
    elif FreeQ(u, x):
        tmp = []
        for i in alst:
            if i.args[1] == u:
                tmp.append(i)
                break
        if tmp == []:
            return u
        return tmp[0][0]
    return u

def KernelSubst(u, x, alst):
    # (* KernelSubst[u,x,alst] returns u with the gensymed names in alst replaced by kernels free of x. *)
    if AtomQ(u):
        tmp = []
        for i in alst:
            if i.args[0] == u:
                tmp.append(i)
                break
        if tmp == []:
            return u
        return tmp[0][1]
    elif IntegerPowerQ(u):
        tmp = KernelSubst(u.base, x, alst)
        if u.args[1] < 0 and ZeroQ(tmp):
            return 'Indeterminate'
        return tmp**u.exp
    elif ProductQ(u) or SumQ(u):
        return u.func(*[KernelSubst(i, x, alst) for i in u.args])
    return u

def ExpandExpression(u, x):
    if AlgebraicFunctionQ(u, x) and Not(RationalFunctionQ(u, x)):
        v = ExpandAlgebraicFunction(u, x)
    else:
        v = S(0)
    if SumQ(v):
        return ExpandCleanup(v, x)
    v = SmartApart(u, x)
    if SumQ(v):
        return ExpandCleanup(v, x)
    v = SmartApart(RationalFunctionFactors(u, x), x, x)
    if SumQ(v):
        w = NonrationalFunctionFactors(u, x)
        return ExpandCleanup(v.func(*[i*w for i in v.args]), x)
    v = Expand(u)
    if SumQ(v):
        return ExpandCleanup(v, x)
    v = Expand(u)
    if SumQ(v):
        return ExpandCleanup(v, x)
    return SimplifyTerm(u, x)

def Apart(u, x):
    if RationalFunctionQ(u, x):
        return apart(u, x)
    return u

def SmartApart(*args):
    if len(args) == 2:
        u, x = args
        alst = MakeAssocList(u, x)
        tmp = KernelSubst(Apart(GensymSubst(u, x, alst), x), x, alst)
        if tmp == 'Indeterminate':
            return u
        return tmp

    u, v, x = args
    alst = MakeAssocList(u, x)
    tmp = KernelSubst(Apart(GensymSubst(u, x, alst), x), x, alst)
    if tmp == 'Indeterminate':
        return u
    return tmp

def MatchQ(expr, pattern, *var):
    # returns the matched arguments after matching pattern with expression
    match = expr.match(pattern)
    if match:
        return tuple(match[i] for i in var)
    else:
        return None

def PolynomialQuotientRemainder(p, q, x):
    return [quo(p, q), rem(p, q)]

def FreeFactors(u, x):
    # returns the product of the factors of u free of x.
    if ProductQ(u):
        result = 1
        for i in u.args:
            if FreeQ(i, x):
                result *= i
        return result
    elif FreeQ(u, x):
        return u
    else:
        return 1

def NonfreeFactors(u, x):
    # returns the product of the factors of u not free of x.
    if ProductQ(u):
        result = 1
        for i in u.args:
            if not FreeQ(i, x):
                result *= i
        return result
    elif FreeQ(u, x):
        return 1
    else:
        return u

def RemoveContentAux(expr, x):
    if SumQ(expr):
        w_ = Wild('w')
        p_ = Wild('p')
        n_ = Wild('n', exclude=[x])
        m_ = Wild('m', exclude=[x])
        u_ = Wild('u')
        v_ = Wild('v')
        b_ = Wild('b', exclude=[x])
        a_ = Wild('a', exclude=[x])

        pattern = a_**m_*u_ + b_*v_
        match = expr.match(pattern)
        if match:
            keys = [a_, m_, u_, b_, v_]
            if len(keys) == len(match):
                a, m, u, b, v = tuple([match[i] for i in keys])
                if IntegersQ(a, b) & (a + b == 0) & RationalQ(m):
                    if m > 0:
                        return RemoveContentAux(a**(m - 1)*u - v, x)
                    else:
                        return RemoveContentAux(u - a**(1 - m)*v, x)


        pattern = a_**m_*u_ + a_**n_*v_
        match = expr.match(pattern)
        if match:
            keys = [a_, m_, u_, n_, v_]
            if len(keys) == len(match):
                a, m, u, n, v = tuple([match[i] for i in keys])
                if FreeQ(a, x) & RationalQ(m, n) & (n - m >= 0) & (m != 0):
                    return RemoveContentAux(u + a**(n - m)*v, x)

        pattern = a_**m_*u_ + a_**n_*v_ + a_**p_*w_
        match = expr.match(pattern)
        if match:
            keys = [a_, m_, u_, n_, v_, p_, w_]
            if len(keys) == len(match):
                a, m, u, n, v, p, w = tuple([match[i] for i in keys])
                if RationalQ(m, n, p) & (n - m >= 0) & (p - m >= 0):
                    return RemoveContentAux(u + a**(n - m)*v + a**(p - m)*w, x)

        if NegQ(First(expr)):
            return -expr

    return expr


def RemoveContent(u, x):
    v = NonfreeFactors(u, x)
    w = Together(v)

    if EqQ(FreeFactors(w, x), 1):
        return RemoveContentAux(v, x)
    else:
        return RemoveContentAux(NonfreeFactors(w, x), x)


def FreeTerms(u, x):
    # returns the sum of the terms of u free of x.
    if SumQ(u):
        result = 0
        for i in u.args:
            if FreeQ(i, x):
                result += i
        return result
    elif FreeQ(u, x):
        return u
    else:
        return 0

def NonfreeTerms(u, x):
    # returns the sum of the terms of u free of x.
    if SumQ(u):
        result = S(0)
        for i in u.args:
            if not FreeQ(i, x):
                result += i
        return result
    elif not FreeQ(u, x):
        return u
    else:
        return S(0)

def ExpandAlgebraicFunction(expr, x):
    if ProductQ(expr):
        u_ = Wild('u', exclude=[x])
        n_ = Wild('n', exclude=[x])
        v_ = Wild('v')
        pattern = u_*v_
        match = expr.match(pattern)
        if match:
            keys = [u_, v_]
            if len(keys) == len(match):
                u, v = tuple([match[i] for i in keys])
                if SumQ(v):
                    u, v = v, u
                if not FreeQ(u, x) and SumQ(u):
                    result = 0
                    for i in u.args:
                        result += i*v
                    return result

        pattern = u_**n_*v_
        match = expr.match(pattern)
        if match:
            keys = [u_, n_, v_]
            if len(keys) == len(match):
                u, n, v = tuple([match[i] for i in keys])
                if PositiveIntegerQ(n) and SumQ(u):
                    w = Expand(u**n)
                    result = 0
                    for i in w.args:
                        result += i*v
                    return result

    return expr

def CollectReciprocals(expr, x):
    # Basis: e/(a+b x)+f/(c+d x)==(c e+a f+(d e+b f) x)/(a c+(b c+a d) x+b d x^2)
    if SumQ(expr):
        u_ = Wild('u')
        a_ = Wild('a', exclude=[x])
        b_ = Wild('b', exclude=[x])
        c_ = Wild('c', exclude=[x])
        d_ = Wild('d', exclude=[x])
        e_ = Wild('e', exclude=[x])
        f_ = Wild('f', exclude=[x])
        pattern = u_ + e_/(a_ + b_*x) + f_/(c_+d_*x)
        match = expr.match(pattern)
        if match:
            try: # .match() does not work peoperly always
                keys = [u_, a_, b_, c_, d_, e_, f_]
                u, a, b, c, d, e, f = tuple([match[i] for i in keys])
                if ZeroQ(b*c + a*d) & ZeroQ(d*e + b*f):
                    return CollectReciprocals(u + (c*e + a*f)/(a*c + b*d*x**2),x)
                elif ZeroQ(b*c + a*d) & ZeroQ(c*e + a*f):
                    return CollectReciprocals(u + (d*e + b*f)*x/(a*c + b*d*x**2),x)
            except:
                pass
    return expr

def ExpandCleanup(u, x):
    v = CollectReciprocals(u, x)
    if SumQ(v):
        res = 0
        for i in v.args:
            res += SimplifyTerm(i, x)
        v = res
        if SumQ(v):
            return UnifySum(v, x)
        else:
            return v
    else:
        return v

def AlgebraicFunctionQ(u, x, flag=False):
    if ListQ(u):
        if u == []:
            return True
        elif AlgebraicFunctionQ(First(u), x, flag):
            return AlgebraicFunctionQ(Rest(u), x, flag)
        else:
            return False
    elif AtomQ(u) or FreeQ(u, x):
        return True
    elif PowerQ(u):
        if RationalQ(u.args[1]) | flag & FreeQ(u.args[1], x):
            return AlgebraicFunctionQ(u.args[1], x, flag)
    elif ProductQ(u) | SumQ(u):
        for i in u.args:
            if not AlgebraicFunctionQ(i, x, flag):
                return False
        return True
    return False

def Coeff(expr, form, n=1):
    if n == 1:
        return Coefficient(Together(expr), form, n)
    else:
        coef1 = Coefficient(expr, form, n)
        coef2 = Coefficient(Together(expr), form, n)
        if Simplify(coef1 - coef2) == 0:
            return coef1
        else:
            return coef2

def LeadTerm(u):
    if SumQ(u):
        return First(u)
    return u

def RemainingTerms(u):
    if SumQ(u):
        return Rest(u)
    return u

def LeadFactor(u):
    # returns the leading factor of u.
    if ComplexNumberQ(u) and Re(u) == 0:
        if Im(u) == S(1):
            return u
        else:
            return LeadFactor(Im(u))
    elif ProductQ(u):
            return LeadFactor(First(u))
    return u

def RemainingFactors(u):
    # returns the remaining factors of u.
    if ComplexNumberQ(u) and Re(u) == 0:
        if Im(u) == 1:
            return S(1)
        else:
            return I*RemainingFactors(Im(u))
    elif ProductQ(u):
        return RemainingFactors(First(u))*Rest(u)
    return S(1)

def LeadBase(u):
    # returns the base of the leading factor of u.
    v = LeadFactor(u)
    if PowerQ(v):
        return v.base
    return v

def LeadDegree(u):
    # returns the degree of the leading factor of u.
    v = LeadFactor(u)
    if PowerQ(v):
        return v.exp
    return v

def Numer(expr):
    # returns the numerator of u.
    if PowerQ(expr):
        if expr.args[1] < 0:
            return 1
    if ProductQ(expr):
        return Mul(*[Numer(i) for i in expr.args])
    return Numerator(expr)

def Denom(u):
    # returns the denominator of u
    if PowerQ(u):
        if u.args[1] < 0:
            return u.args[0]**(-u.args[1])
    elif ProductQ(u):
        return Mul(*[Denom(i) for i in u.args])
    return Denominator(u)

def hypergeom(n, d, z):
    return hyper(n, d, z)

def Expon(expr, form, h=None):
    if h:
        return Exponent(Together(expr), form, h)
    else:
        return Exponent(Together(expr), form)

def MergeMonomials(expr, x):
    u_ = Wild('u')
    p_ = Wild('p', exclude=[x, 1, 0])
    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x, 0])
    c_ = Wild('c', exclude=[x])
    d_ = Wild('d', exclude=[x, 0])
    n_ = Wild('n', exclude=[x])
    m_ = Wild('m', exclude=[x])

    # Basis: If  m/n\[Element]\[DoubleStruckCapitalZ], then z^m (c z^n)^p==(c z^n)^(m/n+p)/c^(m/n)
    pattern = u_*(a_ + b_*x)**m_*(c_*(a_ + b_*x)**n_)**p_
    match = expr.match(pattern)
    if match:
        keys = [u_, a_, b_, m_, c_, n_, p_]
        if len(keys) == len(match):
            u, a, b, m, c, n, p = tuple([match[i] for i in keys])
            if IntegerQ(m/n):
                return u*(c*(a + b*x)**n)**(m/n + p)/c**(m/n)

    # Basis: If  m\[Element]\[DoubleStruckCapitalZ] \[And] b c-a d==0, then (a+b z)^m==b^m/d^m (c+d z)^m
    pattern = u_*(a_ + b_*x)**m_*(c_ + d_*x)**n_
    match = expr.match(pattern)
    if match:
        keys = [u_, a_, b_, m_, c_, d_, n_]
        if len(keys) == len(match):
            u, a, b, m, c, d, n = tuple([match[i] for i in keys])
            if IntegerQ(m) and ZeroQ(b*c - a*d):
                return u*b**m/d**m*(c + d*x)**(m + n)
    return expr

def PolynomialDivide(p_, q_, x):
    p = poly(p_, x)
    q = poly(q_, x)
    quotient = quo(p, q).as_expr()
    remainder = rem(p, q).as_expr()
    result = quotient
    if SumQ(remainder):
        for i in remainder.args:
            result += i/q_
    else:
        result += remainder/q_
    return result

def BinomialQ(u, x, n=None):
    if ListQ(u):
        for i in u.args:
            if Not(BinomialQ(i, x, n)):
                return False
        return True
    return ListQ(BinomialParts(u, x))

def TrinomialQ(u, x):
    # * If u is equivalent to an expression of the form a+b*x^n+c*x^(2*n) where n, b and c are not 0, *)
    # (* TrinomialQ[u,x] returns True; else it returns False. *)
    if ListQ(u):
        for i in u.args:
            if Not(TrinomialQ(i, x, n)):
                return False
        return True

    check = False
    if PowerQ(u):
        if u.exp == 2 and BinomialQ(u.base, x):
            check = True

    return ListQ(TrinomialParts(u,x)) and Not(QuadraticQ(u, x)) and Not(check)

def GeneralizedBinomialQ(u, x):
    # (* If u is equivalent to an expression of the form a*x^q+b*x^n where n, q and b are not 0, *)
    # (* GeneralizedBinomialQ[u,x] returns True; else it returns False. *)
    if ListQ(u):
        return all(GeneralizedBinomialQ(i, x) for i in u)
    return ListQ(GeneralizedBinomialParts(u, x))

def GeneralizedTrinomialQ(u, x):
    # (* If u is equivalent to an expression of the form a*x^q+b*x^n+c*x^(2*n-q) where n, q, b and c are not 0, *)
    # (* GeneralizedTrinomialQ[u,x] returns True; else it returns False. *)
    if ListQ(u):
        return all(GeneralizedTrinomialQ(i, x) for i in u)
    return ListQ(GeneralizedTrinomialParts(u, x))

def FactorSquareFreeList(poly):
    r = sqf_list(poly)
    result = [[1, 1]]
    for i in r[1]:
        result.append(list(i))
    return result

def PerfectPowerTest(u, x):
    # If u (x) is equivalent to a polynomial raised to an integer power greater than 1,
    # PerfectPowerTest[u,x] returns u (x) as an expanded polynomial raised to the power;
    # else it returns False.
    if PolynomialQ(u, x):
        lst = FactorSquareFreeList(u)
        gcd = 0
        v = 1
        if lst[0] == [1, 1]:
            lst = Rest(lst)
        for i in lst:
            gcd = GCD(gcd, i[1])
        if gcd > 1:
            for i in lst:
                v = v*i[0]**(i[1]/gcd)
            return Expand(v)**gcd
        else:
            return False
    return False

def SquareFreeFactorTest(u, x):
    # If u (x) can be square free factored, SquareFreeFactorTest[u,x] returns u (x) in
    # factored form; else it returns False.
    if PolynomialQ(u, x):
        v = FactorSquareFree(u)
        if PowerQ(v) or ProductQ(v):
            return v
        return False
    return False

def RationalFunctionQ(u, x):
    # If u is a rational function of x, RationalFunctionQ[u,x] returns True; else it returns False.
    if AtomQ(u) or FreeQ(u, x):
        return True
    elif IntegerPowerQ(u):
        return RationalFunctionQ(u.args[0], x)
    elif ProductQ(u) or SumQ(u):
        for i in u.args:
            if Not(RationalFunctionQ(i, x)):
                return False
        return True
    return False

def RationalFunctionFactors(u, x):
    # RationalFunctionFactors[u,x] returns the product of the factors of u that are rational functions of x.
    if ProductQ(u):
        res = 1
        for i in u.args:
            if RationalFunctionQ(i, x):
                res *= i
        return res
    elif RationalFunctionQ(u, x):
        return u
    return S(1)

def NonrationalFunctionFactors(u, x):
    if ProductQ(u):
        res = 1
        for i in u.args:
            if not RationalFunctionQ(i, x):
                res *= i
        return res
    elif RationalFunctionQ(u, x):
        return S(1)
    return u

def Reverse(u):
    if isinstance(u, list):
        return list(reversed(u))
    else:
        l = list(u.args)
        return u.func(*list(reversed(l)))

def RationalFunctionExponents(u, x):
    # (* u is a polynomial or rational function of x. *)
    # (* RationalFunctionExponents[u,x] returns a list of the exponent of the *)
    # (* numerator of u and the exponent of the denominator of u. *)
    if PolynomialQ(u, x):
        return [Exponent(u, x), 0]
    elif IntegerPowerQ(u):
        if PositiveQ(u.args[1]):
            return u.args[1]*RationalFunctionExponents(u.args[0], x)
        return  (-u.args[1])*Reverse(RationalFunctionExponents(u.args[0], x))
    elif ProductQ(u):
        lst1 = RationalFunctionExponents(First(u), x)
        lst2 = RationalFunctionExponents(Rest(u), x)
        return [lst1[0] + lst2[0], lst1[1] + lst2[1]]
    elif SumQ(u):
        v = Together(u)
        if SumQ(v):
            lst1 = RationalFunctionExponents(First(u), x)
            lst2 = RationalFunctionExponents(Rest(u), x)
            return [Max(lst1[0] + lst2[1], lst2[0] + lst1[1]), lst1[1] + lst2[1]]
        else:
            return RationalFunctionExponents(v, x)
    return [0, 0]

def RationalFunctionExpand(expr, x):
    # (* u is a polynomial or rational function of x. *)
    # (* RationalFunctionExpand[u,x] returns the expansion of the factors of u that are rational functions times the other factors. *)
    u_ = Wild('u', exclude=[1, 0])
    v_ = Wild('v', exclude=[1, 0])
    n_ = Wild('n', exclude=[x, 1, -1, 0])
    pattern = u_*v_**n_
    match = expr.match(pattern)
    if match:
        keys = [u_, v_, n_]
        if len(keys) == len(match):
            u, v, n = tuple([match[i] for i in keys])
            if FractionQ(n) and v != x:
                w = RationalFunctionExpand(u, x)
                if SumQ(w):
                    return Add(*[i*v**n for i in w.args])
                else:
                    return w*v**n

    v = ExpandIntegrand(expr, x)
    t = False

    a_ = Wild('a', exclude=[x, 0])
    b_ = Wild('b', exclude=[x, 0])
    c_ = Wild('c', exclude=[x, 0])
    d_ = Wild('d', exclude=[x, 0])
    m_ = Wild('m', exclude=[x])
    p_ = Wild('p', exclude=[x, 1])

    pattern = x**m_*(c_ + d_*x)**p_/(a_ + b_*x**n_)
    match = expr.match(pattern)

    if match:
        keys = [m_, c_, d_, p_, a_, b_, n_]
        if len(keys) == len(match):
            m, c, d, p, a, b, n = tuple([match[i] for i in keys])
            if IntegersQ(m, n) and m == n-1:
                t = True

    u = expr
    if v != u and t:
        return v

    v = ExpandIntegrand(RationalFunctionFactors(u, x), x)
    w = NonrationalFunctionFactors(u, x)
    if SumQ(v):
        return Add(*[i*w for i in v.args])
    return v*w

def ExpandIntegrand(expr, x, extra=None):
    if extra:
        w = ExpandIntegrand(extra, x)
        r = NonfreeTerms(w, x)
        if SumQ(r):
            result = [expr*FreeTerms(w, x)]
            for i in r.args:
                result.append(MergeMonomials(expr*i, x))
            return r.func(*result)
        else:
            return expr*FreeTerms(w, x) + MergeMonomials(expr*r, x)


    u_ =  Wild('c', exclude=[0, 1])
    a_ =  Wild('c', exclude=[x])
    b_ =  Wild('c', exclude=[x, 0])
    m_ =  Wild('c', exclude=[x, 0])
    f_ =  Wild('c', exclude=[x, 0, 1])
    e_ =  Wild('c', exclude=[x, 0])
    c_ =  Wild('c', exclude=[x])
    d_ =  Wild('c', exclude=[x, 0])
    n_ =  Wild('c', exclude=[x, 0])
    pattern = u_*(a_ + b_*x)**m_*f_**(e_*(c_ + d_*x)**n_)
    match = expr.match(pattern)
    if match:
        keys = [u_, a_, b_, m_, f_, e_, c_, d_, n_]
        if len(keys) == len(match):
            u, a, b, m, f, e, c, d, n = tuple([match[i] for i in keys])
            if PolynomialQ(u, x):
                v = ExpandIntegrand(u*(a + b*x)**m)
                if SumQ(v):
                    return Add(*[f**(e*(c + d*x)**n)*i for i in v.args])

    m_ =  Wild('m', exclude=[x, 0])
    e_ =  Wild('e', exclude=[x, 0])
    f_ =  Wild('f', exclude=[x, 0])
    p_ =  Wild('p', exclude=[x, 0])
    F_ =  Wild('F', exclude=[x, 0, 1])
    a_ =  Wild('a', exclude=[x])
    b_ =  Wild('b', exclude=[x, 0])
    c_ =  Wild('c', exclude=[x])
    d_ =  Wild('d', exclude=[x, 0])
    n_ =  Wild('n', exclude=[x, 0])
    pattern = x**m_*(e_ + f_*x)**p_*F_**(a_ + b_*(c_ + d_*x)**n_)
    match = expr.match(pattern)
    if match:
        keys = [m_, e_, f_, p_, F_, a_, b_, c_, d_, n_]
        if len(keys) == len(match):
            m, e, f, p, F, a, b, c, d, n = tuple([match[i] for i in keys])
            if PositiveIntegerQ(m, p) and m<=p and (EqQ(n, 1) or ZeroQ(d*e - c*f)):
                return ExpandLinearProduct((e + f*x)**p*F**(a + b*(c + d*x)**n), x**m, e, f, x)
            elif PositiveIntegerQ(p):
                v = Expand((e + f*x)**p)
                if SumQ(v):
                    return Add(*[x**m*F**(a + b*(c + d*x)**n)*i for i in v.args])
                else:
                    return x**m*F**(a + b*(c + d*x)**n)*v
            return ExpandIntegrand(F**(a + b*(c + d*x)**n), x**m*(e + f*x)**p, x)


    F_ =  Wild('F', exclude=[x, 0, 1])
    b_ =  Wild('b', exclude=[x, 0])
    c_ =  Wild('c', exclude=[x])
    d_ =  Wild('d', exclude=[x, 0])
    e_ =  Wild('e', exclude=[x, 0])
    f_ =  Wild('f', exclude=[x, 0])
    m_ =  Wild('m', exclude=[x, 0])
    n_ =  Wild('n', exclude=[x, 0])
    p_ =  Wild('p', exclude=[x, 0])
    pattern = x**m_*(e_ + f_*x)**p_*F_**(b_*(c_ + d_*x)**n_)
    match = expr.match(pattern)
    if match:
        keys = [m_, e_, f_, p_, F_, b_, c_, d_, n_]
        if len(match) == len(keys):
            m, e, f, p, F, b, c, d, n = tuple([match[i] for i in keys])
            if PositiveIntegerQ(m, p) and m<=p and (EqQ(n, 1) or ZeroQ(d*e - c*f)):
                return ExpandLinearProduct((e + f*x)**p*F**(b*(c + d*x)**n), x**m, e, f, x)
            elif PositiveIntegerQ(p):
                u = ((e + f*x)**p).expand()
                if SumQ(u):
                    return Add(*[x**m*F**(b*(c + d*x)**n)*i for i in u.args])
                else:
                    return x**m*F**(b*(c + d*x)**n)*u
            else:
                return Expand(F**(b*(c + d*x)**n), x, x**m*(e + f*x)**p)

    k, q, i = symbols('k q i')
    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x, 0])
    c_ = Wild('c', exclude=[x])
    d_ = Wild('d', exclude=[x, 0])
    e_ = Wild('e', exclude=[x, 0])
    f_ = Wild('f', exclude=[x])
    g_ = Wild('g', exclude=[x])
    h_ = Wild('h', exclude=[x, 0])
    n_ = Wild('n', exclude=[x, 0])
    m_ = Wild('m', exclude=[x, 0])
    # Basis: (a+b x)^m/(c+d x)==(b (a+b x)^(m-1))/d+((a d-b c) (a+b x)^(m-1))/(d (c+d x))
    pattern = (a_ + b_*x)**m_*f_**(e_*(c_ + d_*x)**n_)/(g_+h_*x)
    match = expr.match(pattern)
    if match:
        keys = [a_, b_, c_, d_, e_, f_, g_, h_, m_, n_]
        if len(keys) == len(match):
            a, b, c, d, e, f, g, h, m, n = tuple([match[i] for i in keys])
            if PositiveIntegerQ(m) and ZeroQ(b*c - a*d):
                tmp = a*h - b*g
                return SimplifyTerm(tmp**m/h**m, x)*f**(e*(c + d*x)**n)/(g + h*x) + Sum(SimplifyTerm(b*tmp**(k-1)/h**k, x)*f**(e*(c+d*x)**n)*(a + b*x)**(m-k), (k, 1, m)).doit()

    a_ = Wild('a', exclude=[x, 0])
    b_ = Wild('b', exclude=[x, 0])
    c_ = Wild('c', exclude=[x, 0])
    d_ = Wild('d', exclude=[x, 0])
    n_ = Wild('n', exclude=[x, 1, 0])
    m_ = Wild('m', exclude=[x, 0])
    F_ = Wild('F', exclude=[x, 1, 0])
    v_ = Wild('v', exclude=[0])
    u_ = Wild('u', exclude=[0])
    pattern = u_*(a_ + b_*F_**v_)**m_*(c_ + d_*F_**v_)**n_
    match = expr.match(pattern)
    if match:
        keys = [u_, a_, b_, F_, v_, m_, c_, d_, n_]
        if len(keys) == len(match):
            u, a, b, F, v, m, c, d, n = tuple([match[i] for i in keys])
            if IntegersQ(m, n) and NegativeQ(n):
                w = ReplaceAll(ExpandIntegrand((a + b*x)**m*(c + d*x)**n, x), {x: F**v})
                result = []
                for i in w.args:
                    result.append(i*u)
                return w.func(*result)

    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x, 0])
    c_ = Wild('c', exclude=[x, 0])
    d_ = Wild('d', exclude=[x])
    e_ = Wild('e', exclude=[x, 0])
    n_ = Wild('n', exclude=[x, 0])
    m_ = Wild('m', exclude=[x, 0])
    p_ = Wild('p', exclude=[x, 0])
    u_ = Wild('u', exclude=[0, 1])
    pattern = u_*(a_ + b_*x)**m_*Log(c_*(d_ + e_*x**n_)**p_)
    match = expr.match(pattern)
    if match:
        keys = [u_, a_, b_, m_, c_, d_, e_, n_, p_]
        if len(keys) == len(match):
            u, a, b, m, c, d, e, n, p = tuple([match[i] for i in keys])
            if PolynomialQ(u, x):
                return ExpandIntegrand(Log(c*(d + e*x**n)**p), x, u*(a + b*x)**m)

    c_ = Wild('c', exclude=[x])
    d_ = Wild('d', exclude=[x, 0])
    e_ = Wild('e', exclude=[x, 0])
    f_ = Wild('f', exclude=[x, 0, 1])
    u_ = Wild('u', exclude=[0, 1])
    n_ = Wild('n', exclude=[x, 0])
    pattern = u_*f_**(e_*(c_ + d_*x)**n_)
    match = expr.match(pattern)
    if match:
        keys = [u_, f_, e_, c_, d_, n_]
        if len(keys) == len(match):
            u, f, e, c, d, n = tuple([match[i] for i in keys])
            if PolynomialQ(u,x):
                if EqQ(n, 1):
                    return ExpandIntegrand(f**(e*(c + d*x)**n), x, u)
                else:
                    return ExpandLinearProduct(f**(e*(c + d*x)**n), u, c, d, x)

    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x, 0])
    c_ = Wild('c', exclude=[x, 0])
    d_ = Wild('d', exclude=[x, 0])
    e_ = Wild('e', exclude=[x])
    f_ = Wild('f', exclude=[x, 0])
    u_ = Wild('u', exclude=[0, 1])
    p_ = Wild('p', exclude=[x, 0])
    q_ = Wild('q', exclude=[x, 0])
    n_ = Wild('n', exclude=[x, 0])
    pattern = u_*(a_ + b_*Log(c_*(d_*(e_ + f_*x)**p_)**q_))**n_
    match = expr.match(pattern)
    if match:
        keys = [u_, a_, b_, c_, d_, e_, f_, p_, q_, n_]
        if len(keys) == len(match):
            u, a, b, c, d, e, f, p, q, n = tuple([match[i] for i in keys])
            if PolynomialQ(u, x):
                return ExpandLinearProduct((a + b*Log(c*(d*(e + f*x)**p)**q))**n, u, e, f, x)

    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x, 0])
    c_ = Wild('c', exclude=[x, 0])
    u_ = Wild('u', exclude=[0, 1])
    n_ = Wild('n', exclude=[0])
    j_ = Wild('j', exclude=[0])
    p_ = Wild('p', exclude=[0])
    pattern = (a_ + b_*u_**n_ + c_*u_**j_)**p_
    match = expr.match(pattern)
    if match:
        keys = [a_, b_, u_, n_, c_, j_, p_]
        if len(keys) == len(match):
            a, b, u, n, c, j, p = tuple([match[i] for i in keys])
            if IntegerQ(n) and ZeroQ(j - 2*n) and NegativeIntegerQ(p) and NonzeroQ(b**2 - 4*a*c):
                ReplaceAll(ExpandIntegrand(S(1)/(4**p*c**p), x, (b - q + 2*c*x)**p*(b + q + 2*c*x)**p), {q: Rt(b**2-4*a*c,2), x: u**n})

    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x, 0])
    c_ = Wild('c', exclude=[x, 0])
    n_ = Wild('n', exclude=[0])
    m_ = Wild('m', exclude=[0])
    j_ = Wild('j', exclude=[0])
    p_ = Wild('p', exclude=[1])
    u_ = Wild('u', exclude=[0, 1])
    pattern = u_**m_*(a_ + b_*u_**n_ + c_*u_**j_)**p_
    match = expr.match(pattern)
    if match:
        keys = [u_, m_, a_, b_, n_, c_, j_, p_]
        if len(keys) == len(match):
            u, m, a, b, n, c, j, p = tuple([match[i] for i in keys])
            if IntegersQ(m, n, j) and ZeroQ(j - 2*n) and NegativeIntegerQ(p) and 0<m<2*n and Not(m == n and p == -1) and NonzeroQ(b**2 - 4*a*c):
                return ReplaceAll(ExpandIntegrand(S(1)/(4**p*c**p), x, x**m*(b - q + 2*c*x**n)**p*(b + q+ 2*c*x**n)**p), {q: Rt(b**2 - 4*a*c, 2),x: u})

    a_ = Wild('a', exclude=[x, 0])
    c_ = Wild('c', exclude=[x, 0])
    u_ = Wild('u', exclude=[1, 0])
    n_ = Wild('n', exclude=[x, 0])
    p_ = Wild('p', exclude=[0, 1])
    # Basis: If  q=Sqrt[-a c], then a+c z^2==((-q+c z)(q+c z))/c
    pattern = (a_ + c_*u_**n_)**p_
    match = expr.match(pattern)
    if match:
        keys = [a_, c_, u_, n_, p_]
        if len(keys) == len(match):
            a, c, u, n, p = tuple([match[i] for i in keys])
            if IntegerQ(n/2) and NegativeIntegerQ(p):
                return ReplaceAll(ExpandIntegrand(S(1)/c**p, x, (-q + c*x)**p*(q + c*x)**p), {q: Rt(-a*c,2),x: u**(n/2)})

    u_ = Wild('u', exclude=[0, 1])
    m_ = Wild('m', exclude=[x, 0])
    a_ = Wild('a', exclude=[x])
    c_ = Wild('c', exclude=[x, 1])
    n_ = Wild('n', exclude=[x, 0])
    p_ = Wild('p', exclude=[x, 1, 0])
    pattern = u_**m_*(a_ + c_*u_**n_)**p_
    match = expr.match(pattern)
    if match:
        keys = [u_, m_, a_, c_, n_, p_]
        if len(keys) == len(match):
            u, m, a, c, n, p = tuple([match[i] for i in keys])
            if IntegersQ(m, n/2) and NegativeIntegerQ(p) and 0 < m < n and (m != n/2):
                return ReplaceAll(ExpandIntegrand(S(1)/c**p, x, x**m*(-q + c*x**(n/2))**p*(q + c*x**(n/2))**p),{q: Rt(-a*c, 2), x: u})

    u_ = Wild('u', exclude=[0])
    a_ = Wild('a', exclude=[x, 0])
    b_ = Wild('b', exclude=[x, 0])
    c_ = Wild('c', exclude=[x, 0])
    d_ = Wild('d', exclude=[x, 0])
    n_ = Wild('n', exclude=[x, 0, 1])
    j_ = Wild('j', exclude=[x, 0, 1])
    # Basis: 1/(a x^n+b Sqrt[c+d x^(2 n)])==(a x^n-b Sqrt[c+d x^(2 n)])/(-b^2 c+(a^2-b^2 d) x^(2 n))
    pattern = u_/(a_*x**n_ + b_*Sqrt(c_ + d_*x**j_))
    match = expr.match(pattern)
    if match:
        keys = [u_, a_, n_, b_, c_, d_, j_]
        if len(keys) == len(match):
            u, a, n, b, c, d, j = tuple([match[i] for i in keys])
            if ZeroQ(j - 2*n):
                return ExpandIntegrand(u*(a*x**n - b*Sqrt(c + d*x**(2*n)))/(-b**2*c + (a**2 - b**2*d)*x**(2*n)), x)

    d_ = Wild('d', exclude=[x])
    e_ = Wild('e', exclude=[x, 0])
    f_ = Wild('f', exclude=[x])
    g_ = Wild('g', exclude=[x, 0])
    u_ = Wild('u', exclude=[0])
    n_ = Wild('n', exclude=[x, 0])
    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x, 0])
    c_ = Wild('c', exclude=[x, 0])
    j_ = Wild('j', exclude=[x, 0])
    # Basis: If  q=Sqrt[b^2-4a c] and r=(2 c d-b e)/q, then (d+e z)/(a+b z+c z^2)==(e+r)/(b-q+2 c z)+(e-r)/(b+q+2 c z)*)
    pattern = (d_ + e_*(f_ + g_*u_**n_))/(a_ + b_*u_**n_ + c_*u_**j_)
    match = expr.match(pattern)
    if match:
        keys = [d_, e_, f_, g_, u_, n_, a_, b_, c_, j_]
        if len(keys) == len(match):
            d, e, f, g, u, n, a, b, c, j = tuple([match[i] for i in keys])
            if ZeroQ(j - 2*n) and NonzeroQ(b*2 - 4*a*c):
                q = Rt(b**2 - 4*a*c, 2)
                r = TogetherSimplify((2*c*(d + e*f) - b*e*g)/q)
                return (e*g + r)/(b - q + 2*c*u**n) + (e*g - r)/(b + q + 2*c*u**n)

    a_ = Wild('a', exclude=[x, 0])
    b_ = Wild('b', exclude=[x, 0])
    c_ = Wild('c', exclude=[x, 0])
    d_ = Wild('d', exclude=[x, 0])
    A_ = Wild('A', exclude=[x, 0])
    B_ = Wild('B', exclude=[x, 0])
    m_ = Wild('m', exclude=[x, 0])
    pattern = (a_ + b_*x)**m_*(A_ + B_*x)/(c_ + d_*x)
    match = expr.match(pattern)
    if match:
        keys = [a_, b_, m_, A_, B_, c_, d_]
        if len(match) == len(keys):
            a, b, m, A, B, c, d = tuple([match[i] for i in keys])
            if PositiveIntegerQ(m):
                if RationalQ(a, b, c, d, A, B):
                    return ExpandExpression((a + b*x)**m*(A + B*x)/(c + d*x), x)
                else:
                    tmp1 = (A*d - B*c)/d
                    tmp2 = ExpandIntegrand((a + b*x)**m/(c + d*x), x)
                    if SumQ(tmp2):
                        tmp2 = Add(*[SimplifyTerm(tmp1*i, x) for i in tmp2.args])
                    else:
                        tmp2 = SimplifyTerm(tmp1*tmp2, x)
                    return SimplifyTerm(B/d, x)*(a + b*x)**m + tmp2

    c_ = Wild('c', exclude=[x])
    d_ = Wild('d', exclude=[x, 0])
    u_ = Wild('u', exclude=[0])
    m_ = Wild('m', exclude=[x, 0])
    n_ = Wild('n', exclude=[x, 0, 1])
    e_ = Wild('e', exclude=[x, 0])
    p_ = Wild('p', exclude=[x, 0, 1])
    f_ = Wild('f', exclude=[x, 0])
    q_ = Wild('q', exclude=[x, 0, 1])
    a_ = Wild('a', exclude=[x, 0])
    b_ = Wild('b', exclude=[x, 0])
    # Basis: If (m|n,p,q)\[Element]\[DoubleStruckCapitalZ] \[And] 0<=m<p<q<n, let r/s=(-(a/b))^(1/n), then  (c + d*z^m + e*z^p + f*z^q)/(a + b*z^n) == (r*Sum[(c + (d*(r/s)^m)/(-1)^(2*k*(m/n)) + (e*(r/s)^p)/(-1)^(2*k*(p/n)) + (f*(r/s)^q)/(-1)^(2*k*(q/n)))/(r - (-1)^(2*(k/n))*s*z), {k, 1, n}])/(a*n)
    pattern = (c_ + d_*u_**m_ + e_*u_**p_ + f_*u_**q_)/(a_ + b_*u_**n_)
    match = expr.match(pattern)
    if match:
        keys = [c_, d_, u_, m_, e_, p_, f_, q_, a_, b_, n_]
        if len(keys) == len(match):
            c, d, u, m, e, p, f, q, a, b, n = tuple([match[i] for i in keys])
            if IntegersQ(m, n, p, q) & 0<m<p<q<n:
                r = Numerator(Rt(-a/b, n))
                s = Denominator(Rt(-a/b, n))
                return Sum((r*c + r*d*(r/s)**m*(-1)**(-2*k*m/n) + r*e*(r/s)**p*(-1)**(-2*k*p/n) + r*f*(r/s)**q*(-1)**(-2*k*q/n))/(a*n*(r - (-1)**(2*k/n)*s*u)),(k,1,n)).doit()

    c_ = Wild('c', exclude=[x])
    d_ = Wild('d', exclude=[x, 0])
    u_ = Wild('u', exclude=[0])
    m_ = Wild('m', exclude=[x, 0])
    e_ = Wild('e', exclude=[x, 0])
    p_ = Wild('p', exclude=[x, 0, 1])
    a_ = Wild('a', exclude=[x, 0])
    b_ = Wild('b', exclude=[x, 0])
    n_ = Wild('n', exclude=[x, 0, 1])
    # Basis: If (m|n,p)\[Element]\[DoubleStruckCapitalZ] \[And] 0<=m<p<n, let r/s=(-(a/b))^(1/n), then  (c + d*z^m + e*z^p)/(a + b*z^n) == (r*Sum[(c + (d*(r/s)^m)/(-1)^(2*k*(m/n)) + (e*(r/s)^p)/(-1)^(2*k*(p/n)))/(r - (-1)^(2*(k/n))*s*z), {k, 1, n}])/(a*n)
    pattern = (c_ + d_*u_**m_ + e_*u_**p_)/(a_ + b_*u_**n_)
    match = expr.match(pattern)
    if match:
        keys = [c_, d_, u_, m_, e_, p_, a_, b_, n_]
        if len(keys) == len(match):
            c, d, u, m, e, p, a, b, n = tuple([match[i] for i in keys])
            if IntegersQ(m, n, p) & 0<m<p<n:
                r = Numerator(Rt(-a/b, n))
                s = Denominator(Rt(-a/b, n))
                return Sum((r*c + r*d*(r/s)**m*(-1)**(-2*k*m/n) + r*e*(r/s)**p*(-1)**(-2*k*p/n))/(a*n*(r - (-1)**(2*k/n)*s*u)),(k, 1, n)).doit()

    F_ = Wild('F', exclude=[0, 1])
    m_ = Wild('m', exclude=[x, 0])
    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x, 0])
    G_ = Wild('G', exclude=[0, 1])
    n_ = Wild('n', exclude=[x, 0])
    pattern = F_**m_*(a_ + b_*G_)**n_
    match = expr.match(pattern)
    if match:
        keys = [F_, m_, a_, b_, G_, n_]
        if len(match) == len(keys):
            F, m, a, b, G, n = tuple([match[i] for i in keys])
            if ((1/F).args == G.args and F*G ==1 and IntegersQ(m, n)):
                return ReplaceAll(ExpandIntegrand((a + b*x)**n/x**m, x), {x: G})


    a_ = Wild('a', exclude=[x, 0])
    b_ = Wild('b', exclude=[x, 0])
    m_ = Wild('m', exclude=[x, 0, 1, -1])
    c_ = Wild('c', exclude=[x, 0])
    d_ = Wild('d', exclude=[x, 0])
    # Basis: (a+b x)^m/(c+d x)==(b (a+b x)^(m-1))/d+((a d-b c) (a+b x)^(m-1))/(d (c+d x))
    pattern = (a_ + b_*x)**m_/(c_ + d_*x)
    match = expr.match(pattern)
    if match:
        keys = [a_, b_, c_, d_, m_]
        if len(keys) == len(match):
            a, b, c, d, m = tuple([match[i] for i in keys])
            if PositiveIntegerQ(m):
                if RationalQ(a, b, c, d):
                    return ExpandExpression((a + b*x)**m/(c + d*x), x)
                else:
                    tmp = a*d - b*c
                    result = SimplifyTerm(tmp**m / d**m, x)/(c + d*x)
                    for k in range(1, m + 1):
                        result += SimplifyTerm(b*tmp**(k - 1)/d**k, x)*(a + b*x)**(m - k)
                    return result

    c_ = Wild('c', exclude=[x])
    d_ = Wild('d', exclude=[x, 0])
    u_ = Wild('u', exclude=[0, 1])
    m_ = Wild('m', exclude=[x, 0])
    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x, 0])
    n_ = Wild('n', exclude=[x, 0])
    pattern = (c_ + d_*u_**m_)/(a_ + b_*u_**n_)
    match = expr.match(pattern)
    if match:
        keys = [c_, d_, u_, m_, a_, b_, n_]
        if len(keys) == len(match):
            c, d, u, m, a, b, n = tuple([match[i] for i in keys])
            if IntegersQ(m,n) & 0<m<n:
                # *Basis: If (m|n)\[Element]\[DoubleStruckCapitalZ] \[And] 0<=m<n, let r/s=(-(a/b))^(1/n), then  (c + d*z^m)/(a + b*z^n) == (r*Sum[(c + (d*(r/s)^m)/(-1)^(2*k*(m/n)))/(r - (-1)^(2*(k/n))*s*z), {k, 1, n}])/(a*n)
                r = Numerator(Rt(-a/b, n))
                s = Denominator(Rt(-a/b, n))
                return Sum((r*c + r*d*(r/s)**m*(-1)**(-2*k*m/n))/(a*n*(r-(-1)**(2*k/n)*s*u)),(k,1,n)).doit()
            elif ZeroQ(n - 2*m):
                # q=Sqrt[-(a/b)], then (c+d z)/(a+b z^2)==-((c-d q)/(2 b q(q+z)))-(c+d q)/(2 b q(q-z))
                j = n
                n = m
                q = Rt(-a/b, 2)
                return -(c - d*q)/(2*b*q*(q + u**n)) - (c + d*q)/(2*b*q*(q - u**n))

    a_ = Wild('a', exclude=[x, 0])
    b_ = Wild('b', exclude=[x, 0])
    u_ = Wild('u', exclude=[0, 1])
    n_ = Wild('n', exclude=[x, 0, -1])
    pattern = 1/(a_ + b_*u_**n_)
    match = expr.match(pattern)
    if match:
        keys = [a_, b_, u_, n_]
        if len(keys) == len(match):
            a, b, u, n = tuple([match[i] for i in keys])
            if n == 3:
                # Basis: Let r/s=(-(a/b))^(1/3), then  1/(a+b z^3)==r/(3a(r-s z))+(r(2 r+s z))/(3a(r^2+r s z+s^2 z^2))
                r = Numerator(Rt(-a/b, 3))
                s = Denominator(Rt(-a/b, 3))
                return r/(S(3)*a*(r - s*u)) + r*(2*r + s*u)/(3*a*(r**2 + r*s*u + s**2*u**2))
            elif PositiveIntegerQ(n/4):
                # Let r/s=Sqrt[-(a/b)], then  1/(a+b z^2)==r/(2a(r-s z))+r/(2a(r+s z))
                r = Numerator(Rt(-a/b,2))
                s = Denominator(Rt(-a/b, 2))
                return r/(2*a*(r - s*u**(n/2))) + r/(2*a*(r + s*u**(n/2)))
            elif IntegerQ(n) & PositiveQ(n):
                # Basis: If  n\[Element]SuperPlus[\[DoubleStruckCapitalZ]], let r/s=(-(a/b))^(1/n), then  1/(a + b*z^n) == (r*Sum[1/(r - (-1)^(2*(k/n))*s*z), {k, 1, n}])/(a*n)
                r = Numerator(Rt(-a/b, n))
                s = Denominator(Rt(-a/b, n))
                return Sum(r/(a*n*(r - (-1)**(2*k/n)*s*u)),(k, 1, n)).doit()

    u_ = Wild('u', exclude=[0, 1])
    m_ = Wild('m', exclude=[x])
    a_ = Wild('a', exclude=[x, 0])
    b_ = Wild('b', exclude=[x, 0])
    n_ = Wild('n', exclude=[x, 0, 1])
    # Basis: If  (m|(n-1)/2)\[Element]\[DoubleStruckCapitalZ] \[And] 0<=m<n, let r/s=(a/b)^(1/n), then z^m/(a + b*z^n) == (r*(-(r/s))^m*Sum[1/((-1)^(2*k*(m/n))*(r + (-1)^(2*(k/n))*s*z)), {k, 1, n}])/(a*n) == (r*(-(r/s))^m*Sum[(-1)^(2*k*((m + 1)/n))/((-1)^(2*(k/n))*r + s*z), {k, 1, n}])/(a*n)
    pattern = u_**m_/(a_ + b_*u_**n_)
    match = expr.match(pattern)
    if match:
        keys = [u_, m_, a_, b_, n_]
        if len(keys) == len(match):
            u, m, a, b, n = tuple([match[i] for i in keys])
            if IntegersQ(m, n) and 0 < m < n and OddQ(n/GCD(m, n)) and PosQ(a/b):
                g = GCD(m, n)
                r = Numerator(Rt(a/b, n/GCD(m, n)))
                s = Denominator(Rt(a/b, n/GCD(m, n)))
                if CoprimeQ(m + g, n):
                    return Sum(r*(-r/s)**(m/g)*(-1)**(-2*k*m/n)/(a*n*(r + (-1)**(2*k*g/n)*s*u**g)),(k, 1, n/g)).doit()
                else:
                    return Sum(r*(-r/s)**(m/g)*(-1)**(2*k*(m+g)/n)/(a*n*((-1)**(2*k*g/n)*r + s*u**g)),(k, 1, n/g)).doit()
            elif IntegersQ(m, n) and 0 < m < n:
                g = GCD(m, n)
                r = Numerator(Rt(-a/b, n/GCD(m, n)))
                s = Denominator(Rt(-a/b, n/GCD(m, n)))
                if n/g == 2:
                    return s/(2*b*(r + s*u**g)) - s/(2*b*(r - s*u**g))
                else:
                    if CoprimeQ[m+g,n]:
                        return Sum(r*(r/s)**(m/g)*(-1)**(-2*k*m/n)/(a*n*(r - (-1)**(2*k*g/n)*s*u**g)),(k,1,n/g)).doit()
                    else:
                        return Sum(r*(r/s)**(m/g)*(-1)**(2*k*(m+g)/n)/(a*n*((-1)**(2*k*g/n)*r - s*u**g)),(k,1,n/g)).doit()

    u_ = Wild('u', exclude=[0, 1])
    a_ = Wild('a', exclude=[x, 0])
    b_ = Wild('b', exclude=[x, 0])
    m_ = Wild('m', exclude=[x, -1, 0])
    # If u is a polynomial in x, ExpandIntegrand[u*(a+b*x)^m,x] expand u*(a+b*x)^m into a sum of terms of the form A*(a+b*x)^n.
    pattern = u_*(a_ + b_*x)**m_
    match = expr.match(pattern)
    if match:
        keys = [u_, a_, b_, m_]
        if len(keys) == len(match):
            u, a, b, m = tuple([match[i] for i in keys])
            w_ = Wild('w', exclude=[0])
            c_ = Wild('c', exclude=[x, 0])
            d_ = Wild('d', exclude=[x, 0])
            p_ = Wild('p', exclude=[x, 0, 1])

            if PolynomialQ(u, x):
                if PositiveIntegerQ(m):
                    return (u*(a + b*x)**m).expand()

            match = u.match(w_*(c_+d_*x)**p_)
            if match:
                if IntegerQ(match[p_]) and match[p_] >= m:
                    res = True
                else:
                    res = False
            else:
                res = False

            if PolynomialQ(u, x) and Not(PositiveIntegerQ(m) and res) and (u != 1):
                tmp1 = ExpandLinearProduct((a+b*x)**m, u, a, b, x)
                if not IntegerQ(m):
                    return tmp1
                else:
                    tmp2 = ExpandExpression(u*(a+b*x)**m, x)
                    if SumQ(tmp2) and (LeafCount(tmp2) <= LeafCount(tmp1)+2):
                        return tmp2
                    else:
                        return tmp1

    u_ = Wild('u', exclude=[0, 1])
    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x, 0])
    F_ = Wild('F', exclude=[0])
    c_ = Wild('c', exclude=[x])
    d_ = Wild('d', exclude=[x, 0])
    n_ = Wild('n', exclude=[0, 1])
    pattern = u_*(a_ + b_*F_)**n_
    match = expr.match(pattern)
    if match:
        if MemberQ([asin, acos, asinh, acosh], match[F_].func):
            keys = [u_, a_, b_, F_, n_]
            if len(match) == len(keys):
                u, a, b, F, n = tuple([match[i] for i in keys])
                match = F.args[0].match(c_ + d_*x)
                if match:
                    keys = c_, d_
                    if len(keys) == len(match):
                        c, d = tuple([match[i] for i in keys])
                        if PolynomialQ(u, x):
                            F = F.func
                            return ExpandLinearProduct((a + b*F(c + d*x))**n, u, c, d, x)

    u_ = Wild('u', exclude=[0, 1])
    v_ = Wild('v', exclude=[0, 1])
    n_ = Wild('n', exclude=[x, 1, 0])
    a_ = Wild('a', exclude=[x, 0])
    b_ = Wild('b', exclude=[x, 0])
    m_ = Wild('m', exclude=[x, 0, 1])
    # (* If u is a polynomial in x, ExpandIntegrand[u*(a+b*x)^m,x] expand u*(a+b*x)^m into a sum of terms of the form A*(a+b*x)^n. *)
    pattern = u_*v_**n_*(a_ + b_*x)**m_
    match = expr.match(pattern)
    if match:
        keys = [u_, v_, n_, a_, b_, m_]
        if len(keys) == len(match):
            u, v, n, a, b, m = tuple([match[i] for i in keys])
            if NegativeIntegerQ(n) & Not(IntegerQ(m)) & PolynomialQ(u, x) & PolynomialQ(v, x) & RationalQ(m) & (m < -1) & (Exponent(u, x) >= -(n+IntegerPart(m))*Exponent(v, x)):
                pr = PolynomialQuotientRemainder(u, v**(-n)*(a + b*x)**(-IntegerPart(m)), x)
                return ExpandIntegrand(pr[0]*(a + b*x)**FractionalPart(m), x) + ExpandIntegrand(pr[1]*v**n*(a + b*x)**m, x)
            elif NegativeIntegerQ(n) & Not(IntegerQ(m)) & PolynomialQ(u, x) & PolynomialQ(v, x) & (Exponent(u, x) >= -n*Exponent(v, x)):
                pr = PolynomialQuotientRemainder(u, v**(-n),x)
                return ExpandIntegrand(pr[0]*(a + b*x)**m, x) + ExpandIntegrand(pr[1]*v**n*(a + b*x)**m, x)

    u_ = Wild('u', exclude=[0, 1])
    v_ = Wild('v', exclude=[0, 1])
    p_ = Wild('p', exclude=[0, 1])
    pattern = u_*v_**p_
    match = expr.match(pattern)
    if match:
        keys = [u_, v_, p_]
        if len(match) == len(keys):
            u, v, p = tuple([match[i] for i in keys])
            if Not(IntegerQ(p)):
                if PolynomialQ(u, x) and FreeQ(v/x, x):
                    return ExpandToSum((v)**p, u, x)
                else:
                    return ExpandIntegrand(NormalizeIntegrand(v**p, x), x, u)

    u_ = Wild('u', exclude=[0, 1])
    v_ = Wild('v', exclude=[0, 1])
    pattern = u_/v_
    match = expr.match(pattern)
    if match:
        keys = [u_, v_]
        if len(keys) == len(match):
            u, v = Numerator(expr), Denominator(expr)
            if PolynomialQ(u, x) and PolynomialQ(v, x) and BinomialQ(v,x) and (Exponent(u, x) == Exponent(v, x)-1 >= 2):
                lst = CoefficientList(u, x)
                result = lst[-1]*x**Exponent(u,x)/v
                for i in range(0, Exponent(u,x) + 1):
                    result += lst[i-1]*x**S(i-1)/v
            elif PolynomialQ(u, x) and PolynomialQ(v, x) and Exponent(u, x) >= Exponent(v, x):
                return PolynomialDivide(u, v, x)

    return ExpandExpression(expr, x)

def SimplerQ(u, v):
    # If u is simpler than v, SimplerQ(u, v) returns True, else it returns False.  SimplerQ(u, u) returns False
    if IntegerQ(u):
        if IntegerQ(v):
            if Abs(u)==Abs(v):
                return v<0
            else:
                return Abs(u)<Abs(v)
        else:
            return True
    if IntegerQ(v):
        return False
    if FractionQ(u):
        if FractionQ(v):
            if Denominator(u)==Denominator(v):
                return SimplerQ(Numerator(u), Numerator(v))
            else:
                return Denominator(u)<Denominator(v)
        else:
            return True
    if FractionQ(v):
        return False
    if (Re(u)==0 or Re(u)==0.0) and (Re(v)==0 or Re(v)==0.0):
        return SimplerQ(Im(u), Im(v))
    if ComplexNumberQ(u):
        if ComplexNumberQ(v):
            if Re(u)==Re(v):
                return SimplerQ(Im(u), Im(v))
            else:
                return SimplerQ(Re(u),Re(v))
        else:
            return False
    if NumberQ(u):
        if NumberQ(v):
            return OrderedQ([u,v])
        else:
            return True
    if NumberQ(v):
        return False
    if AtomQ(u):
        if AtomQ(v):
            return OrderedQ([u,v])
        else:
            return True
    if AtomQ(v):
        return False
    if Head(u)==Head(v):
        if Length(u)==Length(v):
            for i in range(len(u)):
                if not u[i]==v[i]:
                    return SimplerQ(u[i], v[i])
                else:
                    return False
        return Length(u)<Length(v)
    if LeafCount(u)<LeafCount(v):
        return True
    if LeafCount(v)<LeafCount(u):
        return False
    else:
        return Not(OrderedQ([v,u]))

def SimplerSqrtQ(u, v):
    # If Rt(u, 2) is simpler than Rt(v, 2), SimplerSqrtQ(u, v) returns True, else it returns False.  SimplerSqrtQ(u, u) returns False
    if NegativeQ(v) and Not(NegativeQ(u)):
        return True
    if NegativeQ(u) and Not(NegativeQ(v)):
        return False
    sqrtu = Rt(u, 2)
    sqrtv = Rt(v, 2)
    if IntegerQ(sqrtu):
        if IntegerQ(sqrtv):
            return sqrtu<sqrtv
        else:
            return True
    if IntegerQ(sqrtv):
        return False
    if RationalQ(sqrtu):
        if RationalQ(sqrtv):
            return sqrtu<sqrtv
        else:
            return True
    if RationalQ(sqrtv):
        return False
    if PosQ(u):
        if PosQ(v):
            return LeafCount(sqrtu)<LeafCount(sqrtv)
        else:
            return True
    if PosQ(v):
        return False
    if LeafCount(sqrtu)<LeafCount(sqrtv):
        return True
    if LeafCount(sqrtv)<LeafCount(sqrtu):
        return False
    else:
        return Not(OrderedQ([v, u]))

def SumSimplerQ(u, v):
    if RationalQ(u, v):
        if v == S(0):
            return False
        elif v > S(0):
            return u < -S(1)
        else:
            return u >= -v
    else:
        return SumSimplerAuxQ(Expand(u), Expand(v))

def SumSimplerAuxQ(u, v):
    if SumQ(v):
        return (RationalQ(First(v)) or SumSimplerAuxQ(u, First(v))) and (RationalQ(Rest(v)) or SumSimplerAuxQ(u, Rest(v)))
    elif SumQ(u):
        return SumSimplerAuxQ(First(u), v) or SumSimplerAuxQ(Rest(u), v)
    else:
        return v!=S(0) and NonnumericFactors(u)==NonnumericFactors(v) and (NumericFactor(u)/NumericFactor(v)< -S(1)/2 or NumericFactor(u)/NumericFactor(v)== -S(1)/2 and NumericFactor(u)<S(0))

def BinomialDegree(u, x):
    # if u is a binomial. BinomialDegree[u,x] returns the degree of x in u.
    return BinomialParts(u, x)[2]

def TrinomialDegree(u, x):
    # If u is equivalent to a trinomial of the form a + b*x^n + c*x^(2*n) where n!=0, b!=0 and c!=0, TrinomialDegree[u,x] returns n
    t = TrinomialParts(u, x)
    if t:
        return t[3]
    return t

def CancelCommonFactors(u, v):
    # CancelCommonFactors[u,v] returns {u',v'} are the noncommon factors of u and v respectively.
    com_fac = gcd(u, v)
    return [cancel(u/com_fac), cancel(v/com_fac)]

def SimplerIntegrandQ(u, v, x):
    lst = CancelCommonFactors(u, v)
    u1 = lst[0]
    v1 = lst[1]
    if Head(u1) == Head(v1) and Length(u1) == 1 and Length(v1) == 1:
        return SimplerIntegrandQ(u1.args[0], v1.args[0], x)
    if LeafCount(u1)<3/4*LeafCount(v1):
        return True
    if RationalFunctionQ(u1, x):
        if RationalFunctionQ(v1, x):
            t1 = 0
            t2 = 0
            for i in RationalFunctionExponents(u1, x):
                t1 += i
            for i in RationalFunctionExponents(v1, x):
                t2 += i
            return t1 < t2
        else:
            return True
    else:
        return False

def GeneralizedBinomialDegree(u, x):
    b = GeneralizedBinomialParts(u, x)
    if b:
        return b[2] - b[3]

def GeneralizedBinomialParts(expr, x):
    expr = Expand(expr)
    if GeneralizedBinomialMatchQ(expr, x):
        a = Wild('a', exclude=[x])
        b = Wild('b', exclude=[x])
        n = Wild('n', exclude=[x])
        m = Wild('m', exclude=[x])
        q = Wild('q', exclude=[x])
        u = Wild('u')
        Match = expr.match(a*x**q + b*x**n)
        if Match and PosQ(Match[q] - Match[n]):
            return [Match[b], Match[a], Match[q], Match[n]]
    else:
        return False

def GeneralizedTrinomialDegree(u, x):
    t = GeneralizedTrinomialParts(u, x)
    if t:
        return t[3] - t[4]

def GeneralizedTrinomialParts(expr, x):
    expr = Expand(expr)
    if GeneralizedTrinomialMatchQ(expr, x):
        a = Wild('a', exclude=[x])
        b = Wild('b', exclude=[x])
        c = Wild('c', exclude=[x])
        n = Wild('n', exclude=[x])
        q = Wild('q', exclude=[x])
        Match = expr.match(a*x**q+b*x**n+c*x**(2*n-q))
        if Match and expr.is_Add:
            return [Match[c], Match[b], Match[a], Match[n], 2*Match[n]-Match[q]]
    else:
        return False

def MonomialQ(u, x):
    # If u is of the form a*x^n where n!=0 and a!=0, MonomialQ[u,x] returns True; else False
    if isinstance(u, list):
        return all(MonomialQ(i) for i in u)
    else:
        a = Wild('a', exclude=[x])
        b = Wild('b', exclude=[x])
        re = u.match(a*x**b)
        if re:
            return True

def MonomialSumQ(u, x):
    # if u(x) is a sum and each term is free of x or an expression of the form a*x^n, MonomialSumQ(u, x) returns True; else it returns False
    if SumQ(u):
        for i in u.args:
            if Not(FreeQ(i, x) or MonomialQ(i, x)):
                return False
        return True

def MinimumMonomialExponent(u, x):
    # u is sum whose terms are monomials.  MinimumMonomialExponent(u, x) returns the exponent of the term having the smallest exponent
    lst = []
    for i in u.args:
        lst = lst + [MonomialExponent(i, x)]
    return min(lst)

def MonomialExponent(u, x):
    # u is a monomial. MonomialExponent(u, x) returns the exponent of x in u
    a = Wild('a', exclude=[x])
    b = Wild('b', exclude=[x])
    re = u.match(a*x**b)
    if re:
        return re[b]

def LinearMatchQ(u, x):
    # LinearMatchQ(u, x) returns True iff u matches patterns of the form a+b*x where a and b are free of x
    if isinstance(u, list):
        return all(LinearMatchQ(i, x) for i in u)
    else:
        a = Wild('a', exclude=[x])
        b = Wild('b', exclude=[x])
        re = u.match(a + b*x)
        if re:
            return True

def PowerOfLinearMatchQ(u, x):
    if isinstance(u, list):
        for i in u:
            if not PowerOfLinearMatchQ(i, x):
                return False
        return True
    else:
        a = Wild('a', exclude=[x])
        b = Wild('b', exclude=[x])
        m = Wild('m', exclude=[x])
        Match = u.match((a + b*x)**m)
        if Match and Match[a] and Match[b] and Match[m]:
            try:
                return True
            except KeyError:
                return False
        else:
            return False

def QuadraticMatchQ(u, x):
    return QuadraticQ(u, x)

def CubicMatchQ(u, x):
    if isinstance(u, list):
        for i in u:
            if not CubicMatchQ(i, x):
                return False
        return True
    else:
        a = Wild('a', exclude=[x])
        b = Wild('b', exclude=[x])
        c = Wild('c', exclude=[x])
        d = Wild('d', exclude=[x])
        Match = Expand(u).match(a + b*x + c*x**2 + d*x**3)
        if Match and Match[a] and Match[d]:
            return True
        else:
            return False

def BinomialMatchQ(u, x):
    if isinstance(u, list):
        for i in u:
            if not BinomialMatchQ(i, x):
                return False
        return True
    else:
        a = Wild('a', exclude=[x])
        b = Wild('b', exclude=[x])
        n = Wild('n', exclude=[x])
        Match = u.match(a + b*x**n)
        if Match and Match[a] and Match[b] and Match[n]:
            return True
        else:
            return False

def TrinomialMatchQ(u, x):
    if isinstance(u, list):
        for i in u:
            if not TrinomialMatchQ(i, x):
                return False
        return True
    else:
        a = Wild('a', exclude=[x])
        b = Wild('b', exclude=[x])
        n = Wild('n', exclude=[x])
        c = Wild('c', exclude=[x])
        Match = Expand(u).match(a + b*x**n + c*x**(2*n))
        if Match and Match[a] and Match[b] and Match[n] and Match[c]:
            return True
        else:
            return False

def GeneralizedBinomialMatchQ(u, x):
    if isinstance(u, list):
        for i in u:
            if not GeneralizedBinomialMatchQ(i, x):
                return False
        return True
    else:
        a = Wild('a', exclude=[x])
        b = Wild('b', exclude=[x])
        n = Wild('n', exclude=[x])
        q = Wild('q', exclude=[x])
        Match = u.match(a*x**q + b*x**n)
        if Match and len(Match) == 4 and Match[q] != 0 and Match[n] != 0:
            return True
        else:
            return False

def GeneralizedTrinomialMatchQ(u, x):
    if isinstance(u, list):
        for i in u:
            if not GeneralizedTrinomialMatchQ(i, x):
                return False
        return True
    else:
        a = Wild('a', exclude=[x])
        b = Wild('b', exclude=[x])
        n = Wild('n', exclude=[x])
        c = Wild('c', exclude=[x])
        q = Wild('q', exclude=[x])
        Match = u.match(a*x**q + b*x**n + c*x**(2*n - q))
        if Match and len(Match) == 5 and 2*Match[n] - Match[q] != 0 and Match[n] != 0:
            return True
        else:
            return False

def QuotientOfLinearsMatchQ(u, x):
    if isinstance(u, list):
        for i in u:
            if not QuotientOfLinearsMatchQ(i, x):
                return False
        return True
    else:
        a = Wild('a', exclude=[x])
        b = Wild('b', exclude=[x])
        d = Wild('d', exclude=[x])
        c = Wild('c', exclude=[x])
        e = Wild('e')
        Match = u.match(e*(a + b*x)/(c + d*x))
        if Match and len(Match) == 5:
            return True
        else:
            return False

def PolynomialTermQ(u, x):
    a = Wild('a', exclude=[x])
    n = Wild('n', exclude=[x])
    Match = u.match(a*x**n)
    if Match and IntegerQ(Match[n]) and Match[n]>0:
        return True
    else:
        return False

def PolynomialTerms(u, x):
    s = 0
    for i in u.args:
        if PolynomialTermQ(i, x):
            s = s + i
    return s

def NonpolynomialTerms(u, x):
    s = 0
    for i in u.args:
        if not PolynomialTermQ(i, x):
            s = s + i
    return s

def PseudoBinomialParts(u, x):
    a = Wild('a', exclude=[x])
    b = Wild('b', exclude=[x])
    d = Wild('d', exclude=[x])
    c = Wild('c', exclude=[x])
    n = Wild('n', exclude=[x])
    Match = u.match(a + b*(c + d*x)**n)
    if Match and Match[n]>2 and Match[a] and Match[b] and Match[c] and Match[d] and Match[n]:
        return [Match[a], Match[b], Match[c], Match[d], Match[n]]
    else:
        return False

def NormalizePseudoBinomial(u, x):
    lst = PseudoBinomialParts(u, x)
    if lst:
        return (lst[0] + lst[1]*(lst[2] + lst[3]*x)**lst[4])

def PseudoBinomialPairQ(u, v, x):
    lst1 = PseudoBinomialParts(u, x)
    if not lst1:
        return False
    lst2 = PseudoBinomialParts(v, x)
    if not lst2:
        return False
    else:
        return Drop(lst1, 2) == Drop(lst2, 2)

def PseudoBinomialQ(u, x):
    lst = PseudoBinomialParts(u, x)
    if lst:
        return True
    else:
        return False

def PolynomialGCD(f, g):
    return gcd(f, g)

def PolyGCD(u, v, x):
    # (* u and v are polynomials in x. *)
    # (* PolyGCD[u,v,x] returns the factors of the gcd of u and v dependent on x. *)
    return NonfreeFactors(PolynomialGCD(u, v), x)

def AlgebraicFunctionFactors(u, x, flag=False):
    # (* AlgebraicFunctionFactors[u,x] returns the product of the factors of u that are algebraic functions of x. *)
    if ProductQ(u):
        result = 1
        for i in u.args:
            if AlgebraicFunctionQ(i, x, flag):
                result *= i
        return result
    if AlgebraicFunctionQ(u, x, flag):
        return u
    return 1

def NonalgebraicFunctionFactors(u, x):
    # (* NonalgebraicFunctionFactors[u,x] returns the product of the factors of u that are not algebraic functions of x. *)
    if ProductQ(u):
        result = 1
        for i in u.args:
            if not AlgebraicFunctionQ(i, x):
                result *= i
        return result
    if AlgebraicFunctionQ(u, x):
        return 1
    return u

def QuotientOfLinearsP(u, x):
    if LinearQ(u, x):
        return True
    elif SumQ(u):
        if FreeQ(u.args[0], x):
            return QuotientOfLinearsP(Rest(u), x)
    elif LinearQ(Numerator(u), x) and LinearQ(Denominator(u), x):
        return True
    elif ProductQ(u):
        if FreeQ(First(u), x):
            return QuotientOfLinearsP(Rest(u), x)
    elif Numerator(u) == 1 and PowerQ(u):
        return QuotientOfLinearsP(Denominator(u), x)
    return u == x or FreeQ(u, x)

def QuotientOfLinearsParts(u, x):
    # If u is equivalent to an expression of the form (a+b*x)/(c+d*x), QuotientOfLinearsParts[u,x]
    #	returns the list {a, b, c, d}.
    if LinearQ(u, x):
        return [Coefficient(u, x, 0), Coefficient(u, x, 1), 1, 0]
    elif PowerQ(u):
        if Numerator(u) == 1:
            u = Denominator(u)
            r = QuotientOfLinearsParts(u, x)
            return [r[2], r[3], r[0], r[1]]
    elif SumQ(u):
        a = First(u)
        if FreeQ(a, x):
            u = Rest(u)
            r = QuotientOfLinearsParts(u, x)
            return [r[0] + a*r[2], r[1] + a*r[3], r[2], r[3]]
    elif ProductQ(u):
        a = First(u)
        if FreeQ(a, x):
            r = QuotientOfLinearsParts(Rest(u), x)
            return [a*r[0], a*r[1], r[2], r[3]]
        a = Numerator(u)
        d = Denominator(u)
        if LinearQ(a, x) and LinearQ(d, x):
            return [Coefficient(a, x, 0), Coefficient(a, x, 1), Coefficient(d, x, 0), Coefficient(d, x, 1)]
    elif u == x:
        return [0, 1, 1, 0]
    elif FreeQ(u, x):
        return [u, 0, 1, 0]
    return [u, 0, 1, 0]

def QuotientOfLinearsQ(u, x):
    # (*QuotientOfLinearsQ[u,x] returns True iff u is equivalent to an expression of the form (a+b x)/(c+d x) where b!=0 and d!=0.*)
    if ListQ(u):
        for i in u:
            if not QuotientOfLinearsQ(i, x):
                return False
        return True
    q = QuotientOfLinearsParts(u, x)
    return QuotientOfLinearsP(u, x) and NonzeroQ(q[1]) and NonzeroQ(q[3])

def Flatten(l):
    return flatten(l)

def Sort(u, r=False):
    return sorted(u, key=lambda x: x.sort_key(), reverse=r)

# (*Definition: A number is absurd if it is a rational number, a positive rational number raised to a fractional power, or a product of absurd numbers.*)
def AbsurdNumberQ(u):
    # (* AbsurdNumberQ[u] returns True if u is an absurd number, else it returns False. *)
    if PowerQ(u):
        v = u.exp
        u = u.base
        return RationalQ(u) and u > 0 and FractionQ(v)
    elif ProductQ(u):
        return all(AbsurdNumberQ(i) for i in u.args)
    return RationalQ(u)

def AbsurdNumberFactors(u):
    # (* AbsurdNumberFactors[u] returns the product of the factors of u that are absurd numbers. *)
    if AbsurdNumberQ(u):
        return u
    elif ProductQ(u):
        result = S(1)
        for i in u.args:
            if AbsurdNumberQ(i):
                result *= i
        return result
    return NumericFactor(u)

def NonabsurdNumberFactors(u):
    # (* NonabsurdNumberFactors[u] returns the product of the factors of u that are not absurd numbers. *)
    if AbsurdNumberQ(u):
        return S(1)
    elif ProductQ(u):
        result = 1
        for i in u.args:
            result *= NonabsurdNumberFactors(i)
        return result
    return NonnumericFactors(u)

def SumSimplerAuxQ(u, v):
    if SumQ(v):
        return (RationalQ(First(v)) or SumSimplerAuxQ(u,First(v))) and (RationalQ(Rest(v)) or SumSimplerAuxQ(u,Rest(v)))
    elif SumQ(u):
        return SumSimplerAuxQ(First(u), v) or SumSimplerAuxQ(Rest(u), v)
    else:
        return v!=0 and NonnumericFactors(u)==NonnumericFactors(v) and (NumericFactor(u)/NumericFactor(v)<-1/2 or NumericFactor(u)/NumericFactor(v)==-1/2 and NumericFactor(u)<0)

def SumSimplerQ(u, v):
    # (* If u+v is simpler than u, SumSimplerQ[u,v] returns True, else it returns False. *)
    # (* If for every term w of v there is a term of u equal to n*w where n<-1/2, u+v will be simpler than u. *)
    if RationalQ(u, v):
        if v==0:
            return False
        elif v > 0:
            return u < -1
        else:
            return u >= -v
    else:
        return SumSimplerAuxQ(Expand(u), Expand(v))

def Prepend(l1, l2):
    if not isinstance(l2, list):
        return [l2] + l1
    return l2 + l1

def Drop(lst, n):
    if isinstance(lst, list):
        if isinstance(n, list):
            lst = lst[:(n[0]-1)] + lst[n[1]:]
        elif n > 0:
            lst = lst[n:]
        elif n < 0:
            lst = lst[:-n]
        else:
            return lst
        return lst
    return lst.func(*[i for i in Drop(list(lst.args), n)])

def CombineExponents(lst):
    if Length(lst) < 2:
        return lst
    elif lst[0][0] == lst[1][0]:
        return CombineExponents(Prepend(Drop(lst,2),[lst[0][0], lst[0][1] + lst[1][1]]))
    return Prepend(CombineExponents(Rest(lst)), First(lst))

def FactorInteger(n, l=None):
    return sorted(factorint(n, limit=l).items())

def FactorAbsurdNumber(m):
    # (* m must be an absurd number.  FactorAbsurdNumber[m] returns the prime factorization of m *)
    # (* as list of base-degree pairs where the bases are prime numbers and the degrees are rational. *)
    if RationalQ(m):
        return FactorInteger(m)
    elif PowerQ(m):
        r = FactorInteger(m.base)
        return [r[0], r[1]*m.exp]
    #CombineExponents[Sort[Flatten[Map[FactorAbsurdNumber,Apply[List,m]],1], Function[i1[[1]]<i2[[1]]]]]]]
    return CombineExponents

def SubstForInverseFunction(*args):
    # (* SubstForInverseFunction[u,v,w,x] returns u with subexpressions equal to v replaced by x
    # and x replaced by w. *)
    if len(args) == 3:
        u, v, x = args[0], args[1], args[2]
        return SubstForInverseFunction(u, v, (-Coefficient(v.args[0], x, 0) + InverseFunction(Head(v))(x))/Coefficient(v.args[0], x, 1), x)
    elif len(args) == 4:
        u, v, w, x = args[0], args[1], args[2], args[3]
        if AtomQ(u):
            if u == x:
                return w
            return u
        elif Head(u) == Head(v) and ZeroQ(u.args[0] - v.args[0]):
            return x
        res = [SubstForInverseFunction(i, v, w, x) for i in u.args]
        return u.func(*res)

def SubstForFractionalPower(u, v, n, w, x):
    # (* SubstForFractionalPower[u,v,n,w,x] returns u with subexpressions equal to v^(m/n) replaced
    # by x^m and x replaced by w. *)
    if AtomQ(u):
        if u == x:
            return w
        return u
    elif FractionalPowerQ(u) and ZeroQ(u.args[0] - v):
        return x**(n*u.args[1])
    res = [SubstForFractionalPower(i, v, n, w, x) for i in u.args]
    return u.func(*res)

def SubstForFractionalPowerOfQuotientOfLinears(u, x):
    # (* If u has a subexpression of the form ((a+b*x)/(c+d*x))^(m/n) where m and n>1 are integers,
    # SubstForFractionalPowerOfQuotientOfLinears[u,x] returns the list {v,n,(a+b*x)/(c+d*x),b*c-a*d} where v is u
    # with subexpressions of the form ((a+b*x)/(c+d*x))^(m/n) replaced by x^m and x replaced
    lst = FractionalPowerOfQuotientOfLinears(u, 1, False, x)
    if AtomQ(lst) or AtomQ(lst[1]):
        return False
    n = lst[0]
    tmp = lst[1]
    lst = QuotientOfLinearsParts(tmp, x)
    a, b, c, d = lst[0], lst[1], lst[2], lst[3]
    if ZeroQ(d):
        return False
    lst = Simplify(x**(n - 1)*SubstForFractionalPower(u, tmp, n, (-a + c*x**n)/(b - d*x**n), x)/(b - d*x**n)**2)
    return [NonfreeFactors(lst, x), n, tmp, FreeFactors(lst, x)*(b*c - a*d)]

def FractionalPowerOfQuotientOfLinears(u, n, v, x):
    # (* If u has a subexpression of the form ((a+b*x)/(c+d*x))^(m/n),
    # FractionalPowerOfQuotientOfLinears[u,1,False,x] returns {n,(a+b*x)/(c+d*x)}; else it returns False. *)
    if AtomQ(u) or FreeQ(u, x):
        return [n, v]
    elif CalculusQ(u):
        return False
    elif FractionalPowerQ(u) and QuotientOfLinearsQ(u.args[0], x) and Not(LinearQ(u.args[0], x)) and (FalseQ(v) or ZeroQ(u.args[0] - v)):
        return [LCM(Denominator(u.exp), n), u.base]
    lst = [n, v]
    for i in u.args:
        lst = FractionalPowerOfQuotientOfLinears(i, lst[0], lst[1],x)
        if AtomQ(lst):
            return False
    return lst

def SubstForFractionalPowerQ(u, v, x):
    # (* If the substitution x=v^(1/n) will not complicate algebraic subexpressions of u,
    # SubstForFractionalPowerQ[u,v,x] returns True; else it returns False. *)
    if AtomQ(u) or FreeQ(u, x):
        return True
    elif FractionalPowerQ(u):
        return SubstForFractionalPowerAuxQ(u, v, x)
    return all(SubstForFractionalPowerQ(i, v, x) for i in u.args)

def SubstForFractionalPowerAuxQ(u, v, x):
    if AtomQ(u):
        return False
    elif FractionalPowerQ(u) and ZeroQ(u.args[0] - v):
        return True
    return any(SubstForFractionalPowerAuxQ(i, v, x) for i in u.args)

def FractionalPowerOfSquareQ(u):
    # (* If a subexpression of u is of the form ((v+w)^2)^n where n is a fraction, *)
    # (* FractionalPowerOfSquareQ[u] returns (v+w)^2; else it returns False. *)
    if AtomQ(u):
        return False
    elif FractionalPowerQ(u):
        a_ = Wild('a')
        b_ = Wild('b')
        c_ = Wild('c')
        match = u.match(a_*(b_ + c_)**S(2))
        if match:
            keys = [a_, b_, c_]
            if len(keys) == len(match):
                a, b, c = tuple(match[i] for i in keys)
                if NonsumQ(a):
                    return (b + c)**S(2)
    for i in u.args:
        tmp = FractionalPowerOfSquareQ(i)
        if Not(FalseQ(tmp)):
            return tmp
    return False

def FractionalPowerSubexpressionQ(u, v, w):
    # (* If a subexpression of u is of the form w^n where n is a fraction but not equal to v, *)
    # (* FractionalPowerSubexpressionQ[u,v,w] returns True; else it returns False. *)
    if AtomQ(u):
        return False
    elif FractionalPowerQ(u) and PositiveQ(u.args[0]/w):
        return Not(u.args[0] == v) and LeafCount(w) < 3*LeafCount(v)
    for i in u.args:
        if FractionalPowerSubexpressionQ(i, v, w):
            return True
    return False

def Apply(f, lst):
    return f(*lst)

def FactorNumericGcd(u):
    # (* FactorNumericGcd[u] returns u with the gcd of the numeric coefficients of terms of sums factored out. *)
    if PowerQ(u):
        if RationalQ(u.exp):
            return FactorNumericGcd(u.base)**u.exp
    elif ProductQ(u):
        res = [FactorNumericGcd(i) for i in u.args]
        return Mul(*res)
    elif SumQ(u):
        g = GCD(*[NumericFactor(i) for i in u.args])
        r = Add(*[i/g for i in u.args])
        return g*r
    return u

def MergeableFactorQ(bas, deg, v):
    # (* MergeableFactorQ[bas,deg,v] returns True iff bas equals the base of a factor of v or bas is a factor of every term of v. *)
    if bas == v:
        return RationalQ(deg + S(1)) and (deg + 1>=0 or RationalQ(deg) and deg>0)
    elif PowerQ(v):
        if bas == v.base:
            return RationalQ(deg+v.exp) and (deg+v.exp>=0 or RationalQ(deg) and deg>0)
        return SumQ(v.base) and IntegerQ(v.exp) and (Not(IntegerQ(deg) or IntegerQ(deg/v.exp))) and MergeableFactorQ(bas, deg/v.exp, v.base)
    elif ProductQ(v):
        return MergeableFactorQ(base, deg, First(v)) or MergeableFactorQ(bas, deg, Rest(v))
    return SumQ(v) and MergeableFactorQ(bas, deg, First(v)) and MergeableFactorQ(bas, deg, Rest(v))

def MergeFactor(bas, deg, v):
    # (* If MergeableFactorQ[bas,deg,v], MergeFactor[bas,deg,v] return the product of bas^deg and v,
    # but with bas^deg merged into the factor of v whose base equals bas. *)
    if bas == v:
        return bas**(deg + 1)
    elif PowerQ(v):
        if bas == v.base:
            return bas**(deg + b.exp)
        return MergeFactor(bas, deg/b.exp, v.base**v.exp)
    elif ProductQ(v):
        if MergeableFactorQ(bas, deg, First(v)):
            return MergeFactor(bas, deg, First(v))*Rest(v)
        return First(v)*MergeFactor(bas, deg, Rest(v))
    return MergeFactor(bas, deg, First(v) + MergeFactor(bas, deg, Rest(v)))

def MergeFactors(u, v):
    # (* MergeFactors[u,v] returns the product of u and v, but with the mergeable factors of u merged into v. *)
    if ProductQ(u):
        return MergeFactors(Rest(u), MergeFactors(First(u), v))
    elif PowerQ(u):
        if MergeableFactorQ(u.base, u.exp, v):
            return MergeFactor(u.base, u.exp, v)
        elif RationalQ(u.exp) and u.exp < -1 and MergeableFactorQ(u.base, -S(1), v):
            return MergeFactors(u.base**(u.exp + 1), MergeFactor(u.base, -S(1), v))
        return u*v
    elif MergeableFactorQ(u, S(1), v):
        return MergeFactor(u, S(1), v)
    return u*v

def TrigSimplifyQ(u):
    return u == TrigSimplify(u)

def TrigSimplify(u):
    return trigsimp(u)

def Order(expr1, expr2):
    if expr1 == expr2:
        return 0
    elif expr1.sort_key() > expr2.sort_key():
        return -1
    return 1

def FactorOrder(u, v):
    if u == 1:
        if v == 1:
            return 0
        return -1
    elif v == 1:
        return 1
    return Order(u, v)

def Smallest(num1, num2=None):
    if num2 == None:
        lst = num1
        num = lst[0]
        for i in Rest(lst):
            num = Smallest(num, i)
        return num
    return Min(num1, num2)


def OrderedQ(l):
    return l == Sort(l)

def MinimumDegree(deg1, deg2):
    if RationalQ(deg1):
        if RationalQ(deg2):
            return Min(deg1, deg2)
        return deg1
    elif RationalQ(deg2):
        return deg2

    deg = Simplify(deg1- deg2)

    if RationalQ(deg):
        if deg > 0:
            return deg2
        return deg1
    elif OrderedQ([deg1, deg2]):
        return deg1
    return deg2

def PositiveFactors(u):
    # (* PositiveFactors[u] returns the positive factors of u *)
    if ZeroQ(u):
        return S(1)
    elif RationalQ(u):
        return Abs(u)
    elif PositiveQ(u):
        return u
    elif ProductQ(u):
        res = 1
        for i in u.args:
            res *= PositiveFactors(i)
        return res
    return 1

def Sign(u):
    return sign(u)

def NonpositiveFactors(u):
    # (* NonpositiveFactors[u] returns the nonpositive factors of u *)
    if ZeroQ(u):
        return u
    elif RationalQ(u):
        return Sign(u)
    elif PositiveQ(u):
        return S(1)
    elif ProductQ(u):
        res = S(1)
        for i in u.args:
            res *= NonpositiveFactors(i)
        return res
    return u

def PolynomialInAuxQ(u, v, x):
    if u == v:
        return True
    elif AtomQ(u):
        return u != x
    elif PowerQ(u):
        if PowerQ(v):
            if u.base == v.base:
                return PositiveIntegerQ(u.exp/v.exp)
        return PositiveIntegerQ(u.exp) and PolynomialInAuxQ(u.base, v, x)
    elif SumQ(u) or ProductQ(u):
        for i in u.args:
            if Not(PolynomialInAuxQ(i, v, x)):
                return False
        return True
    return False

def PolynomialInQ(u, v, x):
    # If u is a polynomial in v[x], PolynomialInQ[u,v,x] returns True; else it returns False.
    return PolynomialInAuxQ(u, NonfreeFactors(NonfreeTerms(v, x), x), x)

def ExponentInAux(u, v, x):
    if u == v:
        return S(1)
    elif AtomQ(u):
        return S(0)
    elif PowerQ(u):
        if PowerQ(v):
            if u.base == v.base:
                return u.exp/v.exp
        return u.exp*ExponentInAux(u.base, v, x)
    elif ProductQ(u):
        return Add(*[ExponentInAux(i, v, x) for i in u.args])
    return Max(*[ExponentInAux(i, v, x) for i in u.args])

def ExponentIn(u, v, x):
    return ExponentInAux(u, NonfreeFactors(NonfreeTerms(v, x), x), x)

def PolynomialInSubstAux(u, v, x):
    if u == v:
        return x
    elif AtomQ(u):
        return u
    elif PowerQ(u):
        if PowerQ(v):
            if u.base == v.base:
                return x**(u.exp/v.exp)
        return PolynomialInSubstAux(u.base, v, x)**u.exp
    return u.func(*[PolynomialInSubstAux(i, v, x) for i in u.args])

def PolynomialInSubst(u, v, x):
    # If u is a polynomial in v[x], PolynomialInSubst[u,v,x] returns the polynomial u in x.
    w = NonfreeTerms(v, x)
    return ReplaceAll(PolynomialInSubstAux(u, NonfreeFactors(w, x), x), {x: x - FreeTerms(v, x)/FreeFactors(w, x)})

def Distrib(u, v):
    # Distrib[u,v] returns the sum of u times each term of v.
    if SumQ(v):
        return Add(*[u*i for i in v.args])
    return u*v

def DistributeDegree(u, m):
    # DistributeDegree[u,m] returns the product of the factors of u each raised to the mth degree.
    if AtomQ(u):
        return u**m
    elif PowerQ(u):
        return u.base**(u.exp*m)
    elif ProductQ(u):
        return Mul(*[DistributeDegree(i, m) for i in u.args])
    return u**m

def FunctionOfPower(*args):
    # FunctionOfPower[u,x] returns the gcd of the integer degrees of x in u.
    if len(args) == 2:
        return FunctionOfPower(args[0], None, args[1])

    u, n, x = args

    if FreeQ(u, x):
        return n
    elif u == x:
        return S(1)
    elif PowerQ(u):
        if u.base == x and IntegerQ(u.exp):
            if n == None:
                return u.exp
            return GCD(n, u.exp)
    tmp = n
    for i in u.args:
        tmp = FunctionOfPower(i, tmp, x)
    return tmp

def DivideDegreesOfFactors(u, n):
    # DivideDegreesOfFactors[u,n] returns the product of the base of the factors of u raised to
    # the degree of the factors divided by n.
    if ProductQ(u):
        return Mul(*[LeadBase(i)**(LeadDegree(i)/n) for i in u.args])
    return LeadBase(u)**(LeadDegree(u)/n)

def MonomialFactor(u, x):
    # MonomialFactor[u,x] returns the list {n,v} where x^n*v==u and n is free of x.
    if AtomQ(u):
        if u == x:
            return [S(1), S(1)]
        return [S(0), u]
    elif PowerQ(u):
        if IntegerQ(u.exp):
            lst = MonomialFactor(u.base, x)
            return [lst[0]*u.exp, lst[1]**u.exp]
        elif u.base == x and FreeQ(u.exp, x):
            return [u.exp, S(1)]
        return [S(0), u]
    elif ProductQ(u):
        lst1 = MonomialFactor(First(u), x)
        lst2 = MonomialFactor(Rest(u), x)
        return [lst1[0] + lst2[0], lst1[1]*lst2[1]]
    elif SumQ(u):
        lst = [MonomialFactor(i, x) for i in u.args]
        deg = lst[0][0]
        for i in Rest(lst):
            deg = MinimumDegree(deg, i[0])
        if ZeroQ(deg) or RationalQ(deg) and deg < 0:
            return [S(0), u]
        return [deg, Add(*[x**(i[0] - deg)*i[1] for i in lst])]
    return [S(0), u]

def FullSimplify(expr):
    return simplify(expr)

def FunctionOfLinearSubst(u, a, b, x):
    if FreeQ(u, x):
        return u
    elif LinearQ(u, x):
        tmp = Coefficient(u, x, 1)
        if tmp == b:
            tmp = S(1)
        else:
            tmp = tmp/b
        return Coefficient(u, x, S(0)) - a*tmp + tmp*x
    elif PowerQ(u) and FreeQ(u.base[0], x):
        return E**(FullSimplify(FunctionOfLinearSubst(Log(u.base*u.exp, a, b, x))))
    lst = MonomialFactor(u, x)
    if ProductQ(u) and NonzeroQ(lst[0]):
        if RationalQ(LeadFactor(lst[1])) and LeadFactor(lst[1]) < 0:
            return  -FunctionOfLinearSubst(DivideDegreesOfFactors(-lst[1], lst[0])*x, a, b, x)**lst[0]
        return FunctionOfLinearSubst(DivideDegreesOfFactors(lst[1], lst[0])*x, a, b, x)**lst[0]
    return u.func(*[FunctionOfLinearSubst(i, a, b, x) for i in u.args])

def FunctionOfLinear(*args):
    # (* If u (x) is equivalent to an expression of the form f (a+b*x) and not the case that a==0 and
    # b==1, FunctionOfLinear[u,x] returns the list {f (x),a,b}; else it returns False. *)
    if len(args) == 2:
        u, x = args
        lst = FunctionOfLinear(u, False, False, x, False)
        if AtomQ(lst) or FalseQ(lst[0]) or (lst[0] == 0 and lst[1] == 1):
            return False
        return [FunctionOfLinearSubst(u, lst[0], lst[1], x), lst[0], lst[1]]
    u, a, b, x, flag = args
    if FreeQ(u, x):
        return [a, b]
    elif CalculusQ(u):
        return False
    elif LinearQ(u, x):
        if FalseQ(a):
            return [Coefficient(u, x, 0), Coefficient(u, x, 1)]
        lst = CommonFactors([b, Coefficient(u, x, 1)])
        if ZeroQ(Coefficient(u, x, 0)) and Not(flag):
            return [0, lst[0]]
        elif ZeroQ(b*Coefficient(u, x, 0) - a*Coefficient(u, x, 1)):
            return [a/lst[1], lst[0]]
        return [0, 1]
    elif PowerQ(u) and FreeQ(u.args[0], x):
        return FunctionOfLinear(log(u.base)*u.exp, a, b, x, False)
    lst = MonomialFactor(u, x)
    if ProductQ(u) and NonzeroQ(lst[0]):
        if False and IntegerQ(lst[0]) and lst[0] != -1 and FreeQ(lst[1], x):
            if RationalQ(LeadFactor(lst[1])) and LeadFactor(lst[1]) < 0:
                return FunctionOfLinear(DivideDegreesOfFactors(-lst[1], lst[0])*x, a, b, x, False)
            return FunctionOfLinear(DivideDegreesOfFactors(lst[1], lst[0])*x, a, b, x, False)
        return False
    lst = [a, b]
    for i in u.args:
        lst = FunctionOfLinear(i, lst[0], lst[1], x, SumQ(u))
        if AtomQ(lst):
            return False
    return lst

def NormalizeIntegrand(u, x):
    v = NormalizeLeadTermSigns(NormalizeIntegrandAux(u, x))
    if v == NormalizeLeadTermSigns(u):
        return u
    else:
        return v

def NormalizeIntegrandAux(u, x):
    if SumQ(u):
        l = 0
        for i in u.args:
            l += NormalizeIntegrandAux(i, x)
        return l
    if ProductQ(MergeMonomials(u, x)):
        l = 1
        for i in MergeMonomials(u, x).args:
            l *= NormalizeIntegrandFactor(i, x)
        return l
    else:
        return NormalizeIntegrandFactor(MergeMonomials(u, x), x)

def NormalizeIntegrandFactor(u, x):
    if PowerQ(u) and FreeQ(u.args[1], x):
        bas = NormalizeIntegrandFactorBase(u.args[0], x)
        deg = u.args[1]
        if IntegerQ(deg) and SumQ(bas):
            if all(MonomialQ(i, x) for i in bas.args):
                mi = MinimumMonomialExponent(bas, x)
                q = 0
                for i in bas.args:
                    q += Simplify(i/x**mi)
                return x**(mi*deg)*q**deg
            else:
                return bas**deg
        else:
            return bas**deg
    if PowerQ(u) and FreeQ(u.args[0], x):
        return u.args[0]**NormalizeIntegrandFactorBase(u.args[1], x)
    bas = NormalizeIntegrandFactorBase(u, x)
    if SumQ(bas):
        if all(MonomialQ(i, x) for i in bas.args):
            mi = MinimumMonomialExponent(bas, x)
            z = 0
            for j in bas.args:
                z += j/x**mi
            return x**mi*z
        else:
            return bas
    else:
        return bas

def NormalizeIntegrandFactorBase(expr, x):
    m = Wild('m', exclude=[x])
    u = Wild('u')
    match = expr.match(x**m*u)
    if match and SumQ(u):
        l = 0
        for i in u.args:
            l += NormalizeIntegrandFactorBase((x**m*i), x)
        return l
    if BinomialQ(expr, x):
        if BinomialMatchQ(expr, x):
            return expr
        else:
            return ExpandToSum(expr, x)
    elif TrinomialQ(expr, x):
        if TrinomialMatchQ(expr, x):
            return expr
        else:
            return ExpandToSum(expr, x)
    elif ProductQ(expr):
        l = 1
        for i in expr.args:
            l *= NormalizeIntegrandFactor(i, x)
        return l
    elif PolynomialQ(expr, x) and Exponent(expr, x)<=4:
        return ExpandToSum(expr, x)
    elif SumQ(expr):
        w = Wild('w')
        m = Wild('m', exclude=[x])
        v = TogetherSimplify(expr)
        if SumQ(v) or v.match(x**m*w) and SumQ(w) or LeafCount(v)>LeafCount(expr)+2:
            return UnifySum(expr, x)
        else:
            return NormalizeIntegrandFactorBase(v, x)
    else:
        return expr

def NormalizeTogether(u):
    return NormalizeLeadTermSigns(Together(u))

def NormalizeLeadTermSigns(u):
    if ProductQ(u):
        t = 1
        for i in u.args:
            lst = SignOfFactor(i)
            if lst[0] == 1:
                t *= lst[1]
            else:
                t *= AbsorbMinusSign(lst[1])
        return t
    else:
        lst = SignOfFactor(u)
    if lst[0] == 1:
        return lst[1]
    else:
        return AbsorbMinusSign(lst[1])

def AbsorbMinusSign(expr, *x):
    m = Wild('m', exclude=[x])
    u = Wild('u')
    v = Wild('v')
    match = expr.match(u*v**m)
    if match:
        if len(match) == 3:
            if SumQ(match[v]) and OddQ(match[m]):
                return match[u]*(-match[v])**match[m]

    return -expr
# yet todo
def NormalizeSumFactors(u):
    if AtomQ(u):
        return u
    if ProductQ(u):
        k = 1
        for i in u.args:
            k *= NormalizeSumFactors(i)
        return SignOfFactor(k)[0]*SignOfFactor(k)[1]
    else:
        k = 0
        for i in u.args:
            k += NormalizeSumFactors(i)
        return k

def SignOfFactor(u):
    if RationalQ(u) and u<0 or SumQ(u) and NumericFactor(First(u))<0:
        return [-1, -u]
    elif IntegerPowerQ(u) and SumQ(u.args[0]) and NumericFactor(First(u.args[0]))<0:
        return [(-1)**u.args[1], (-u.args[0])**u.args[1]]
    elif ProductQ(u):
        k = 1
        h = 1
        for i in u.args:
            k *= SignOfFactor(i)[0]
            h *= SignOfFactor(i)[1]
        return [k, h]
    else:
        return [1, u]

def NormalizePowerOfLinear(u, x):
    v = FactorSquareFree(u)
    if PowerQ(v) and LinearQ(v.args[0], x) and FreeQ(v.args[1], x):
        return ExpandToSum(v.args[0], x)**v.args[1]
    else:
        return ExpandToSum(v, x)

def SimplifyIntegrand(u, x):
    v = NormalizeLeadTermSigns(NormalizeIntegrandAux(Simplify(u), x))
    if LeafCount(v) < 4/5*LeafCount(u):
        return v
    if v != NormalizeLeadTermSigns(u):
        return v
    else:
        return u

def SimplifyTerm(u, x):
    v = Simplify(u)
    w = Together(v)
    if LeafCount(v) < LeafCount(w):
        return NormalizeIntegrand(v, x)
    else:
        return NormalizeIntegrand(w, x)

def TogetherSimplify(u):
    v = Together(Simplify(Together(u)))
    return FixSimplify(v)

def SmartSimplify(u):
    v = Simplify(u)
    w = factor(v)
    if LeafCount(w)<LeafCount(v):
        v = w
    if Not(FalseQ(w= FractionalPowerOfSquareQ(v))) and FractionalPowerSubexpressionQ(u, w, Expand(w)):
        v = SubstForExpn(v, w, Expand(w))
    else:
        v = FactorNumericGcd(v)

def SubstForExpn(u, v, w):
    if u == v:
        return w
    if AtomQ(u):
        return u
    else:
        k = 0
        for i in u.args:
            k +=  SubstForExpn(i, v, w)
        return k

def ExpandToSum(u, *x):
    if len(x) == 1:
        x = x[0]
        expr = 0
        if S(u).is_polynomial(x):
            for t in Exponent(u, x, List):
                expr += Coeff(u, x, t)*x**t
            return expr
        if BinomialQ(u, x):
            i = BinomialParts(u, x)
            expr += i[0] + i[1]*x**i[2]
            return expr
        if TrinomialQ(u, x):
            i = TrinomialParts(u, x)
            expr += i[0] + i[1]*x**i[3] + i[2]*x**(2*i[3])
            return expr
        if GeneralizedBinomialMatchQ(u, x):
            i = GeneralizedBinomialParts(u, x)
            expr += i[0]*x**i[3] + i[1]*x**i[2]
            return expr
        if GeneralizedTrinomialMatchQ(u, x):
            i = GeneralizedTrinomialParts(u, x)
            expr += i[0]*x**i[4] + i[1]*x**i[3] + i[2]*x**(2*i[3]-i[4])
            return expr
        else:
            #print("Warning: Unrecognized expression for expansion")
            return Expand(u)
    else:
        v = x[0]
        x = x[1]
        w = ExpandToSum(v, x)
        r = NonfreeTerms(w, x)
        if SumQ(r):
            k = u*FreeTerms(w, x)
            for i in r.args:
                k += MergeMonomials(u*i, x)
            return k
        else:
            return u*FreeTerms(w, x) + MergeMonomials(u*r, x)

def UnifySum(u, x):
    if SumQ(u):
        t = 0
        lst = []
        for i in u.args:
            lst += [i]
        for j in UnifyTerms(lst, x):
            t += j
        return t
    else:
        return SimplifyTerm(u, x)

def UnifyTerms(lst, x):
    if lst==[]:
        return lst
    else:
        return UnifyTerm(First(lst), UnifyTerms(Rest(lst), x), x)

def UnifyTerm(term, lst, x):
    if lst==[]:
        return [term]
    tmp = Simplify(First(lst)/term)
    if FreeQ(tmp, x):
        return Prepend(Rest(lst), [(1+tmp)*term])
    else:
        return Prepend(UnifyTerm(term, Rest(lst), x), [First(lst)])

def CalculusQ(u):
    return False

def FunctionOfInverseLinear(*args):
    # (* If u is a function of an inverse linear binomial of the form 1/(a+b*x),
    # FunctionOfInverseLinear[u,x] returns the list {a,b}; else it returns False. *)
    if len(args) == 2:
        u, x = args
        return FunctionOfInverseLinear(u, None, x)
    u, lst, x = args

    if FreeQ(u, x):
        return lst
    elif u == x:
        return False
    elif QuotientOfLinearsQ(u, x):
        tmp = Drop(QuotientOfLinearsParts(u, x), 2)
        if tmp[1] == 0:
            return False
        elif lst == None:
            return tmp
        elif ZeroQ(lst[0]*tmp[1] - lst[1]*tmp[0]):
            return lst
        return False
    elif CalculusQ(u):
        return False
    tmp = lst
    for i in u.args:
        tmp = FunctionOfInverseLinear(i, tmp, x)
        if AtomQ(tmp):
            return False
    return tmp

def PureFunctionOfSinhQ(u, v, x):
    # (* If u is a pure function of Sinh[v] and/or Csch[v], PureFunctionOfSinhQ[u,v,x] returns True;
    # else it returns False. *)
    if AtomQ(u):
        return u != x
    elif CalculusQ(u):
        return False
    elif HyperbolicQ(u) and ZeroQ(u.args[0] - v):
        return SinhQ(u) or CschQ(u)
    for i in u.args:
        if Not(PureFunctionOfSinhQ(i, v, x)):
            return False
    return True

def PureFunctionOfTanhQ(u, v , x):
    # (* If u is a pure function of Tanh[v] and/or Coth[v], PureFunctionOfTanhQ[u,v,x] returns True;
    # else it returns False. *)
    if AtomQ(u):
        return u != x
    elif CalculusQ(u):
        return False
    elif HyperbolicQ(u) and ZeroQ(u.args[0] - v):
        return TanhQ(u) or CothQ(u)
    for i in u.args:
        if Not(PureFunctionOfTanhQ(i, v, x)):
            return False
    return True

def PureFunctionOfCoshQ(u, v, x):
    # (* If u is a pure function of Cosh[v] and/or Sech[v], PureFunctionOfCoshQ[u,v,x] returns True;
    # else it returns False. *)
    if AtomQ(u):
        return u != x
    elif CalculusQ(u):
        return False
    elif HyperbolicQ(u) and ZeroQ(u.args[0] - v):
        return CoshQ(u) or SechQ(u)
    for i in u.args:
        if Not(PureFunctionOfCoshQ(i, v, x)):
            return False
    return True

def IntegerQuotientQ(u, v):
    # (* If u/v is an integer, IntegerQuotientQ[u,v] returns True; else it returns False. *)
    return IntegerQ(Simplify(u/v))

def OddQuotientQ(u, v):
    # (* If u/v is odd, OddQuotientQ[u,v] returns True; else it returns False. *)
    return OddQ(Simplify(u/v))

def EvenQuotientQ(u, v):
    # (* If u/v is even, EvenQuotientQ[u,v] returns True; else it returns False. *)
    return EvenQ(Simplify(u/v))

def FindTrigFactor(func1, func2, u, v, flag):
    # (* If func[w]^m is a factor of u where m is odd and w is an integer multiple of v,
    # FindTrigFactor[func1,func2,u,v,True] returns the list {w,u/func[w]^n}; else it returns False. *)
    # (* If func[w]^m is a factor of u where m is odd and w is an integer multiple of v not equal to v,
    # FindTrigFactor[func1,func2,u,v,False] returns the list {w,u/func[w]^n}; else it returns False. *)
    if u == 1:
        return False
    elif (Head(LeadBase(u)) == func1 or Head(LeadBase(u)) == func2) and OddQ(LeadDegree(u)) and IntegerQuotientQ(LeadBase(u).args[0], v) and (flag or NonzeroQ(LeadBase(u).args[0] - v)):
        return [LeadBase[u].args[0], RemainingFactors(u)]
    lst = FindTrigFactor(func1, func2, RemainingFactors(u), v, flag)
    if AtomQ(lst):
        return False
    return [lst[0], LeadFactor(u)*lst[1]]

def FunctionOfSinhQ(u, v, x):
    # (* If u is a function of Sinh[v], FunctionOfSinhQ[u,v,x] returns True; else it returns False. *)
    if AtomQ(u):
        return u != x
    elif CalculusQ(u):
        return False
    elif HyperbolicQ(u) and IntegerQuotientQ(u.args[0], v):
        if OddQuotientQ(u.args[0], v):
            # (* Basis: If m odd, Sinh[m*v]^n is a function of Sinh[v]. *)
            return SinhQ(u) or CschQ(u)
        # (* Basis: If m even, Cos[m*v]^n is a function of Sinh[v]. *)
        return CoshQ(u) or SechQ(u)
    elif IntegerPowerQ(u) and HyperbolicQ(u.args[0]) and IntegerQuotientQ(u.args[0].args[0], v):
        if EvenQ(u.args[1]):
            # (* Basis: If m integer and n even, Hyper[m*v]^n is a function of Sinh[v]. *)
            return True
        return FunctionOfSinhQ(u.args[0], v, x)
    elif ProductQ(u):
        if CoshQ(u.args[0]) and SinhQ(u.args[1]) and ZeroQ(u.args[0].args[0] - v/2) and ZeroQ(u.args[1].args[0] - v/2):
            return FunctionOfSinhQ(Drop(u, 2), v, x)
        lst = FindTrigFactor(Sinh, Csch, u, v, False)
        if ListQ(lst) and EvenQuotientQ(lst[0], v):
            # (* Basis: If m even and n odd, Sinh[m*v]^n == Cosh[v]*u where u is a function of Sinh[v]. *)
            return FunctionOfSinhQ(Cosh(v)*lst[1], v, x)
        lst = FindTrigFactor(Cosh, Sech, u, v, False)
        if ListQ(lst) and OddQuotientQ(lst[0], v):
            # (* Basis: If m odd and n odd, Cosh[m*v]^n == Cosh[v]*u where u is a function of Sinh[v]. *)
            return FunctionOfSinhQ(Cosh(v)*lst[1], v, x)
        lst = FindTrigFactor(Tanh, Coth, u, v, True)
        if ListQ(lst):
            # (* Basis: If m integer and n odd, Tanh[m*v]^n == Cosh[v]*u where u is a function of Sinh[v]. *)
            return FunctionOfSinhQ(Cosh(v)*lst[1], v, x)
        return all(FunctionOfSinhQ(i, v, x) for i in u.args)
    return all(FunctionOfSinhQ(i, v, x) for i in u.args)

def FunctionOfCoshQ(u, v, x):
    #(* If u is a function of Cosh[v], FunctionOfCoshQ[u,v,x] returns True; else it returns False. *)
    if AtomQ(u):
        return u != x
    elif CalculusQ(u):
        return False
    elif HyperbolicQ(u) and IntegerQuotientQ(u.args[0], v):
        # (* Basis: If m integer, Cosh[m*v]^n is a function of Cosh[v]. *)
        return CoshQ(u) or SechQ(u)
    elif IntegerPowerQ(u) and HyperbolicQ(u.args[0]) and IntegerQuotientQ(u.args[0].args[0], v):
        if EvenQ(u.args[1]):
            # (* Basis: If m integer and n even, Hyper[m*v]^n is a function of Cosh[v]. *)
            return True
        return FunctionOfCoshQ(u.args[0], v, x)
    elif ProductQ(u):
        lst = FindTrigFactor(Sinh, Csch, u, v, False)
        if ListQ(lst):
            # (* Basis: If m integer and n odd, Sinh[m*v]^n == Sinh[v]*u where u is a function of Cosh[v]. *)
            return FunctionOfCoshQ(Sinh(v)*lst[1], v, x)
        lst = FindTrigFactor(Tanh, Coth, u, v, True)
        if ListQ(lst):
            # (* Basis: If m integer and n odd, Tanh[m*v]^n == Sinh[v]*u where u is a function of Cosh[v]. *)
            return FunctionOfCoshQ(Sinh(v)*lst[1], v, x)
        return all(FunctionOfCoshQ(i, v, x) for i in u.args)
    return all(FunctionOfCoshQ(i, v, x) for i in u.args)

def OddHyperbolicPowerQ(u, v, x):
    if SinhQ(u) or CoshQ(u) or SechQ(u) or CschQ(u):
        return OddQuotientQ(u.args[0], v)
    if PowerQ(u):
        return Odd(u.args[1]) and OddHyperbolicPowerQ(u.base, v, x)
    if ProductQ(u):
        if Not(Eq(FreeFactors(u, x), 1)):
            return OddHyperbolicPowerQ(NonfreeFactors(u, x), v, x)
        lst = []
        for i in u.args:
            if Not(FunctionOfTanhQ(i, v, x)):
                lst.append(i)
        if lst == []:
            return True
        return Length(lst)==1 and OddHyperbolicPowerQ(lst[0], v, x)
    if SumQ(u):
        return All(OddHyperbolicPowerQ(i, v, x) for i in u.args)
    return False

def FunctionOfTanhQ(u, v, x):
    #(* If u is a function of the form f[Tanh[v],Coth[v]] where f is independent of x,
    # FunctionOfTanhQ[u,v,x] returns True; else it returns False. *)
    if AtomQ(u):
        return u != x
    elif CalculusQ(u):
        return False
    elif HyperbolicQ(u) and IntegerQuotientQ(u.args[0], v):
        return TanhQ(u) or CothQ(u) or EvenQuotientQ(u.args[0], v)
    elif PowerQ(u):
        if EvenQ(u.args[1]) and HyperbolicQ(u.args[0]) and IntegerQuotientQ(u.args[0].args[0], v):
            return True
        elif EvenQ(u.args[1]) and SumQ(u.args[0]):
            return FunctionOfTanhQ(Expand(u.args[0]**2, v, x))
    if ProductQ(u):
        lst = []
        for i in u.args:
            if Not(FunctionOfTanhQ(i, v, x)):
                lst.append(i)
        if lst == []:
            return True
        return Length(lst)==2 and OddHyperbolicPowerQ(lst[0], v, x) and OddHyperbolicPowerQ(lst[1], v, x)
    return all(FunctionOfTanhQ(i, v, x) for i in u.args)

def FunctionOfTanhWeight(u, v, x):
    # (* u is a function of the form f[Tanh[v],Coth[v]] where f is independent of x.
    # FunctionOfTanhWeight[u,v,x] returns a nonnegative number if u is best considered a function
    # of Tanh[v]; else it returns a negative number. *)
    if AtomQ(u):
        return S(0)
    elif CalculusQ(u):
        return S(0)
    elif HyperbolicQ(u) and IntegerQuotientQ(u.args[0], v):
        if TanhQ(u) and ZeroQ(u.args[0] - v):
            return S(1)
        elif CothQ(u) and ZeroQ(u.args[0] - v):
            return S(-1)
        return S(0)
    elif PowerQ(u):
        if EvenQ(u.exp) and HyperbolicQ(u.base) and IntegerQuotientQ(u.base.args[0], v):
            if TanhQ(u.base) or CoshQ(u.base) or SechQ(u.base):
                return S(1)
            return S(-1)
    if ProductQ(u):
        if all(FunctionOfTanhQ(i, v, x) for i in u.args):
            return Add(*[FunctionOfTanhWeight(i, v, x) for i in u.args])
        return S(0)
    return Add(*[FunctionOfTanhWeight(i, v, x) for i in u.args])

def FunctionOfHyperbolicQ(u, v, x):
    # (* If u (x) is equivalent to a function of the form f (Sinh[v],Cosh[v],Tanh[v],Coth[v],Sech[v],Csch[v])
    # where f is independent of x, FunctionOfHyperbolicQ[u,v,x] returns True; else it returns False. *)
    if AtomQ(u):
        return u != x
    elif CalculusQ(u):
        return False
    elif HyperbolicQ(u) and IntegerQuotientQ(u.args[0], v):
        return True
    return all(FunctionOfHyperbolicQ(i, v, x) for i in u.args)

def SmartNumerator(expr):
    if PowerQ(expr):
        n = expr.exp
        u = expr.base
        if RationalQ(n) and n < 0:
            return SmartDenominator(u**(-n))
    elif ProductQ(expr):
        return Mul(*[SmartNumerator(i) for i in expr.args])
    return Numerator(expr)

def SmartDenominator(expr):
    if PowerQ(expr):
        u = expr.base
        n = expr.exp
        if RationalQ(n) and n < 0:
            return SmartNumerator(u**(-n))
    elif ProductQ(expr):
        return Mul(*[SmartDenominator(i) for i in expr.args])
    return Denominator(expr)

def SubstForAux(u, v, x):
    # (* u is a function of v.  SubstForAux[u,v,x] returns u with v replaced by x. *)
    if u == v:
        return x
    elif AtomQ(u):
        if PowerQ(v):
            if FreeQ(v.exp, x) and ZeroQ(u - v.base):
                return x**(Simplify(1/v.exp))
        return u
    elif PowerQ(u):
        if FreeQ(u.exp, x):
            if ZeroQ(u.base - v):
                return x**(u.exp)
            elif PowerQ(v):
                if FreeQ(v.exp, x) and ZeroQ(u.base - v.base):
                    return x**Simplify(u.exp/v.exp)
            return SubstForAux(u.base, v, x)**u.exp
    if ProductQ(u) and Not(EqQ(FreeFactors(u, x), 1)):
        return FreeFactors(u, x)*SubstForAux(NonfreeFactors(u, x), v, x)
    elif ProductQ(u) and ProductQ(v):
        return SubstForAux(First(u), First(v), x)
    return u.func(*[SubstForAux(i, v, x) for i in u.args])

def ActivateTrig(u):
    return u

def ExpandTrig(*args):
    if len(args) == 2:
        u, x = args
        return ActivateTrig(ExpandIntegrand(u, x))
    u, v, x = args
    w = ExpandTrig(v, x)
    z = ActivateTrig(u)
    if SumQ(w):
        return w.func(*[z*i for i in w.args])
    return z*w

def TrigExpand(u):
    return expand_trig(u)

def SubstForTrig(u, sin , cos, v, x):
    # (* u (v) is an expression of the form f (Sin[v],Cos[v],Tan[v],Cot[v],Sec[v],Csc[v]). *)
    # (* SubstForTrig[u,sin,cos,v,x] returns the expression f (sin,cos,sin/cos,cos/sin,1/cos,1/sin). *)
    if AtomQ(u):
        return u
    elif TrigQ(u) and IntegerQuotientQ(u.args[0], v):
        if u.args[0] == v or ZeroQ(u.args[0] - v):
            if SinQ(u):
                return sin(x)
            elif CosQ(u):
                return cos(x)
            elif TanQ(u):
                return sin(x)/cos(x)
            elif CotQ(u):
                return cos(x)/sin(x)
            elif SecQ(u):
                return 1/cos(x)
            return 1/sin(x)
        r = ReplaceAll(TrigExpand(Head(u)(Simplify(u.args[0]/v*x))), {x: v})
        return r.func(*[SubstForTrig(i, sin, cos, v, x) for i in r.args])
    if ProductQ(u) and CosQ(u.args[0]) and SinQ(u.args[1]) and ZeroQ(u.args[0].args[0] - v/2) and ZeroQ(u.args[1].args[0] - v/2):
        return sin(x)/2*SubstForTrig(Drop(u, 2), sin, cos, v, x)
    return u.func(*[SubstForTrig(i, sin, cos, v, x) for i in u.args])

def SubstForHyperbolic(u, sinh, cosh, v, x):
    # (* u (v) is an expression of the form f (Sinh[v],Cosh[v],Tanh[v],Coth[v],Sech[v],Csch[v]). *)
    # (* SubstForHyperbolic[u,sinh,cosh,v,x] returns the expression
    # f (sinh,cosh,sinh/cosh,cosh/sinh,1/cosh,1/sinh). *)
    if AtomQ(u):
        return u
    elif HyperbolicQ(u) and IntegerQuotientQ(u.args[0], v):
        if u.args[0] == v or ZeroQ(u.args[0] - v):
            if SinhQ(u):
                return sinh(x)
            elif CoshQ(u):
                return cosh(x)
            elif TanhQ(u):
                return sinh(x)/cosh(x)
            elif CothQ(u):
                return cosh(x)/sinh(x)
            if SechQ(u):
                return 1/cosh(x)
            return 1/sinh(x)
        r = ReplaceAll(TrigExpand(Head(u)(Simplify(u.args[0]/v)*x)), {x: v})
        return r.func(*[SubstForHyperbolic(i, sinh, cosh, v, x) for i in r.args])
    elif ProductQ(u) and CoshQ(u.args[0]) and SinhQ(u.args[1]) and ZeroQ(u.args[0].args[0] - v/2) and ZeroQ(u.args[1].args[0] - v/2):
        return sinh(x)/2*SubstForHyperbolic(Drop(u, 2), sinh, cosh, v, x)
    return u.func(*[SubstForHyperbolic(i, sinh, cosh, v, x) for i in u.args])

def InertTrigFreeQ(u):
    return FreeQ(u, sin) and FreeQ(u, cos) and FreeQ(u, tan) and FreeQ(u, cot) and FreeQ(u, sec) and FreeQ(u, csc)

def LCM(a, b):
    return lcm(a, b)

def SubstForFractionalPowerOfLinear(u, x):
    # (* If u has a subexpression of the form (a+b*x)^(m/n) where m and n>1 are integers,
    # SubstForFractionalPowerOfLinear[u,x] returns the list {v,n,a+b*x,1/b} where v is u
    # with subexpressions of the form (a+b*x)^(m/n) replaced by x^m and x replaced
    # by -a/b+x^n/b, and all times x^(n-1); else it returns False. *)
    lst = FractionalPowerOfLinear(u, S(1), False, x)
    if AtomQ(lst) or FalseQ(lst[1]):
        return False
    n = lst[0]
    a = Coefficient(lst[1], x, 0)
    b = Coefficient(lst[1], x, 1)
    tmp = Simplify(x**(n-1)*SubstForFractionalPower(u, lst[1], n, -a/b + x**n/b, x))
    return [NonfreeFactors(tmp, x), n, lst[1], FreeFactors(tmp, x)/b]

def FractionalPowerOfLinear(u, n, v, x):
    # (* If u has a subexpression of the form (a+b*x)^(m/n),
    # FractionalPowerOfLinear[u,1,False,x] returns {n,a+b*x}; else it returns False. *)
    if AtomQ(u) or FreeQ(u, x):
        return [n, v]
    elif CalculusQ(u):
        return False
    elif FractionalPowerQ(u) and LinearQ(u.args[0], x) and (FalseQ(v) or ZeroQ(u.args[0] - v)):
        return [LCM(Denominator(u.exp), n), u.base]
    lst = [n, v]
    for i in u.args:
        lst = FractionalPowerOfLinear(i, lst[0], lst[1], x)
        if AtomQ(lst):
            return False
    return lst

def InverseFunctionOfLinear(u, x):
    # (* If u has a subexpression of the form g[a+b*x] where g is an inverse function,
    # InverseFunctionOfLinear[u,x] returns g[a+b*x]; else it returns False. *)
    if AtomQ(u) or CalculusQ(u) or FreeQ(u, x):
        return False
    elif InverseFunctionQ(u) and LinearQ(u.args[0], x):
        return u
    for i in u.args:
        tmp = InverseFunctionOfLinear(i, x)
        if Not(AtomQ(tmp)):
            return tmp
    return False

def InertTrigQ(*args):
    if len(args) == 1:
        f = args[0]
        l = [sin,cos,tan,cot,sec,csc]
        return any(Head(f) == i for i in l)
    elif len(args) == 2:
        f, g = args
        if f == g:
            return InertTrigQ(f)
        return InertReciprocalQ(f, g) or InertReciprocalQ(g, f)
    else:
        f, g, h = args
        return InertTrigQ(g, f) and InertTrigQ(g, h)

def InertReciprocalQ(f, g):
    return (f.func == sin and g.func == csc) or (f.func == cos and g.func == sec) or (f.func == tan and g.func == cot)

def ActivateTrig(u):
    return u

def DeactivateTrig(u, x):
    # (* u is a function of trig functions of a linear function of x. *)
    # (* DeactivateTrig[u,x] returns u with the trig functions replaced with inert trig functions. *)
    return FixInertTrigFunction(DeactivateTrigAux(u, x), x)

def FixInertTrigFunction(u, x):
    return u

def DeactivateTrigAux(u, x):
    if AtomQ(u):
        return u
    elif TrigQ(u) and LinearQ(u.args[0], x):
        v = ExpandToSum(u.args[0], x)

def DeactivateTrigAux(u, x):
    if AtomQ(u):
        return u
    elif TrigQ(u) and LinearQ(u.args[0], x):
        v = ExpandToSum(u.args[0], x)
        if SinQ(u):
            return sin(v)
        elif CosQ(u):
            return cos(v)
        elif TanQ(u):
            return tan(u)
        elif CotQ(u):
            return cot(v)
        elif SecQ(u):
            return csc(v)
        return csc(v)
    elif HyperbolicQ(u) and LinearQ(u.args[0], x):
        v = ExpandToSum(I*u.args[0], x)
        if SinhQ(u):
            return -I*sin(v)
        elif CoshQ(u):
            return cos(v)
        elif TanhQ(u):
            return -I*tan(v)
        elif CothQ(u):
            I*cot(v)
        elif SechQ(u):
            return sec(v)
        return I*csc(v)
    return u.func(*[DeactivateTrigAux(i, x) for i in u.args])

def PowerOfInertTrigSumQ(u, func, x):
    p_ = Wild('p', exclude=[x])
    q_ = Wild('q', exclude=[x])
    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x])
    c_ = Wild('c', exclude=[x])
    d_ = Wild('d', exclude=[x])
    n_ = Wild('n', exclude=[x])
    w_ = Wild('w')

    pattern = (a_ + b_*(c_*func(w_))**p_)**n_
    match = u.match(pattern)
    if match:
        keys = [a_, b_, c_, n_, p_, w_]
        if len(keys) == len(match):
            return True

    pattern = (a_ + b_*(d_*func(w_))**p_ + c_*(d_*func(w_))**q_)**n_
    match = u.match(pattern)
    if match:
        keys = [a_, b_, c_, d_, n_, p_, q_, w_]
        if len(keys) == len(match):
            return True
    return False

def PiecewiseLinearQ(*args):
    # (* If the derivative of u wrt x is a constant wrt x, PiecewiseLinearQ[u,x] returns True;
    # else it returns False. *)
    if len(args) == 3:
        u, v, x = args
        return PiecewiseLinearQ(u, x) and PiecewiseLinearQ(v, x)

    u, x = args
    if LinearQ(u, x):
        return True

    c_ = Wild('c', exclude=[x])
    F_ = Wild('F', exclude=[x])
    v_ = Wild('v')
    match = u.match(log(c_*F_**v_))
    if match:
        if len(match) == 3:
            if LinearQ(match[v_], x):
                return True

    try:
        F = type(u)
        G = type(u.args[0])
        v = u.args[0].args[0]
        if LinearQ(v, x):
            if MemberQ([[atanh, tanh], [atanh, coth], [acoth, coth], [acoth, tanh], [atan, tan], [atan, cot], [acot, cot], [acot, tan]], [F, G]):
                return True
    except:
        pass
    return False

def KnownTrigIntegrandQ(lst, u, x):
    if u == 1:
        return True
    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x])
    func_ = Wild('func')
    m_ = Wild('m', exclude=[x])
    A_ = Wild('A', exclude=[x])
    B_ = Wild('B', exclude=[x])
    C_ = Wild('C', exclude=[x])
    F_ = Wild('F')

    match = u.match((a_ + b_*func_)**m_)
    if match:
        func = match[func_]
        if LinearQ(func.args[0], x) and MemberQ(lst, func.func):
            return True

    match = u.match((a_ + b_*func_)**m_*(A_ + B_*func_))
    if match:
        func = match[func_]
        if LinearQ(func.args[0], x) and MemberQ(lst, func.func):
            return True

    match = u.match(A_ + C_*func_**2)
    if match:
        func = match[func_]
        if LinearQ(func.args[0], x) and MemberQ(lst, func.func):
            return True

    match = u.match(A_ + B_*func_ + C_*func_**2)
    if match:
        func = match[func_]
        if LinearQ(func.args[0], x) and MemberQ(lst, func.func):
            return True

    match = u.match((a_ + b_*func_)**m_*(A_ + C_*func_**2))
    if match:
        func = match[func_]
        if LinearQ(func.args[0], x) and MemberQ(lst, func.func):
            return True

    match = u.match((a_ + b_*func_)**m_*(A_ + B_*func_ + C_*func_**2))
    if match:
        func = match[func_]
        if LinearQ(func.args[0], x) and MemberQ(lst, func.func):
            return True

    return False

def KnownSineIntegrandQ(u, x):
    return KnownTrigIntegrandQ([sin, cos], u, x)

def KnownTangentIntegrandQ(u, x):
    return KnownTrigIntegrandQ([tan], u, x)

def KnownCotangentIntegrandQ(u, x):
    return KnownTrigIntegrandQ([cot], u, x)

def KnownSecantIntegrandQ(u, x):
    return KnownTrigIntegrandQ([sec, csc], u, x)

def ExpandTrigExpand(u, F, v, m, n, x):
    w = Expand(TrigExpand(F.subs(x, n*x))**m, x).subs(x, v)
    if SumQ(w):
        t = 0
        for i in w.args:
            t += u*i
        return t
    else:
        return u*w

def ExpandTrigReduce(u, v, x):
    w = ExpandTrigReduce(v, x)
    if SumQ(w):
        t = 0
        for i in w.args:
            t += u*i
        return t
    else:
            return u*w

def TryPureTanSubst(u, x):
    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x])
    c_ = Wild('c', exclude=[x])
    G_ = Wild('G')

    F = u.func
    try:
        if MemberQ([atan, acot, atanh, acoth], F):
            match = u.args[0].match(c_*(a_ + b_*G_))
            if match:
                if len(match) == 4:
                    G = match[G_]
                    if MemberQ([tan, cot, tanh, coth], G.func):
                        if LinearQ(G.args[0], x):
                            return True
    except:
        pass

    return False

def TryTanhSubst(u, x):
    if u.func == log:
        return False
    elif not FalseQ(FunctionOfLinear(u, x)):
        return False

    a_ = Wild('a', exclude=[x])
    m_ = Wild('m', exclude=[x])
    p_ = Wild('p', exclude=[x])
    r_, s_, t_, n_, b_, f_, g_ = map(Wild, 'rstnbfg')

    match = u.match(r_*(s_ + t_)**n_)
    if match:
        if len(match) == 4:
            r, s, t, n = [match[i] for i in [r_, s_, t_, n_]]
            if IntegerQ(n) and PositiveQ(n):
                return False

    match = u.match(1/(a_ + b_*f_**n_))
    if match:
        if len(match) == 4:
            a, b, f, n = [match[i] for i in [a_, b_, f_, n_]]
            if SinhCoshQ(f) and IntegerQ(n) and n > 2:
                return False

    match = u.match(f_*g_)
    if match:
        if len(match) == 2:
            f, g = match[f_], match[g_]
            if SinhCoshQ(f) and SinhCoshQ(g):
                if IntegersQ(f.args[0]/x, g.args[0]/x):
                    return False

    match = u.match(r_*(a_*s_**m_)**p_)
    if match:
        if len(match) == 5:
            r, a, s, m, p = [match[i] for i in [r_, a_, s_, m_, p_]]
            if Not(m==2 and (s == Sech(x) or s == Csch(x))):
                return False

    if u != ExpandIntegrand(u, x):
        return False

    return True

def TryPureTanhSubst(u, x):
    F = u.func
    a_ = Wild('a', exclude=[x])
    G_ = Wild('G')

    if F == log:
        return False

    match = u.args[0].match(a_*G_)
    if match and len(match) == 2:
        G = match[G_].func
        if MemberQ([atanh, acoth], F) and MemberQ([tanh, coth], G):
            return False

    if u != ExpandIntegrand(u, x):
        return False

    return True

def AbsurdNumberGCD(*seq):
    # (* m, n, ... must be absurd numbers.  AbsurdNumberGCD[m,n,...] returns the gcd of m, n, ... *)
    lst = list(seq)
    if Length(lst) == 1:
        return First(lst)
    return AbsurdNumberGCDList(FactorAbsurdNumber(First(lst)), FactorAbsurdNumber(AbsurdNumberGCD(*Rest(lst))))

def AbsurdNumberGCDList(lst1, lst2):
    # (* lst1 and lst2 must be absurd number prime factorization lists. *)
    # (* AbsurdNumberGCDList[lst1,lst2] returns the gcd of the absurd numbers represented by lst1 and lst2. *)
    if lst1 == []:
        return Mul(*[i[0]**Min(i[1],0) for i in lst2])
    elif lst2 == []:
        return Mul(*[i[0]**Min(i[1],0) for i in lst1])
    elif lst1[0][0] == lst2[0][0]:
        if lst1[0][1] <= lst2[0][1]:
            return lst1[0][0]**lst1[0][1]*AbsurdNumberGCDList(Rest(lst1), Rest(lst2))
        return lst1[0][0]**lst2[0][1]*AbsurdNumberGCDList(Rest(lst1), Rest(lst2))
    elif lst1[0][0] < lst2[0][0]:
        if lst1[0][1] < 0:
            return lst1[0][0]**lst1[0][1]*AbsurdNumberGCDList(Rest(lst1), lst2)
        return AbsurdNumberGCDList(Rest(lst1), lst2)
    elif lst2[0][1] < 0:
        return lst2[0][0]**lst2[0][1]*AbsurdNumberGCDList(lst1, Rest(lst2))
    return AbsurdNumberGCDList(lst1, Rest(lst2))

def Map2(f, lst1, lst2):
    result = []
    for i in range(0, len(lst1)):
        result.append(f(lst1[i], lst2[i]))
    return result

def ConstantFactor(u, x):
    # (* ConstantFactor[u,x] returns a 2-element list of the factors of u[x] free of x and the
    # factors not free of u[x].  Common constant factors of the terms of sums are also collected. *)
    if FreeQ(u, x):
        return [u, S(1)]
    elif AtomQ(u):
        return [S(1), u]
    elif PowerQ(u):
        if FreeQ(u.exp, x):
            lst = ConstantFactor(u.base, x)
            if IntegerQ(u.exp):
                return [lst[0]**u.exp, lst[1]**u.exp]
            tmp = PositiveFactors(lst[0])
            if tmp == 1:
                return [S(1), u]
            return [tmp**u.exp, (NonpositiveFactors(lst[0])*lst[1])**u.exp]
    elif ProductQ(u):
        lst = [ConstantFactor(i, x) for i in u.args]
        return [Mul(*[First(i) for i in lst]), Mul(*[i[1] for i in lst])]
    elif SumQ(u):
        lst1 = [ConstantFactor(i, x) for i in u.args]
        if SameQ(*[i[1] for i in lst1]):
            return [Add(*[i[0] for i in lst]), lst1[0][1]]
        lst2 = CommonFactors([First(i) for i in lst1])
        return [First(lst2), Add(*Map2(Mul, Rest(lst2), [i[1] for i in lst1]))]
    return [S(1), u]

def SameQ(*args):
    for i in range(0, len(args) - 1):
        if args[i] != args[i+1]:
            return False
    return True

def ReplacePart(lst, a, b):
    lst[b] = a
    return lst

def CommonFactors(lst):
    # (* If lst is a list of n terms, CommonFactors[lst] returns a n+1-element list whose first
    # element is the product of the factors common to all terms of lst, and whose remaining
    # elements are quotients of each term divided by the common factor. *)
    lst1 = [NonabsurdNumberFactors(i) for i in lst]
    lst2 = [AbsurdNumberFactors(i) for i in lst]
    num = AbsurdNumberGCD(*lst2)
    common = num
    lst2 = [i/num for i in lst2]
    while (True):
        lst3 = [LeadFactor(i) for i in lst1]

        if SameQ(*lst3):
            common = common*lst3[0]
            lst1 = [RemainingFactors(i) for i in lst1]
        elif (all((LogQ(i) and IntegerQ(First(i)) and First(i) > 0) for i in lst3) and
            all(RationalQ(i) for i in [FullSimplify(j/First(lst3)) for j in lst3])):
            lst4 = [FullSimplify(j/First(lst3)) for j in lst3]
            num = GCD(*lst4)
            common = common*Log((First(lst3)[0])**num)
            lst2 = [lst2[i]*lst4[i]/num for i in range(0, len(lst2))]
            lst1 = [RemainingFactors(i) for i in lst1]
        lst4 = [LeadDegree(i) for i in lst1]
        if SameQ(*[LeadBase(i) for i in lst1]) and RationalQ(*lst4):
            num = Smallest(lst4)
            base = LeadBase(lst1[0])
            if num != 0:
                common = common*base**num
            lst2 = [lst2[i]*base**(lst4[i] - num) for i in range(0, len(lst2))]
            lst1 = [RemainingFactors(i) for i in lst1]
        elif (Length(lst1) == 2 and ZeroQ(LeadBase(lst1[0]) + LeadBase(lst1[1])) and
            NonzeroQ(lst1[0] - 1) and IntegerQ(lst4[0]) and FractionQ(lst4[1])):
            num = Min(lst4)
            base = LeadBase(lst1[1])
            if num != 0:
                common = common*base**num
            lst2 = [lst2[0]*(-1)**lst4[0], lst2[1]]
            lst2 = [lst2[i]*base**(lst4[i] - num) for i in range(0, len(lst2))]
            lst1 = [RemainingFactors(i) for i in lst1]
        elif (Length(lst1) == 2 and ZeroQ(lst1[0] + LeadBase(lst1[1])) and
            NonzeroQ(lst1[1] - 1) and IntegerQ(lst1[1]) and FractionQ(lst4[0])):
            num = Min(lst4)
            base = LeadBase(lst1[0])
            if num != 0:
                common = common*base**num
            lst2 = [lst2[0], lst2[1]*(-1)**lst4[1]]
            lst2 = [lst2[i]*base**(lst4[i] - num) for i in range(0, len(lst2))]
            lst1 = [RemainingFactors(i) for i in lst1]
        else:
            num = MostMainFactorPosition(lst3)
            lst2 = ReplacePart(lst2, lst3[num]*lst2[num], num)
            lst1 = ReplacePart(lst1, RemainingFactors(lst1[num]), num)
        if all(i==1 for i in lst1):
            return Prepend(lst2, common)

def MostMainFactorPosition(lst):
    factor = S(1)
    num = 0
    for i in range(0, Length(lst)):
        if FactorOrder(lst[i], factor) > 0:
            factor = lst[i]
            num = i
    return num

SbaseS, SexponS = None, None
SexponFlagS = False
def FunctionOfExponentialQ(u, x):
    # (* FunctionOfExponentialQ[u,x] returns True iff u is a function of F^v where F is a constant and v is linear in x, *)
    # (* and such an exponential explicitly occurs in u (i.e. not just implicitly in hyperbolic functions). *)
    global SbaseS, SexponS, SexponFlagS
    SbaseS, SexponS = None, None
    SexponFlagS = False
    return FunctionOfExponentialTest(u, x) and SexponFlagS

def FunctionOfExponential(u, x):
    global SbaseS, SexponS, SexponFlagS
    # (* u is a function of F^v where v is linear in x.  FunctionOfExponential[u,x] returns F^v. *)
    SbaseS, SexponS = None, None
    SexponFlagS = False
    FunctionOfExponentialTest(u, x)
    return SbaseS**SexponS

def FunctionOfExponentialFunction(u, x):
    global SbaseS, SexponS, SexponFlagS
    # (* u is a function of F^v where v is linear in x.  FunctionOfExponentialFunction[u,x] returns u with F^v replaced by x. *)
    SbaseS, SexponS = None, None
    SexponFlagS = False
    FunctionOfExponentialTest(u, x)
    return SimplifyIntegrand(FunctionOfExponentialFunctionAux(u, x), x)

def FunctionOfExponentialFunctionAux(u, x):
    # (* u is a function of F^v where v is linear in x, and the fluid variables $base$=F and $expon$=v. *)
    # (* FunctionOfExponentialFunctionAux[u,x] returns u with F^v replaced by x. *)
    global SbaseS, SexponS, SexponFlagS
    if AtomQ(u):
        return u
    elif PowerQ(u) and FreeQ(u.args[0], x) and LinearQ(u.args[1], x):
        if ZeroQ(Coefficient(SexponS, x, 0)):
            return u.base**Coefficient(u.exp, x, 0)*x**FullSimplify(log(u.base)*Coefficient(u.exp, x, 1)/(Log(SbaseS)*Coefficient(SexponS, x, 1)))
        return x**FullSimplify(log(u.base)*Coefficient(u.exp, x, 1)/(Log(SbaseS)*Coefficient(SexponS, x, 1)))
    elif HyperbolicQ(u) and LinearQ(u.args[0], x):
        tmp = x**FullSimplify(Coefficient(u.args[0], x, 1)/(Log(SbaseS)*Coefficient(SexponS, x, 1)))
        if SinhQ(u):
            return tmp/2 - 1/(2*tmp)
        elif CoshQ(u):
            return tmp/2 + 1/(2*tmp)
        elif TanhQ(u):
            return (tmp - 1/tmp)/(tmp + 1/tmp)
        elif CothQ(u):
            return (tmp + 1/tmp)/(tmp - 1/tmp)
        elif SechQ(u):
            return 2/(tmp + 1/tmp)
        return 2/(tmp - 1/tmp)
    elif PowerQ(u) and FreeQ(u.args[0], x) and SumQ(u.args[1]):
        return FunctionOfExponentialFunctionAux(u.base**First(u.exp), x)*FunctionOfExponentialFunctionAux(u.base*Rest(u.exp), x)
    return u.func(*[FunctionOfExponentialFunctionAux(i, x) for i in u.args])

def FunctionOfExponentialTest(u, x):
    # (* FunctionOfExponentialTest[u,x] returns True iff u is a function of F^v where F is a constant and v is linear in x. *)
    # (* Before it is called, the fluid variables $base$ and $expon$ should be set to Null and $exponFlag$ to False. *)
    # (* If u is a function of F^v, $base$ and $expon$ are set to F and v, respectively. *)
    # (* If an explicit exponential occurs in u, $exponFlag$ is set to True. *)
    global SbaseS, SexponS, SexponFlagS
    if FreeQ(u, x):
        return True
    elif u == x or CalculusQ(u):
        return False
    elif PowerQ(u) and FreeQ(u.args[0], x) and LinearQ(u.args[1], x):
        SexponFlagS = True
        return FunctionOfExponentialTestAux(u.base, u.exp, x)
    elif HyperbolicQ(u) and LinearQ(u.args[0], x):
        return FunctionOfExponentialTestAux(E, u.args[0], x)
    elif PowerQ(u) and FreeQ(u.args[0], x) and SumQ(u.args[1]):
        return FunctionOfExponentialTest(u.base**First(u.exp), x) and FunctionOfExponentialTest(u.base**Rest(u.exp), x)
    return all(FunctionOfExponentialTest(i, x) for i in u.args)

def FunctionOfExponentialTestAux(base, expon, x):
    global SbaseS, SexponS, SexponFlagS
    if SbaseS == None:
        SbaseS = base
        SexponS = expon
        return True
    tmp = FullSimplify(Log(base)*Coefficient(expon, x, 1)/(Log(SbaseS)*Coefficient(SexponS, x, 1)))
    if Not(RationalQ(tmp)):
        return False
    elif ZeroQ(Coefficient(SexponS, x, 0)) or NonzeroQ(tmp - FullSimplify(Log(base)*Coefficient(expon, x, 0)/(Log(SbaseS)*Coefficient(SexponS, x, 0)))):
        if PositiveIntegerQ(base, SbaseS) and base<SbaseS:
            SbaseS = base
            SexponS = expon
            tmp = 1/tmp
        SexponS = Coefficient(SexponS, x, 1)*x/Denominator(tmp)
        if tmp < 0 and NegQ(Coefficient(SexponS, x, 1)):
            SexponS = -SexponS
            return True
        else:
            return True
        if PositiveIntegerQ(base, SbaseS) and base < SbaseS:
            SbaseS = base
            SexponS = expon
            tmp = 1/tmp
    SexponS = SexponS/Denominator(tmp)
    if tmp < 0 and NegQ(Coefficient(SexponS, x, 1)):
        SexponS = -SexponS
        return True
    return True
