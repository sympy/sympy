'''
Utility functions for Rubi integration
'''

from sympy.functions.elementary.integers import floor, frac
from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf)
from sympy.functions.elementary.hyperbolic import acosh, asinh, atanh, acoth, acsch, acsch, cosh, sinh, tanh, coth, sech, csch
from sympy.functions.elementary.trigonometric import atan, acsc, asin, acot, acos, asec
from sympy.polys.polytools import degree, Poly, quo, rem
from sympy.simplify.simplify import fraction, simplify, count_ops, factor
from sympy.core.expr import UnevaluatedExpr
from sympy.utilities.iterables import postorder_traversal
from sympy.core.expr import UnevaluatedExpr
from sympy.functions.elementary.complexes import im, re, Abs
from sympy.simplify.simplify import nthroot
from sympy.core.exprtools import factor_terms
from sympy import (exp, polylog, N, Wild, factor, gcd, Sum, S, I, Mul, hyper,
    Symbol, symbols, sqf_list, sqf, Max, gcd, hyperexpand)
from mpmath import ellippi, ellipe, ellipf, appellf1
from sympy.utilities.iterables import flatten

def Int(expr, var):
    from .rubi import rubi_integrate
    if expr == None:
        return None
    return rubi_integrate(expr, var)

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

def PositiveIntegerQ(var):
    return var.is_Integer and PositiveQ(var)

def NegativeIntegerQ(var):
    return var.is_Integer and NegativeQ(var)

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

def PosQ(var):
    return PositiveQ(var)

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
    try:
        if isinstance(expr, list):
            return expr[0]
        else:
            return expr.args[0]
    except:
        return d

def Rest(expr):
    if isinstance(expr, list):
        return expr[1:]
    else:
        if SumQ(expr) or ProductQ(expr):
            return expr.func(*expr.args[1:])
        else:
            return expr.args[1]

def SqrtNumberQ(expr):
    # SqrtNumberQ[u] returns True if u^2 is a rational number; else it returns False.
    if expr.is_Pow:
        m = expr.base
        n = expr.exp
        return IntegerQ(n) & SqrtNumberQ(m) | IntegerQ(n-1/2) & RationalQ(m)
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
    return hyperexpand(hyper([a, b], [c], z))

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
    if isinstance(expr, list):
        for e in expr:
            if not e.is_Atom:
                return False
        return True
    else:
        return expr.is_Atom

def ExpQ(u):
    return Head(u) == exp

def LogQ(u):
    return u.func == log

def Head(u):
    return u.func

def MemberQ(l, u):
    return u in l

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
    return nthroot(val, n)

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
    [u, m] = expr.args
    if FreeQ(m, x) and PolynomialQ(u, x):
        if IntegerQ(m):
            FactorSquareFree(u).match(w**n)
            return FreeQ(n, x) and LinearQ(w, x)
        else:
            return LinearQ(u, x)

def Exponent(expr, x, *k):
    if expr.is_number or (not expr.has(x)):
        return 0
    if not k:
        if expr.is_polynomial(x):
            return degree(expr, gen = x)
        else:
            return 0
    else:
        lst=[]
        if expr.is_Add:
            for t in expr.args:
                if t.is_Pow:
                    lst = lst + [t.args[1]]
                if t.is_Mul:
                    if t.args[1].is_Pow:
                        lst = lst + [(t.args[1]).args[1]]
                    if AtomQ(t.args[1]):
                        lst = lst + [1]
                if t.is_Number:
                    lst = lst + [0]
                if AtomQ(t) and not t.is_Number:
                    lst = lst + [1]
            return lst
        else:
            return [degree(expr, gen = x)]


def QuadraticQ(u, x):
    if ListQ(u):
        for expr in u:
            if Not(PolyQ(expr, x, 2) and Not(Coefficient(expr, x, 0) == 0 and Coefficient(expr, x, 1) == 0)):
                return False
        return True
    else:
        return PolyQ(u, x, 2) and Not(Coefficient(u, x, 0) == 0 and Coefficient(u, x, 1) == 0)

def LinearPairQ(u, v, x):
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
    if PolynomialQ(u, x):
        lst = CoefficientList(u, x)
        if len(lst)<3 or EvenQ(len(lst)) or ZeroQ((len(lst)+1)/2):
            return False
        #Catch(
         #   Scan(Function(if ZeroQ(lst), Null, Throw(False), Drop(Drop(Drop(lst, [(len(lst)+1)/2]), 1), -1];
          #  [First(lst), lst[(len(lst)+1)/2], Last(lst), (len(lst)-1)/2]):
    if PowerQ(u):
        if EqQ(u.args[1], 2):
            lst = BinomialParts(u.args[0], x)
            if AtomQ(lst) or ZeroQ(lst[0]):
                return False
            else:
                return [lst[0]**2, 2*lst[0]*lst[1], lst[1]**2, lst[2]]
        else:
            return False
    if ProductQ(u):
        if FreeQ(First(u), x):
            lst2 = TrinomialParts(Rest(u), x)
            if AtomQ(lst2):
                return False
            else:
                return [First(u)*lst2[0], First(u)*lst2[1], First(u)*lst2[2], lst2[3]]
        if FreeQ(Rest(u), x):
            lst1 = TrinomialParts(First(u), x)
            if AtomQ(lst1):
                return False
            else:
                return [Rest(u)*lst1[0], Rest(u)*lst1[1], Rest(u)*lst1[2], lst1[3]]
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
        if EqQ(m, n) and NonzeroQ(a*d+b*c):
            return [a*c, a*d + b*c, b*d, m]
        else:
            return False
    if SumQ(u):
        if FreeQ(First(u), x):
            lst2 = TrinomialParts(Rest(u), x)
            if AtomQ(lst2):
                return False
            else:
                return [First(u)+lst2[0], lst2[1], lst2[2], lst2[3]]
        if FreeQ(Rest(u), x):
            lst1 = TrinomialParts(First(u), x)
            if AtomQ(lst1):
                return False
            else:
                return [Rest(u)+lst1[0], lst1[1], lst1[2], lst1[3]]
        lst1 = TrinomialParts(First(u), x)
        if AtomQ(lst1):
            lst3 = BinomialParts(First(u), x)
            if AtomQ(lst3):
                return False
            lst2 = TrinomialParts(Rest(u), x)
            if AtomQ(lst2):
                lst4 = BinomialParts(Rest(u), x)
                if AtomQ(lst4):
                    return False
                if EqQ(lst3[3], 2*lst4[3]):
                    return [lst3[0]+lst4[0], lst4[1], lst3[1], lst4[2]]
                if EqQ(lst4[3], 2*lst3[3]):
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
            if AtomQ(lst4):
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

def NormalizeIntegrand(u, x):
    # returns u in a standard form recognizable by integration rules.
    return u

def SimplifyTerm(u, x):
    v = Simplify(u)
    w = Together(v)
    if LeafCount(v) < LeafCount(w):
        return NormalizeIntegrand(v, x)
    else:
        return NormalizeIntegrand(w, x)

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

def ExpandExpression(expr, x):
    return expr.expand()

def MatchQ(expr, pattern, *var):
    # returns the matched arguments after matching pattern with expression
    '''
    Example:
    >>> a_ = Wild('a', exclude=[x])
    >>> b_ = Wild('b', exclude=[x])
    >>> c_ = Wild('c', exclude=[x])
    >>> MatchQ(a + b, a_ + b_, a_, b_)
    (a, b) # or {a_: a, b_: b}
    '''
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
        result = 0
        for i in u.args:
            if not FreeQ(i, x):
                result += i
        return result
    elif not FreeQ(u, x):
        return u
    else:
        return 0


def ExpandAlgebraicFunction(expr, x):
    if ProductQ(expr):
        u_ = Wild('u', exclude=[x])
        n_ = Wild('n', exclude=[x])
        v_ = Wild('v')
        pattern = u_*v_
        match = expr.match(pattern)
        if match:
            keys = [u_, v_]
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

def UnifySum(term, x):
    # returns u with terms having indentical nonfree factors of x collected into a single term.
    return term

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
        return v.args[0]
    return v

def LeadDegree(u):
    # returns the degree of the leading factor of u.
    v = LeadFactor(u)
    if PowerQ(v):
        return v.args[1]
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
    p_ = Wild('p', exclude=[x])
    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x])
    c_ = Wild('c', exclude=[x])
    d_ = Wild('d', exclude=[x])
    n_ = Wild('n', exclude=[x])
    m_ = Wild('m', exclude=[x])

    # Basis: If  m\[Element]\[DoubleStruckCapitalZ] \[And] b c-a d==0, then (a+b z)^m==b^m/d^m (c+d z)^m
    pattern = u_*(a_ + b_*x)**m_*(c_ + d_*x)**n_
    match = expr.match(pattern)
    if match:
        keys = [u_, a_, b_, m_, c_, d_, n_]
        if len(keys) == len(match):
            u, a, b, m, c, d, n = tuple([match[i] for i in keys])
            if IntegerQ(m) and ZeroQ(b*c - a*d):
                return u*b**m/d**m*(c + d*x)**(m + n)

    # Basis: If  m/n\[Element]\[DoubleStruckCapitalZ], then z^m (c z^n)^p==(c z^n)^(m/n+p)/c^(m/n)
    pattern = u_*(a_ + b_*x)**m_*(c_*(a_ + b_*x)**n_)**p_
    match = expr.match(pattern)
    if match:
        keys = [u_, a_, b_, m_, c_, n_, p_]
        if len(keys) == len(match):
            u, a, b, m, c, n, p = tuple([match[i] for i in keys])
            if IntegerQ(m/n):
                return u*(c*(a + b*x)**n)**(m/n + p)/c**(m/n)
    return expr

def PolynomialDivide(p, q, x):
    return quo(p, q) + rem(p, q)/q

def BinomialQ(u, x, n=None):
    if ListQ(u):
        for i in u.args:
            if Not(BinomialQ(i, x, n)):
                return False
        return True
    return ListQ(BinomialParts(u, x))

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
    u_ = Wild('u')
    v_ = Wild('v')
    n_ = Wild('n', exclude=[x])
    pattern = u_*v_**n_
    match = expr.match(pattern)
    if match:
        print(match)
        keys = [u_, v_, n_]
        if len(keys) == len(match):
            u, v, n = tuple([match[i] for i in keys])
            if FractionQ(n) and v != x:
                w = RationalFunctionExpand(u, x)
                if SumQ(w):
                    result = 0
                    for i in w.args:
                        result += i*v**n
                    return result
                else:
                    return w*v**n

    v = ExpandIntegrand(expr, x)
    t = False

    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x])
    c_ = Wild('c', exclude=[x])
    d_ = Wild('d', exclude=[x])
    m_ = Wild('m', exclude=[x])
    p_ = Wild('p', exclude=[x])

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
    else:
        v = ExpandIntegrand(RationalFunctionFactors(u, x), x)
        w = NonrationalFunctionFactors(u, x)
        if SumQ(v):
            result = 0
            for i in v:
                result += i*w
            return result

def ExpandIntegrand(expr, x, extra=None):
    w_ = Wild('w')
    p_ = Wild('p', exclude=[x])
    q_ = Wild('q', exclude=[x])
    u_ = Wild('u')
    v_ = Wild('v')
    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x])
    c_ = Wild('c', exclude=[x])
    d_ = Wild('d', exclude=[x])
    e_ = Wild('e', exclude=[x])
    f_ = Wild('f', exclude=[x])
    g_ = Wild('g', exclude=[x])
    h_ = Wild('h', exclude=[x])
    j_ = Wild('j', exclude=[x])
    n_ = Wild('n', exclude=[x])
    m_ = Wild('m', exclude=[x])
    F_ = Wild('F', exclude=[x])
    k, q, i = symbols('k q i')

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

    pattern = u_*(a_ + b_*F_**v_)**m_*(c_ + d_*F_**v_)**n_
    match = expr.match(pattern)
    if match:
        keys = [u_, a_, b_, F_, v_, m_, c_, d_, n_]
        if len(keys) == len(match):
            u, a, b, F, v, m, c, d, n = tuple([match[i] for i in keys])
            if IntegersQ(m, n) and NegativeQ(n):
                w = ReplaceAll(ExpandIntegrand((a+b*x)**m*(c + d*x)**n, x), {x: F**v})
                result = []
                for i in w.args:
                    result.append(i*u)
                return w.func(*result)

    pattern = u_*(a_ + b_*x)**m_*Log(c_*(d_ + e_*x**n_)**p_)
    match = expr.match(pattern)
    if match:
        keys = [u_, a_, b_, m_, c_, d_, e_, n_, p_]
        if len(keys) == len(match):
            u, a, b, m, c, d, e, n, p = tuple([match[i] for i in keys])
            if PolynomialQ(u, x):
                return ExpandIntegrand(Log(c*(d + e*x**n)**p), x, u*(a + b*x)**m)

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

    pattern = u_*(a_ + b_*Log(c_*(d_*(e_ + f_*x)**p_)**q_))**n_
    match = expr.match(pattern)
    if match:
        keys = [u_, a_, b_, c_, d_, e_, f_, p_, q_, n_]
        if len(keys) == len(match):
            u, a, b, c, d, e, f, p, q, n = tuple([match[i] for i in keys])
            if PolynomialQ(u, x):
                return ExpandLinearProduct((a + b*Log(c*(d*(e + f*x)**p)**q))**n, u, e, f, x)

    pattern = (a_ + b_*u_**n_ + c_*u_**j_)**p_
    match = expr.match(pattern)
    if match:
        keys = [a_, b_, u_, n_, c_, j_, p_]
        if len(keys) == len(match):
            a, b, u, n, c, j, p = tuple([match[i] for i in keys])
            if IntegerQ(n) and ZeroQ(j - 2*n) and NegativeIntegerQ(p) and NonzeroQ(b**2 - 4*a*c):
                ReplaceAll(ExpandIntegrand(1/(4**p*c**p), x, (b - q + 2*c*x)**p*(b + q + 2*c*x)**p), {q: Rt(b**2-4*a*c,2),x: u**n})

    pattern = u_**m_*(a_ + b_*u_**n_ + c_*u_**j_)**p_
    match = expr.match(pattern)
    if match:
        keys = [u_, m_, a_, b_, n_, c_, j_, p_]
        if len(keys) == len(match):
            u, m, a, b, n, c, j, p = tuple([match[i] for i in keys])
            if IntegersQ(m, n, j) and ZeroQ(j - 2*n) and NegativeIntegerQ(p) and 0<m<2*n and Not(m == n and p == -1) and NonzeroQ(b**2 - 4*a*c):
                return ReplaceAll(ExpandIntegrand(1/(4**p*c**p), x, x**m*(b - q + 2*c*x**n)**p*(b + q+ 2*c*x**n)**p), {q: Rt(b**2 - 4*a*c, 2),x: u})

    # Basis: If  q=Sqrt[-a c], then a+c z^2==((-q+c z)(q+c z))/c
    pattern = (a_ + c_*u_**n_)**p_
    match = expr.match(pattern)
    if match:
        keys = [a_, c_, u_, n_, p_]
        if len(keys) == len(match):
            a, c, u, n, p = tuple([match[i] for i in keys])
            if IntegerQ(n/2) and NegativeIntegerQ(p):
                return ReplaceAll(ExpandIntegrand(1/c**p, x, (-q + c*x)**p*(q + c*x)**p), {q: Rt(-a*c,2),x: u**(n/2)})

    pattern = u_**m_*(a_ + c_*u_**n_)**p_
    match = expr.match(pattern)
    if match:
        keys = [u_, m_, a_, c_, n_, p_]
        if len(keys) == len(match):
            u, m, a, c, n, p = tuple([match[i] for i in keys])
            if IntegersQ(m, n/2) and NegativeIntegerQ(p) and 0<m<n and (m != n/2):
                return ReplaceAll(ExpandIntegrand(1/c**p, x, x**m*(-q + c*x**(n/2))**p*(q + c*x**(n/2))**p),{q: Rt(-a*c, 2), x: u})

    # Basis: 1/(a x^n+b Sqrt[c+d x^(2 n)])==(a x^n-b Sqrt[c+d x^(2 n)])/(-b^2 c+(a^2-b^2 d) x^(2 n))
    pattern = u_/(a_*x**n_ + b_*Sqrt(c_ + d_*x**j_))
    match = expr.match(pattern)
    if match:
        keys = [u_, a_, n_, b_, c_, d_, j_]
        if len(keys) == len(match):
            u, a, n, b, c, d, j = tuple([match[i] for i in keys])
            if ZeroQ(j - 2*n):
                return ExpandIntegrand(u*(a*x**n - b*Sqrt(c + d*x**(2*n)))/(-b**2*c + (a**2 - b**2*d)*x**(2*n)), x)

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
                return r/(3*a*(r - s*u)) + r*(2*r + s*u)/(3*a*(r**2 + r*s*u + s**2*u**2))
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

    # If u is a polynomial in x, ExpandIntegrand[u*(a+b*x)^m,x] expand u*(a+b*x)^m into a sum of terms of the form A*(a+b*x)^n.
    pattern = u_*(a_ + b_*x)**m_
    match = expr.match(pattern)
    if match:
        keys = [u_, a_, b_, m_]
        if len(keys) == len(match):
            u, a, b, m = tuple([match[i] for i in keys])
            try:
                w, c, d, p = MatchQ(u, w_*(c_+d_*x)**p_, w_, c_, d_, p_)
                res = IntegerQ(p) & p > m
            except: # if not matched
                res = False
            if PolynomialQ(u, x) and Not(PositiveIntegerQ(m) and res):
                tmp1 = ExpandLinearProduct((a+b*x)**m, u, a, b, x)
                if not IntegerQ(m):
                    return tmp1
                else:
                    tmp2 = ExpandExpression(u*(a+b*x)**m, x)
                    if SumQ(tmp2) and (LeafCount(tmp2) <= LeafCount(tmp1)+2):
                        return tmp2
                    else:
                        return tmp1

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

    # Basis: If  (m|(n-1)/2)\[Element]\[DoubleStruckCapitalZ] \[And] 0<=m<n, let r/s=(a/b)^(1/n), then z^m/(a + b*z^n) == (r*(-(r/s))^m*Sum[1/((-1)^(2*k*(m/n))*(r + (-1)^(2*(k/n))*s*z)), {k, 1, n}])/(a*n) == (r*(-(r/s))^m*Sum[(-1)^(2*k*((m + 1)/n))/((-1)^(2*(k/n))*r + s*z), {k, 1, n}])/(a*n)
    pattern = u_**m_/(a_+b_*u_**n_)
    match = expr.match(pattern)
    if match:
        keys = [u_, m_, a_, b_, n_]
        if len(keys) == len(match):
            u, m, a, b, n = tuple([match[i] for i in keys])
            if IntegersQ(m, n) & 0 < m < n & OddQ(n/GCD(m, n)) & PosQ(a/b):
                g = GCD(m, n)
                r = Numerator(Rt(a/b, n/GCD(m, n)))
                s = Denominator(Rt(a/b, n/GCD(m, n)))
                if CoprimeQ(m + g, n):
                    return Sum(r*(-r/s)**(m/g)*(-1)**(-2*k*m/n)/(a*n*(r + (-1)**(2*k*g/n)*s*u**g)),(k, 1, n/g)).doit()
                else:
                    return Sum(r*(-r/s)**(m/g)*(-1)**(2*k*(m+g)/n)/(a*n*((-1)**(2*k*g/n)*r + s*u**g)),(k, 1, n/g)).doit()
            elif IntegersQ(m, n) & 0<m<n:
                g = GCD(m, n)
                r = Numerator(Rt(-a/b, n/GCD(m, n)))
                s = Denominator(Rt(-a/b, n/GCD(m, n)))
                if n/g == 2:
                    return s/(2*b*(r+s*u^g)) - s/(2*b*(r-s*u^g))
                else:
                    if CoprimeQ[m+g,n]:
                        return Sum(r*(r/s)**(m/g)*(-1)**(-2*k*m/n)/(a*n*(r - (-1)**(2*k*g/n)*s*u**g)),(k,1,n/g)).doit()
                    else:
                        return Sum(r*(r/s)**(m/g)*(-1)**(2*k*(m+g)/n)/(a*n*((-1)**(2*k*g/n)*r - s*u**g)),(k,1,n/g)).doit()

    pattern = u_/v_
    match = expr.match(pattern)
    if match:
        keys = [u_, v_]
        if len(keys) == len(match):
            u, v = tuple([match[i] for i in keys])
            if PolynomialQ(u, x) and PolynomialQ(v, x) and BinomialQ(v,x) and (Exponent(u, x) == Exponent(v, x)-1 >= 2):
                lst = CoefficientList(u, x)
                return lst[-1]*x**Exponent[u,x]/v + Sum(lst[i]*x**(i-1),(i, 1, Exponent[u,x])).doit()/v
            elif PolynomialQ(u, x) and PolynomialQ(v, x) and Exponent(u, x)>=Exponent(v, x):
                return PolynomialDivide(u, v, x)

    #return None
    return expr.expand()

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
        return QuotientOfLinearsP(Denominator(u))
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
        result = 1
        for i in u:
            if AbsurdNumberQ(i):
                result *= i
        return result
    return NumericFactor(u)

def NonabsurdNumberFactors(u):
    # (* NonabsurdNumberFactors[u] returns the product of the factors of u that are not absurd numbers. *)
    if AbsurdNumberQ(u):
        return 1
    elif ProductQ(u):
        result = 1
        for i in u.args:
            result *= NonabsurdNumberFactors(i)
        return result
    return NonnumericFactors(u)
