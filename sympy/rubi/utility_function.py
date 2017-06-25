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
from sympy import exp, polylog, N, Wild, factor, gcd, Sum
from mpmath import hyp2f1, ellippi, ellipe, ellipf, appellf1, nthroot

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
    return u < 0

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

def PositiveIntegerQ(var):
    return var.is_Integer and var > 0

def NegativeIntegerQ(var):
    return var.is_Integer and var < 0

def PositiveQ(var):
    return var > 0

def IntegerQ(var):
    if isinstance(var, int):
        return True
    else:
        return var.is_Integer

def IntegersQ(*var):
    return all(IntegerQ(i) for i in var)

def ComplexNumberQ(*var):
    return all(i.is_complex for i in var)

def RealNumericQ(u):
    return u.is_real

def PositiveOrZeroQ(u):
    return u.is_real and u >= 0

def NegativeOrZeroQ(u):
    return u.is_real and u <= 0

def FractionOrNegativeQ(u):
    return FractionQ(u) or u < 0

def PosQ(var):
    return var > 0

def NegQ(var):
    return var < 0

def Equal(a, b):
    return a == b

def Unequal(a, b):
    return a != b

def IntPart(u):
    # IntPart[u] returns the sum of the integer terms of u.
    if ProductQ(u):
        if IntegerQ(First(u)):
            return First(u)*IntPart(Rest(u))

    if IntegerQ(u):
        return u
    elif FractionQ(u):
        return IntegerPart(u)
    elif SumQ(u):
        res = 0
        for i in u.args:
            res += IntPart(u)
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
            return expr.func(expr.args[0])
    except:
        return d

def Rest(expr):
    if isinstance(expr, list):
        return expr[1:]
    else:
        return expr.func(*expr.args[1:])

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
    if degree(expr, gen=x) == 1:
        return True
    else:
        return False

def Sqrt(a):
    return sqrt(a)

def ArcCosh(a):
    return acosh(a)

def Coefficient(expr, var, n):
    a = Poly(expr, var)
    if (degree(a) - n) < 0:
        return 0
    else:
        return a.all_coeffs()[degree(a) - n]

def Denominator(var):
    return fraction(var)[1]

def Hypergeometric2F1(a, b, c, z):
    return hyp2f1(a, b, c, z)

def ArcTan(a):
    return atan(a)

def Not(var):
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
        if args[i] > args[i + 1]:
            return False
    return True

def Less(*args):
    for i in range(0, len(args) - 1):
        if args[i] >= args[i + 1]:
            return False
    return True

def Greater(*args):
    for i in range(0, len(args) - 1):
        if args[i] <= args[i + 1]:
            return False
    return True

def GreaterEqual(*args):
    for i in range(0, len(args) - 1):
        if args[i] < args[i + 1]:
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
    return PowerQ(u) and IntegerQ(u.args[1]) and u.args[1]>0

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
    return LogQ(u) or InverseTrigQ(u) and Length(u)<=1 or InverseHyperbolicQ(u) or u.func == polylog

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
        if InverseFunctionQ(u) | CalculusQ(u) | u.func == hyp2f1 | u.func == appellf1:
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
    e = factor(u).args
    mul_fac = 1
    for expr in e:
        if expr.is_Pow:
            mul_fac = mul_fac*expr
    return mul_fac*simplify(u/mul_fac)

def PowerOfLinearQ(expr, x):
    [u, m] = expr.args
    if FreeQ(m, x) and PolynomialQ(u, x):
        if IntegerQ(m):
            FactorSquareFree(u).match(w**n)
            return FreeQ(n, x) and LinearQ(w, x)
        else:
            return LinearQ(u, x)

def Exponent(expr, x, *k):
    if not k:
        return degree(expr, gen = x)
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
        if Exponent(u, x)>0:
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
    return u
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

def ReplaceAll(expr, x, substitution):
    return expr.subs(x, substitution)

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
        lst = CoefficientList(ReplaceAll(u, x, (x - a)/b), x)
        lst = [SimplifyTerm(i, x) for i in lst]
        res = 0
        for k in range(1, len(lst)+1):
            res = res + v*lst[k-1]*(a + b*x)**(k - 1)
        return res

def GCD(*args):
    return gcd(*args)

def ContentFactor(expn):
    return ContentFactorAux(expn)

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
        if RationalQ(u.args[0]) and RationalQ(u.args[1]):
            if u.args[1] > 0:
                return 1/Denominator(u.args[1])
            else:
                return 1/(1/Denominator(u.args[1]))
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
    else:
        return S(1)

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

def RemoveContent(u, x):
    # returns u with the content free of x removed.
    return None

def ExpandIntegrand(expr, x):

    w_ = Wild('w')
    p_ = Wild('p')
    u_ = Wild('u')
    v_ = Wild('v')
    a_ = Wild('a', exclude=[x])
    b_ = Wild('b', exclude=[x])
    c_ = Wild('c', exclude=[x])
    d_ = Wild('d', exclude=[x])
    n_ = Wild('n', exclude=[x])
    m_ = Wild('m', exclude=[x])

    # Basis: (a+b x)^m/(c+d x)==(b (a+b x)^(m-1))/d+((a d-b c) (a+b x)^(m-1))/(d (c+d x))
    pattern = (a_ + b_*x)**m_/(c_ + d_*x)
    match = expr.match(pattern)
    if match:
        keys = [a_, b_, c_, d_, m_]
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

    # If u is a polynomial in x, ExpandIntegrand[u*(a+b*x)^m,x] expand u*(a+b*x)^m into a sum of terms of the form A*(a+b*x)^n.
    pattern = u_*(a_ + b_*x)**m_
    match = expr.match(pattern)
    if match:
        keys = [u_, a_, b_, m_]
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

    # *Basis: If (m|n)\[Element]\[DoubleStruckCapitalZ] \[And] 0<=m<n, let r/s=(-(a/b))^(1/n), then  (c + d*z^m)/(a + b*z^n) == (r*Sum[(c + (d*(r/s)^m)/(-1)^(2*k*(m/n)))/(r - (-1)^(2*(k/n))*s*z), {k, 1, n}])/(a*n)
    pattern = (c_ + d_*u_**m_)/(a_ + b_*u_**n_)
    match = expr.match(pattern)
    if match:
        keys = [c_, d_, u_, m_, a_, b_, u_, n_]
        c, d, u, m, a, b, u, n = tuple([match[i] for i in keys])
        if IntegersQ[m,n] & 0<m<n:
            r = Numerator(Rt(-a/b, n))
            s = Denominator(Rt(-a/b, n))
            k = 1
            return Sum((r*c + r*d*(r/s)**m*(-1)**(-2*k*m/n))/(a*n*(r-(-1)**(2*k/n)*s*u)),(k,1,n))


    pattern = u_*v_**n_*(a_ + b_*x)**m_
    match = expr.match(pattern)
    if match:
        keys = [u_, v_, n_, a_, b_, m_]
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
        if IntegerQ(m, n) & 0<m<n & OddQ(n/GCD(m, n)) & PosQ(a/b):
            g = GCD(m, n)
            r = Numerator(Rt(a/b, n/GCD(m, n)))
            s = Denominator(Rt(a/b, n/GCD(m, n)))
            k = 1
            if CoprimeQ(m + g, n):
                return Sum(r*(-r/s)**(m/g)*(-1)**(-2*k*m/n)/(a*n*(r + (-1)**(2*k*g/n)*s*u**g)),(k, 1, n/g))
            else:
                return Sum(r*(-r/s)**(m/g)*(-1)**(2*k*(m+g)/n)/(a*n*((-1)**(2*k*g/n)*r + s*u**g)),(k, 1, n/g))
        elif IntegersQ(m, n) & 0<m<n:
            g = GCD(m, n)
            r = Numerator(Rt(-a/b, n/GCD(m, n)))
            s = Denominator(Rt(-a/b, n/GCD(m, n)))
            if n/g == 2:
                return s/(2*b*(r+s*u^g)) - s/(2*b*(r-s*u^g))
            else:
                if CoprimeQ[m+g,n]:
                    return Sum(r*(r/s)**(m/g)*(-1)**(-2*k*m/n)/(a*n*(r - (-1)**(2*k*g/n)*s*u**g)),(k,1,n/g))
                else:
                    return Sum(r*(r/s)**(m/g)*(-1)**(2*k*(m+g)/n)/(a*n*((-1)**(2*k*g/n)*r - s*u**g)),(k,1,n/g))

    #return None
    return expr.expand()
