'''
Utility functions for Constraints in Rubi
'''

from sympy.functions.elementary.integers import floor, frac
from sympy.functions import (log, sin, cos, tan, cot, csc, sec, sqrt, erf)
from sympy.functions.elementary.integers import floor, frac
from sympy.functions.elementary.hyperbolic import acosh, asinh, atanh, acoth, acsch, acsch, cosh, sinh, tanh, coth, sech, csch
from sympy.functions.elementary.trigonometric import atan, acsc, asin, acot, acos, asec
from sympy.polys.polytools import degree, Poly
from sympy.simplify.simplify import fraction, simplify, count_ops, factor
from sympy.core.expr import UnevaluatedExpr
from sympy.utilities.iterables import postorder_traversal
from sympy.core.expr import UnevaluatedExpr
from sympy.functions.elementary.complexes import im, re, Abs
from sympy import exp, polylog, N

from mpmath import hyp2f1, ellippi, ellipe, ellipf, appellf1, nthroot

#from .rubi import rubi_integrate

#def Int(expr, var):
#    return rubi_integrate(expr, var)

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

def FracPart(var):
    return frac(var)

def IntPart(var):
    return floor(var)

def RationalQ(*nodes):
    return all(var.is_Rational for var in nodes)

def ProductQ(expr):
    return expr.is_Mul

def SumQ(expr):
    return expr.is_Add

def NonsumQ(expr):
    return not SumQ(expr)

def Subst(a, x, y):
    return a.subs(x, y)

def First(expr, d=None):
    # gives the first element if it exists, or d otherwise.
    try:
        return expr[0]
    except:
        return d

def Rest(l):
    return l[1:]

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
    return FracPart(a)

def IntegerPart(a):
    return IntPart(a)

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
    return u.is_Number

def Length(expr):
    # returns number of elements in the experssion
    return len(expr.args)

def ListQ(u):
    return isinstance(u, list)

def Im(u):
    return im(u)

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

def First(nodes):
    if isinstance(nodes, dict):
        return nodes[list(nodes.keys())[0]]
    elif isinstance(nodes, list):
        return nodes[0]
    elif ExpQ(nodes):
        return exp
    else:
        return nodes.args[0]

def Rest(nodes):
    if isinstance(nodes, dict):
        del nodes[list(nodes.keys())[0]]
        return nodes
    elif isinstance(nodes, list):
        del nodes[0]
        return nodes
    elif ExpQ(nodes):
        return nodes.args[0]
    elif nodes.is_Add:
        return nodes - nodes.args[0]
    elif nodes.is_Pow:
        return nodes.args[1]
    elif nodes.is_Mul:
        return nodes/nodes.args[0]

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

def PolyQ(u, x, n):
    if ListQ(u):
        for expr in u:
            if Not(PolynomialQ(expr, x) and Exponent(expr, x) == n and Coefficient(expr, x, n) != 0):
                return False
        return True
    else:
        return PolynomialQ(u, x) and Exponent(u, x) == n and Coefficient(u, x, n) != 0


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

def MonomialQ(u, x):
    if ListQ(u):
        for expr in u:
            if not MonomialQ(expr, x):
                return False
        return True 
    [a, n] = u.args

def MonomialSumQ(u, x):
    if SumQ(u):
        for expr in u.args:
            if Not(FreeQ(expr, x) or MonomialQ(expr, x)):
                return False
        return True

#def MonomialExponent(a_.*x_^n_., x):
 # n /; 
#FreeQ[{a,n},x]

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

def CoefficientList(u, x):
    a = Poly(u, x)
    return a.all_coeffs()

def EvenQ(u):
    # gives True if expr is an even integer, and False otherwise.
    return u%2 == 0

def OddQ(u):
    # gives True if expr is an odd integer, and False otherwise.
    return u%2 == 1

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
                return [lst[0]^2, 2*lst[0]*lst[1], lst[1]^2, lst[2]]
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

'''
def CancelCommonFactors(u v)
    if ProductQ(u):
        if ProductQ(v):
            if MemberQ(v, First(u)):
                return CancelCommonFactors(Rest(u), DeleteCases(v, First(u), 1, 1))
            return Function({First(u) ((1))((2))})(CancelCommonFactors(Rest(u)v))
        if MemberQ(u, v):
            return [DeleteCases(u, v, 1, 1), 1]
        return [u, v]
    if ProductQ(v)
        if MemberQ(v, u)
            return [1, DeleteCases(v, u, 1, 1)]
        return [uv]
    return [uv]

def SimplerQ(u, v):
    if IntegerQ(u):
        if IntegerQ(v):
            if Abs(u) == Abs(v):
                return v<0
            else:
                return Abs(u)<Abs(v)
        else:
            return True
    elif not IntegerQ(u) and IntegerQ(v):
        return False
    if FractionQ(u):
        if FractionQ(v):
            if Denominator(u) == Denominator(v):
                return SimplerQ(Numerator(u), Numerator(v))
            else:
                return Denominator(u)<Denominator(v):
        else:
            return True
    elif not FractionQ(u) and FractionQ(v)
        return False
    if (re(u) == 0 or re(u) == 0.0) and (re(v) == 0 or re(v) == 0.0):
        return SimplerQ(im(u), im(v))
    if ComplexNumberQ(u):
        if ComplexNumberQ(v):
            if re(u) == re(v):
                return SimplerQ(im(u), im(v))
            else:
                return SimplerQ(re(u), re(v))
        else:
            return False
    if NumberQ(u):
        if NumberQ(v):
            return OrderedQ((u, v))
        else:
            return True
    if not NumberQ(u) and NumberQ(v):
        return False
    if AtomQ(u):
        if AtomQ(v):
            return OrderedQ((u, v))
        else:
            return True
    if not AtomQ(u) and AtomQ(v):
        return False
    if Head(u) == Head(v):
        if Length(u) == Length(v):
            for i in range(0, len(u)):
                if not u.args[i] == v.args[i]:
                    return SimplerQ(u.args[i], v.args[i])
        else:
         return Length(u)<Length(v)
    if LeafCount(u)<LeafCount(v):
       return True
    if LeafCount(v)<LeafCount(u):
       return False
    else:
        return Not(OrderedQ(v, u))

def SimplerSqrtQ(u, v):
    if NegativeQ(v) and Not(NegativeQ(u)):
        return True
    if NegativeQ(u) and Not(NegativeQ(v)):
        return False
    sqrtu = Rt(u, 2)
    sqrtv = Rt(v,2)
    if IntegerQ(sqrtu):
        if IntegerQ(sqrtv):
            return sqrtu<sqrtv
        else:
            return True
    elif IntegerQ(sqrtv):
        return False
    if RationalQ(sqrtu):
        if RationalQ(sqrtv):
            return sqrtu<sqrtv
        else:
            return True
    elif RationalQ(sqrtv):
        return False
    if PosQ(u):
        if PosQ(v):
            return LeafCount(sqrtu)<LeafCount(sqrtv)
        else:
            return True
    elif PosQ(v):
        return False
    if LeafCount(sqrtu)<LeafCount(sqrtv):
        return True
    if LeafCount(sqrtv)<LeafCount(sqrtu):
        return False
    return Not(OrderedQ(v,u))

def SumSimplerQ(u, v):
    if RationalQ(u, v):
        if v == 0:
            return False
        elif v>0:
            return u<-1
        else:
            return u>=-v
    else:
        return SumSimplerAuxQ(Expand(u), Expand(v))


def SumSimplerAuxQ(u, v):
    if (RationalQ(First(v)) or SumSimplerAuxQ(u, First(v))) and 
    (RationalQ(Rest(v)) or SumSimplerAuxQ(u, Rest(v)))
        return SumQ(v)
    if SumSimplerAuxQ(First(u),v) or SumSimplerAuxQ(Rest(u),v)
        return SumQ(u)
    if not v == 0 and NonnumericFactors(u) == NonnumericFactors(v) and 
    (NumericFactor(u)/NumericFactor(v)<-1/2 or NumericFactor(u)/NumericFactor(v)==-1/2 and NumericFactor(u)<0)

def NonnumericFactors(u):
    if NumberQ(u):
        if ZeroQ(Im(u)):
            return 1        
        elif ZeroQ(Re(u)):
            return I
        else:
            return u
    if PowerQ(u):
        if RationalQ(u.args[0]) && FractionQ(u.args[1]):
            return u/NumericFactor(u)
        else:
            return u
    if ProductQ(u):
        return Map(NonnumericFactors,u)
    if SumQ(u)
        if LeafCount(u)<50,                         (* Eliminate this kludge! *)
            Function(if SumQ(#), u, NonnumericFactors(#)))(ContentFactor(u))
        With({n=NumericFactor(u)}
        Map(Function(#/n),u)))
    u))))

'''