from sympy.utilities import all, any
from sympy.functions import sqrt, sign
from sympy.core import Pow, S, Number, Wild, Rational, sympify

def sqrtdenest (expr):
    """
    Denests an expression that contains nested square roots.
    This algorithm is based on <http://www.almaden.ibm.com/cs/people/fagin/symb85.pdf>.
    """
    expr = sympify(expr)
    if expr.is_Pow and expr.exp is S.Half: #If expr is a square root
        return denester([expr])[0]
    return expr

def denester (nested):
    """
    Denests a list of expressions that contain nested square roots.
    This method should not be called directly - use 'denest' instead.
    This algorithm is based on <http://www.almaden.ibm.com/cs/people/fagin/symb85.pdf>.

    It is assumed that all of the elements of 'nested' share the same
    bottom-level radicand. (This is stated in the paper, on page 177, in
    the paragraph immediately preceding the algorithm.)

    When evaluating all of the arguments in parallel, the bottom-level
    radicand only needs to be denested once. This means that calling
    denester with x arguments results in a recursive invocation with x+1
    arguments; hence denester has polynomial complexity.

    However, if the arguments were evaluated separately, each call would
    result in two recursive invocations, and the algorithm would have
    exponential complexity.

    This is discussed in the paper in the middle paragraph of page 179.
    """
    if all((n**2).is_Number for n in nested): #If none of the arguments are nested
        for f in subsets(len(nested)): #Test subset 'f' of nested
            p = prod(nested[i]**2 for i in range(len(f)) if f[i]).expand()
            if 1 in f and f.count(1) > 1 and f[-1]: p = -p
            if sqrt(p).is_Number: return sqrt(p), f #If we got a perfect square, return its square root.
        return nested[-1], [0]*len(nested) #Otherwise, return the radicand from the previous invocation.
    else:
        a, b, r, R = Wild('a'), Wild('b'), Wild('r'), None
        values = [expr.match(sqrt(a + b * sqrt(r))) for expr in nested]
        for v in values:
            if r in v: #Since if b=0, r is not defined
                if R is not None: assert R == v[r] #All the 'r's should be the same.
                else: R = v[r]
        d, f = denester([sqrt((v[a]**2).expand()-(R*v[b]**2).expand()) for v in values] + [sqrt(R)])
        if not any([f[i] for i in range(len(nested))]): #If f[i]=0 for all i < len(nested)
            v = values[-1]
            return sqrt(v[a] + v[b]*d), f
        else:
            v = prod(nested[i]**2 for i in range(len(nested)) if f[i]).expand().match(a+b*sqrt(r))
            if 1 in f and f.index(1) < len(nested) - 1 and f[len(nested)-1]:
                v[a] = -1 * v[a]
                v[b] = -1 * v[b]
            if not f[len(nested)]: #Solution denests with square roots
                return (sqrt((v[a]+d).expand()/2)+sign(v[b])*sqrt((v[b]**2*R/(2*(v[a]+d))).expand())).expand(), f
            else: #Solution requires a fourth root
                FR, s = (R.expand()**Rational(1,4)), sqrt((v[b]*R).expand()+d)
                return (s/(sqrt(2)*FR) + v[a]*FR/(sqrt(2)*s)).expand(), f

def subsets(n):
    """
    Returns all possible subsets of the set (0, 1, ..., n-1) except the empty set.
    """
    binary = lambda x: x>0 and binary(x>>1) + [x&1] or []
    pad = lambda l: [0]*(n-len(l)) + l #Always returns a list of length 'n'
    return [pad(binary(i)) for i in range(1, 2**n)]

def prod(n):
    """
    Returns the product of all elements of n, as a Rational.
    """
    product = S.One
    for i in n:
        product = product * i
    return product
