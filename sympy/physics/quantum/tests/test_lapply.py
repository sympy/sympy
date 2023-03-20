from sympy.core.mul import Mul
from sympy.core.numbers import Rational
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.core import Pow

from sympy.physics.quantum.lapply import lapply, c_nc_ncef #------------ adapt path to final location of lapply

def test_expansion_NCSym():
    # check expansion behaviour on options mul and power_base,
    # power_exp on commutative factors. lapply() applies hint
    # functions for options  mul, power_base and power_exp of expand()

    A1, A2, A3 = symbols("A1 A2 A3", commutative=False)

    p = (2*(A1 + A2) + 2*(A1 - A2))*A3  # distribute factor A3 over sums
    assert lapply(p, mul=False) == 2*(A1*A3 - A2*A3) + 2*(A1*A3 + A2*A3)
    # distribute commutative factors over sums, too  (default)
    assert lapply(p, mul=True) == 4*A1*A3

    u, v = symbols("u v", commutative=True)
    w, z = symbols("w z", integer=True, positive=True)

    p = (u*(A1 + A2) + u*(A1 - A2))*A3
    assert lapply(p, mul=False) == u*(A1*A3 - A2*A3) + u*(A1*A3 + A2*A3)
    assert lapply(p, mul=True) == 2*u*A1*A3 # mul=True is default

    p = (u*(A1 + A2) + v*(A1 - A2))*A3
    assert lapply(p, mul=False) == u*(A1*A3 + A2*A3) + v*(A1*A3 - A2*A3)
    assert lapply(p, mul=True) == u*A1*A3 + u*A2*A3 + v*A1*A3 - v*A2*A3

    p = (u*(A1 + A2) + u*(A1 - A2))**2 * A3
    r = u*(u*(A1**2*A3 - A2*A1*A3) + u*(A1**2*A3 + A2*A1*A3) + \
        Mul(-1, (u*(A1*A2*A3 - A2**2*A3) + u*(A1*A2*A3 + A2**2*A3)), evaluate=False)) +\
        u*(u*(A1**2*A3 - A2*A1*A3) + u*(A1**2*A3 + A2*A1*A3) +  \
        u*(A1*A2*A3 - A2**2*A3) + u*(A1*A2*A3 + A2**2*A3) )
    assert lapply(p, mul=False) == r
    assert lapply(p) == 4 * u**2 * A1**2 * A3 # mul=True is default

    # Check  effects of expand options power_base and power_exp on
    # communicative factors. Note that sqrt(u*(u+1) + w*(w+1)) is unmodified.
    p = u*(w + sqrt(u*(u + 1) + w*(w + 1)))**(w + z)
    # commutative factors left untouched
    assert lapply(p, mul=False, power_base=False, power_exp=False) == p
    # mul/add expansion has no effect on p
    assert lapply(p, mul=True, power_base=False, power_exp=False) == p
    # but power_exp has
    assert lapply(p, mul=False, power_base=False, power_exp=True) == \
        u*(w + sqrt(u*(u + 1) + w*(w + 1)))**w*(w + sqrt(u*(u + 1) + w*(w + 1)))**z
    assert lapply(p, mul=True, power_base=False, power_exp=True) == \
        u*(w + sqrt(u*(u + 1) + w*(w + 1)))**w*(w + sqrt(u*(u + 1) + w*(w + 1)))**z

    p = (u*(w + sqrt(u*(u + 1) + w*(w + 1))))**(w + z)
    # commutative factors left untouched
    assert lapply(p, mul=False, power_base=False, power_exp=False) == p
    # power_base has effect
    assert lapply(p, mul=False, power_base=True, power_exp=False) == \
        u**(w + z) * (w + sqrt(u*(u + 1) + w*(w + 1)))**(w + z)
    # mul has effect in basis
    assert lapply(p, mul=True, power_base=False, power_exp=False) == \
        (u*w + u*sqrt(u*(u + 1) + w*(w + 1)))**(w + z)
    assert lapply(p, mul=False, power_base=True, power_exp=True) == \
        u**w * u**z * (w + sqrt(u*(u + 1) + w*(w + 1)))**w * \
        (w + sqrt(u*(u + 1) + w*(w + 1)))**z # power_xxx has effect


# Various test cases for lapply and matrices
from sympy.matrices import (Matrix, ImmutableMatrix, MatrixSymbol,
            MatPow, MatMul, DiagMatrix, Identity)
from sympy import simplify
from sympy.core.expr import unchanged

def test_powers2023_01_Matrix():
    # Use MatrixSymbols and Identity as generic nc linear operators. Note
    # that non-commutative factors are not implemented in MatMul, so we
    # use MatrixSymbol throughout.
    n = symbols("n", integer=True, positive=True)
    A = MatrixSymbol('A', n, n)
    B = MatrixSymbol('B', n, n) # generic n-dimensional linear operator
    In = Identity(n)
    c, d = symbols("c d", commutative=True)
    f, g = symbols("f g", commutative=True, nonnegative=True)
    o = symbols("o", commutative=True, integer=True, positive=True)
    r = symbols("r", real=True)

    # Results verified using Pow() resp. expand().
    # Also shows that lapply does expansion like expand() without deep=True
    p = Pow(c*d*f*In, c)  # all commutative, f positive, unknown exponent
    # Identity In is reduced to 1 by lapply though c, d are scalars. f is pulled out.
    assert lapply(p) == f**c*(c*d)**c
    p = Pow(c*d*A*In, -o)  # negative integer exponent. Fully expands.
    assert lapply(p) == c**(-o)*d**(-o)*A**(-o)
    p = Pow(c*d*A*In, Rational(1/2)) # fractional exponent. No extraction.
    assert lapply(p) == sqrt(c*d*A)
    p = Pow(c*d*f*A*In, Rational(1/2))
    # fractional exponent, f non-negative. expands() same
    assert lapply(p) == sqrt(f)*sqrt(c*d*A)
    p = Pow(c*d*f*A*In, c)  # arbitrary exponent. expand() same
    assert lapply(p) == f**c*(c*d*A)**c

    # numeric positive integer exponents; Pow() extracts same
    p = Pow(c*d*f*A*In, 4)*Pow(d*g*A*In, 3)
    assert lapply(p) == c**4*d**7*f**4*g**3*A**7
    # symbolic positive integer exponents; expand() extracts same
    p = Pow(c*d*f*A*In, o+4)*Pow(d*g*In*A, o)
    assert lapply(p) == c**(o + 4)*d**o*d**(o + 4)* \
                                                f**(o + 4)*g**o*A**(2*o + 4)
    # symbolic positive and negative integer exponents; expand extracts same
    p = Pow(c*d*f*A*In, o+4)*Pow(d*g*In*A, o-3)
    assert lapply(p) == c**(o + 4)*d**(o - 3)*d**(o + 4)* \
                                          f**(o + 4)*g**(o - 3)*A**(2*o + 1)
    # symbolic real exponents non-integer; expand extracts same
    p = Pow(d*f*A*In, o+r)*Pow(d*g*In*A, o-r)
    assert lapply(p) == d**(2*o) * f**(o + r) * g**(o - r) * A**(2*o) #2023-01-31
    # symbolic real exponents; expand extracts same
    p = Pow(d*f*A*In, r)*Pow(d*g*A*In, -r)
    assert lapply(p) == f**r*g**(-r)
    # numeric rational exp; expand extracts same
    p = Pow(c*A*In, 1+Rational(1,2))*Pow(d*A, Rational(1,2))
    assert lapply(p) == c*sqrt(c*A)*A*sqrt(d*A)

    # positive integer exponents, nested pos integer exp; expand extracts same
    p = Pow(Pow(c*f*A*In, o+3)*Pow(d*g*B*In, 3), 7)
    assert lapply(p) == c**(7*o + 21)*d**21*f**(7*o + 21)*g**21* \
                                                        (A**(o + 3)*B**3)**7
    # positive integer exponents, nested pos integer exp
    p = Pow(Pow(c*f*A*In, o+3)*Pow(d*g*A, 3), 7)
    assert lapply(p) == c**(7*o + 21)*d**21*f**(7*o + 21)*g**21*A**(7*o + 42)
    # pos int exp, unknown int exp, nested pos int exp: expand extracts same
    p = Pow(Pow(c*f*A*In, o+3)*Pow(d*g*A, 3-o), o+2)
    assert lapply(p) == c**(o**2 + 5*o + 6)*d**(-o**2 + o + 6)* \
                        f**(o**2 + 5*o + 6)*g**(-o**2 + o + 6)*A**(6*o + 12)
    assert lapply(p, nested_exp=False) == c**((o + 2)*(o + 3))* \
      d**((3 - o)*(o + 2))*f**((o + 2)*(o + 3))*g**((3 - o)*(o + 2))*A**(6*o + 12)

    # pos int exp, real exp, nested pos int exp: expand extracts same
    p = Pow(Pow(c*f*A*In, o+3)*Pow(d*g*A, r), o+2)
    assert lapply(p) == c**(o**2 + 5*o + 6)*f**(o**2 + 5*o + 6)* \
                        g**(o*r + 2*r)*(A**(o + 3)*(d*A)**r)**(o + 2)
    assert lapply(p, nested_exp=False) == c**((o + 2)*(o + 3))* \
        f**((o + 2)*(o + 3))*g**(r*(o + 2))*(A**(o + 3)*(d*A)**r)**(o + 2)

    # pos int exp, unknown int exp, nested real exp: expand extracts same
    p = Pow(Pow(c*f*A*In, o+3)*Pow(d*g*A, 3-o), r)
    assert lapply(p) == f**(o*r + 3*r)*(c**(o + 3)*d**(3 - o)*g**(3 - o)*A**6)**r
    assert lapply(p, nested_exp=False) == f**(r*(o + 3))* \
                                 (c**(o + 3)*d**(3 - o)*g**(3 - o)*A**6)**r

    # pos exp, unknown int exp, nested unknown exp: expand extracts same
    p = Pow(Pow(c*f*A*In, f)*Pow(d*g*A, 3-o), c)
    assert lapply(p) == f**(c*f)*(d**(3 - o)*g**(3 - o)*(c*A)**f*A**(3 - o))**c
    # unknown exp, unknown integer exp, nested unknown exp: no extraction
    p = Pow(Pow(c*f*A*In, d)*Pow(d*g*A, 3-o), c)
    assert lapply(p) == (d**(3 - o)*f**d*g**(3 - o)*(c*A)**d*A**(3 - o))**c
    # lapply won't expand the exponent, expand does
    p = Pow(c*f*A*In, f*(f+c)*(d-g))
    assert lapply(p) == f**(f*(c + f)*(d - g))*(c*A)**(f*(c + f)*(d - g))

    from sympy import evaluate, expand
    with evaluate(False): # creates objects in non-canonical form
        p0 = Pow(Pow(Pow(Pow(2*A, c), d), c), d)
        p1 = Pow(c, A) * Pow(c, B) # Pow of atrix "not implemented"
        p2 = Pow(c, A) * Pow(c, d) # c**A * c**d unsorted
        #p3 = Pow(p2, d*(A+B))
        #p4 = Pow(c, A*(A*c + d*In)*B)

    # prove coincidence with expand(). Using str() ignores
    # type differences Mul--MatMul, Pow--MatPow
    assert str(lapply(p0)) == str(expand(p0))
    assert lapply(p1) == expand(p1) # c**A * c**B, not a**(A+B)
    assert lapply(p2) == expand(p2.doit())
    #assert lapply(p3) == p3.doit() # Pow of matrix "not implemented"
    #assert lapply(p3, apply_exp=True) == expand(p3.doit())
    #assert lapply(p4, apply_exp=True) == c**(c*A**2*B + d*A*B)


def test_2023_02_11_MatSym():
    # Cases with Powers of MatrixSymbols and generic nc scalars
    A = symbols('A',commutative=False)
    c, d = symbols("c d", commutative=True)
    # expansion behaviour of x**(a+b) changed with 1.12, PR #23962
    f, g = symbols("f g", commutative=True, positive=True)
    o, u = symbols("o u", commutative=True, integer=True, positive=True)
    r = symbols("r", real=True)

    Ud = MatrixSymbol("Ud", 2, 2)
    Ud4 = Ud ** 4

    # Test expansion of exp ((c+d)*(f+g)) when apply_exp=True
    p = Pow(Ud, ((c+d)*(f+g))) * Ud**o * Ud**3
    q = lapply(p, apply_exp=True)
    assert q == Ud**(c*f + c*g + d*f + d*g + o + 3)

    p = Pow(A, ((c+d)*(f+g)), evaluate=False) * A**o * A**3
    q = lapply(p, apply_exp=True)
    assert q == A**(o + c*f + c*g + d*f + d*g + 3)

    # test nested powers
    p = Pow(Ud, o+1) * Pow(Pow(Ud, o+1), o)
    q = lapply(p, nexted_exp=True, apply_exp=True)
    assert q == Ud**(o**2 + 2*o - 3) * Ud4
    p = Pow(Ud, r) * Pow(Pow(Ud, r), o)
    q = lapply(p, nexted_exp=True)
    assert q == Ud**(o*r +r)

    # test apply_exp option
    p = Pow(Ud, 5*(o+1)*(u+1))
    q = lapply(p, apply_exp=True)
    assert q == Ud**(5*o*u+5*o+5*u+5)

    p = Pow(Ud, 5*(o+1)*(u+1))
    q = lapply(p) # exp remains (5*o+5)*(u+1)
    assert q == Ud**(5*(o+1)*(u+1)) # no expansion of exp

    from sympy import expand
    # some checks for commutative terms and power_xxx options:
    e = Pow(Pow(o*f*g, o+3)*Pow(f*g*o, 3-o), o+2)
    ep = lapply(e, power_exp=True, power_base=True)
    # lapply hands power_xxx options to expand(power_xxx=True)
    assert expand(e, mul=False, power_exp=True, power_base=True) == ep
    # but expand(power_exp=True) doesn't execute power_exp option here:
    # so the result is lengthy
    assert ep == f**((3 - o)*(o + 2))*f**((o + 2)*(o + 3))*g**((3 - o)*(o + 2))*\
           g**((o + 2)*(o + 3))*o**((3 - o)*(o + 2))*o**((o + 2)*(o + 3))

    e = Pow(Pow(o*f*g, o+3)*Pow(f*g*o, 3-o)+17, o+2)
    ep = lapply(e, power_exp=True, power_base=True)
    # lapply hands power_xxx options to expand(power_xxx=True)
    epx = expand(e, mul=False, power_exp=True, power_base=True)
    # and this time expand() executes power_exp option:
    # expansion behaviour of x**(a+b) changed with 1.12, PR #23962
    assert epx == (f**6*g**6*o**6 + 17)**o * \
                  (f**12*g**12*o**12 + 34*f**6*g**6*o**6 + 289)
    assert ep == (f**6*g**6*o**6 + 17)**2*(f**6*g**6*o**6 + 17)**o


def test_2023_03_matrix():
    ## real test cases
    a, b = symbols("a b")
    M1 = ImmutableMatrix([[a, 2], [3, 4]])
    A = MatrixSymbol("A", 2, 2)
    B = MatrixSymbol("B", 2, 2)
    Idb = ImmutableMatrix([[b, 0], [0, b]]) # instead of b*Identity(2)

    M2 = ImmutableMatrix([[0, sqrt(a*(a+1))], [sqrt(a*(a-1)), 0]])
    M3 = ImmutableMatrix([[0, a], [b, 0]]) # M3*M3=a*b*Identity(2)

    # lapply tries to identify multiples of the Identity and
    # returns a scalar in that case
    z = (2*M2.I.doit()) * M2 * A
    assert z.simplify() == z.expand() == Matrix([[2, 0], [0, 2]])*A
    assert lapply(z) == 2*A

    # lapply adopts usance to treat a*Identity as scalar a
    Da = DiagMatrix(Matrix([a, a])) # type: DiagMatrix
    assert lapply(Da) == a

    z = MatPow(M1 * M3, 10) # keeps unevaluated
    # gapply relies on native MatPow.doit() to evaluate powers, so
    # avoids unrolling of powers where possible
    assert z.doit().expand() == lapply(z).expand()

    assert M2.expand() != M2 # expands sqrt(a*(a+1)) to sqrt(a**2+a)
    assert lapply(M2) == M2  # no unnecessary expands()

    # lapply applies .doit(), .simplify() at intermediate stages of
    # computation, so may help to compute and simplify expressions
    # that become complicated at intermediate stages.
    z = M2.T * (A * M2.I).T
    assert z.doit() != A.T and z.expand().simplify().doit() == A.T
    assert lapply(z) == A.T  # one-stop-shop

    z = M1.T * (A*M1.I + B*M1.I).T # final result expected: A.T + B.T
    assert z.doit() != A.T + B.T and z.expand().simplify() != A.T + B.T
    assert z.expand().simplify().doit() == A.T + B.T
    assert lapply(z) == A.T + B.T # one-stop-shop

    z = (M2.T - Idb).I * (Idb - M2.T) + (M2)**2 + b*Idb
    assert not z.is_diagonal() and not z.doit().is_diagonal()
    d = z.simplify()
    assert d.is_diagonal() and simplify(d[0,0] - d[1,1]) == 0
    # lapply returns scalar instead of k*Identity(2)
    assert simplify(d[0,0] - lapply(z)) == 0

    z = MatPow(M2, 12)
    zz = z.doit().simplify()
    assert zz.is_diagonal()
    assert simplify(zz[0,0] - z[1,1]) == 0
    assert simplify(zz[0,0] - lapply(z)) == 0

    # lapply relies on capabilities of MatPow
    o = symbols("o", integer=True, nonnegative=True)
    z = MatPow(M1, o + 2) # unevaluated
    zz = M1**(o + 2) # computes via diagonalisation of M1!
    assert isinstance(zz, ImmutableMatrix) # complicated term
    assert isinstance(lapply(z), ImmutableMatrix) # complicated term

    # lapply tries to handle powers with MatrixSymbols in it,
    # e.g. does evaluation of integer parts of symbolic exponents
    z = MatPow((M1 + A), o + 2) # unevaluated
    assert unchanged(MatPow, (M1 + A), o + 2) # verified
    zz = (M1 + A)**(o + 2)      # unevaluated
    assert zz == z # verified
    F = (M1+A)**o
    q = lapply(z)
    assert q == F*M1**2 + F*A*M1 + F*M1*A + F*A**2

    # lapply uses shortcuts on some powers if base permits,
    # e.g. involutoric factors in base
    z = MatPow(M3 * A * M3, o+10) # unevaluated
    assert unchanged(MatPow, M3 * A * M3, o+10) # verified
    zz = (M3 * A * M3)**(o + 10)    # doesn't compute
    assert z.expand().simplify().doit() == z == zz # verified
    assert lapply(z) == MatMul(a**(o+9) * b**(o+9), M3, A**(o+10), M3)

    z = MatPow(M3 * A * M3, 1000) # unevaluated
    assert unchanged(MatPow, M3 * A * M3, 1000) # verified
    zz = (M3 * A * M3)**1000      # doesn't compute
    assert z.doit() == z.expand() == z.simplify() == z == zz
    assert lapply(z) == MatMul(a**999 * b**999, M3, A**1000, M3)

    # lapply uses shortcuts on some powers if base permits,
    # e.g. base change W * A * W**-1
    z = MatPow(M1 * A * M1.I * Idb, 1000)    # unevaluated
    assert unchanged(MatPow, M1 * A * M1.I * Idb, 1000) # verified
    zz = (M1 * A * M1.I * Idb)**1000         # doesn't compute
    assert z.doit() == z.expand() == z == zz # verified
    assert (M1 * A * M1.I.simplify() * Idb)**1000 == zz.simplify()
    assert lapply(z) == MatMul(b**999, M1, A**1000, (b*M1.I.simplify()))

    # Example: lapply rolls offs the power because M2*Idb*M1
    # computes into one matrix F, while MatPow keeps it unevaluated
    z = MatPow(M1 * A * M2 * Idb, 4) # unevaluated
    assert unchanged(MatPow, M1 * A * M2 * Idb, 4) # verified
    zz = (M1 * A * M2 * Idb)**4      # doesn't compute
    assert (M1 * A * (M2*Idb).expand())**4 == zz.expand()
    assert z.doit() == zz.simplify() == z == zz # verfied
    F = M2 * Idb * M1
    assert lapply(z) == M1 * A*F* A*F* A*F* A * M2*Idb

# doc examples for lapply
def test_lapply_doc():

    a, b, c, d, x = symbols("a b c d x", commutative=True)
    M1 = ImmutableMatrix([[a, 2], [3, 4]])
    M2 = ImmutableMatrix([[0, sqrt(a*(a+1))], [sqrt(a*(a-1)), 0]])
    M3 = ImmutableMatrix([[0, a], [b, 0]]) # M3*M3=a*b*Identity(2)

    A = MatrixSymbol("A", 2, 2)
    B = MatrixSymbol("B", 2, 2)
    Idb = ImmutableMatrix([[b, 0], [0, b]]) # instead of b*Identity(2)
    o = symbols("o", integer=True, nonnegative=True, even=True)

    # Replace multiples of the Identity by scalar
    z = (2*M2.I.doit()) * M2 * A
    assert z.simplify() == Matrix([[2, 0], [0, 2]]) * A
    assert str(lapply(z)) == "2*A"
    # Evaluate term 1
    z = M2.T * (A * M2.I.doit()).T
    assert len(str(z.doit())) > 120 # complex term
    assert z.expand().doit() == A.T
    assert str(lapply(z)) == "A.T"
    # Evaluate term 2  IN DOCS
    Z = M1.T * (A*M1.I + B*M1.I).T # final result expected: A.T + B.T
    assert len(str(Z.expand())) > 250 # complex term
    assert len(str(Z.expand().simplify())) > 50 # still complex
    assert str(Z.expand().simplify().doit()) == "A.T + B.T"
    assert str(lapply(Z)) == "A.T + B.T"

    # Compute a large term without manual simplification steps IN DOCS
    # (expand.simplify.doit can't handle as they don't simplify multiples
    # of Identity!)
    # Example may not use powers as these are computed by MatPow within
    # lapply, and not just varying factors as these are computed in
    # lapply1_types(Mul). It must be factors that Mul resp. MatMul leave
    # unevaluated. This can only be MatrixSymbols.
    # Note that Z - A.T - B.T == 0 as seen above. But lapply, expand,
    # simplify etc. just expand the term. Key for simplification is reducing
    # multiples of the identity to scalars:
    p = Mul(*[(k * Identity(2) + Z - A.T - B.T) for k in range(1,3)])
    assert str(lapply(p*A)) == "2*A"
    pp = (p*A).expand().simplify().doit() # complex term with k*Identity
    assert str(lapply(pp)) == "2*A"

    # involutoric short cut; multiples of Identity become scalar IN DOCS
    assert str(lapply(M3**(o+10000))) == "a**(o/2 + 5000)*b**(o/2 + 5000)"

    # Exploit involutoric factors in base IN DOCS
    Z = MatPow(M3 * A * M3, o+10000)   # unevaluated
    zz = (M3 * A * M3)**(o + 10000)    # doesn't compute
    assert unchanged(MatPow, M3 * A * M3, o+10000)
    assert Z == zz == Z.expand() == Z.simplify() == Z.doit()
    assert lapply(Z) == \
        MatMul(a**(o+9999) * b**(o+9999), M3, A**(o+10000), M3)

    # uses shortcuts on powers if base permits, in DOCS
    # e.g. base change W * A * W**-1
    Z = MatPow(M1 * A * M1.I * Idb, 10000)    # unevaluated
    zz = (M1 * A * M1.I * Idb)**10000         # doesn't compute
    assert unchanged(MatPow, M1 * A * M1.I * Idb, 10000)
    assert Z == zz == Z.expand() == Z.doit()
    #assert (M1 * A * M1.I.simplify() * Idb)**1000 == zz.simplify() # no real progress
    assert lapply(Z) == \
            MatMul(b**9999, M1, A**10000, (b*M1.I.simplify()))

    # Example: lapply rolls offs the power because M2*Idb*M1
    # computes into one matrix F, while MatPow keeps it unevaluated
    Z = MatPow(M1 * A * M2 * Idb, 4) # unevaluated
    zz = (M1 * A * M2 * Idb)**4      # doesn't compute
    assert unchanged(MatPow, M1 * A * M2 * Idb, 4)
    assert (M1 * A * (M2*Idb).expand())**4 == zz.expand() #no real progress
    assert Z == zz == zz.simplify() == Z.doit()
    F = M2 * Idb * M1
    assert lapply(Z) == M1 * A*F* A*F* A*F* A * M2*Idb

    # Matrix representation of Quaternion(a,b,c,d)
    Mq = ImmutableMatrix([[a, -b, -c, -d], [b, a, -d, c],
                          [c,  d,  a, -b], [d, -c, b, a]])
    cp = Mq.charpoly().eval(x)
    # any matrix annuls its characteristic polynomial
    assert str(lapply(cp.subs(x, Mq))) == "0"

    # As an example, lapply may help to work with Quaternions:
    # Note: Don't use complex numbers, as they are flagged .is_commutative=True
    # but they don't commute with Quaternions over the Reals. And don't use
    # Quaternions over complex numbers with lapply as they are not associative!
    from sympy.algebras import Quaternion
    q = Quaternion(a, b, c, d)

    quatAdd = (lambda *summands, evaluate=True: sum(summands).expand())
    assert str(lapply(cp.subs(x, q), Add=quatAdd)) == "0 + 0*i + 0*j + 0*k"

    # Example for c_nc_ncef: identify real quaternions as commutative
    c_nc_ncef.add((Quaternion,), (lambda q, to_L, **options:
        ([q.scalar_part()], [], S.One) if q.vector_part().is_zero_quaternion()
         else ([],[], q) ))
    assert str(lapply(q * q.conjugate())) == "a**2 + b**2 + c**2 + d**2"
