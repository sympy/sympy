from sympy.core.compatibility import as_int, is_sequence
from sympy.core.numbers import oo, mod_inverse
from sympy.core.relational import Eq
from sympy.core.symbol import symbols
from sympy.polys.domains import FiniteField, QQ, RationalField, FF
from sympy.solvers.solvers import solve
from sympy.polys import poly
from .factor_ import divisors
from .residue_ntheory import polynomial_congruence



    
class EllipticCurve:
    """
    Create the following Elliptic Curve over domain.

    `y^{2} + a_{1} x y + a_{3} y = x^{3} + a_{2} x^{2} + a_{4} x + a_{6}`

    The default domain is ``QQ``. If no coefficient ``a1``, ``a2``, ``a3``,
    it create curve as following form.

    `y^{2} = x^{3} + a_{4} x + a_{6}`

    Examples
    ========

    References
    ==========

    [1] J. Silverman "A Friendly Introduction to Number Theory" Third Edition
    [2] http://mathworld.wolfram.com/EllipticDiscriminant.html
    [3] G. Hardy, E. Wright "An Introduction to the Theory of Numbers" Sixth Edition

    """

    def __init__(self, a4, a6, a1=0, a2=0, a3=0, modulus = 0):
        if modulus == 0:
            domain = QQ
        else:
            domain = FF(modulus)
        a1, a2, a3, a4, a6 = map(domain.convert, (a1, a2, a3, a4, a6))
        self._domain = domain
        self.modulus = modulus
        # Calculate discriminant
        b2 = a1**2 + 4 * a2
        b4 = 2 * a4 + a1 * a3
        b6 = a3**2 + 4 * a6
        b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
        self._b2, self._b4, self._b6, self._b8 = b2, b4, b6, b8
        self._discrim = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
        self._a1 = a1
        self._a2 = a2
        self._a3 = a3
        self._a4 = a4
        self._a6 = a6
        x, y, z = symbols('x y z')
        self.x, self.y, self.z = x, y, z
        self._eq = Eq(y**2*z + a1*x*y*z + a3*y*z**2, x**3 + a2*x**2*z + a4*x*z**2 + a6*z**3)
        if isinstance(self._domain, FiniteField):
            self._rank = 0
        elif isinstance(self._domain, RationalField):
            self._rank = None

        #For Division polynomials
        X = symbols('X')
        Y = symbols('Y')
        div = {}
        if domain != QQ:
            self.X = X
            self.Y = Y
            div[-1] = poly(-1, X, Y, modulus=modulus)
            div[0] = poly(0, X, Y, modulus=modulus)
            div[1] = poly(1, X, Y, modulus=modulus)
            div[2] = poly(2*Y, X, Y, modulus=modulus)
            div[3] = poly(3*X**4 + 6*a4*X**2 + 12*a6*X - a4*a4, X, Y, modulus=modulus)
            div[4] = poly(4*Y*(X**6 + 5*a4*X**4 + 20*a6*X**3 - 5*a4*a4*X**2 - 4*a4*a6*X - 8*a6**2 - a4**3), modulus=modulus)
            self.div = div

    def __call__(self, x, y, z=1):
        return EllipticCurvePoint(x, y, z, self)

    def __contains__(self, point):
        if is_sequence(point):
            if len(point) == 2:
                z1 = 1
            else:
                z1 = point[2]
            x1, y1 = point[:2]
        elif isinstance(point, EllipticCurvePoint):
            x1, y1, z1 = point.x, point.y, point.z
        else:
            raise ValueError('Invalid point.')
        if self.characteristic == 0 and z1 == 0:
            return True
        return self._eq.subs({self.x: x1, self.y: y1, self.z: z1})

    def __repr__(self):
        return 'E({}): {}'.format(self._domain, self._eq)

    def minimal(self):
        """
        Return minimal Weierstrass equation.

        Examples
        ========

        >>> from sympy.ntheory.elliptic_curve import EllipticCurve

        >>> e1 = EllipticCurve(-10, -20, 0, -1, 1)
        >>> e1.minimal()
        E(QQ): Eq(y**2*z, x**3 - 13392*x*z**2 - 1080432*z**3)

        """
        char = self.characteristic
        if char == 2:
            return self
        if char == 3:
            return EllipticCurve(self._b4/2, self._b6/4, a2=self._b2/4, modulus=self.modulus)
        c4 = self._b2**2 - 24*self._b4
        c6 = -self._b2**3 + 36*self._b2*self._b4 - 216*self._b6
        return EllipticCurve(-27*c4, -54*c6, modulus=self.modulus)

    def points(self):
        """
        Return points of curve over Finite Field.

        Examples
        ========

        >>> from sympy.ntheory.elliptic_curve import EllipticCurve
        >>> e2 = EllipticCurve(1, 1, 1, 1, 1, modulus=5)
        >>> e2.points()
        {(0, 2), (1, 4), (2, 0), (2, 2), (3, 0), (3, 1), (4, 0)}

        """

        char = self.characteristic
        all_pt = set()
        if char >= 1:
            for i in range(char):
                congruence_eq = ((self._eq.lhs - self._eq.rhs).subs({self.x: i, self.z: 1}))
                sol = polynomial_congruence(congruence_eq, char)
                for num in sol:
                    all_pt.add((i, num))
            return all_pt
        else:
            raise ValueError("Infinitely many points")

    def points_x(self, x):
        "Returns points on with curve where xcoordinate = x"
        pt = []
        if self._domain == QQ:
            for y in solve(self._eq.subs(self.x, x)):
                    pt.append((x, y))
        congruence_eq = ((self._eq.lhs - self._eq.rhs).subs({self.x: x, self.z: 1}))
        for y in polynomial_congruence(congruence_eq, self.characteristic):
            pt.append((x, y))
        return pt

    def torsion_points(self):
        """
        Return torsion points of curve over Rational number.

        Return point objects those are finite order.
        According to Nagell-Lutz theorem, torsion point p(x, y)
        x and y are integers, either y = 0 or y**2 is divisor
        of discriminent. According to Mazur's theorem, there are
        at most 15 points in torsion collection.

        Examples
        ========

        >>> from sympy.ntheory.elliptic_curve import EllipticCurve
        >>> e2 = EllipticCurve(-43, 166)
        >>> sorted(e2.torsion_points())
        [(-5, -16), (-5, 16), O, (3, -8), (3, 8), (11, -32), (11, 32)]

        """
        if self.characteristic > 0:
            raise ValueError("No torsion point for Finite Field.")
        l = [EllipticCurvePoint.point_at_infinity(self)]
        for xx in solve(self._eq.subs({self.y: 0, self.z: 1})):
            if xx.is_rational:
                l.append(self(xx, 0))
        for i in divisors(self.discriminant, generator=True):
            j = int(i**.5)
            if j**2 == i:
                for xx in solve(self._eq.subs({self.y: j, self.z: 1})):
                    if not xx.is_rational:
                        continue
                    p = self(xx, j)
                    if p.order() != oo:
                        l.extend([p, -p])
        return l

    @property
    def characteristic(self):
        """
        Return domain characteristic.

        Examples
        ========

        >>> from sympy.ntheory.elliptic_curve import EllipticCurve
        >>> e2 = EllipticCurve(-43, 166)
        >>> e2.characteristic
        0

        """
        return self._domain.characteristic()

    @property
    def discriminant(self):
        """
        Return curve discriminant.

        Examples
        ========

        >>> from sympy.ntheory.elliptic_curve import EllipticCurve
        >>> e2 = EllipticCurve(0, 17)
        >>> e2.discriminant
        -124848

        """
        return int(self._discrim)

    @property
    def is_singular(self):
        """
        Return True if curve discriminant is equal to zero.
        """
        return self.discriminant == 0

    @property
    def j_invariant(self):
        """
        Return curve j-invariant.

        Examples
        ========

        >>> from sympy.ntheory.elliptic_curve import EllipticCurve
        >>> e1 = EllipticCurve(-2, 0, 0, 1, 1)
        >>> e1.j_invariant
        1404928/389

        """
        c4 = self._b2**2 - 24*self._b4
        return self._domain.to_sympy(c4**3 / self._discrim)

    @property
    def order(self):
        """
        Number of points in Finite field.

        Examples
        ========

        >>> from sympy.polys.domains import FF
        >>> from sympy.ntheory.elliptic_curve import EllipticCurve
        >>> e2 = EllipticCurve(1, 0, modulus=19)
        >>> e2.order
        19

        """
        if self.characteristic == 0:
            raise NotImplementedError("Still not implemented")
        return len(list(self.points()))

    @property
    def rank(self):
        """
        Number of independent points of infinite order.

        For Finite field, it must be 0.
        """
        if self._rank is not None:
            return self._rank
        raise NotImplementedError("Still not implemented")

    def div_poly(self, n):
        """
        Division Polynomials of an elliptic curve over finite fields.

        """
        from sympy import simplify
        if self._domain == QQ:
            raise ValueError("Only defined for curves over finite fields")
        if n in self.div:
            return self.div[n]
        f = lambda a : self.div_poly(a)
        g = lambda a : poly(a.subs(self.Y**2, self.X**3 + self._a4*self.X + self._a6), self.X, self.Y, modulus=self.modulus)
        if n % 2 == 0:
            n //= 2
            self.div[2*n] = f(n) * (f(n + 2)*f(n - 1)**2 - f(n - 2)*f(n + 1)**2)
            self.div[2*n] = simplify((self.div[2*n] * mod_inverse(2, self.modulus)) / self.Y)
            self.div[2*n] = g(self.div[2*n])
            return self.div[2*n]
        n = (n - 1) // 2
        self.div[2*n + 1] = poly(g(f(n + 2)*g(g(f(n)**2)*f(n))) - g(g(g(f(n + 1)**2)*f(n + 1))*f(n - 1)), self.X, self.Y, modulus=self.modulus)
        self.div[2*n + 1] = g(self.div[2*n + 1])
        return self.div[2*n + 1]

    def schoof(self):
        """
        Computes the order of the curve using Schoof's Algorithm.

        References
        ==========

        [1] http://jjmcgee.asp.radford.edu/SchoofsAlgorithm06.pdf

        """
        from sympy.ntheory import prime
        from sympy import sqrt
        from sympy.polys.polytools import gcd
        l = 1
        i = 1
        p = self.modulus
        X = self.X
        Y = self.Y
        L = []
        m = 1
        f2 = {}
        y2 = poly(X**3 + self._a4*X+self._a6, X, modulus=p)
        while(m <= 4*sqrt(p)):
            l = prime(i)
            if p % l > 1 or l == 2:
                m *= l
                L.append(l)
            i += 1

        def exp(a, e, mod):
            """
            Binary Exponentiation

            """
            if e == 1:
                return a % mod
            if e == 0:
                return poly(1, X, Y, modulus = p) % mod
            if e % 2 == 0:
                return exp(a, e // 2 , mod)**2 % mod
            return exp(a, e - 1, mod)*a % mod

        def eq16(l, short_com):
            k = p % l
            p1 = short_com["p3"]
            if k not in f2:
                f2[k] = f(k)**2
            if k % 2 == 0:
                p_16 = (p1 - poly(X, X, modulus=p))*f2[k]*y2 + short_com["p1"]
            else:
                p_16 = (p1 - poly(X,X,modulus=p))*f2[k] + short_com["p1"]*y2
            return  p_16 % f(l)


        def f(k):
            if k % 2 == 0:
                return poly((self.div_poly(k).subs(Y, 1)), X, Y, modulus=p)
            return self.div_poly(k)

        def eq17(l, w):
            fl = f(l)
            p1 = exp(poly(X, X, modulus=p), p, fl)
            if w not in f2:
                f2[w] = f(w)**2
            p2 = (p1 - poly(X, X, modulus=p)) * f2[w] % fl
            p3 = f(w - 1)*f(w + 1)
            if w % 2 == 0:
                p17 = p2*y2 + p3
            else:
                p17 = p2 + p3*y2
            return p17 % fl

        def eq18(l, w):
            fl = f(l)
            if w % 2 == 0:
                k = (p + 3) // 2
            else:
                k = (p - 1) // 2
            p1 = exp(y2, k, fl)
            if w not in f2:
                f2[w] = f(w)**2
            p2 = f(w) * (f2[w] % fl) % fl
            if w-1 not in f2:
                f2[w - 1] = f(w - 1)**2
            if w + 1 not in f2:
                f2[w + 1] = f(w + 1)**2
            p3 = f(w + 2) * (f2[w - 1] % fl) % fl
            p4 = f(w - 2) * (f2[w + 1] % fl) % fl
            p18y = (poly(4, X, modulus=p)*p1*p2 - p3 + p4) % fl
            return p18y

        def alpha_cal(l, short_com):
            k = p % l
            fl = f(l)
            y4 = y2**2 % fl
            if k - 1 not in f2:
                f2[k - 1] = f(k - 1)**2
            p1 = f2[k - 1]%fl
            p1 = f(k + 2)*p1 % fl
            if k + 1 not in f2:
                f2[k + 1] = f(k + 1)**2
            p2 = f2[k + 1]%fl
            p2 = f(k - 2)*p2 % fl
            p3 = short_com["y2_p2"]
            p4 = short_com["p5"]*f(k) % fl
            if k % 2 == 0:
                p5 = poly(4, X, modulus=p)*y4*p3*p4 % fl
                alpha = (p1 - p2 - p5)%fl
            else:
                p5 = poly(4, X, modulus=p)*y2*p3*p4 % fl
                alpha = (y2*(p1 - p2) - p5) % fl
            return alpha

        def beta_cal(l, short_com):
            k = p % l
            fl = f(l)
            p1 = short_com["p3"]
            p1 = poly(X, X, modulus=p) - p1
            p2 = short_com["p1"]
            p3 = short_com["p5"]
            if k % 2 == 0:
                p4 = (p1*y2*p3 - p2 ) % fl
                beta = poly(4, X, modulus=p)*y2*f(k)*p4 % fl
            else:
                p4 = (p1*p3 - y2*p2 ) % fl
                beta = poly(4, X, modulus=p)*f(k)*p4 % fl
            return beta

        def com19x(l, tau, alpha, beta, short_com):
            fl = f(l)
            k = p % l
            p1 = short_com["p1"]
            p5 = short_com["p5"]
            p6 = f(tau - 1)*f(tau + 1) % fl
            p6 = exp(p6, p, fl)
            if k % 2 == 0:
                alpha2 = y2*short_com["alpha2"] % fl
                beta2 = short_com["beta2"] % fl
                p7 = (p1- y2*short_com["p5p4"]) % fl
                p7 = (p7*beta2 + y2*p5*alpha2) % fl
            else:
                alpha2 = short_com["alpha2"] % fl
                beta2 = y2*short_com["beta2"] % fl
                p7 = (y2*p1 - short_com["p5p4"]) % fl
                p7 = (p7*beta2 + p5*alpha2) % fl
            if tau not in f2:
                f2[tau] = f(tau)**2
            p8 = exp(f2[tau] % fl, p, fl)
            if k % 2 == 0 and tau % 2 == 0:
                p9 = short_com["p9"]
                p19 = (p7*p9*p8 + p6*beta2*y2*p5) % fl
            if k % 2 ==0 and tau % 2 == 1:
                p9 = short_com["p9"]*y2 %fl
                p19 = (p7*p8 + p9*p6*beta2*p5) % fl
            if k % 2 == 1 and tau % 2 == 0:
                p9 = short_com["p9"]
                p19 = (p7*p8*p9 + p6*beta2*p5) % fl
            if k % 2 == 1 and tau % 2 == 1:
                p9 = short_com["p9"]
                p19 = (p7*p8 + p9*p6*beta2*p5 ) % fl
            return p19

        def com19y(l, tau, alpha, beta, short_com):
            fl = f(l)
            k = p % l
            a3 = alpha*short_com["alpha2"] % fl
            b2 = short_com["beta2"] % fl
            b3 = short_com["beta2"]*beta % fl
            p1 = alpha*b2*short_com["p1"] % fl
            p2 = short_com["p5"]
            p3 = short_com["p3"]
            p3 = (poly(2, X, modulus=p)*p3 + poly(X, X, modulus=p)) % fl
            p4 =  poly(4, X, modulus=p) * exp(f(tau), 3*p, fl) % fl
            if tau - 1 not in f2:
                f2[tau - 1] = f(tau - 1)**2
            p5 = f(tau + 2) * (f2[tau - 1] % fl) % fl
            if tau + 1 not in f2:
                f2[tau + 1] = f(tau + 1)**2
            p5 = (p5 - f(tau - 2) * (f2[tau + 1]%fl)) % fl
            p5 = exp(p5, p, fl)
            p6 = short_com["y2_p2"]
            if k % 2 == 0:
                p19 = (p3*alpha*b2 - b3*p6 - y2*a3) % fl
                p19 = p4*(p19*y2*p2 - p1) % fl
                if tau % 2 == 0:
                    p19 *= exp(y2, (3*p + 1) // 2, fl)
                    p19 = (p19 - b3*y2*p2*p5) % fl
                else:
                    yq = exp(y2, (p + 1) // 2, fl)
                    p19 = (p19 - b3*p2*p5*yq) % fl
            else:
                y4 = y2**2 % fl
                p7 = b3*y4*p6 % fl
                p19 = (p3*alpha*b2*y2 - p7- a3) % fl
                p19 = (p4*(p19*p2 - y4*p1)) % fl
                if tau % 2 == 0:
                    p19 *= exp(y2, (3*p - 3) // 2, fl)
                    p19 = (p19 - b3*p2*p5) % fl
                else:
                    yq = exp(y2, (p + 3) // 2, fl)
                    p19 = (p19 - b3*p2*p5*yq) % fl
            return p19

        def case1(l):
            from sympy.ntheory import is_quad_residue, sqrt_mod
            if not is_quad_residue(p, l):
                return 0
            w = sqrt_mod(p, l)
            fl = f(l)
            p17 = eq17(l, w)
            g = gcd(p17, fl)
            if g == 1:
                return 0
            p18 = eq18(l, w)
            g = gcd(p18, fl)
            if g == 1:
                t1 = -2*w % l
            else:
                t1 = 2*w % l
            return t1

        def tmod2(p):
            xp = exp(poly(X, X, modulus=p), p, y2)
            g = gcd(xp-poly(X, X, modulus=p) , y2)
            if g == 1:
                return 1
            return 0
        def _schoof():
            from sympy.ntheory.modular import solve_congruence
            t1 = tmod2(p)
            t = {2 : t1}
            for i in range(1, len(L)):
                l = L[i]
                fl= f(l)
                short_com = {}
                k = p % l
                short_com["p1"] = f(k - 1)*f(k + 1) % fl
                short_com["p2"] = exp(poly(X, X, modulus=p), p, fl)
                short_com["p3"] = exp(short_com["p2"], p, fl)

                p16 = eq16(l, short_com)
                g = gcd(p16, fl)
                if g != 1:
                    t1 = case1(l)
                    t[l] = t1
                else:
                    if k not in f2:
                        f2[k] = f2(k)**2
                    short_com["p5"] = f2[k] % fl
                    short_com["y2_p2"] = exp(y2, (p**2 - 1) // 2, fl)
                    alpha = alpha_cal(l, short_com)
                    beta = beta_cal(l, short_com)
                    short_com["p4"] = short_com["p3"] + short_com["p2"] + poly(X, X, modulus=p)

                    short_com["p9"] = exp(y2, p, fl)
                    short_com["alpha2"] = alpha**2 % fl
                    short_com["beta2"] = beta**2 % fl
                    short_com["p5p4"] = short_com["p4"]*short_com["p5"] % fl
                    for tau in range(1, (l + 1) // 2):
                        p19x = com19x(l, tau, alpha, beta, short_com)
                        g = gcd(p19x, fl)
                        if g != 1:
                            p19y = com19y(l, tau ,alpha, beta, short_com)
                            gy = gcd(p19y, fl)
                            if gy != 1 :
                                t1 = tau
                            else:
                                t1 = l - tau
                            t[l] = t1
                            break

            a = []
            b = []
            for i in t:
                a.append(i)
                b.append(t[i])
            ans = solve_congruence(*zip(b, a))
            fans = ans[0]
            if  ans[0] > 2*sqrt(p):
                fans -= ans[1]
            return p + 1 - fans

        return _schoof()


class EllipticCurvePoint:
    """
    Point of Elliptic Curve

    Examples
    ========

    >>> from sympy.ntheory.elliptic_curve import EllipticCurve
    >>> e1 = EllipticCurve(-17, 16)
    >>> p1 = e1(0, -4, 1)
    >>> p2 = e1(1, 0)
    >>> p1 + p2
    (15, -56)
    >>> e3 = EllipticCurve(-1, 9)
    >>> e3(1, -3) * 3
    (664/169, 17811/2197)
    >>> (e3(1, -3) * 3).order()
    oo
    >>> e2 = EllipticCurve(-2, 0, 0, 1, 1)
    >>> p = e2(-1,1)
    >>> q = e2(0, -1)
    >>> p+q
    (4, 8)
    >>> p-q
    (1, 0)
    >>> 3*p-5*q
    (328/361, -2800/6859)
    """

    @staticmethod
    def point_at_infinity(curve):
        return EllipticCurvePoint(0, 1, 0, curve)

    def __init__(self, x, y, z, curve):
        dom = curve._domain.convert
        self.x = dom(x)
        self.y = dom(y)
        self.z = dom(z)
        self._curve = curve
        self._domain = self._curve._domain
        if not self._curve.__contains__(self):
            raise ValueError("The curve does not contain this point")

    def __add__(self, p):
        if self.z == 0:
            return p
        if p.z == 0:
            return self
        x1, y1 = self.x/self.z, self.y/self.z
        x2, y2 = p.x/p.z, p.y/p.z
        a1 = self._curve._a1
        a2 = self._curve._a2
        a3 = self._curve._a3
        a4 = self._curve._a4
        a6 = self._curve._a6
        if x1 != x2:
            slope = (y1 - y2) / (x1 - x2)
            yint = (y1 * x2 - y2 * x1) / (x2 - x1)
        else:
            if (y1 + y2) == 0:
                return self.point_at_infinity(self._curve)
            slope = (3 * x1**2 + 2*a2*x1 + a4 - a1*y1) / (a1 * x1 + a3 + 2 * y1)
            yint = (-x1**3 + a4*x1 + 2*a6 - a3*y1) / (a1*x1 + a3 + 2*y1)
        x3 = slope**2 + a1*slope - a2 - x1 - x2
        y3 = -(slope + a1) * x3 - yint - a3
        return self._curve(x3, y3, 1)

    def __lt__(self, other):
        return (self.x, self.y, self.z) < (other.x, other.y, other.z)

    def __mul__(self, n):
        n = as_int(n)
        r = self.point_at_infinity(self._curve)
        if n == 0:
            return r
        if n < 0:
            return -self * -n
        p = self
        while n:
            if n & 1:
                r = r + p
            n >>= 1
            p = p + p
        return r

    def __rmul__(self, n):
        return self * n

    def __neg__(self):
        return EllipticCurvePoint(self.x, -self.y - self._curve._a1*self.x - self._curve._a3, self.z, self._curve)

    def __repr__(self):
        if self.z == 0:
            return 'O'
        dom = self._curve._domain
        try:
            return '({}, {})'.format(dom.to_sympy(self.x), dom.to_sympy(self.y))
        except TypeError:
            pass
        return '({}, {})'.format(self.x, self.y)

    def __sub__(self, other):
        return self.__add__(-other)

    def order(self):
        """
        Return point order n where nP = 0.

        """
        if self.z == 0:
            return 1
        if self.y == 0:  # P = -P
            return 2
        p = self * 2
        if p.y == -self.y:  # 2P = -P
            return 3
        i = 2
        if self._domain != QQ:
            while int(p.x) == p.x and int(p.y) == p.y:
                p = self + p
                i += 1
                if p.z == 0:
                    return i
            return oo
        while p.x.numerator == p.x and p.y.numerator == p.y:
            p = self + p
            i += 1
            if i > 12:
                return oo
            if p.z == 0:
                return i
        return oo
