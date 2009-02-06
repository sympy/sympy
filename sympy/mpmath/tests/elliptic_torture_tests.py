#!/usr/bin/env python
"""
    Tests of the elliptic functions module.  If all unittests are run
    at 100 digit precision, it takes ~ 12 m on a Athlon 64 X2 3800+.

    Author: M.T. Taschuk

    References:

    [1] Abramowitz & Stegun. 'Handbook of Mathematical Functions, 9th Ed.',
        (Dover duplicate of 1972 edition)
    [2] Whittaker 'A Course of Modern Analysis, 4th Ed.', 1946,
        Cambridge Univeristy Press

"""
__version__ = '$Id:$'

import unittest
#import mpmath.mptypes
import mpmath
import random

from sympy.mpmath.elliptic import *

class precisemathTests(unittest.TestCase):

    def testCalculateNome(self):
        mpmath.mpf.dps = 100
        mpmath.mp.dps = 100
        testlimit = mpmath.mpf('10')**(-1*(mpmath.mpf.dps - 2))
        #print >> sys.stderr, testlimit

        q = calculate_nome(mpmath.mpf('0'))
        self.assertEquals(mpmath.mpf('0'), mpmath.mpf('0'))

        mathematica = [ (0.1,   0.00658465),
                        (0.2,   0.0139429),
                        (0.3,   0.0222774),
                        (0.4,   0.0318833),
                        (0.5,   0.0432139),
                        (0.6,   0.0570203),
                        (0.7,   0.0746899),
                        (0.8,   0.0992737),
                        (0.9,   0.140173),
                        (0.99,  0.262196)]

        for i in mathematica:
            m = mpmath.mpf(i[0])
            value = calculate_nome(m.sqrt())
            self.assertEquals(round(i[1], 6), round(value, 6))

    def testJacobiTheta1(self):
        mpmath.mpf.dps = 100
        mpmath.mp.dps = 100
        testlimit = mpmath.mpf('10')**(-1*(mpmath.mpf.dps - 2))
        #print >> sys.stderr, testlimit

        z = mpmath.mpf('0')
        q = mpmath.mpf('0')
        value = jacobi_theta_1(z, q)
        self.assertEquals(value, mpmath.mpf('0'))

        z = mpmath.mpf('0')
        m = mpmath.pi
        self.assertRaises(ValueError, jacobi_theta_1, z, m)
        m = mpmath.mpf('1')
        self.assertRaises(ValueError, jacobi_theta_1, z, m)

        # Mathematica value for v1(u = 0.1, q = 0.1) = 0.108958
        # q = 0.1, therefore m = 0.802403, according to Mathematica
        z = mpmath.mpf('0.1')
        m = mpmath.mpf('0.802403')

        result = jacobi_theta_1(z, m)
        self.assertEquals(round(result, 6), 0.108958)
        self.assertTrue(isinstance(result, mpmath.mpf))

        z = mpmath.pi                           # test for sin zeros
        for i in range(10):
            qstring = str(random.random())
            q = mpmath.mpf(qstring)

            result = jacobi_theta_1(z, q)
            if result < testlimit:
                self.assertEquals(True, True)
            else:
                print >> sys.stderr, 'Result (~ 0): %e' % result
                self.assertEquals(False, True)

    def testJacobiTheta2(self):
        mpmath.mpf.dps = 100
        mpmath.mp.dps = 100
        testlimit = mpmath.mpf('10')**(-1*(mpmath.mpf.dps - 2))
        #print >> sys.stderr, testlimit

        z = mpmath.mpf('0')
        q = mpmath.mpf('0')
        value = jacobi_theta_2(z, q)
        self.assertEquals(value, mpmath.mpf('0'))

        # Mathematica value for v2(z = 0.1, q = 0.1) = 1.12981
        # q = 0.1, therefore m = 0.802403, according to Mathematica
        z = mpmath.mpf('0.1')
        m = mpmath.mpf('0.802403')

        result = jacobi_theta_2(z, m)     # verbosity on
        self.assertEquals(round(result, 5), 1.12981)

        z = (mpmath.pi)/(mpmath.mpf('2'))       # test for cos zeros
        for i in range(10):
            qstring = str(random.random())
            q = mpmath.mpf(qstring)

            result = jacobi_theta_2(z, q)
            if result < testlimit:
                self.assertEquals(True, True)
            else:
                print >> sys.stderr, result
                self.assertEquals(False, True)

        z = mpmath.mpf('0')
        q = mpmath.pi / mpmath.mpf('2')
        self.assertRaises(ValueError, jacobi_theta_1, z, q)
        q = mpmath.mpf('1')
        self.assertRaises(ValueError, jacobi_theta_1, z, q)

        z = mpmath.mpf('0.1')
        q = mpmath.mpf('0.1')
        value = jacobi_theta_2(z, q)
        self.assertTrue(isinstance(value, mpmath.mpf))

    def testJacobiTheta3(self):
        mpmath.mpf.dps = 100
        mpmath.mp.dps = 100
        testlimit = mpmath.mpf('10')**(-1*(mpmath.mpf.dps - 2))
        #print >> sys.stderr, testlimit

        one = mpmath.mpf('1')

        z = mpmath.mpf('0')
        q = mpmath.mpf('0')
        value = jacobi_theta_3(z, q)
        self.assertEquals(mpmath.mpf('1'), value)
        self.assertTrue(isinstance(value, mpmath.mpf))

        # Mathematica value for v3(z = 0.1, q = 0.1) = 1.1962
        # q = 0.1, therefore m = 0.802403, according to Mathematica
        z = mpmath.mpf('0.1')
        m = mpmath.mpf('0.802403')

        result = jacobi_theta_3(z, m)
        self.assertEquals(round(result, 4), 1.1962)

        mpmath.mpf.dps = 2
        z = mpmath.mpf('0')
        q = mpmath.pi / mpmath.mpf('2')
        self.assertRaises(ValueError, jacobi_theta_1, z, q)
        q = mpmath.mpf('1')
        self.assertRaises(ValueError, jacobi_theta_1, z, q)

        z = mpmath.mpf('0.1')
        q = mpmath.mpf('0.1')
        value = jacobi_theta_2(z, q)
        self.assertTrue(isinstance(value, mpmath.mpf))

    def testJacobiTheta4(self):
        mpmath.mpf.dps = 100
        mpmath.mp.dps = 100
        testlimit = mpmath.mpf('10')**(-1*(mpmath.mpf.dps - 2))
        #print >> sys.stderr, testlimit

        z = mpmath.mpf('0')
        q = mpmath.mpf('0')
        value = jacobi_theta_4(z, q)
        self.assertEquals(value, mpmath.mpf('1.0'))

        # Mathematica value for v4(z = 0.1, q = 0.1) = 0.804171
        # q = 0.1, therefore m = 0.802403, according to Mathematica
        z = mpmath.mpf('0.1')
        m = mpmath.mpf('0.802403')

        result = jacobi_theta_4(z, m)
        self.assertEquals(round(result, 6), 0.804171)

        z = mpmath.mpf('0')
        q = mpmath.pi / mpmath.mpf('2')
        self.assertRaises(ValueError, jacobi_theta_1, z, q)
        q = mpmath.mpf('1')
        self.assertRaises(ValueError, jacobi_theta_1, z, q)

        z = mpmath.mpf('0.1')
        q = mpmath.mpf('0.1')
        value = jacobi_theta_4(z, q)
        self.assertTrue(isinstance(value, mpmath.mpf))

    def testJacobiThetaIdentidies(self):
        """
        Tests the some of the jacobi identidies found in Abramowitz,
        Sec. 16.28, Pg. 576.  The identies are tested to 1 part in 10^98.

        Warning: running these tests takes a long time.
        """
        mpmath.mpf.dps = 110        # offset to satisfy target of 1 in 10^98
        mpmath.mp.dps = 110
        #mpmath.mpf.dps = 30        # testing version
        #mpmath.mp.dps = 30
        testlimit = mpmath.mpf('10')**(-1*(mpmath.mpf.dps - 10))
        #print >> sys.stderr, testlimit

        zero = mpmath.mpf('0')

        # Abramowitz 16.28.1
        #
        # v_1(z, q)**2 * v_4(0, q)**2 =   v_3(z, q)**2 * v_2(0, q)**2
        #                               - v_2(z, q)**2 * v_3(0, q)**2

        for i in range(10):
            qstring = str(random.random())
            q = mpmath.mpf(qstring)

            zstring = str(10*random.random())
            z = mpmath.mpf(zstring)

            term1 = (jacobi_theta_1(z, q)**2) * (jacobi_theta_4(zero, q)**2)
            term2 = (jacobi_theta_3(z, q)**2) * (jacobi_theta_2(zero, q)**2)
            term3 = (jacobi_theta_2(z, q)**2) * (jacobi_theta_3(zero, q)**2)

            equality = term1 - term2 + term3

            if equality < testlimit:
                self.assertEquals(True, True)
            else:
                print >> sys.stderr, 'Equality (~ zero): %e' % equality
                self.assertEquals(False, True)

        # Abramowitz 16.28.2
        #
        # v_2(z, q)**2 * v_4(0, q)**2 =   v_4(z, q)**2 * v_2(0, q)**2
        #                               - v_1(z, q)**2 * v_3(0, q)**2

        for i in range(10):
            qstring = str(random.random())
            q = mpmath.mpf(qstring)

            zstring = str(100*random.random())
            z = mpmath.mpf(zstring)

            term1 = (jacobi_theta_2(z, q)**2) * (jacobi_theta_4(zero, q)**2)
            term2 = (jacobi_theta_4(z, q)**2) * (jacobi_theta_2(zero, q)**2)
            term3 = (jacobi_theta_1(z, q)**2) * (jacobi_theta_3(zero, q)**2)

            equality = term1 - term2 + term3

            if equality < testlimit:
                self.assertEquals(True, True)
            else:
                print >> sys.stderr, 'Equality (~ zero): %e' % equality
                self.assertEquals(False, True)

        # Abramowitz 16.28.3
        #
        # v_3(z, q)**2 * v_4(0, q)**2 =   v_4(z, q)**2 * v_3(0, q)**2
        #                               - v_1(z, q)**2 * v_2(0, q)**2

        for i in range(10):
            qstring = str(random.random())
            q = mpmath.mpf(qstring)

            zstring = str(100*random.random())
            z = mpmath.mpf(zstring)

            term1 = (jacobi_theta_3(z, q)**2) * (jacobi_theta_4(zero, q)**2)
            term2 = (jacobi_theta_4(z, q)**2) * (jacobi_theta_3(zero, q)**2)
            term3 = (jacobi_theta_1(z, q)**2) * (jacobi_theta_2(zero, q)**2)

            equality = term1 - term2 + term3

            if equality < testlimit:
                self.assertEquals(True, True)
            else:
                print >> sys.stderr, 'Equality (~ zero): %e' % equality
                self.assertEquals(False, True)

        # Abramowitz 16.28.4
        #
        # v_4(z, q)**2 * v_4(0, q)**2 =   v_3(z, q)**2 * v_3(0, q)**2
        #                               - v_2(z, q)**2 * v_2(0, q)**2

        for i in range(10):
            qstring = str(random.random())
            q = mpmath.mpf(qstring)

            zstring = str(100*random.random())
            z = mpmath.mpf(zstring)

            term1 = (jacobi_theta_4(z, q)**2) * (jacobi_theta_4(zero, q)**2)
            term2 = (jacobi_theta_3(z, q)**2) * (jacobi_theta_3(zero, q)**2)
            term3 = (jacobi_theta_2(z, q)**2) * (jacobi_theta_2(zero, q)**2)

            equality = term1 - term2 + term3

            if equality < testlimit:
                self.assertEquals(True, True)
            else:
                print >> sys.stderr, 'Equality (~ zero): %e' % equality
                self.assertEquals(False, True)

        # Abramowitz 16.28.5
        #
        # v_2(0, q)**4 + v_4(0, q)**4 == v_3(0, q)**4

        for i in range(10):
            qstring = str(random.random())
            q = mpmath.mpf(qstring)

            term1 = (jacobi_theta_2(zero, q))**4
            term2 = (jacobi_theta_4(zero, q))**4
            term3 = (jacobi_theta_3(zero, q))**4

            equality = term1 + term2 - term3

            if equality < testlimit:
                self.assertEquals(True, True)
            else:
                print >> sys.stderr, 'Equality (~ zero): %e' % equality
                self.assertEquals(False, True)

    def testJacobiEllipticSn(self):
        """
        Test some special cases of the sn(z, q) function.

        This is an intensive test, so precision turned down during
        development.
        """
        mpmath.mpf.dps = 100
        mpmath.mp.dps = 100
        #mpmath.mpf.dps = 20             # testing version
        #mpmath.mp.dps = 20
        testlimit = mpmath.mpf('10')**(-1*(mpmath.mpf.dps - 4))
        #print >> sys.stderr, testlimit

        zero = mpmath.mpf('0')
        one = mpmath.mpf('1')

        # trival case
        result = jacobi_elliptic_sn(zero, zero)
        self.assertEquals(result, zero)

        # Abramowitz Table 16.5
        #
        # sn(0, m) = 0

        for i in range(10):
            qstring = str(random.random())
            q = mpmath.mpf(qstring)

            equality = jacobi_elliptic_sn(zero, q)

            if equality < testlimit:
                self.assertEquals(True, True)
            else:
                print >> sys.stderr, 'Sn, sn(0, m) == 0: %e' % equality
                self.assertEquals(False, True)

        # Abramowitz Table 16.5
        #
        # sn(K, m) = 1; K is K(k), first complete elliptic integral
        for i in range(10):
            mstring = str(random.random())
            m = mpmath.mpf(mstring)
            k = m.sqrt()

            K = mpmath.ellipk(k**2)

            equality = abs(one - jacobi_elliptic_sn(K, m))

            #print >> sys.stderr, '\nk:', k,
            #print >> sys.stderr, '\tK(k):', K,
            #print >> sys.stderr, '\tsn(K(k), m):', (equality + 1)

            if equality < testlimit:
                self.assertEquals(True, True)
            else:
                print >> sys.stderr, 'Sn, sn(K, m) - 1 == 0: %e' % equality
                equality = jacobi_elliptic_sn(K, m, True)   # verbose
                self.assertEquals(False, True)

        # Abramowitz Table 16.6.1
        #
        # sn(z, 0) = sin(z), m == 0
        #
        # sn(z, 1) = tanh(z), m == 1
        #
        # It would be nice to test these, but I find that they run
        # in to numerical trouble.  I'm currently treating as a boundary
        # case for sn function.

        # Mathematica value for sn(z = 0.1, m = 0.1) = 0.0998169
        arg = mpmath.mpf('0.1')
        result = jacobi_elliptic_sn(arg, arg)
        self.assertEquals(round(result, 7), 0.0998169)

    def testJacobiEllipticCn(self):
        """
        Test some special cases of the cn(z, q) function.

        This is an intensive test, so precision turned down during
        development.
        """
        mpmath.mpf.dps = 100
        mpmath.mp.dps = 100
        #mpmath.mpf.dps = 20             # testing version
        #mpmath.mp.dps = 20
        testlimit = mpmath.mpf('10')**(-1*(mpmath.mpf.dps - 4))
        #print >> sys.stderr, testlimit

        zero = mpmath.mpf('0')
        one = mpmath.mpf('1')

        # Abramowitz Table 16.5
        #
        # cn(0, q) = 1

        for i in range(10):
            qstring = str(random.random())
            q = mpmath.mpf(qstring)

            cn = jacobi_elliptic_cn(zero, q)
            equality = one - cn

            if equality < testlimit:
                self.assertEquals(True, True)
            else:
                print >> sys.stderr, 'Cn (~ 1): %e' % cn
                print >> sys.stderr, 'Equality (~ zero): %e' % equality
                self.assertEquals(False, True)

        # Abramowitz Table 16.5
        #
        # cn(K, q) = 0; K is K(k), first complete elliptic integral

        for i in range(10):
            mstring = str(random.random())
            m = mpmath.mpf(mstring)
            k = m.sqrt()

            K = mpmath.ellipk(k**2)

            equality = jacobi_elliptic_cn(K, m)

            if equality < testlimit:
                self.assertEquals(True, True)
            else:
                print >> sys.stderr, '\n**** Cn failure ****'
                print >> sys.stderr, '\nK: %e' % K,
                print >> sys.stderr, '\tm: %f' % m,
                print >> sys.stderr, '\tcn: %e' % equality
                equality = jacobi_elliptic_cn(K, k, True)
                self.assertEquals(False, True)

        # Abramowitz Table 16.6.2
        #
        # cn(u, 0) = cos(u), m == 0
        #
        # cn(u, 1) = sech(z), m == 1
        #
        # It would be nice to test these, but I find that they run
        # in to numerical trouble.  I'm currently treating as a boundary
        # case for cn function.

    def testJacobiEllipticDn(self):
        """
        Test some special cases of the dn(z, q) function.
        """
        mpmath.mpf.dps = 100
        mpmath.mp.dps = 100
        #mpmath.mpf.dps = 20             # testing version
        #mpmath.mp.dps = 20
        testlimit = mpmath.mpf('10')**(-1*(mpmath.mpf.dps - 4))
        #print >> sys.stderr, testlimit

        zero = mpmath.mpf('0')
        one = mpmath.mpf('1')

        # Abramowitz Table 16.5
        #
        # dn(0, q) = 1

        for i in range(10):
            mstring = str(random.random())
            m = mpmath.mpf(mstring)

            dn = jacobi_elliptic_dn(zero, m)
            equality = one - dn

            if equality < testlimit:
                self.assertEquals(True, True)
            else:
                print >> sys.stderr, '\n**** Dn failure ****'
                print >> sys.stderr, '\tm: %f' % m,
                print >> sys.stderr, '\tdn: %e' % dn,
                print >> sys.stderr, '\tEquality: %e' % equality
                equality = jacobi_elliptic_dn(zero, m, True)
                self.assertEquals(False, True)

        # Abramowitz Table 16.6.3
        #
        # dn(z, 0) = 1, m == 0

        for i in range(10):
            zstring = str(random.random())
            z = mpmath.mpf(zstring)

            value = jacobi_elliptic_dn(z, zero)

            equality = value - one

            if equality < testlimit:
                self.assertEquals(True, True)
            else:
                print >> sys.stderr, 'Equality (~ zero): %e' % equality
                self.assertEquals(False, True)

    def testJacobiEllipticIdentidies(self):
        """
        Tests the some of the jacobi elliptic function identidies found
        on Mathworld.  Havne't found in Abramowitz.  The identies are
        tested to 1 part in 10^98.
        """
        mpmath.mpf.dps = 100
        mpmath.mp.dps = 100
        #mpmath.mpf.dps = 20             # testing version
        #mpmath.mp.dps = 20
        testlimit = mpmath.mpf('10')**(-1*(mpmath.mpf.dps - 4))
        #print >> sys.stderr, testlimit

        zero = mpmath.mpf('0')
        one = mpmath.mpf('1')

        # MathWorld
        #
        # sn(z, q)**2 + cn(z, q)**2 == 1

        for i in range(10):
            qstring = str(random.random())
            q = mpmath.mpf(qstring)

            zstring = str(100*random.random())
            z = mpmath.mpf(zstring)

            term1 = jacobi_elliptic_sn(z, q)**2
            term2 = jacobi_elliptic_cn(z, q)**2

            equality = one - term1 - term2

            if equality < testlimit:
                self.assertEquals(True, True)
            else:
                print >> sys.stderr, 'Equality (~ zero): %e' % equality
                self.assertEquals(False, True)

        # MathWorld
        #
        # k**2 * sn(z, m)**2 + dn(z, m)**2 == 1

        for i in range(10):
            mstring = str(random.random())
            m = mpmath.mpf(qstring)
            k = m.sqrt()

            zstring = str(10*random.random())
            z = mpmath.mpf(zstring)

            term1 = k**2 * jacobi_elliptic_sn(z, m)**2
            term2 = jacobi_elliptic_dn(z, m)**2

            equality = one - term1 - term2

            if equality < testlimit:
                self.assertEquals(True, True)
            else:
                print >> sys.stderr, 'Equality (~ zero): %e' % equality
                self.assertEquals(False, True)

        # MathWorld
        #
        # k**2 * cn(z, m)**2 + (1 - k**2) = dn(z, m)**2

        for i in range(10):
            mstring = str(random.random())
            m = mpmath.mpf(qstring)
            k = m.sqrt()

            zstring = str(random.random())
            z = mpmath.mpf(zstring)

            term1 = k**2 * jacobi_elliptic_cn(z, m)**2
            term2 = 1 - k**2
            term3 = jacobi_elliptic_dn(z, m)**2

            equality = term3 - term1 - term2

            if equality < testlimit:
                self.assertEquals(True, True)
            else:
                print >> sys.stderr, 'Equality (~ zero): %e' % equality
                self.assertEquals(False, True)


if __name__ == '__main__':          # if run as script, run tests
    unittest.main()


