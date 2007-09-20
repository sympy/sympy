from sympy.sandbox.core import *

def test_arithmetic():
    x = Symbol('x')
    y = Symbol('y')
    assert x+1 == 1+x
    assert x+y == y+x
    assert (x+y)+3 == x+(y+3)
    assert 2*x == x*2
    assert 0*x == 0
    assert 1*x == x
    assert 2*x - x == x

def test_integer_arithmetic():
    i = Integer(2)
    assert i is Integer(2)
    j = Integer(-5)
    assert +i==2
    assert -i==-2
    assert i+j==-3
    assert i-j==7
    assert i*j==-10
    assert i/j==Fraction(-2,5)
    assert j**i==25
    assert i**j==Fraction(1,32)

    i = Integer(2)
    j = -5
    assert i+j==-3
    assert i-j==7
    assert i*j==-10
    assert i/j==Fraction(-2,5)
    assert j**i==25
    assert i**j==Fraction(1,32)

    i = 2
    j = Integer(-5)
    assert i+j==-3
    assert i-j==7
    assert i*j==-10
    assert i/j==Fraction(-2,5)
    assert j**i==25
    assert i**j==Fraction(1,32)

def test_fraction_arithmetic():

    i = Fraction(2,3)
    assert i is Fraction(2,3)
    j = Fraction(-5,4)
    assert +i==Fraction(2,3)
    assert -i==Fraction(-2,3)
    assert i+j==Fraction(-7,12)
    assert i-j==Fraction(23,12)
    assert i*j==Fraction(-5,6)
    assert i/j==Fraction(-8,15)

    i = Fraction(2,3)
    j = Integer(5)
    assert i+j==Fraction(17,3)
    assert i-j==Fraction(-13,3)
    assert i*j==Fraction(10,3)
    assert i/j==Fraction(2,15)
    assert i**j==Fraction(32,243)

    j = Fraction(2,3)
    i = Integer(5)
    assert i+j==Fraction(17,3)
    assert i-j==Fraction(13,3)
    assert i*j==Fraction(10,3)
    assert i/j==Fraction(15,2)

    i = Fraction(2,3)
    j = 5
    assert i+j==Fraction(17,3)
    assert i-j==Fraction(-13,3)
    assert i*j==Fraction(10,3)
    assert i/j==Fraction(2,15)
    assert i**j==Fraction(32,243)

    j = Fraction(2,3)
    i = 5
    assert i+j==Fraction(17,3)
    assert i-j==Fraction(13,3)
    assert i*j==Fraction(10,3)
    assert i/j==Fraction(15,2)

def test_float_arithmetic():

    tofloat = lambda n: repr(float())[:12]
    
    i = Float(1.2)
    assert i is Float(1.2)
    j = Float(-3.4)
    assert +i==1.2
    assert float(-i)==-1.2
    assert tofloat(i+j)==tofloat(-2.2)
    assert float(i-j)==4.6
    assert float(i*j)==-4.08
    assert float(i/j)==1.2/(-3.4)
    assert tofloat(i**j)==tofloat(1.2**(-3.4))

    i = Float(1.2)
    j = -3.4
    assert tofloat(i+j)==tofloat(-2.2)
    assert tofloat(i-j)==tofloat(4.6)
    assert tofloat(i*j)==tofloat(-4.08)
    assert float(i/j)==1.2/(-3.4)
    assert tofloat(i**j)==tofloat(1.2**(-3.4))

    i = 1.2
    j = Float(-3.4)
    assert tofloat(i+j)==tofloat(-2.2)
    assert tofloat(i-j)==tofloat(4.6)
    assert tofloat(i*j)==tofloat(-4.08)
    assert float(i/j)==1.2/(-3.4)
    assert tofloat(i**j)==tofloat(1.2**(-3.4))

    i = Float(1.2)
    j = Integer(2)
    assert tofloat(i+j)==tofloat(3.2)
    assert tofloat(i-j)==tofloat(-0.8)
    assert tofloat(i*j)==tofloat(2.4)
    assert float(i/j)==1.2/2
    assert float(i**j)==1.2**2

    j = Float(1.2)
    i = Integer(2)
    assert tofloat(i+j)==tofloat(3.2)
    assert tofloat(i-j)==tofloat(0.8)
    assert tofloat(i*j)==tofloat(2.4)
    assert float(i/j)==2/1.2
    assert tofloat(i**j)==tofloat(2/1.2)

    i = Float(1.2)
    j = 2
    assert tofloat(i+j)==tofloat(3.2)
    assert tofloat(i-j)==tofloat(-0.8)
    assert tofloat(i*j)==tofloat(2.4)
    assert float(i/j)==1.2/2
    assert float(i**j)==1.2**2

    j = Float(1.2)
    i = 2
    assert tofloat(i+j)==tofloat(3.2)
    assert tofloat(i-j)==tofloat(0.8)
    assert tofloat(i*j)==tofloat(2.4)
    assert float(i/j)==2/1.2
    assert tofloat(i**j)==tofloat(2**1.2)

    i = Float(1.2)
    j = Fraction(-17,5)
    assert tofloat(i+j)==tofloat(-2.2)
    assert tofloat(i-j)==tofloat(4.6)
    assert tofloat(i*j)==tofloat(-4.08)
    assert tofloat(i/j)==tofloat(1.2/(-3.4))
    assert tofloat(i**j)==tofloat(1.2**(-3.4))

    i = Fraction(6,5)
    j = Float(-3.4)
    assert tofloat(i+j)==tofloat(-2.2)
    assert tofloat(i-j)==tofloat(4.6)
    assert tofloat(i*j)==tofloat(-4.08)
    assert float(i/j)==(6.0/5)/(-3.4)
    assert tofloat(i**j)==tofloat(1.2**(-3.4))

def timedtest():
    from time import clock
    i = 20
    t1 = clock()
    while i:
        i -= 1
        test_integer_arithmetic()
        test_fraction_arithmetic()
        test_float_arithmetic()
    t2 = clock()
    return t2-t1

if __name__=='__main__':
    print timedtest()
