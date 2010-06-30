from sympy.physics.hilbert import l2, L2, FockSpace, TensorProductHilbertSpace, DirectSumHilbertSpace, DirectPowerHilbertSpace

from sympy import Interval, oo, Symbol

def test_l2_int():
    s1 = l2(2)
    s2 = l2(2)
    assert isinstance(s1, l2)
    assert s1.dimension == 2
    assert s1 == s2
    assert s1.subs(2, 42) == l2(42)

def test_l2_oo():
    s2 = l2(oo)
    s3 = l2(oo)
    assert isinstance(s2, l2)
    assert s2.dimension == oo
    assert s2 == s3
    assert s2.subs(oo, 42) == l2(42)

def test_l2_symbol():
    x = Symbol('x')
    s3 = l2(x)
    s4 = l2(x)
    assert isinstance(s3, l2)
    assert s3.dimension == x
    assert s3 == s4
    assert s3.subs(x, 42) == l2(42)

def test_L2_int():
    b1 = L2(Interval(-42, 42))
    b2 = L2(Interval(-42, 42))
    assert isinstance(b1, L2)
    assert b1.dimension == oo
    assert b1.interval == Interval(-42, 42)
    assert b1 == b2
    assert b1.subs(-42, 2) == L2(Interval(2, 42))

def test_L2_oo():
    b2 = L2(Interval(-oo, oo))
    b3 = L2(Interval(-oo, oo))
    assert isinstance(b2, L2)
    assert b2.dimension == oo
    assert b2.interval == Interval(-oo, oo)
    assert b2 == b3
    assert b2.subs(-oo, -42) == L2(Interval(-42, oo, True))

def test_L2_symbol():
    x = Symbol('x', real = True)
    y = Symbol('y', real = True)
    b3 = L2(Interval(x, y))
    b4 = L2(Interval(x, y))
    assert isinstance(b3, L2)
    assert b3.dimension == oo
    assert b3.interval == Interval(x, y)
    assert b3 == b4
    assert b3.subs(x, -y) == L2(Interval(-y, y))

def test_FockSpace():
    f1 = FockSpace()
    f2 = FockSpace()
    assert isinstance(f1, FockSpace)
    assert f1.dimension == oo
    assert f1 == f2

def test_TensorProductHilbertSpace_l2_int():
    s1 = l2(2)
    s2 = l2(42)
    s3 = s1*s2
    assert isinstance(s3, TensorProductHilbertSpace)
    assert s3.dimension == 84   #(2*42)
    assert s3.subs(42, 3) == l2(2)*l2(3)
    assert s3.spaces == (l2(2), l2(42))

def test_TensorProductHilbertSpace_l2_oo():
    s1 = l2(oo)
    s2 = l2(42)
    s3 = s1*s2
    assert isinstance(s3, TensorProductHilbertSpace)
    assert s3.dimension == oo
    assert s3.subs(oo, 2) == l2(2)*l2(42)
    assert s3.spaces == (l2(oo), l2(42))

def test_TensorProductHilbertSpace_l2_symbol():
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')
    s1 = l2(x)
    s2 = l2(y)
    s3 = s1*s2
    assert isinstance(s3, TensorProductHilbertSpace)
    assert s3.dimension == x*y
    assert s3.subs(y, z) == l2(x)*l2(z)
    assert s3.spaces == (l2(x), l2(y))

def test_TensorProductHilbertSpace_L2_int():
    b1 = L2(Interval(-42, 42))
    b2 = L2(Interval(-21, 21))
    b3 = b1*b2
    assert isinstance(b3, TensorProductHilbertSpace)
    assert b3.dimension == oo
    assert b3.subs(-21, -42) == L2(Interval(-42, 42))*L2(Interval(-42, 21))
    assert b3.spaces == (L2(Interval(-42, 42)), L2(Interval(-21, 21)))

def test_TensorProductHilbertSpace_L2_oo():
    b1 = L2(Interval(-42, oo))
    b2 = L2(Interval(-oo, 42))
    b3 = b1*b2
    assert isinstance(b3, TensorProductHilbertSpace)
    assert b3.dimension == oo
    assert b3.subs(oo, 42) == L2(Interval(-42, 42, False, True))*L2(Interval(-oo, 42))
    assert b3.spaces == (L2(Interval(-42, oo)), L2(Interval(-oo, 42)))

def test_TensorProductHilbertSpace_L2_symbol():
    x = Symbol('x', real = True)
    y = Symbol('y', real = True)
    q = Symbol('q', real = True)
    p = Symbol('p', real = True)
    b1 = L2(Interval(x, y))
    b2 = L2(Interval(q, p))
    b3 = b1*b2
    assert isinstance(b3, TensorProductHilbertSpace)
    assert b3.dimension == oo
    assert b3.subs(q, x) == L2(Interval(x, y))*L2(Interval(x, p))
    assert b3.spaces == (L2(Interval(x, y)), L2(Interval(q, p)))

def test_TensorProductHilbertSpace_FockSpace():
    f1 = FockSpace()
    f2 = FockSpace()
    f3 = f1*f2
    # will never be a tensor product instance due to power combining (Fock spaces can't be different as of now)
    assert isinstance(f3, DirectPowerHilbertSpace)
    assert f3.dimension == oo

def test_TensorProductHilbertSpace_mixed():
    s1 = l2(2)
    s2 = l2(42)
    s3 = s1*s2
    s4 = l2(oo)
    s5 = l2(oo)
    s6 = s4*s5
    x = Symbol('x', real = True)
    y = Symbol('y', real = True)
    s7 = l2(x)
    s8 = l2(x)
    s9 = s7*s8
    b1 = L2(Interval(-42, 42))
    b2 = L2(Interval(-21, 21))
    b3 = b1*b2
    b4 = L2(Interval(-oo, oo))
    b5 = L2(Interval(-oo, oo))
    b6 = b4*b5
    q = Symbol('q', real = True)
    p = Symbol('p', real = True)
    b7 = L2(Interval(x, y))
    b8 = L2(Interval(q, p))
    b9 = b7*b8
    f1 = FockSpace()
    f2 = FockSpace()
    f3 = f1*f2
    true_test = s3*s6*s9*b3*b6*b9*f3
    assert isinstance(true_test, TensorProductHilbertSpace)
    assert true_test.dimension == oo
    assert true_test.spaces == (l2(2), l2(42), l2(oo)**2, l2(x)**2, L2(Interval(-42, 42)), L2(Interval(-21, 21)), L2(Interval(-oo, oo))**2, L2(Interval(x, y)), L2(Interval(q, p)), FockSpace()**2)

def test_DirectSumHilbertSpace_l2_int():
    s1 = l2(2)
    s2 = l2(42)
    s3 = s1+s2
    assert isinstance(s3, DirectSumHilbertSpace)
    assert s3.dimension == 44   #(2+42)
    assert s3.subs(2, 21) == l2(21)+l2(42)
    assert s3.spaces == (l2(2), l2(42))

def test_DirectSumHilbertSpace_l2_oo():
    s1 = l2(oo)
    s2 = l2(oo)
    s3 = s1+s2
    assert isinstance(s3, DirectSumHilbertSpace)
    assert s3.dimension == oo
    assert s3.subs(oo, 42) == l2(42)+l2(42)
    assert s3.spaces == (l2(oo), l2(oo))

def test_DirectSumHilbertSpace_l2_symbol():
    x = Symbol('x')
    y = Symbol('y')
    s1 = l2(x)
    s2 = l2(y)
    s3 = s1+s2
    assert isinstance(s3, DirectSumHilbertSpace)
    assert s3.dimension == x+y
    assert s3.subs(y, x) == l2(x)+l2(x)
    assert s3.spaces == (l2(x), l2(y))

def test_DirectSumHilbertSpace_L2_int():
    b1 = L2(Interval(-42, 42))
    b2 = L2(Interval(-21, 21))
    b3 = b1+b2
    assert isinstance(b3, DirectSumHilbertSpace)
    assert b3.dimension == oo
    assert b3.subs(-21, -42) == L2(Interval(-42, 42))+L2(Interval(-42, 21))
    assert b3.spaces == (L2(Interval(-42, 42)), L2(Interval(-21, 21)))

def test_DirectSumHilbertSpace_L2_oo():
    b1 = L2(Interval(-oo, oo))
    b2 = L2(Interval(-oo, oo))
    b3 = b1+b2
    assert isinstance(b3, DirectSumHilbertSpace)
    assert b3.dimension == oo
    assert b3.subs(oo, 42) == L2(Interval(-oo, 42, False, True))+L2(Interval(-oo, 42, False, True))
    assert b3.spaces == (L2(Interval(-oo, oo)), L2(Interval(-oo, oo)))

def test_DirectSumHilbertSpace_L2_symbol():
    x = Symbol('x', real = True)
    y = Symbol('y', real = True)
    q = Symbol('q', real = True)
    p = Symbol('p', real = True)
    b1 = L2(Interval(x, y))
    b2 = L2(Interval(q, p))
    b3 = b1+b2
    assert isinstance(b3, DirectSumHilbertSpace)
    assert b3.dimension == oo
    assert b3.subs(x, q) == L2(Interval(q, y))+L2(Interval(q, p))
    assert b3.spaces == (L2(Interval(x, y)), L2(Interval(q, p)))

def test_DirectSumHilbertSpace_FockSpace():
    f1 = FockSpace()
    f2 = FockSpace()
    f3 = f1+f2
    assert isinstance(f3, DirectSumHilbertSpace)
    assert f3.dimension == oo

def test_DirectSumHilbertSpace_mixed():
    s1 = l2(2)
    s2 = l2(42)
    s3 = s1+s2
    s4 = l2(oo)
    s5 = l2(oo)
    s6 = s4+s5
    x = Symbol('x', real = True)
    y = Symbol('y', real = True)
    s7 = l2(x)
    s8 = l2(x)
    s9 = s7+s8
    b1 = L2(Interval(-42, 42))
    b2 = L2(Interval(-21, 21))
    b3 = b1+b2
    b4 = L2(Interval(-oo, oo))
    b5 = L2(Interval(-oo, oo))
    b6 = b4+b5
    q = Symbol('q', real = True)
    p = Symbol('p', real = True)
    b7 = L2(Interval(x, y))
    b8 = L2(Interval(q, p))
    b9 = b7+b8
    f1 = FockSpace()
    f2 = FockSpace()
    f3 = f1+f2
    true_test = s3+s6+s9+b3+b6+b9+f3
    assert isinstance(true_test, DirectSumHilbertSpace)
    assert true_test.dimension == oo
    assert true_test.subs(-oo, -42) == s3+s6+s9+b3+L2(Interval(-42, oo, True))+L2(Interval(-42, oo, True))+b9+f3

def test_DirectPowerHilbertSpace_l2_int():
    x = Symbol('x')
    s1 = l2(5)
    s2 = s1*s1
    s3 = s1*s1*s1
    s4 = s1**2
    s5 = s1**3
    s6 = s3*s1**x
    assert isinstance(s2, DirectPowerHilbertSpace)
    assert isinstance(s3, DirectPowerHilbertSpace)
    assert isinstance(s4, DirectPowerHilbertSpace)
    assert isinstance(s5, DirectPowerHilbertSpace)
    assert isinstance(s6, DirectPowerHilbertSpace)
    assert s2 == s4
    assert s3 == s5
    assert s4.dimension == s2.dimension == 25  #(5**2)
    assert s5.dimension == s3.dimension == 125 #(5**3)
    assert s6.dimension == 5**(3+x)
    assert s4.base == s2.base == l2(5)
    assert s4.exp == s4.exp == 2
    assert s5.base == s3.base == l2(5)
    assert s5.exp == s3.exp == 3
    assert s6.base == l2(5)
    assert s6.exp == 3+x
    assert s5.subs(3, 4) == s1**4

def test_DirectPowerHilbertSpace_l2_oo():
    x = Symbol('x')
    s1 = l2(oo)
    s2 = s1*s1
    s3 = s1*s1*s1
    s4 = s1**2
    s5 = s1**3
    s6 = s3*s1**x
    assert isinstance(s2, DirectPowerHilbertSpace)
    assert isinstance(s3, DirectPowerHilbertSpace)
    assert isinstance(s4, DirectPowerHilbertSpace)
    assert isinstance(s5, DirectPowerHilbertSpace)
    assert isinstance(s6, DirectPowerHilbertSpace)
    assert s2 == s4
    assert s3 == s5
    assert s4.dimension == s2.dimension == oo
    assert s5.dimension == s3.dimension == oo
    assert s6.dimension == oo
    assert s4.base == s2.base == l2(oo)
    assert s4.exp == s4.exp == 2
    assert s5.base == s3.base == l2(oo)
    assert s5.exp == s3.exp == 3
    assert s6.base == l2(oo)
    assert s6.exp == 3+x
    assert s5.subs(3, 4) == s1**4

def test_DirectPowerHilbertSpace_l2_symbol():
    x = Symbol('x')
    y = Symbol('y')
    s1 = l2(y)
    s2 = s1*s1
    s3 = s1*s1*s1
    s4 = s1**2
    s5 = s1**3
    s6 = s3*s1**x
    assert isinstance(s2, DirectPowerHilbertSpace)
    assert isinstance(s3, DirectPowerHilbertSpace)
    assert isinstance(s4, DirectPowerHilbertSpace)
    assert isinstance(s5, DirectPowerHilbertSpace)
    assert isinstance(s6, DirectPowerHilbertSpace)
    assert s2 == s4
    assert s3 == s5
    assert s4.dimension == s2.dimension == y**2
    assert s5.dimension == s3.dimension == y**3
    assert s6.dimension == y**(3+x)
    assert s4.base == s2.base == l2(y)
    assert s4.exp == s4.exp == 2
    assert s5.base == s3.base == l2(y)
    assert s5.exp == s3.exp == 3
    assert s6.base == l2(y)
    assert s6.exp == 3+x
    assert s6.subs(x, y) == s1**(3+y)

def test_DirectPowerHilbertSpace_L2_int():
    b1 = L2(Interval(-42, 42))
    b3 = b1*b1
    b4 = b3*b1**2*b1
    assert isinstance(b3, DirectPowerHilbertSpace)
    assert isinstance(b4, DirectPowerHilbertSpace)
    assert b3 == b1**2
    assert b4 == b1**5
    assert b3.dimension == b4.dimension == oo
    assert b3.base == b4.base == L2(Interval(-42, 42))
    assert b4.exp == 5
    assert b3.exp == 2

def test_DirectPowerHilbertSpace_L2_oo():
    b1 = L2(Interval(-oo, oo))
    b3 = b1*b1
    b4 = b3*b1**2*b1
    assert isinstance(b3, DirectPowerHilbertSpace)
    assert isinstance(b4, DirectPowerHilbertSpace)
    assert b3 == b1**2
    assert b4 == b1**5
    assert b3.dimension == b4.dimension == oo
    assert b3.base == b4.base == L2(Interval(-oo, oo))
    assert b4.exp == 5
    assert b3.exp == 2

def test_DirectPowerHilbertSpace_L2_symbol():
    x = Symbol('x', real = True)
    y = Symbol('y', real = True)
    b1 = L2(Interval(x, y))
    b3 = b1*b1
    b4 = b3*b1**2*b1
    b5 = b4*b1**x
    assert isinstance(b3, DirectPowerHilbertSpace)
    assert isinstance(b4, DirectPowerHilbertSpace)
    assert isinstance(b5, DirectPowerHilbertSpace)
    assert b3 == b1**2
    assert b4 == b1**5
    assert b5 == b1**(5+x)
    assert b3.dimension == b4.dimension == b5.dimension == oo
    assert b3.base == b4.base == b5.base == L2(Interval(x, y))
    assert b4.exp == 5
    assert b3.exp == 2
    assert b5.exp == (5+x)

def test_DirectPowerHilbertSpace_FockSpace():
    x = Symbol('x')
    y = Symbol('y')
    f1 = FockSpace()
    f2 = f1*f1*f1*f1**4*f1*f1**x*f1**2*f1**y
    assert isinstance(f2, DirectPowerHilbertSpace)
    assert f2.dimension == oo
    assert f2 == f1**(10+x+y)
    assert f2.base == FockSpace()
    assert f2.exp == (10+x+y)

def test_DirectPowerHilbertSpace_mixed():
    x = Symbol('x', real = True)
    y = Symbol('y', real = True)
    s1 = l2(5)
    s2 = s1*s1
    s3 = s1*s1*s1
    s6 = s3*s1**x
    s12 = l2(oo)
    s22 = s12*s12
    s32 = s12*s12*s12
    s62 = s32*s12**x
    s13 = l2(y)
    s23 = s13*s13
    s33 = s13*s13*s13
    s63 = s33*s13**x
    b1 = L2(Interval(-42, 42))
    b3 = b1*b1
    b4 = b3*b1**2*b1
    b12 = L2(Interval(-oo, oo))
    b32 = b12*b12
    b42 = b32*b12**2*b12
    b13 = L2(Interval(x, y))
    b33 = b13*b13
    b43 = b33*b13**2*b13
    b53 = b43*b13**x
    f1 = FockSpace()
    f2 = f1*f1*f1*f1**4*f1*f1**x*f1**2*f1**y
    true_test = s2*s3*s6*s22*s32*s62*s23*s33*s63*b3*b4*b32*b42*b33*b43*b53*f2
    assert isinstance(true_test, TensorProductHilbertSpace)
    assert true_test.dimension == oo
    assert true_test == l2(5)**(8+x)*l2(oo)**(8+x)*l2(y)**(8+x)*L2(Interval(-42,42))**7*L2(Interval(-oo,oo))**7*L2(Interval(x, y))**(12+x)*FockSpace()**(10+x+y)
    assert true_test.spaces == (l2(5)**(8+x), l2(oo)**(8+x), l2(y)**(8+x), L2(Interval(-42,42))**7, L2(Interval(-oo,oo))**7, L2(Interval(x, y))**(12+x), FockSpace()**(10+x+y))
