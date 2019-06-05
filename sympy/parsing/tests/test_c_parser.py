from sympy.parsing.c.c_parser import convert_c_file
import os

def test_variable():
    c_src1 = (
        'int a;' + '\n' +
        'int b;' + '\n'
    )
    c_src2 = (
        'float a;' + '\n' +
        'float b;' + '\n'
    )
    c_src3 = (
        'int a;' + '\n' +
        'float b;' + '\n' +
        'int c;'
    )
    file1 = open('..test1.h','w')
    file2 = open('..test2.h', 'w')
    file3 = open('..test3.h', 'w')

    file1.write(c_src1)
    file2.write(c_src2)
    file3.write(c_src3)

    file1.close()
    file2.close()
    file3.close()

    res1 = convert_c_file('..test1.h')
    res2 = convert_c_file('..test2.h')
    res3 = convert_c_file('..test3.h')

    cmp1 = ['a = 0', 'b = 0']
    cmp2 = ['a = 0.0', 'b = 0.0']
    cmp3 = ['a = 0', 'b = 0.0', 'c = 0']

    assert res1 == cmp1
    assert res2 == cmp2
    assert res3 == cmp3

    os.remove('..test1.h')
    os.remove('..test2.h')
    os.remove('..test3.h')

def test_int():
    c_src1 = 'int a = 1;'
    c_src2 = (
        'int a = 1;' + '\n' +
        'int b = 2;' + '\n'
    )
    file1 = open('..test1.h','w')
    file2 = open('..test2.h', 'w')

    file1.write(c_src1)
    file2.write(c_src2)

    file1.close()
    file2.close()

    res1 = convert_c_file('..test1.h')
    res2 = convert_c_file('..test2.h')

    cmp1 = ['a = 1']
    cmp2 = ['a = 1', 'b = 2']

    assert res1 == cmp1
    assert res2 == cmp2

    os.remove('..test1.h')
    os.remove('..test2.h')

def test_float():
    c_src1 = 'float a = 1.0;'
    c_src2 = (
        'float a = 1.25;' + '\n' +
        'float b = 2.39;' + '\n'
    )
    file1 = open('..test1.h','w')
    file2 = open('..test2.h', 'w')

    file1.write(c_src1)
    file2.write(c_src2)

    file1.close()
    file2.close()

    res1 = convert_c_file('..test1.h')
    res2 = convert_c_file('..test2.h')

    cmp1 = ['a = 1.0']
    cmp2 = ['a = 1.25', 'b = 2.39']

    assert res1 == cmp1
    assert res2 == cmp2

    os.remove('..test1.h')
    os.remove('..test2.h')


def test_function():
    c_src1 = (
        'void fun1()' + '\n' +
        '{' + '\n' +
        'int a;' + '\n' +
        '}'
    )
    c_src2 = (
        'int fun2()' + '\n' +
        '{'+ '\n' +
        'int a;' + '\n' +
        'return a;' + '\n' +
        '}'
    )
    c_src3 = (
        'float fun3()' + '\n' +
        '{' + '\n' +
        'float b;' + '\n' +
        'return b;' + '\n' +
        '}'
    )
    file1 = open('..test1.h','w')
    file2 = open('..test2.h', 'w')
    file3 = open('..test3.h', 'w')

    file1.write(c_src1)
    file2.write(c_src2)
    file3.write(c_src3)

    file1.close()
    file2.close()
    file3.close()

    res1 = convert_c_file('..test1.h')
    res2 = convert_c_file('..test2.h')
    res3 = convert_c_file('..test3.h')

    str1 = (
        'def fun1():' + '\n' +
        '    ' + 'a = 0'
    )

    str2 = (
        'def fun2():' + '\n' +
        '    ' + 'a = 0' + '\n' +
        '    ' + 'return a'
    )

    str3 = (
        'def fun3():' + '\n' +
        '    ' + 'b = 0.0' +'\n' +
        '    ' + 'return b'
    )
    cmp1 = [str1]
    cmp2 = [str2]
    cmp3 = [str3]

    assert res1 == cmp1
    assert res2 == cmp2
    assert res3 == cmp3

    os.remove('..test1.h')
    os.remove('..test2.h')
    os.remove('..test3.h')

def test_parameters():
    c_src1 = (
        'void fun1( int a)' + '\n' +
        '{' + '\n' +
        'int i;' + '\n' +
        '}'
    )
    c_src2 = (
        'int fun2(float x, float y)' + '\n' +
        '{'+ '\n' +
        'int a;' + '\n' +
        'return a;' + '\n' +
        '}'
    )
    c_src3 = (
        'float fun3(int p, float q, int r)' + '\n' +
        '{' + '\n' +
        'float b;' + '\n' +
        'return b;' + '\n' +
        '}'
    )
    file1 = open('..test1.h','w')
    file2 = open('..test2.h', 'w')
    file3 = open('..test3.h', 'w')

    file1.write(c_src1)
    file2.write(c_src2)
    file3.write(c_src3)

    file1.close()
    file2.close()
    file3.close()

    res1 = convert_c_file('..test1.h')
    res2 = convert_c_file('..test2.h')
    res3 = convert_c_file('..test3.h')

    str1 = (
        'def fun1(a):' + '\n' +
        '    ' + 'i = 0'
    )

    str2 = (
        'def fun2(x, y):' + '\n' +
        '    ' + 'a = 0' + '\n' +
        '    ' + 'return a'
    )

    str3 = (
        'def fun3(p, q, r):' + '\n' +
        '    ' + 'b = 0.0' +'\n' +
        '    ' + 'return b'
    )
    cmp1 = [str1]
    cmp2 = [str2]
    cmp3 = [str3]

    assert res1 == cmp1
    assert res2 == cmp2
    assert res3 == cmp3

    os.remove('..test1.h')
    os.remove('..test2.h')
    os.remove('..test3.h')

def test_function_call():
    c_src1 = 'x = fun1(2);'
    c_src2 = 'y = fun2(2, 3, 4);'
    c_src3 = (
        'int p, q, r;' + '\n' +
        'z = fun3(p, q, r);'
    )
    c_src4 = (
        'float x, y;' + '\n' +
        'int z;' + '\n' +
        'i = fun4(x, y, z)'
    )
    file1 = open('..test1.h','w')
    file2 = open('..test2.h', 'w')
    file3 = open('..test3.h', 'w')
    file4 = open('..test4.h','w')

    file1.write(c_src1)
    file2.write(c_src2)
    file3.write(c_src3)
    file4.write(c_src4)

    file1.close()
    file2.close()
    file3.close()
    file4.close()

    res1 = convert_c_file('..test1.h')
    res2 = convert_c_file('..test2.h')
    res3 = convert_c_file('..test3.h')
    res4 = convert_c_file('..test4.h')

    cmp1 = ['x = fun1(2)']
    cmp2 = ['y = fun2(2, 3, 4)']
    cmp3 = ['p = 0', 'q = 0', 'r = 0', 'z = fun3(p, q, r)']
    cmp4 = ['x = 0.0', 'y = 0.0', 'z = 0', 'i = fun4(x, y, z)']

    assert res1 == cmp1
    assert res2 == cmp2
    assert res3 == cmp3
    assert res4 == cmp4

    os.remove('..test1.h')
    os.remove('..test2.h')
    os.remove('..test3.h')
    os.remove('..test4.h')

test_variable()
test_int()
test_float()
test_function()
test_parameters()
test_function_call()
