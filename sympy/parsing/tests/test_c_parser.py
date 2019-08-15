from sympy.parsing.sym_expr import SymPyExpression
from sympy.utilities.pytest import raises
from sympy.external import import_module

cin = import_module('clang.cindex', __import__kwargs = {'fromlist': ['cindex']})

if cin:
    from sympy.codegen.ast import (Variable, IntBaseType, FloatBaseType, String,
                                   Return, FunctionDefinition, Integer, Float,
                                   Declaration, CodeBlock, FunctionPrototype,
                                   FunctionCall, NoneToken)
    from sympy import Symbol
    import os

    def test_variable():
        c_src1 = (
            'int a;' + '\n' +
            'int b;' + '\n'
        )
        c_src2 = (
            'float a;' + '\n'
            + 'float b;' + '\n'
        )
        c_src3 = (
            'int a;' + '\n' +
            'float b;' + '\n' +
            'int c;'
        )

        res1 = SymPyExpression(c_src1, 'c').return_expr()
        res2 = SymPyExpression(c_src2, 'c').return_expr()
        res3 = SymPyExpression(c_src3, 'c').return_expr()

        assert res1[0] == Declaration(
            Variable(
                Symbol('a'),
                type=IntBaseType(String('integer')),
                value=Integer(0)
            )
        )
        assert res1[1] == Declaration(
            Variable(
                Symbol('b'),
                type=IntBaseType(String('integer')),
                value=Integer(0)
            )
        )

        assert res2[0] == Declaration(
            Variable(
                Symbol('a'),
                type=FloatBaseType(String('real')),
                value=Float('0.0', precision=53)
            )
        )
        assert res2[1] == Declaration(
            Variable(
                Symbol('b'),
                type=FloatBaseType(String('real')),
                value=Float('0.0', precision=53)
            )
        )

        assert res3[0] == Declaration(
            Variable(
                Symbol('a'),
                type=IntBaseType(String('integer')),
                value=Integer(0)
            )
        )

        assert res3[1] == Declaration(
            Variable(
                Symbol('b'),
                type=FloatBaseType(String('real')),
                value=Float('0.0', precision=53)
            )
        )

        assert res3[2] == Declaration(
            Variable(
                Symbol('c'),
                type=IntBaseType(String('integer')),
                value=Integer(0)
            )
        )


    def test_int():
        c_src1 = 'int a = 1;'
        c_src2 = (
            'int a = 1;' + '\n' +
            'int b = 2;' + '\n'
        )

        res1 = SymPyExpression(c_src1, 'c').return_expr()
        res2 = SymPyExpression(c_src2, 'c').return_expr()

        assert res1[0] == Declaration(
            Variable(
                Symbol('a'),
                type=IntBaseType(String('integer')),
                value=Integer(1)
            )
        )

        assert res2[0] == Declaration(
            Variable(
                Symbol('a'),
                type=IntBaseType(String('integer')),
                value=Integer(1)
            )
        )

        assert res2[1] == Declaration(
            Variable(
                Symbol('b'),
                type=IntBaseType(String('integer')),
                value=Integer(2)
            )
        )


    def test_float():
        c_src1 = 'float a = 1.0;'
        c_src2 = (
            'float a = 1.25;' + '\n' +
            'float b = 2.39;' + '\n'
        )

        res1 = SymPyExpression(c_src1, 'c').return_expr()
        res2 = SymPyExpression(c_src2, 'c').return_expr()

        assert res1[0] == Declaration(
            Variable(
                Symbol('a'),
                type=FloatBaseType(String('real')),
                value=Float('1.0', precision=53)
                )
            )

        assert res2[0] == Declaration(
            Variable(
                Symbol('a'),
                type=FloatBaseType(String('real')),
                value=Float('1.25', precision=53)
            )
        )

        assert res2[1] == Declaration(
            Variable(
                Symbol('b'),
                type=FloatBaseType(String('real')),
                value=Float('2.3900000000000001', precision=53)
            )
        )


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
        c_src4 = (
            'float fun4()' + '\n' +
            '{}'
        )

        res1 = SymPyExpression(c_src1, 'c').return_expr()
        res2 = SymPyExpression(c_src2, 'c').return_expr()
        res3 = SymPyExpression(c_src3, 'c').return_expr()
        res4 = SymPyExpression(c_src4, 'c').return_expr()

        assert res1[0] == FunctionDefinition(
            NoneToken(),
            name=String('fun1'),
            parameters=(),
            body=CodeBlock(
                Declaration(
                    Variable(
                        Symbol('a'),
                        type=IntBaseType(String('integer')),
                        value=Integer(0)
                    )
                )
            )
        )

        assert res2[0] == FunctionDefinition(
            IntBaseType(String('integer')),
            name=String('fun2'),
            parameters=(),
            body=CodeBlock(
                Declaration(
                    Variable(
                        Symbol('a'),
                        type=IntBaseType(String('integer')),
                        value=Integer(0)
                    )
                ),
                Return('a')
            )
        )

        assert res3[0] == FunctionDefinition(
            FloatBaseType(String('real')),
            name=String('fun3'),
            parameters=(),
            body=CodeBlock(
                Declaration(
                    Variable(
                        Symbol('b'),
                        type=FloatBaseType(String('real')),
                        value=Float('0.0', precision=53)
                    )
                ),
                Return('b')
            )
        )

        assert res4[0] == FunctionPrototype(
            FloatBaseType(String('real')),
            name=String('fun4'),
            parameters=()
        )


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

        res1 = SymPyExpression(c_src1, 'c').return_expr()
        res2 = SymPyExpression(c_src2, 'c').return_expr()
        res3 = SymPyExpression(c_src3, 'c').return_expr()

        assert res1[0] == FunctionDefinition(
            NoneToken(),
            name=String('fun1'),
            parameters=(
                Variable(
                    Symbol('a'),
                    type=IntBaseType(String('integer')),
                    value=Integer(0)
                ),
            ),
            body=CodeBlock(
                Declaration(
                    Variable(
                        Symbol('i'),
                        type=IntBaseType(String('integer')),
                        value=Integer(0)
                    )
                )
            )
        )

        assert res2[0] == FunctionDefinition(
            IntBaseType(String('integer')),
            name=String('fun2'),
            parameters=(
                Variable(
                    Symbol('x'),
                    type=FloatBaseType(String('real')),
                    value=Float('0.0', precision=53)
                ),
                Variable(
                    Symbol('y'),
                    type=FloatBaseType(String('real')),
                    value=Float('0.0', precision=53)
                )
            ),
            body=CodeBlock(
                Declaration(
                    Variable(
                        Symbol('a'),
                        type=IntBaseType(String('integer')),
                        value=Integer(0)
                    )
                ),
                Return('a')
            )
        )

        assert res3[0] == FunctionDefinition(
            FloatBaseType(String('real')), name=String('fun3'),
            parameters=(
                Variable(
                    Symbol('p'),
                    type=IntBaseType(String('integer')),
                    value=Integer(0)
                ),
                Variable(
                    Symbol('q'),
                    type=FloatBaseType(String('real')),
                    value=Float('0.0', precision=53)
                ),
                Variable(
                    Symbol('r'),
                    type=IntBaseType(String('integer')),
                    value=Integer(0)
                )
            ),
            body=CodeBlock(
                Declaration(
                    Variable(
                        Symbol('b'),
                        type=FloatBaseType(String('real')),
                        value=Float('0.0', precision=53)
                    )
                ),
                Return('b')
            )
        )


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
        c_src5 = 'a = fun()'

        res1 = SymPyExpression(c_src1, 'c').return_expr()
        res2 = SymPyExpression(c_src2, 'c').return_expr()
        res3 = SymPyExpression(c_src3, 'c').return_expr()
        res4 = SymPyExpression(c_src4, 'c').return_expr()
        res5 = SymPyExpression(c_src5, 'c').return_expr()

        assert res1[0] == Declaration(
            Variable(
                Symbol('x'),
                value=FunctionCall(
                    String('fun1'),
                    function_args=([2, ])
                )
            )
        )

        assert res2[0] == Declaration(
            Variable(
                Symbol('y'),
                value=FunctionCall(
                    String('fun2'),
                    function_args=([2, 3, 4])
                )
            )
        )

        assert res3[0] == Declaration(
            Variable(
                Symbol('p'),
                type=IntBaseType(String('integer')),
                value=Integer(0)
            )
        )

        assert res3[1] == Declaration(
            Variable(
                Symbol('q'),
                type=IntBaseType(String('integer')),
                value=Integer(0)
            )
        )

        assert res3[2] == Declaration(
            Variable(
                Symbol('r'),
                type=IntBaseType(String('integer')),
                value=Integer(0)
            )
        )

        assert res3[3] == Declaration(
            Variable(
                Symbol('z'),
                value=FunctionCall(
                    String('fun3'),
                    function_args=([Symbol('p'), Symbol('q'), Symbol('r')])
                )
            )
        )

        assert res4[0] == Declaration(
            Variable(
                Symbol('x'),
                type=FloatBaseType(String('real')),
                value=Float('0.0', precision=53)
            )
        )

        assert res4[1] == Declaration(
            Variable(
                Symbol('y'),
                type=FloatBaseType(String('real')),
                value=Float('0.0', precision=53)
            )
        )

        assert res4[2] == Declaration(
            Variable(
                Symbol('z'),
                type=IntBaseType(String('integer')),
                value=Integer(0)
            )
        )

        assert res4[3] == Declaration(
            Variable(
                Symbol('i'),
                value=FunctionCall(
                    String('fun4'),
                    function_args=([Symbol('x'), Symbol('y'), Symbol('z')])
                )
            )
        )
        assert res5[0] == Declaration(
            Variable(
                Symbol('a'),
                value=FunctionCall(String('fun'), function_args=())
            )
        )


    def test_parse():
        c_src1 = (
            'int a;' + '\n' +
            'int b;' + '\n'
        )
        c_src2 = (
            'void fun1()' + '\n' +
            '{' + '\n' +
            'int a;' + '\n' +
            '}'
        )

        f1 = open('..a.h', 'w')
        f2 = open('..b.h', 'w')

        f1.write(c_src1)
        f2. write(c_src2)

        f1.close()
        f2.close()

        res1 = SymPyExpression('..a.h', 'c').return_expr()
        res2 = SymPyExpression('..b.h', 'c').return_expr()

        os.remove('..a.h')
        os.remove('..b.h')

        assert res1[0] == Declaration(
            Variable(
                Symbol('a'),
                type=IntBaseType(String('integer')),
                value=Integer(0)
            )
        )
        assert res1[1] == Declaration(
            Variable(
                Symbol('b'),
                type=IntBaseType(String('integer')),
                value=Integer(0)
            )
        )
        assert res2[0] == FunctionDefinition(
            NoneToken(),
            name=String('fun1'),
            parameters=(),
            body=CodeBlock(
                Declaration(
                    Variable(
                        Symbol('a'),
                        type=IntBaseType(String('integer')),
                        value=Integer(0)
                    )
                )
            )
        )

else:
    def test_raise():
        from sympy.parsing.c.c_parser import CCodeConverter
        raises(ImportError, lambda: CCodeConverter())
        raises(ImportError, lambda: SymPyExpression(' ', mode = 'c'))
