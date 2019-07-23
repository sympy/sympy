from sympy.parsing.sym_expr import SymPyExpression
import re

expr1 = SymPyExpression()
expr2 = SymPyExpression()
src = """\
integer :: a, b, c, d
real :: p, q, r, s
"""

def test_sym_expr():
    src1 = (
        src +
        """\
        d = a + b -c
        """
    )
    expr3 = SymPyExpression(src,'f')
    expr4 = SymPyExpression(src1,'f')
    for iter in expr3.return_expr():
        assert re.match(r'Declaration\(.*?\)', str(iter))
    ls = expr4.return_expr()
    for i in range(0,7):
        assert re.match(r'Declaration\(.*?\)', str(ls[i]))
    assert re.match(r'Assignment\(.*?\)', str(ls[8]))

def test_assignment():
    src1 = (
        src +
        """\
        a = b
        c = d
        p = q
        r = s
        """
    )
    expr1.convert_to_expr(src1, 'f')
    ls1 = expr1.return_expr()
    for iter in range(0, 12):
        if (iter < 8):
            pass
        else:
            assert re.match(
                r'Assignment\((?P<var>Variable)\(\w+\),\s(?P=var)\(\w+\)\)',
                str(ls1[iter])
            )

def test_binop_add():
    src1 = (
        src +
        """\
        c = a + b
        d = a + c
        s = p + q + r
        """
    )
    expr1.convert_to_expr(src1, 'f')
    ls1 = expr1.return_expr()
    for iter in range(0, 11):
        if(iter<8):
            pass
        else:
            assert re.match(
                r'Assignment\(Variable\(\w+\),\s\w+(\s\+\s\w+)+\)',
                str(ls1[iter])
            )

def test_binop_sub():
    src1 = (
        src +
        """\
        c = a - b
        d = a - c
        s = p - q - r
        """
    )
    expr1.convert_to_expr(src1, 'f')
    ls1 = expr1.return_expr()
    for iter in range(0, 11):
        if(iter<8):
            pass
        else:
            assert re.match(
                r'Assignment\(Variable\(\w+\),\s\w+(\s-\s\w+)+\)',
                str(ls1[iter])
            )

def test_binop_mul():
    src1 = (
        src +
        """\
        c = a * b
        d = a * c
        s = p * q * r
        """
    )
    expr1.convert_to_expr(src1, 'f')
    ls1 = expr1.return_expr()
    for iter in range(0, 11):
        if(iter<8):
            pass
        else:
            assert re.match(
                r'Assignment\(Variable\(\w+\),\s\w+(\*\w+)+\)',
                str(ls1[iter])
            )

def test_binop_div():
    src1 = (
        src +
        """\
        c = a / b
        d = a / c
        s = p / q
        r = q / p
        """
    )
    expr1.convert_to_expr(src1, 'f')
    ls1 = expr1.return_expr()
    for iter in range(0, 12):
        if(iter<8):
            pass
        else:
            assert re.match(
                r'Assignment\(Variable\(\w+\),\s\w+(/\w+)+\)',
                str(ls1[iter])
            )

def test_mul_binop():
    src1 = (
        src +
        """\
        d = a + b - c
        c = a * b + d
        s = p * q / r
        r = p * s + q / p
        """
    )
    expr1.convert_to_expr(src1, 'f')
    ls1 = expr1.return_expr()
    for iter in range(0, 12):
        if(iter<8):
            pass
        else:
            assert re.match(
                r'Assignment\(Variable\(\w+\),\s\w+(\s?[*/+-]\s?\w+)+\)',
                str(ls1[iter])
            )


def test_function():
    src1 = """\
    integer function f(a,b)
    integer :: x, y
    f = x + y
    end function
    """
    expr1.convert_to_expr(src1, 'f')
    for iter in expr1.return_expr():
        assert re.match(r'FunctionDefinition\(.*\)', str(iter))
        #assert re.match(
        #    r"FunctionDefinition\(.*,\sname\=String\(\'\w+\'\),\sparameters\=(.*?),\sbody\=CodeBlock\(.*\),\\n\s*Return\(\.*\)\)",
        #    str(iter)
        #)

def test_var():
    expr1.convert_to_expr(src, 'f')
    for iter in expr1.return_expr():
        assert re.match(
            r'Declaration\(Variable\(\w+,\s*type\=(integer|real),\s*value\=.*\)\)',
            str(iter)
        )

def test_convert_py():
    expr1.convert_to_expr(src, 'f')
    for iter in expr1.convert_to_python():
        assert re.match(
            r'\w+\s\=\s\d',
            str(iter)
        )

def test_convert_fort():
    expr1.convert_to_expr(src, 'f')
    for iter in expr1.convert_to_fortran():
        assert re.match(
            r'\s*(integer|real)\*\d\s*.*',
            str(iter)
        )

def test_convert_c():
    expr1.convert_to_expr(src, 'f')
    for iter in expr1.convert_to_c():
        assert re.match(
            r'\s*(int|double)\s\w+\s\=\s.*',
            str(iter)
        )
