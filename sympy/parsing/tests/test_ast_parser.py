from sympy import symbols, S, Rational, Lambda
from sympy.parsing.ast_parser import parse_expr
from sympy.utilities.pytest import raises
from sympy.core.sympify import SympifyError

def test_parse_expr():
    a, b = symbols('a, b')
    # tests issue_16393
    parse_expr('a + b', {}) == a + b
    raises(SympifyError, lambda: parse_expr('a + ', {}))

    # tests Transform.visit_Num
    parse_expr('1 + 2', {}) == S(3)
    parse_expr('1 + 2.0', {}) == S(3.0)

    # tests Transform.visit_Name
    parse_expr('Rational(1, 2)', {}) == S(1)/2
    parse_expr('a', {'a': a}) == a
    
    #tests Transform.visit_Lambda
    #Currently still causes a bug
    try:
        parse_expr('lambda a: 2',{}) == a
    except:
        """Lambda AST Parsing still has issue with recognizing expression type in compile()
        We can dump a, the transformed node expression and we see the expected:
        
        Expression(
            body=Call(
                func=Name(id='Lambda', ctx=Load()), 
                args=[Tuple(elts=[arg(arg='a', annotation=None, type_comment=None)], ctx=Load()), 
                      Call(func=Name(id='Integer', ctx=Load()), args=[Constant(value=2, kind=None)], keywords=[])
                ]#Args
                , keywords=[]
            )#Body Call
        )#Expression
        
        But when we pass it into compile we get a type error claiming this object is a AST arg object.
        
        expected some sort of expr, but got <_ast.arg object at 0x10075da00>
        
        May be an issue with fix_missing_locations?"""
        print("Lambda AST Parsing Error)
