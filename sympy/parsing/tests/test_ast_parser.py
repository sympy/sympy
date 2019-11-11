from sympy.parsing.ast_parser import parse_expr

def test_parse_expr():
    parse_expr('a + b', {}) == a + b
