#!/usr/bin/env python
from sympy.parsing.latex import parse_latex_lark

# Test cases
test_cases = [
    '[x+y]z',
    '{x+y}z',
    '\\{x+y\\}z',
    'x+y',
]

for test in test_cases:
    try:
        result = parse_latex_lark(test)
        print(f'✓ Test "{test}": {result}')
    except Exception as e:
        print(f'✗ Test "{test}": {type(e).__name__}: {e}')
