import sympy as sp
import re
from sympy.parsing.sympy_parser import parse_expr
from io import StringIO
from tokenize import (generate_tokens, untokenize, TokenError, NUMBER, STRING, NAME, OP, ENDMARKER, ERRORTOKEN, NEWLINE)
from typing import Tuple as tTuple, Dict as tDict, Any, Callable, \
    List, Optional, Union as tUnion
def subfinder(tokens, pattern):
    text = ''.join([tok.string for tok in pattern])
    lpat = len(pattern)
    lstr = pattern[-1].end[1]-pattern[0].start[1]
    spat = [(toknum,tokval) for toknum, tokval, _, _, _ in pattern]
    stok = [(toknum,tokval) for toknum, tokval, _, _, _ in tokens]
    res  = stok.copy()
    for i in range(len(tokens)):
        if stok[i] == spat[0] and stok[i:i+lpat] == spat and (tokens[i+lpat-1].end[1]-tokens[i].start[1] == lstr):
            res[i] = (1,text)
            del res[i+1:i+len(pattern)]
    return res


def test_generate_tokens(n):
    t = []
    input_code = StringIO(n.strip())
    for tok in generate_tokens(input_code.readline):
        t.append(tok)
    return(t)



TOKEN = tTuple[int, str]
DICT = tDict[str, Any]
TRANS = Callable[[List[TOKEN], DICT, DICT], List[TOKEN]]

def stringify_expr(s: str, local_dict: DICT, global_dict: DICT,
        transformations: tTuple[TRANS, ...]) -> str:
    """
    Converts the string ``s`` to Python code, in ``local_dict``

    Generally, ``parse_expr`` should be used.
    """

    tokens = []
    input_code = StringIO(s.strip())
    for toknum, tokval, _, _, _ in generate_tokens(input_code.readline):
        tokens.append((toknum, tokval))

    for transform in transformations:
        tokens = transform(tokens, local_dict, global_dict)

    return untokenize(tokens)


def stringify_expr_new(s: str, local_dict: DICT, global_dict: DICT,
        transformations: tTuple[TRANS, ...]) -> str:
    """
    Converts the string ``s`` to Python code, in ``local_dict``

    Generally, ``parse_expr`` should be used.
    """

    tokens = []
    input_code = StringIO(s.strip())
    for tok in generate_tokens(input_code.readline):
        tokens.append(tok)

    tok = []
    if local_dict:
        for symbol in local_dict.values():
            print(f'Symbols : {symbol}')
            print(f'local_dict = {local_dict}')
            t = test_generate_tokens(str(symbol))
            tokens = subfinder(tokens,t[0:-2])
    else:
        tokens = [(toknum,tokval) for toknum, tokval, _, _, _ in tokens]

    for transform in transformations:
        tokens = transform(tokens, local_dict, global_dict)

    return untokenize(tokens)

def stringify_expr_new_v2(s: str, local_dict: DICT, global_dict: DICT,
        transformations: tTuple[TRANS, ...]) -> str:
    """
    Converts the string ``s`` to Python code, in ``local_dict``

    Generally, ``parse_expr`` should be used.
    """

    tokens = []
    input_code = StringIO(s.strip())
    for tok in generate_tokens(input_code.readline):
        tokens.append(tok)

    tok = []
    if local_dict:
        for symbol in local_dict.values():
            print(f'Symbols : {symbol}')
            
            p = re.compile('[0-9-+. x]') #set illegal characters to omit from dictionary
            s = p.sub('', str(symbol))  #replace illegal char
            local_dict[str(symbol)] = s
            print(f'local_dict = {local_dict}')
            
            t = test_generate_tokens(str(symbol))
            tokens = subfinder(tokens,t[0:-2])
            
    else:
        tokens = [(toknum,tokval) for toknum, tokval, _, _, _ in tokens]

    for transform in transformations:
        tokens = transform(tokens, local_dict, global_dict)

    return untokenize(tokens)




def test_parse_expr(n: str) -> str:
    expr = parse_expr(f'{n} + 1',local_dict={n:sp.Symbol(n)})
   
    return(str(expr))

print(f'=== Parse expressions with local_dict Symbols ===')
n = 'product'
print('parse_expr = ', end = '')
print(parse_expr(f'{n} + 1',local_dict={n:sp.Symbol(n)}))

n = 'pro-duct'
print('parse_expr = ', end = '')
print(parse_expr(f'{n} + 1',local_dict={n:sp.Symbol(n)}))

#------------------------------------------------------------------------------
print('\nstringify_expr = ', end = '')
print(stringify_expr(f'{n} + 1',None,None,[]))

print('stringify_expr_new = ', end = '')
print(stringify_expr_new(f'{n} + 1',None,None,[]))

print('stringify_expr_new_v2 = ', end = '')
print(stringify_expr_new_v2(f'{n} + 1',None,None,[]))

#------------------------------------------------------------------------------
print('\nstringify_expr = ', end = '')
print(stringify_expr(f'{n} + 1',{n:sp.Symbol(n)},None, []))

print('stringify_expr_new = ', end = '')
print(stringify_expr_new(f'{n} + 1',{n:sp.Symbol(n)},None, []))

print('stringify_expr_new_v2 = ', end = '')
print(stringify_expr_new_v2(f'{n} + 1',{n:sp.Symbol(n)},None, []))

#------------------------------------------------------------------------------
print('\nstringify_expr = ', end = '')
print(stringify_expr(f'pro - duct + 1',{n:sp.Symbol(n)},None, []))

print('stringify_expr_new = ', end = '')
print(stringify_expr_new(f'pro - duct + 1',{n:sp.Symbol(n)},None, []))

print('stringify_expr_new_v2 = ', end = '')
print(stringify_expr_new_v2(f'pro - duct + 1',{n:sp.Symbol(n)},None, []))

#------------------------------------------------------------------------------
n = 'pro -duct'
print('\nstringify_expr = ', end = '')
print(stringify_expr(f'{n} + 1',{n:sp.Symbol(n)},None, []))

print('stringify_expr_new = ', end = '')
print(stringify_expr_new(f'{n} + 1',{n:sp.Symbol(n)},None, []))

print('stringify_expr_new_v2 = ', end = '')
print(stringify_expr_new_v2(f'{n} + 1',{n:sp.Symbol(n)},None, []))

