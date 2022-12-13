from sympy import *

x = symbols('x')

funcs = [cos, sin, tan, cot , csc , sec ]
#, acos, asin, atan, acot , acsc , asec, cosh, sinh, tanh, coth, csch, sech, acosh, asinh, atanh, acoth, acsch, asech, I*oo, -I*oo, zoo
points = [
    0, 1, -1, I, -I, 1+I, 1-I, -1+I, -1-I, oo, -oo, pi, pi/2
]

bad = 0
total = 0
for f1 in funcs:
    for f2 in funcs:
        total += 1
        expr1 = f1(x)
        expr2 = expr1.rewrite(f2)
        fine = True
        for p in points:
            if expr1.subs(x, p).n() != expr2.subs(x, p).n():
                print()
                print(f1, f2, p)
                print(expr1, ' --> ', expr2)
                print(expr1.subs(x, p), expr2.subs(x, p))
                print(expr1.subs(x, p).n(5), expr2.subs(x, p).n(5))
                fine = False
        if not fine:
            bad += 1

print(bad, 'bad rewrites out of', total)