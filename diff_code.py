from sympy import *
from sympy.series.formal import *
x, n = symbols('x n')
'''
fps_sym = False

symb, res = S.Zero, S.Zero
if isinstance(f, Mul):
    nterm = S.One
    for term in Mul.make_args(f):
        if isinstance(term, Pow):
            # for x**n symbolic terms
            if term.exp.is_symbol:
                symb = term.exp
                fps_sym = True
            # if power is an integer
            elif sympify(term.exp).is_integer:
                    nterm *= term

            # for x**(n-2) or x**(n/2) symbolic terms
            else:
                a, b, subt = S.Zero, S.Zero, S.Zero
                if isinstance(term.exp, Add):
                    a, b = Add.make_args(term.exp)
                    if a.is_symbol or b.is_symbol:
                        if a.is_symbol:
                            subt = a
                        else:
                            subt = b
                        res = term.exp - subt
                        if sympify(res).is_integer:
                            nterm *= x**res
                        symb = subt
                        fps_sym = True

                if isinstance(term.exp, Mul):
                    symb = term.exp
                    fps_sym = True
        else:
            if not (isinstance(term, Pow) and term.exp.is_symbol):
                nterm *= term
    f = nterm
print(f)

# from here on it's x0=0 and dir=1 handling
k = Dummy('k')

result = hyper_algorithm(f, x, k)

if result is None:
    print("0")

ak = sequence(result[0], (k, result[2], oo))
if fps_sym:
    xk = sequence(x**k*x**(symb), (k, 0, oo))
    ind = expand(result[1]*x**(symb))
else:
    xk = sequence(x**k, (k, 0, oo))
    ind = result[1]

print(f)
'''
ind = sqrt(x) + x**(n/2)
for t in Add.make_args(ind):
    if isinstance(t, Mul):
        for term in Mul.make_args(t):
            if (isinstance(term, Pow) and term.exp.is_symbol):
                print(S.One)
    else:
        if isinstance(t, Pow):
            if isinstance(t.exp, Mul) or t.exp.is_symbol:
                print(S.One)