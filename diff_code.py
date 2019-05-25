from sympy import *
from sympy.series.formal import *
x, n = symbols('x n')
f = x**(n-2)*cos(x)

fps_sym = False

nterms = []
symb, res = S.Zero, S.Zero
if isinstance(f, Mul):
    nterm = S.One
    for term in Mul.make_args(f):
        if isinstance(term, Pow):
            if term.exp.is_symbol:
                symb = term.exp
                fps_sym = True
            elif sympify(term.exp).is_integer:
                    nterm *= x**(term.exp)
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
    nterms.append(nterm)
else:
    nterms.append(f)

f = nterms[0]

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

print(ak, xk, ind)

print(f)
'''print(fps(f,x).truncate())'''
