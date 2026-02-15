from sympy import symbols, Q, refine, Abs, sign, arg, re, im, I, conjugate

x = symbols('x')

print("--- INIZIO BATTERIA DI TEST ---")

# TEST 1: Parte Immaginaria di un numero Reale
# Se x è reale, la sua parte immaginaria deve essere 0.
expr1 = im(x)
res1 = refine(expr1, Q.real(x))
print(f"1. im(x) con x reale -> {res1} (Atteso: 0)")

# TEST 2: Coniugato di un numero Reale
# Se x è reale, il coniugato di x è x stesso.
expr2 = conjugate(x)
res2 = refine(expr2, Q.real(x))
print(f"2. conjugate(x) con x reale -> {res2} (Atteso: x)")

# TEST 3: Argomento di un numero negativo
# L'argomento (angolo) di un numero reale negativo è pi greco (180 gradi).
from sympy import pi
expr3 = arg(x)
res3 = refine(expr3, Q.negative(x))
print(f"3. arg(x) con x negativo -> {res3} (Atteso: pi)")

# TEST 4: Modulo quadro (Abs(x)**2) di numero reale
# Se x è reale, |x|^2 è uguale a x^2.
expr4 = Abs(x)**2
res4 = refine(expr4, Q.real(x))
print(f"4. Abs(x)**2 con x reale -> {res4} (Atteso: x**2)")

print("--- FINE ---")
