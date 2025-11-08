from sympy import simplify

def algebraic_rewrite(expr):
    rewritten = expr.rewrite()
    simplified = simplify(rewritten)
    if simplified != expr:
        print(f"Rewritten {expr} → {simplified}")
    return simplified