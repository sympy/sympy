"""Fast algebraic simplification (expand then cancel then combine exponents)."""

from __future__ import print_function, division

from sympy.core import E, S, Mul, Pow, Add, factor_terms
from sympy.functions import exp


def fastalgsimp(expr, multinomial=False):
    if expr.is_Add:
        basess = _basess_from_add(expr, multinomial=multinomial)
    elif expr.is_Mul:
        basess = _basess_from_mul(expr, multinomial=multinomial)
    elif expr.is_Pow or isinstance(expr, exp):
        basess = _basess_from_pow(expr, multinomial=multinomial)
    else:
        return expr

    # for bases in basess:
    #     print ('+')
    #     for base, exps in bases.items ():
    #         print (f'  {base}^{exps}')

    return expr_from_basess(basess)


def _basess_from_add(expr, multinomial=False): # -> [{base: [exp, + exp, + ...], * ...}, + bases, + ...]
    stack  = list(expr.args)
    basess = [] # [bases, + bases, + ...]

    while stack:
        arg = stack.pop()

        if arg.is_Add:
            stack.extend(arg.args)
        elif arg.is_Mul:
            basess.extend(_basess_from_mul(arg, multinomial=multinomial))
        elif arg.is_Pow or isinstance(arg, exp):
            basess.extend(_basess_from_pow(arg, multinomial=multinomial))
        else: # single non-algebraic factor, terminal node for algebraic processing
            basess.append({arg: [S.One]})

    return basess


def _basess_from_mul(expr, multinomial=False):
    stack  = list(expr.args)
    basess = [{}]

    while stack:
        arg = stack.pop()

        if arg.is_Add:
            _basess_prod(basess, _basess_from_add(arg, multinomial=multinomial))
        elif arg.is_Mul:
            stack.extend(arg.args)
        elif arg.is_Pow or isinstance(arg, exp):
            _basess_prod(basess, _basess_from_pow(arg, multinomial=multinomial))
        else: # single non-algebraic factor, terminal node for algebraic processing
            for bases in basess:
                bases.setdefault(arg, []).extend([S.One])

    return basess


def _basess_from_pow(expr, multinomial=False):
    if expr.is_Pow:
        base = expr.args[0]
        exp_ = expr.args[1]
    else: # isinstance(expr, exp)
        base = E
        exp_ = expr.args[0]

    if base.is_Pow:
        exp_ = list(exp_.args) if exp_.is_Mul else [exp_]

        while base.is_Pow:
            e = base.args[1]

            if e.is_Mul:
                exp_.extend(e.args)
            else:
                exp_.append(e)

            base = base.args[0]

        if isinstance(base, exp):
            e = base.args[0]

            if e.is_Mul:
                exp_.extend(e.args)
            else:
                exp_.append(e)

            base = E

        exp_ = Mul(*exp_) if len(exp_) > 1 else exp_[0]

    if exp_.is_zero and not base.is_zero:
        return [{S.One, [S.One]}]

    # if base is sum and exponent is positive int, exponentiate and return bases
    if base.is_Add:
        if multinomial and exp_.is_Integer and exp_.is_nonnegative:
            basess = _basess_from_add(base, multinomial=multinomial)
            pow_   = basess.copy()

            for _ in range(1, int(exp_)):
                _basess_prod(pow_, basess)

            return pow_

    # if base is a product then raise each factor to exponent if is not sum
    elif base.is_Mul:
        basess = _basess_from_mul(base, multinomial=multinomial)

        # stayed a single product, can multiply exponents and return
        if len(basess) == 1:
            _bases_mul_exp(basess[0], exp_)

            return basess

        # turned into sum of products, exponentiate if exp_ is nonnegative int
        if multinomial and exp_.is_Integer and exp_.is_nonnegative:
            pow_ = basess.copy()

            for _ in range(1, int(exp_)):
                _basess_prod(pow_, basess)

            return pow_

    return [{base: [exp_]}]


def _basess_prod(basess1, basess2): # product(a+b)(c+d) ->(ac+bc+ad+bd), modifies basess1 in place
    l1 = len(basess1)
    l2 = len(basess2)

    for i in range(l2 - 1): #(a+b) ->(a+b+a+b)
        for j in range(l1):
            basess1.append(dict((base, exps[:]) for base, exps in basess1[j].items()))

    for i in range(l2): #(a+b+a+b) ->(ac+bc+ad+bd)
        for base, exps in basess2[i].items():
            for j in range(l1):
                basess1[l1*i+j].setdefault(base, []).extend(exps)


def _bases_mul_exp(bases, exp_): # multiply exponent of each base by exp_, modifies bases in place
    for exps in bases.values():
        if len(exps) == 1:
            exps[0] = Mul(exp_, exps[0])

        else:
            prod = Mul(exp_, Add(*exps))
            del exps[1:]
            exps[0] = prod


def expr_from_basess(basess):
    prods = []

    for bases in basess:
        exps = {} # {exp: [base, * base, ...], ...}

        # combine bases with like exponents into power of product of those bases
        for base, exp_ in bases.items():
            exp_ = Add(*exp_) if len(exp_) > 1 else exp_[0]

            exps.setdefault(exp_, []).append(base)

        prod = []

        for exp_, bases in exps.items():
            if exp_.is_zero:
                prod.append(S.One)
                continue

            base = Mul(*bases) if len(bases) > 1 else bases[0]

            if exp_ is S.One:
                prod.append(base)
            elif base is E:
                prod.append(exp(exp_))
            else:
                prod.append(Pow(base, exp_, evaluate=False))

        prods.append(Mul(*prod) if len(prod) > 1 else prod[0])

    return Add(*prods) if len(prods) > 1 else prods[0]


if __name__ == '__main__': # DEBUG!
    from sympy.core import factor_terms

    # expr = S('x + y *(a + b)', evaluate=False)
    # expr = S('(a + b)**3', evaluate=False)
    # expr = S('(a*b)**3', evaluate=False)
    # expr = S('((a*b**-2)**3)**2', evaluate=False)
    # expr = S('x**2*(x*y**2)**2', evaluate=False)
    # expr = S('(exp(3)*exp(-5))**3', evaluate=False)

    # expr = S('(2*a*(x+y)) / (a*(x+y))', evaluate=False)
    expr = S('(2*a*(x+y)) / (x+y)', evaluate=False)

    print(fastalgsimp(expr))
    print(type (expr))
    print(expr.args)
    print(factor_terms(expr))
