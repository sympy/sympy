"""Fast algebraic simplification (expand then cancel then combine exponents)."""

from __future__ import print_function, division

from sympy.core import E, S, Mul, Pow, Add, factor_terms
from sympy.functions import exp



def _mulsimp(expr):
    """A simple simplify function to prevent expression blowup during multiplication."""

    from sympy.core import Add, Mul, Pow, factor_terms
    from sympy.functions.elementary.exponential import ExpBase
    from sympy.polys import together, factor, cancel
    from sympy.simplify.radsimp import radsimp, fraction, _mexpand
    from sympy.simplify.powsimp import powsimp
    from sympy.simplify.simplify import signsimp, bottom_up, count_ops
    from sympy.utilities.iterables import has_variety

    def shorter(*choices):
        return choices[0] if not has_variety(choices) else min(choices, key=count_ops)

    # return cancel(expr)

    # return together(cancel(expr), deep=True)

    expr2 = cancel(expr)
    return shorter(expr, expr2, together(expr, deep=True), together(expr2, deep=True))

    expr2 = cancel(expr)
    expr3 = shorter(together(expr, deep=True), together(expr2, deep=True))
    expr  = shorter(expr, expr2, expr3)
    return expr2

    # expr  = bottom_up(expr, lambda w: getattr(w, 'normal', lambda: w)())
    # expr  = Mul(*powsimp(expr).as_content_primitive())
    _e    = cancel(expr)
    expr1 = shorter(_e, _mexpand(_e).cancel())  # issue 6829
    expr2 = shorter(together(expr, deep=True), together(expr1, deep=True))
    expr  = shorter(expr2, expr1, expr)
    # expr  = shorter(powsimp(expr, combine='exp', deep=True), powsimp(expr), expr)

    return expr


def fastalgsimp(expr):
    expr = expr.expand(power_exp=False, log=False, multinomial=False, basic=False)
    return _mulsimp(expr)

    if expr.is_Add:
        basess = _basess_from_add(expr)
    elif expr.is_Mul:
        basess = [_bases_from_mul(expr)]
    elif expr.is_Pow or isinstance(expr, exp):
        basess = [_bases_from_pow(expr)]
    else:
        return expr

    # for bases in basess:
    #     print ('+')
    #     for base, exps in bases.items ():
    #         print (f'  ({base})^{exps}')

    expanded = True

    while expanded:
        basess, expanded = _basess_expand_adds (basess)

    return _mulsimp(expr_from_basess(basess))


def _basess_from_add(expr): # -> [{base: [exp, + exp, + ...], * ...}, + bases, + ...]
    stack  = list(expr.args)
    basess = [] # [bases, + bases, + ...]

    while stack:
        arg = stack.pop()

        if arg.is_Add:
            stack.extend(arg.args)
        elif arg.is_Mul:
            basess.append(_bases_from_mul(arg))
        elif arg.is_Pow or isinstance(arg, exp):
            basess.append(_bases_from_pow(arg))
        else: # single non-algebraic factor, terminal node for algebraic processing
            basess.append({arg: [S.One]})

    return basess


def _bases_from_mul(expr):
    stack = list(expr.args)
    bases = {}

    while stack:
        arg = stack.pop()

        if arg.is_Mul:
            stack.extend(arg.args)

        elif arg.is_Pow or isinstance(arg, exp):
            pow_ = _bases_from_pow(arg)

            for base, exp_ in pow_.items():
                bases.setdefault(base, []).extend(exp_)

        else: # single non-algebraic factor, terminal node for algebraic processing
            bases.setdefault(arg, []).extend([S.One])

    return bases


def _bases_from_pow(expr):
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
        return {S.One: [S.One]}

    # if base is a product then raise each factor to exponent if is not sum
    if base.is_Mul:
        bases = _bases_from_mul(base)

        for exps in bases.values():
            if len(exps) == 1:
                exps[0] = Mul(exp_, exps[0])

            else:
                prod = Mul(exp_, Add(*exps))
                del exps[1:]
                exps[0] = prod

        return bases

    return {base: [exp_]}


def _basess_mul(basess1, basess2): # product(a+b)(c+d) ->(ac+bc+ad+bd), modifies basess1 in place
    l1 = len(basess1)
    l2 = len(basess2)

    for i in range(l2 - 1): # (a+b) -> (a+b+a+b)
        for j in range(l1):
            basess1.append(dict((base, exps[:]) for base, exps in basess1[j].items()))

    for i in range(l2): # (a+b+a+b) -> (ac+bc+ad+bd)
        for base, exps in basess2[i].items():
            for j in range(l1):
                basess1[l1*i+j].setdefault(base, []).extend(exps)


# def _bases_mul_exp(bases, exp_): # multiply exponent of each base by exp_, modifies bases in place
#     for exps in bases.values():
#         if len(exps) == 1:
#             exps[0] = Mul(exp_, exps[0])

#         else:
#             prod = Mul(exp_, Add(*exps))
#             del exps[1:]
#             exps[0] = prod


def _basess_expand_adds(basess):
    expanded = False
    basess2  = []

    for bases in basess:
        basess3 = [{S.One: [S.One]}]

        for base, exps in bases.items():
            exp_ = Add(*exps) if len(exps) > 1 else exps[0]

            if exp_.is_zero and not base.is_zero:
                for bases3 in basess3:
                    bases3.setdefault(S.One, []).append(S.One)

            elif base.is_Add and exp_ is S.One:
                expanded = True
                basess4  = _basess_from_add(base)

                _basess_mul(basess3, basess4)

            else:
                for bases3 in basess3:
                    bases3.setdefault(base, []).append(exp_)

        basess2.extend(basess3)

    return basess2, expanded


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

    # expr = S('x + y * (a + b)', evaluate=False)
    # expr = S('(a + b)**3', evaluate=False)
    # expr = S('(a*b)**3', evaluate=False)
    # expr = S('((a*b**-2)**3)**2', evaluate=False)
    # expr = S('x**2*(x*y**2)**2', evaluate=False)
    # expr = S('exp(3)*exp(-2)', evaluate=False)
    # expr = S('(exp(3)*exp(-5))**3', evaluate=False)

    expr = S('(2*a*(x+y)) / (a*(x+y))', evaluate=False)
    # expr = S('(2*a*(x+y)) / (x+y)', evaluate=False)
    # expr = S('(2*a*(x+y)) / (x+z)', evaluate=False)

    print(fastalgsimp(expr))
    # print(type (expr))
    # print(expr.args)
    # print(factor_terms(expr))
