'''Fast CSE

Fast CSE supports expressions, Matrix(es), and collections of these.

TODO:
- accept optimizations
- accept postprocess
- CSE of non-commutative Mul terms
- ?add MatrixExprs support?
- unit tests
'''

from sympy.utilities.iterables import numbered_symbols, ordered
from sympy.core.compatibility import iterable
from sympy.core import Basic, Mul, Add, Pow
from sympy.core.singleton import S
from sympy.core.function import _coeff_isneg


_option_extract_minus_sign = 0
_option_order_mul_args = 1
    
def opt_cse(exprs):
    from sympy.matrices import Matrix

    adds = set()
    muls = set()
    
    opt_subs = dict()
    
    ### Find adds and muls #####################
    
    seen_subexp = set()
    def _find_opts(expr):
        
        if isinstance(expr, Basic) and expr.is_Atom:
           return
        
        if iterable(expr):
            map(_find_opts, expr)
            return
            
        if expr in seen_subexp:
            return expr
        seen_subexp.add(expr)
            
        map(_find_opts, expr.args)
        
        if (_option_extract_minus_sign and expr.could_extract_minus_sign()) \
            or (not _option_extract_minus_sign and _coeff_isneg(expr)):
            neg_expr = -expr
            opt_subs[expr] = Mul, (S.NegativeOne, neg_expr)
            seen_subexp.add(neg_expr)
            expr = neg_expr
            
        if expr.is_Mul:
            muls.add(expr)
            
        elif expr.is_Add:
            adds.add(expr)
            
        elif expr.is_Pow:
            if _coeff_isneg(expr.exp):
                opt_subs[expr] = Pow, (Pow(expr.base, -expr.exp), S.NegativeOne)
            elif _option_extract_minus_sign and expr.exp in opt_subs:
                Op, args = opt_subs[expr.exp]
                if Op is Mul and args[0] is S.NegativeOne:
                    opt_subs[expr] = Pow, (Pow(expr.base, Mul(*args[1:])), S.NegativeOne)
        
        
    _find_opts(exprs)
    
    
    ### Process adds and muls #####################
    
    def _match_common_args(Op, ops):
        ops = list(ordered(ops))
        
        op_args = [set(e.args) for e in ops]
        for i in xrange(len(op_args)):
            for j in xrange(i + 1, len(op_args)):
                com_args = op_args[i].intersection(op_args[j])
                if len(com_args) > 1:
                    com_op = Op(*com_args)
                    
                    diff_i = op_args[i].difference(com_args) 
                    op_args[i] = diff_i | set([com_op])
                    if diff_i:
                        opt_subs[ops[i]] = Op, (Op(*diff_i), com_op)
                    
                    diff_j = op_args[j].difference(com_args)
                    op_args[j] = diff_j | set([com_op])
                    opt_subs[ops[j]] = Op, (Op(*diff_j), com_op)
                    
                    for k in xrange(j + 1, len(op_args)):
                        if not com_args.difference(op_args[k]):
                            diff_k = op_args[k].difference(com_args)
                            op_args[k] = diff_k | set([com_op])
                            opt_subs[ops[k]] = Op, (Op(*diff_k), com_op)
                        

    # split muls into commutative                         
    comutative_muls = set()
    for m in muls:
        c, nc = m.args_cnc(cset=True)
        if c:
            c_mul = Mul(*c)
            if nc:
                opt_subs[m] = Mul, (c_mul, Mul(*nc))
            if len(c) > 1:
                comutative_muls.add(c_mul)
        
        
    _match_common_args(Add, adds)
    _match_common_args(Mul, comutative_muls)

    return opt_subs
    


def tree_cse(exprs, symbols=None, opt_subs=None):
    from sympy.matrices import Matrix
    
    if symbols is None:
        symbols = numbered_symbols()
    else:
        # In case we get passed an iterable with an __iter__ method instead of
        # an actual iterator.
        symbols = iter(symbols)

    if opt_subs is None:
        opt_subs = dict()
        
    
    ### Find repeated subexpressions #####################

    to_eliminate = set()

    seen_subexp = set()
    def _find_repeated(expr):
        if isinstance(expr, Basic) and expr.is_Atom:
           return
        
        if iterable(expr):
            # do not cse iterables
            args = expr
            
        else:
            if expr in seen_subexp:
                to_eliminate.add(expr)
                return
            
            seen_subexp.add(expr)
            
            if expr in opt_subs:
                args = opt_subs[expr][1]
            else:
                args = expr.args
            
        map(_find_repeated, args)
   
    _find_repeated(exprs)
    
    
    
    ### Recreate #####################
    
    replacements = []
    
    subs = dict()
    def _recreate(expr):
        
        if isinstance(expr, Basic) and expr.is_Atom:
            return expr
        
        if iterable(expr):
            new_args = [_recreate(arg) for arg in expr]
            return type(expr)(*new_args)
            
        
        if expr in subs:
            return subs[expr]
        
        if expr in opt_subs:
            Op, args = opt_subs[expr]
        else:
            Op, args = type(expr), expr.args
            if _option_order_mul_args and Op is Mul:
                c, nc = expr.args_cnc()
                args = list(ordered(c)) + nc
        
        new_expr = Op(*map(_recreate, args))

        if expr in to_eliminate:
            sym = next(symbols)
            subs[expr] = sym
            replacements.append((sym, new_expr))
            return sym
        
        else:
            return new_expr
        
    
    single = False
    if isinstance(exprs, Basic) or isinstance(exprs, Matrix): # if only one expression or one matrix is passed
        exprs = [exprs]
        single = True
            
    
    reduced_exprs = []
    
    for expr in exprs:
        if isinstance(expr, Matrix):
            reduced_exprs.append(expr.applyfunc(_recreate))
        else:
            reduced_exprs.append(_recreate(expr))
    
    if single:
        reduced_exprs = reduced_exprs[0]
    
    return replacements, reduced_exprs

def fast_cse(exprs, symbols=None):
    opt_subs = opt_cse(exprs)
    return tree_cse(exprs, symbols, opt_subs)






