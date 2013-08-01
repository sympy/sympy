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
from sympy.core import Basic, Mul, Add


def _is_a_leaf(expr):
    from sympy.matrices.expressions.matexpr import MatrixSymbol
    return (isinstance(expr, Basic) and (expr.is_Atom or isinstance(expr, MatrixSymbol))) or isinstance(expr, bool)

def _is_a_collection_of_branches(expr):
    from sympy.matrices.expressions.matexpr import MatrixExpr
    return iterable(expr) and not isinstance(expr, MatrixExpr)


    
def adds_muls_cse(exprs):
    from sympy.matrices import Matrix
    from sympy.matrices.expressions.matexpr import MatrixExpr
    from sympy.matrices.expressions.matadd import MatAdd
    from sympy.matrices.expressions.matmul import MatMul

    adds = set()
    muls = set()
    matadds = set()
    matmuls = set()
    
    ### Find adds and muls #####################
    
    seen_subexp = set()
    def _find_adds_muls(expr):
        
        if _is_a_leaf(expr):
           return
        
        if _is_a_collection_of_branches(expr):
            args = expr
            
        else:
            if expr in seen_subexp:
                return
            seen_subexp.add(expr)
            
            if isinstance(expr, MatrixExpr):
                if expr.is_MatMul:
                    matmuls.add(expr)
                elif expr.is_MatAdd:
                    matadds.add(expr)
            else:
                if expr.is_Mul:
                    muls.add(expr)
                elif expr.is_Add:
                    adds.add(expr)
            
            args = expr.args
    
        map(_find_adds_muls, args)
   
        
    _find_adds_muls(exprs)
    
    
    ### Process adds and muls #####################

    split_args = dict()
    
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
                        split_args[ops[i]] = Op(*diff_i), com_op
                    
                    diff_j = op_args[j].difference(com_args)
                    op_args[j] = diff_j | set([com_op])
                    split_args[ops[j]] = Op(*diff_j), com_op
                    
                    for k in xrange(j + 1, len(op_args)):
                        if not com_args.difference(op_args[k]):
                            diff_k = op_args[k].difference(com_args)
                            op_args[k] = diff_k | set([com_op])
                            split_args[ops[k]] = Op(*diff_k), com_op
                        

    # split muls into commutative                         
    comutative_muls = set()
    for m in muls:
        c, nc = m.args_cnc(cset=True)
        if c:
            c_mul = Mul(*c)
            if nc:
                split_args[m] = c_mul, Mul(*nc)
            if len(c) > 1:
                comutative_muls.add(c_mul)
        
        
    _match_common_args(Add, adds)
    _match_common_args(Mul, comutative_muls)
    _match_common_args(MatAdd, matadds)

    return split_args
    


def tree_cse(exprs, symbols=None, split_args=None):
    from sympy.matrices import Matrix
    from sympy.matrices.expressions.matexpr import MatrixSymbol, MatrixExpr
    
    if symbols is None:
        symbols = numbered_symbols()
    else:
        # In case we get passed an iterable with an __iter__ method instead of
        # an actual iterator.
        symbols = iter(symbols)

    if split_args is None:
        split_args = dict()
        
    
    ### Find repeated subexpressions #####################

    to_eliminate = set()

    seen_subexp = set()
    def _find_repeated(expr):
        if _is_a_leaf(expr):
           return
        
        if _is_a_collection_of_branches(expr):
            # do not cse iterables
            args = expr
            
        else:
            if expr in seen_subexp:
                to_eliminate.add(expr)
                return
            
            seen_subexp.add(expr)
            
            if expr in split_args:
                args = split_args[expr]
            else:
                args = expr.args
            
        map(_find_repeated, args)
   
    _find_repeated(exprs)
    
    
    
    ### Recreate #####################
    
    replacements = []
    
    subs = dict()
    def _recreate(expr):
        if _is_a_leaf(expr):
            return expr
        
        if _is_a_collection_of_branches(expr):
            new_args = [_recreate(arg) for arg in expr]
            return type(expr)(*new_args)
            
        else:
            if expr in subs:
                return subs[expr]
            
            if expr in split_args:
                args = split_args[expr]
            else:
                args = expr.args
            
            new_expr = type(expr)(*map(_recreate, args))

            if expr in to_eliminate:
                sym = next(symbols)
                if isinstance(expr, MatrixExpr):
                    sym = MatrixSymbol(str(sym).upper(), *expr.shape)
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
    split_args = adds_muls_cse(exprs)
    return tree_cse(exprs, symbols, split_args)






