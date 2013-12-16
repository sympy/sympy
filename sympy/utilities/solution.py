"""Make step-by-step solution. """

import gc, sys
from sympy.printing import latex
from sympy.core import sympify

solution_list = []

def find_names(obj):
    frame = sys._getframe()
    for frame in iter(lambda: frame.f_back, None):
        frame.f_locals
    result = []
    for referrer in gc.get_referrers(obj):
        if isinstance(referrer, dict):
            for k, v in referrer.iteritems():
                if v is obj:
                    result.append(k)
    return result

def find_name(obj):
    frame = sys._getframe()
    for frame in iter(lambda: frame.f_back, None):
        frame.f_locals
    for referrer in gc.get_referrers(obj):
        if isinstance(referrer, dict):
            for k, v in referrer.iteritems():
                if v is obj and k != "variable":
                    return k
    return None

def add_comment(cm):
    solution_list.append('_'+cm)
    
def add_step(variable):
    """Add a variable and its value into solution"""
    var = find_name(variable)
    r = repr(variable)
    try:
        r = latex(sympify(r, evaluate=False))
    except:
        print r
    solution_list.append(var + " = " + r)

def add_eq(l, r):
    """Add an equality into solution"""
    l = repr(l)
    try:
        l = latex(sympify(l, evaluate=False))
    except:
        print l
    r = repr(r)
    try:
        r = latex(sympify(r, evaluate=False))
    except:
        print r
    solution_list.append(l + " = " + r)

    
def add_exp(exp):
    """Add an expression into solution"""
    r = repr(exp)
    try:
        r = latex(sympify(r, evaluate=False))
    except:
        print r
    solution_list.append(r, evaluate=False)

def reset_solution():
    """Clear previos solution"""
    print("New solution")
    del solution_list[:]

def start_subroutine(name):
    """Start add soubroutine steps"""
    print("Start subroutine", name)

def cancel_subroutine():
    """Cancel all steps of current subroutine"""
    print("Cancel subroutine")

def commit_subroutine():
    """Finish current subroutine"""
    print("Finish subroutine")

def last_solution():
    return solution_list

def expr_to_str(expr):
    """Get a text representation of the expression"""
    r = repr(expr)
    try:
        r = latex(sympify(r, evaluate=False))
        return r
    except:
        return r

