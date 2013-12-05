"""Make step-by-step solution. """

import gc, sys

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
    print "->", cm
    solution_list.append(cm)
    
def add_step(variable):
    """Add a variable and its value into solution"""
    var = find_name(variable)
    r = repr(variable)
    print "->", var, "=", r
    solution_list.append(var + " = " + r)
    v = str(variable)
    if v != r:
        print "->", var, "=", v
        solution_list.append(var + " = " + v)
    
def add_eq(l, r):
    """Add an equality into solution"""
    rr = repr(l) + " = " + repr(r)
    print "->", rr
    solution_list.append(rr)
    v = str(r)
    # if v != rr:
        # print "->", repr(l), "=", v
        # solution_list.append(repr(l) + " = " + v)

def add_exp(exp):
    """Add an expression into solution"""
    r = repr(exp)
    print "-> ", r
    solution_list.append(r)

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
