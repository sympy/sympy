'''
Utility functions for Constraints in Rubi
'''

def NonzeroQ(expr):
    return expr != 0

def FreeQ(nodes, var):
    return not any(expr.has(var) for expr in nodes)
