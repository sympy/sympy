# Strategies to traverse a SymPy Tree

def top_down(rule):
    def top_down_traversal(expr):
        expr = rule(expr)
        return expr.__class__(*map(rule, expr.args))
    return top_down_traversal

def bottom_up(rule):
    def bottom_up_traversal(expr):
        return rule(expr.__class__(*map(rule, expr.args)))
    return bottom_up
