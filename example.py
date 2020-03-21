from sympy.parsing.sym_expr import SymPyExpression
src = """
void func(int a=1, float b=2.0){
    int d=0;
}
"""
a=SymPyExpression(src, 'c')
print(a.return_expr())