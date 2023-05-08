from sympy import *
from sympy.printing.c import C89CodePrinter
x = Symbol('x')
y = Symbol('y')

myPrinterObj = C89CodePrinter()

# No CodeBlock functions
expr = x
new_expr, code_blocks = myPrinterObj.doblocks(expr)

print(expr)
print()
print('######', 'Final code to print', '######')
for code_block in code_blocks:
  print(myPrinterObj.doprint(code_block))
print(myPrinterObj.doprint(new_expr))

print()
print()
# Multiple CodeBlock functions 
expr = 2*Sum(x, (x,0,1))+2*Sum(x**2, (x,0,1))
new_expr, code_blocks = myPrinterObj.doblocks(expr)

print(expr)
print()
print('######', 'Final code to print', '######')
for code_block in code_blocks:
  print(myPrinterObj.doprint(code_block))
print(myPrinterObj.doprint(new_expr))