from sympy.utilities.solution import reset_solution, last_solution
from sympy import symbols
from sympy import solve
from sympy import Poly

x, y = symbols('x,y')


# 185 - 189 are pages in Vilenkin
# Here are some tests
n = 0 # number of the test
# 0 - 185
# 1 - 186
# 2 - 187 (1)
# 3 - 187 (2)
# 4 - 189 (# 59 (1))

if (n == 0):
	reset_solution()
	res = (solve([2 * x ** 2 + y - 4, x ** 4 + y ** 2 - 16], [x, y], dict = True))
	for i in res:
		for j in i:
			print j, i[j]
		print
	print
	R = last_solution()
	for i in R:
		print i
	print
elif (n == 1):
	reset_solution()
	res = (solve([x ** 3 - y ** 3 - 3 * x ** 2 + 3 * y ** 2 * x + 2, x ** 2 - x ** 2 * y - 1], [x, y], dict = True))
	for i in res:
		for j in i:
			print j, i[j]
		print
	print
	R = last_solution()
	for i in R:
		print i
	print
elif (n == 2):
	reset_solution()
	res = (solve([x + x * y + y - 11, x ** 2 * y + x * y ** 2 - 30], [x, y], dict = True))
	for i in res:
		for j in i:
			print j, i[j]
		print
	print
	R = last_solution()
	for i in R:
		print i
	print
elif (n == 3):
	reset_solution()
	res = (solve([2 * x ** 2 - 3 * x * y + y ** 2, y ** 2 - x ** 2 - 12], [x, y], dict = True))
	for i in res:
		for j in i:
			print j, i[j]
		print
	print
	R = last_solution()
	for i in R:
		print i
	print
elif (n == 4):
	reset_solution()
	res = (solve([x - 2 * y - 6, 5 * x + 1], [x, y], dict = True))
	for i in res:
		for j in i:
			print j, i[j]
		print
	print
	R = last_solution()
	for i in R:
		print i
	print
