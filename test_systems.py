from sympy import reset_solution, last_solution, symbols, solve

x, y, z = symbols('x, y, z')

# 186 - 189 are pages in Vilenkin
# Here are some tests
# just change the value on "n"
n = 5 # number of the test
# 0 - my small e. g.
# 1 - 186
# 2 - 187 (1)
# 3 - 187 (2)
# 4 - 189 (# 59 (1))
# 5 - my e.g. (Gauss solution)
# 6 - your test

if (n == 0):
	reset_solution()
	res = (solve([x + y - 10, x * y + 100], [x, y], dict = True))
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
elif (n == 5):
	reset_solution()
	res = (solve([x + y + z - 6, x - y - z + 4, x + y + 2 * z - 9], [x, y, z], dict = True))
	for i in res:
		for j in i:
			print j, i[j]
		print
	print
	R = last_solution()
	for i in R:
		print i
	print
elif (n == 6):
	reset_solution()
	res = (solve([], [], dict = True))
	for i in res:
		for j in i:
			print j, i[j]
		print
	print
	R = last_solution()
	for i in R:
		print i
	print