"""Fast algebraic simplification (expand then combine)."""

from __future__ import print_function, division

from sympy.core import Function, S, Mul, Pow, Add


def fastalgsimp (expr):
	if expr.is_Add:
		basess = _basess_from_add (expr)
	elif expr.is_Mul:
		basess = _basess_from_mul (expr)
	elif expr.is_Pow:
		basess = _basess_from_pow (expr)
	else:
		return expr

	# for bases in basess:
	# 	print ('+')
	# 	for base, exps in bases.items ():
	# 		print (f'  {base}^{exps}')

	return expr_from_basess (basess)


def _basess_from_add (expr): # -> [{base: [exp, + exp, + ...], * ...}, + bases, + ...]
	stack  = list (expr.args)
	basess = [] # [bases, + bases, + ...]

	while stack:
		arg = stack.pop ()

		if arg.is_Add:
			stack.extend (arg.args)
		elif arg.is_Mul:
			basess.extend (_basess_from_mul (arg))
		elif arg.is_Pow:
			basess.extend (_basess_from_pow (arg))
		else: # single non-algebraic factor, terminal node for algebraic processing
			basess.append ({arg: [S.One]})

	return basess


def _basess_from_mul (expr):
	stack  = list (expr.args)
	basess = [{}]

	while stack:
		arg = stack.pop ()

		if arg.is_Add:
			_basess_prod (basess, _basess_from_add (arg))
		elif arg.is_Mul:
			stack.extend (arg.args)
		elif arg.is_Pow:
			_basess_prod (basess, _basess_from_pow (arg))
		else: # single non-algebraic factor, terminal node for algebraic processing
			for bases in basess:
				bases.setdefault (arg, []).extend ([S.One])

	return basess


def _basess_from_pow (expr):
	base = expr.args [0]
	exp  = expr.args [1]

	if base.is_Pow:
		exp = list (exp.args) if exp.is_Mul else [exp]

		while base.is_Pow:
			e = base.args [1]

			if e.is_Mul:
				exp.extend (e.args)
			else:
				exp.append (e)

			base = base.args [0]

		exp = Mul (*exp) if len (exp) > 1 else exp [0]

	# # if base is sum and exponent is positive int, exponentiate and return bases
	# if base.is_Add:
	# 	if exp.is_Integer and exp.is_nonnegative:
	# 		if exp.is_zero:
	# 			return [{S.One, [S.One]}]

	# 		basess = _basess_from_add (base)
	# 		pow_   = [basess.copy ()]

	# 		for _ in range (1, int (exp)):
	# 			_basess_prod (pow_, basess)

	# 		return pow_

	# if base is a product then raise each factor to exponent if is single product
	if base.is_Mul:
		basess = _basess_from_mul (base)

		# stayed a single product, can multiply exponents and return
		if len (basess) == 1:
			_bases_mul_exp (basess [0], exp)

			return basess

		# # turned into sum of products, exponentiate if exp is nonnegative int
		# if exp.is_Integer and exp.is_nonnegative:
		# 	if exp.is_zero:
		# 		return [{S.One, [S.One]}]

		# 	pow_ = [basess.copy ()]

		# 	for _ in range (1, int (exp)):
		# 		_basess_prod (pow_, basess)

		# 	return pow_

	return [{base: [exp]}]


def _basess_prod (basess1, basess2): # product (a+b)(c+d) -> (ac+bc+ad+bd), modifies basess in place
	l1 = len (basess1)
	l2 = len (basess2)

	for i in range (l2 - 1): # (a+b) -> (a+b+a+b)
		for j in range (l1):
			basess1.append (dict ((base, exps [:]) for base, exps in basess1 [j].items ()))

	for i in range (l2): # (a+b+a+b) -> (ac+bc+ad+bd)
		for base, exps in basess2 [i].items ():
			for j in range (l1):
				basess1 [l1 * i + j].setdefault (base, []).extend (exps)


def _bases_mul_exp (bases, exp): # multiply exponent of each base by exp, modifies bases in place
	for exps in bases.values ():
		if len (exps) == 1:
			exps [0] = Mul (exp, exps [0])

		else:
			prod = Mul (exp, Add (*exps))
			del exps [1:]
			exps [0] = prod


def expr_from_basess (basess):
	prods = []

	for bases in basess:
		exps = {} # {exp: [base, * base, ...], ...}

		# combine bases with like exponents into power of product of those bases
		for base, exp in bases.items ():
			exp = Add (*exp) if len (exp) > 1 else exp [0]

			exps.setdefault (exp, []).append (base)

		prod = []

		for exp, bases in exps.items ():
			if exp.is_zero:
				prod.append (S.One)
				continue

			base = Mul (*bases) if len (bases) > 1 else bases [0]

			if exp is S.One:
				prod.append (base)
			else:
				prod.append (Pow (base, exp, evaluate = False))

		prods.append (Mul (*prod) if len (prod) > 1 else prod [0])

	return Add (*prods) if len (prods) > 1 else prods [0]


# if __name__ == '__main__':
# 	# expr = S ('x + y * (a + b)', evaluate=False)
# 	# expr = S ('(a + b)**3', evaluate=False)
# 	# expr = S ('(a*b)**3', evaluate=False)
# 	# expr = S ('((a*b**-2)**3)**2', evaluate=False)
# 	expr = S ('x**2*(x*y**2)**2', evaluate=False)
# 	print (fastalgsimp (expr))
