is_commutative == False:
========================

Passing commutative = False as argument in some operations, may lead to unexpected answers.

For ex:

	>>> x,y = symbols("x y")
	>>> a,b = symbols("a b",commutative = True)
	>>> c,d = symbols('c d',commutative = False)
	>>> pprint((x**(a+b)).expand())
	    a  b
	   x ⋅x 
	>>> pprint((x**(c+d)).expand())  
	  c + d
	 x
	
So x**(a + b) ! = (x**a*x**b) when a and b are non commutative.
This is because here the expansion is not true in general, so SymPy does not perform it.

Here's a proof which shows:

	>>> pprint(exp(a + b).expand())
		a  b
	          ℯ ⋅ℯ 		

	>>> pprint(exp(a)*exp(b))
		a  b
	          ℯ ⋅ℯ 
	
i.e exp(a+b) = exp(a)*exp(b)

However proof doesn't work for general a,b.
Since if a,b = non-commutative,

	>>> pprint(exp(a + b).expand() )
		  a + b
		 ℯ 

	>>> pprint(exp(a)*exp(b))
	  a  b
	 ℯ ⋅ℯ

i.e exp(a+b) != exp(a)*exp(b)
To explain the above proof we can also consider 2 non-commutating matrices.

	>>> A = Matrix([[0, 1], [1, 0]])
	>>> B = Matrix([[1, 1], [0, 1]])
	>>> exp(A + B)

	Matrix([
	[                 exp(1 - sqrt(2))/2 + exp(1 + sqrt(2))/2, -sqrt(2)*exp(1 - sqrt(2))/2 + sqrt(2)*exp(1 + sqrt(2))/2],
	[-sqrt(2)*exp(1 - sqrt(2))/4 + sqrt(2)*exp(1 + sqrt(2))/4,                  exp(1 - sqrt(2))/2 + exp(1 + sqrt(2))/2]])
	>>> exp(A)*exp(B)
	Matrix([
	[ E*(exp(-1)/2 + E/2), E*(-exp(-1)/2 + E/2) + E*(exp(-1)/2 + E/2)],
	[E*(-exp(-1)/2 + E/2), E*(-exp(-1)/2 + E/2) + E*(exp(-1)/2 + E/2)]])
	>>> A*B
	Matrix([
	[0, 1],
	[1, 1]])
	>>> B*A
	Matrix([
	[1, 1],
	[1, 0]])


See also:
https://www.physicsforums.com/threads/proof-of-commutative-property-in-exponential-matrices-using-power-series.578887/


Similar behaviour can be seen when non commutative symbols are used with powsimp

	>>> pprint(powsimp(x**a*x**b))
	   a + b
	  x

	>>> pprint(powsimp(x**c*x**d))
	   c  d
	  x ⋅x

In some cases it might be impossible to determine the behaviour of a function if non commutative symbols are taken.


See also:

	There's a function defined by sympy "is_commutative" to define if an expression is commutative or not.
	https://github.com/sympy/sympy/blob/master/sympy/functions/elementary/piecewise.py
