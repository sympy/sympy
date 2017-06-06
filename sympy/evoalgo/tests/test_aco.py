'''
 * Author - Yajnavalkya Bandyopadhyay
 * email- 	yajnab@gmail.com
 * Civil Engineering Student
 * Techno India College of Technology
'''
import numpy as np #Import Numpy
from sympy import *
class prob:

	def func(x):
		val = 2*x[0] + 4*x[1]/x[0] + x[0]*x[1]*x[2]  #Calculate the Value
		return(val)
	#Constraints goes here

	#Number of Constraints in the form of g(x) => 0
	def cons(x):
		n_cnstr = 3
		cnstr = np.ones(n_cnstr, dtype=np.float64) #Declare the constraints array
		cnstr[0] = x[0]*x[1] - 1 #1st Constraints Goes Here
		cnstr[1] = 2*x[0] + x[1];#2nd Constraint Goes Here
		cnstr[2] = x[2]-85#3rd Constraint Goes Here
		return(cnstr)

	n_val=3
	maxiter=3
	maxroute=10
	phe=0.2
	bounds = np.ones((n_val,2), dtype=np.float64) 
	print("Enter the Lower and Upper bounds of each variable \n")
	for i in range(n_val):
		for j in range(2):
			bounds[i,j]=input("Enter");
	m = PyACO(func, cons, n_val, bounds, maxiter, maxroute, phe)
	print("Minimum Value is", m[0])
	print("Variables are", m[1])
