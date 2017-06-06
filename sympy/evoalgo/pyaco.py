'''
 * Author - Yajnavalkya Bandyopadhyay
 * email- 	yajnab@gmail.com
 * Civil Engineering Student
 * Techno India College of Technology
'''

import numpy as np
import random

'''
	func = Function which will be evaluated for maximization or Minimization
	cons = Constraint Functions
	n_val = Number of Variables for the Problem to be Optimized
	bounds = Array Storing bounds where bounds (n_val,0) is the lower limit of the variable and (n_val,1) is the upper bound
	maxiter = Maximum number of Iterations
	maxroute = Maximum number of route the ants can make
	phe = Pheromone Evaporation Rate
'''

def PyACO(func, cons, n_val, bounds, maxiter, maxroute, phe):

	''' Class iterinfo is introduced to carry the best optimized value among the iterations
		Contains 2 elements
		1. value - Value
		2. xval - Array of Generated Number
	'''
	class iterinfo:
		def __init__(self, **kwds):
			self.__dict__.update(kwds)

	min_res =None #Variable for Accessing iterinfo class

	gstr = 0 #Counitng for initial value for minimum feed
	#Code for the Solution generation
	for gi in range(maxiter):

		'''
			class routeinfo is introduced to handle each and every successful values generated satisfying the constraints
			which will have the 3 elements
			1. value - will carry the value
			2. xval - Will carry the Array of Number Generated
			3. ph - Pheromone Count

			Class minv is introduced to handle the position, value and the variable value of the minimum instance
			Elements are
			1. value - Contain the Value
			2. ps - Contain the position respect to the cstr , i.e the variable introduced for counter for Validated results
		'''

		class routeinfo:
			def __init__(self, **kwds):
				self.__dict__.update(kwds)
		route = {} #Array for Accessing Class routeinfo

		class minv:
			def __init__(self, **kwds):
				self.__dict__.update(kwds)



		mins = None #Variable for Accessing Class minv

		cstr =int(0) #Counter for Validated first input

		for rc in range(maxroute):
			tval = np.zeros(n_val, dtype=np.float64);

			for i in range(n_val):
				tval[i]=random.uniform(bounds[i,0],bounds[i,1])

			result = calc(func, cons, tval) #Pass the generated Value to the problem

			'''
				result is a 3 column np array which is returned from the function from prob.py which carries the following information in its indexes
				result[0] - Validation Value
				result[1] - Value of the optimized result
				result[2] - tval passed to it
			'''
			if(result[0] == 1):  #Only Validated Result to be printed
				#print("\n",result[0])
				#print("\n",result[1])
				#print("\n",result[2])

				m = float(str(round(result[1], 3)))*1000
				if (cstr==0):
					route[0] = routeinfo(value=m, xval=result[2], ph=1)
					mins = minv(value=m, ps=cstr)
					cstr += 1
					if(gstr==0):
						min_res = iterinfo(value=m, xval=result[2])
						gstr +=1
				else:
					l1=0
					for lc in range(cstr):
						if(float(str(round(route[lc].value, 3)))*1000==m):
							l1+=1
							break
					if (l1==0):
						route[cstr]=routeinfo(value=m, xval=result[2], ph=1)

						if(route[mins.ps].value<route[cstr].value): #If the minimum Value of the iteration is smaller than Generated value then update the Minimum Value
							phn=route[mins.ps].ph*(1-phe)+1;
							route[mins.ps]=routeinfo(value=m, xval=result[2], ph=phn)
						if(route[mins.ps].value>route[cstr].value):
							route[cstr].ph=route[mins.ps].ph*(1-phe)+1
							route[mins.ps].ph=1;
							mins.ps=cstr;
							mins.value=route[cstr].value;
					cstr+=1
		#print("route[mins.ps].value",route[mins.ps].value);
		if(mins.value<min_res.value):
			min_res=iterinfo(value=mins.value, xval=route[mins.ps].xval)
	#print("min_res.value",min_res.value)
	#print("min_res.xval",min_res.xval)
	a = [min_res.value, min_res.xval]
	return a

def calc(func, cons , tval):

	val = func(tval)
	validate = int(1)
	cnstr=cons(tval)
	n_cnstr = cnstr.size;
	for i in range(n_cnstr):
		if (cnstr[i] < 0):
			validate = 0

	result = [validate, val, tval]
	return result
