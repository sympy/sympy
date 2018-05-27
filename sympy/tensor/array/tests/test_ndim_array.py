import numpy as np
import sympy as sp

def test_array_boundry_checks():
	test_array = sp.Array(np.array([[1,2,3,4,5],[6,7,8,9,10]]))
	assert test_array[:,-1] == sp.Array([5,10])
	assert test_array[:,-2] == sp.Array([4,9])
	assert test_array[:,-3] == sp.Array([3,8])
	assert test_array[:,-4] == sp.Array([2,7])
	assert test_array[:,-5] == sp.Array([1,6])
	assert test_array[:,1] == sp.Array([1,6])
	assert test_array[:,2] == sp.Array([2,7])
	assert test_array[:,3] == sp.Array([3,8])
	assert test_array[:,4] == sp.Array([4,9])
	assert test_array[:,5] == sp.Array([5,10])

