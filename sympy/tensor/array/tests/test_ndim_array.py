from sympy import Array


def test_array_negative_indices():
	test_array = Array([[1,2,3,4,5],
                        [6,7,8,9,10]])
	assert test_array[:, -1] == Array([5, 10])
	assert test_array[:, -2] == Array([4, 9])
	assert test_array[:, -3] == Array([3, 8])
	assert test_array[:, -4] == Array([2, 7])
	assert test_array[:, -5] == Array([1, 6])
	assert test_array[:, 1] == Array([1, 6])
	assert test_array[:, 2] == Array([2, 7])
	assert test_array[:, 3] == Array([3, 8])
	assert test_array[:, 4] == Array([4, 9])
	assert test_array[:, 5] == Array([5, 10])
