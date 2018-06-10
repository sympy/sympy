from sympy import Array


def test_array_negative_indices():
 test_array = Array([[1,2,3,4,5],[6,7,8,9,10]])
 assert test_array[:, -1] == Array([10, 5])
 assert test_array[:, -2] == Array([9, 4])
 assert test_array[:, -3] == Array([8, 3])
 assert test_array[:, -4] == Array([7, 2])
 assert test_array[:, -5] == Array([6, 1])
 assert test_array[:, 0] == Array([1, 6])
 assert test_array[:, 1] == Array([2, 7])
 assert test_array[:, 2] == Array([3, 8])
 assert test_array[:, 3] == Array([4, 9])
 assert test_array[:, 4] == Array([5, 10])
