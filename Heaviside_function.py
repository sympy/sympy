""" heaviside(x) function returns the value 0 for x < 0, 1 for x > 0, and 1/2 for x = 0. """
def heaviside(x):
    if (x < 0):
        return 0
    elif (x > 0):
        return 1
    else:
        return 1/2

"""print(heaviside(0))
   print(heaviside(6))
   print(heaviside(-9)) """