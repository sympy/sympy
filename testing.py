import numpy as np
from sympy.physics.wigner import wigner_3j

y = np.float64(6.0)
x = wigner_3j(2, y, 4.0, 0, 0, 0)
