import numpy as np
from scipy.interpolate import interp1d

a = np.linspace(1,10,10)
b = np.sin(a)
inter = interp1d(a,b)

s0 = inter(a[0])
print(a)
print(b)
print(s0)
