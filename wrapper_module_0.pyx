import numpy as np
cimport numpy as np

cdef extern from 'wrapped_code_0.h':
    void test(double x, double y, double *z)

def test_c(double x, double y):

    cdef double z = 0
    test(x, y, &z)
    return z