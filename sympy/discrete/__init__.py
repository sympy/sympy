"""A module containing discrete functions.

Transforms - fft, ifft, ntt, intt, fwht, ifwht, fzt, ifzt, fmt, ifmt
Convolution - convolution, convolution_fft, convolution_ntt,
    convolution_xor, convolution_and, convolution_or, convolution_subset
Recurrence Evaluation - reval_hcc
"""


from .transforms import (fft, ifft, ntt, intt)
from .convolution import (convolution, convolution_fft, convolution_ntt)
