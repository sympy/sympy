"""A module containing discrete functions.

Transforms - fft, ifft, ntt, intt, fwht, ifwht, fzt, ifzt, fmt, ifmt
Convolution - convolution, convolution_fft, convolution_ntt, convolution_fwht,
    convolution_subset, covering_product, intersecting_product
Recurrence - linrec
"""

from .transforms import (fft, ifft, ntt, intt, fwht, ifwht)
from .convolution import convolution
