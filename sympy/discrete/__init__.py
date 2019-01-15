"""This module contains functions which operate on discrete sequences.

Transforms - ``fft``, ``ifft``, ``ntt``, ``intt``, ``fwht``, ``ifwht``,
            ``mobius_transform``, ``inverse_mobius_transform``, ``dft_matrix``,
            ``ìdft_matrix``

Convolutions - ``convolution``, ``convolution_fft``, ``convolution_ntt``,
            ``convolution_fwht``, ``convolution_subset``,
            ``covering_product``, ``intersecting_product``
"""

from .transforms import (fft, ifft, ntt, intt, fwht, ifwht,
    mobius_transform, inverse_mobius_transform, dft_matrix, idft_matrix)
from .convolutions import convolution, covering_product, intersecting_product
