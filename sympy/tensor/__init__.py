#!/usr/bin/env python27
"""A module to manipulate symbolic objects with indices including tensors

"""
from .indexed import IndexedBase, Idx, Indexed
from .index_methods import get_contraction_structure, get_indices
from .arraypy import Arraypy, Tensor, matrix2arraypy, matrix2tensor,\
     list2arraypy, list2tensor
from .tensor_fields import df, grad, curl, diverg, lie_xy, dw, lie_w
from .riemannian_geometry import scal_prod, christoffel_1, christoffel_2,\
     covar_der, covar_der_XY, riemann, ricci, scal_curv, k_sigma
from .tensor_methods import symmetric, asymmetric, tensor_product, wedge,\
     perm_parity, change_basis, lower_index, raise_index, is_symmetric,\
     is_asymmetric
