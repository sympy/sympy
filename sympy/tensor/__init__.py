#!/usr/bin/env python27
"""A module to manipulate symbolic objects with indices including tensors

"""
from .indexed import IndexedBase, Idx, Indexed
from .index_methods import get_contraction_structure, get_indices
from .arraypy import Arraypy, TensorArray, matrix2arraypy, matrix2tensor,\
    list2arraypy, list2tensor
from .tensor_fields import Wedge_array, df, grad, curl, diverg, lie_xy, dw, \
    lie_w, tensor2wedgearray, wedgearray2tensor, inner_product, g_tensor, \
    g_wedge, hodge_star, codiff
from .helper_functions import check_vector_of_arguments, check_metric_tensor, \
    check_the_vector_field, sign_permutations, delete_index_from_list, \
    replace_index_to_k
from .riemannian_geometry import scal_prod, christoffel_1, christoffel_2,\
    covar_der, covar_der_xy, riemann, ricci, scal_curv, k_sigma, nabla,\
    nabla_x, delta, riemann_li, k_sigma_li
from .tensor_methods import symmetric, asymmetric, tensor_product, wedge,\
    perm_parity, change_basis, lower_index, raise_index, is_symmetric,\
    is_asymmetric
