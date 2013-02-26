"""Compatibility interface between dense and sparse polys. """

from sympy.polys.densearith import dup_add_term
from sympy.polys.densearith import dmp_add_term
from sympy.polys.densearith import dup_sub_term
from sympy.polys.densearith import dmp_sub_term
from sympy.polys.densearith import dup_mul_term
from sympy.polys.densearith import dmp_mul_term
from sympy.polys.densearith import dup_add_ground
from sympy.polys.densearith import dmp_add_ground
from sympy.polys.densearith import dup_sub_ground
from sympy.polys.densearith import dmp_sub_ground
from sympy.polys.densearith import dup_mul_ground
from sympy.polys.densearith import dmp_mul_ground
from sympy.polys.densearith import dup_quo_ground
from sympy.polys.densearith import dmp_quo_ground
from sympy.polys.densearith import dup_exquo_ground
from sympy.polys.densearith import dmp_exquo_ground
from sympy.polys.densearith import dup_lshift
from sympy.polys.densearith import dup_rshift
from sympy.polys.densearith import dup_abs
from sympy.polys.densearith import dmp_abs
from sympy.polys.densearith import dup_neg
from sympy.polys.densearith import dmp_neg
from sympy.polys.densearith import dup_add
from sympy.polys.densearith import dmp_add
from sympy.polys.densearith import dup_sub
from sympy.polys.densearith import dmp_sub
from sympy.polys.densearith import dup_add_mul
from sympy.polys.densearith import dmp_add_mul
from sympy.polys.densearith import dup_sub_mul
from sympy.polys.densearith import dmp_sub_mul
from sympy.polys.densearith import dup_mul
from sympy.polys.densearith import dmp_mul
from sympy.polys.densearith import dup_sqr
from sympy.polys.densearith import dmp_sqr
from sympy.polys.densearith import dup_pow
from sympy.polys.densearith import dmp_pow
from sympy.polys.densearith import dup_pdiv
from sympy.polys.densearith import dup_prem
from sympy.polys.densearith import dup_pquo
from sympy.polys.densearith import dup_pexquo
from sympy.polys.densearith import dmp_pdiv
from sympy.polys.densearith import dmp_prem
from sympy.polys.densearith import dmp_pquo
from sympy.polys.densearith import dmp_pexquo
from sympy.polys.densearith import dup_rr_div
from sympy.polys.densearith import dmp_rr_div
from sympy.polys.densearith import dup_ff_div
from sympy.polys.densearith import dmp_ff_div
from sympy.polys.densearith import dup_div
from sympy.polys.densearith import dup_rem
from sympy.polys.densearith import dup_quo
from sympy.polys.densearith import dup_exquo
from sympy.polys.densearith import dmp_div
from sympy.polys.densearith import dmp_rem
from sympy.polys.densearith import dmp_quo
from sympy.polys.densearith import dmp_exquo
from sympy.polys.densearith import dup_max_norm
from sympy.polys.densearith import dmp_max_norm
from sympy.polys.densearith import dup_l1_norm
from sympy.polys.densearith import dmp_l1_norm
from sympy.polys.densearith import dup_expand
from sympy.polys.densearith import dmp_expand
from sympy.polys.densebasic import dup_LC
from sympy.polys.densebasic import dmp_LC
from sympy.polys.densebasic import dup_TC
from sympy.polys.densebasic import dmp_TC
from sympy.polys.densebasic import dmp_ground_LC
from sympy.polys.densebasic import dmp_ground_TC
from sympy.polys.densebasic import dmp_true_LT
from sympy.polys.densebasic import dup_degree
from sympy.polys.densebasic import dmp_degree
from sympy.polys.densebasic import dmp_degree_in
from sympy.polys.densebasic import dmp_degree_list
from sympy.polys.densebasic import dup_strip
from sympy.polys.densebasic import dmp_strip
from sympy.polys.densebasic import dmp_validate
from sympy.polys.densebasic import dup_reverse
from sympy.polys.densebasic import dup_copy
from sympy.polys.densebasic import dmp_copy
from sympy.polys.densebasic import dup_to_tuple
from sympy.polys.densebasic import dmp_to_tuple
from sympy.polys.densebasic import dup_normal
from sympy.polys.densebasic import dmp_normal
from sympy.polys.densebasic import dup_convert
from sympy.polys.densebasic import dmp_convert
from sympy.polys.densebasic import dup_from_sympy
from sympy.polys.densebasic import dmp_from_sympy
from sympy.polys.densebasic import dup_nth
from sympy.polys.densebasic import dmp_nth
from sympy.polys.densebasic import dmp_ground_nth
from sympy.polys.densebasic import dmp_zero_p
from sympy.polys.densebasic import dmp_zero
from sympy.polys.densebasic import dmp_one_p
from sympy.polys.densebasic import dmp_one
from sympy.polys.densebasic import dmp_ground_p
from sympy.polys.densebasic import dmp_ground
from sympy.polys.densebasic import dmp_zeros
from sympy.polys.densebasic import dmp_grounds
from sympy.polys.densebasic import dmp_negative_p
from sympy.polys.densebasic import dmp_positive_p
from sympy.polys.densebasic import dup_from_dict
from sympy.polys.densebasic import dup_from_raw_dict
from sympy.polys.densebasic import dmp_from_dict
from sympy.polys.densebasic import dup_to_dict
from sympy.polys.densebasic import dup_to_raw_dict
from sympy.polys.densebasic import dmp_to_dict
from sympy.polys.densebasic import dmp_swap
from sympy.polys.densebasic import dmp_permute
from sympy.polys.densebasic import dmp_nest
from sympy.polys.densebasic import dmp_raise
from sympy.polys.densebasic import dup_deflate
from sympy.polys.densebasic import dmp_deflate
from sympy.polys.densebasic import dup_multi_deflate
from sympy.polys.densebasic import dmp_multi_deflate
from sympy.polys.densebasic import dup_inflate
from sympy.polys.densebasic import dmp_inflate
from sympy.polys.densebasic import dmp_exclude
from sympy.polys.densebasic import dmp_include
from sympy.polys.densebasic import dmp_inject
from sympy.polys.densebasic import dmp_eject
from sympy.polys.densebasic import dup_terms_gcd
from sympy.polys.densebasic import dmp_terms_gcd
from sympy.polys.densebasic import dmp_list_terms
from sympy.polys.densebasic import dup_apply_pairs
from sympy.polys.densebasic import dmp_apply_pairs
from sympy.polys.densebasic import dup_slice
from sympy.polys.densebasic import dmp_slice
from sympy.polys.densebasic import dmp_slice_in
from sympy.polys.densebasic import dup_random
from sympy.polys.densetools import dup_integrate
from sympy.polys.densetools import dmp_integrate
from sympy.polys.densetools import dmp_integrate_in
from sympy.polys.densetools import dup_diff
from sympy.polys.densetools import dmp_diff
from sympy.polys.densetools import dmp_diff_in
from sympy.polys.densetools import dup_eval
from sympy.polys.densetools import dmp_eval
from sympy.polys.densetools import dmp_eval_in
from sympy.polys.densetools import dmp_eval_tail
from sympy.polys.densetools import dmp_diff_eval_in
from sympy.polys.densetools import dup_trunc
from sympy.polys.densetools import dmp_trunc
from sympy.polys.densetools import dmp_ground_trunc
from sympy.polys.densetools import dup_monic
from sympy.polys.densetools import dmp_ground_monic
from sympy.polys.densetools import dup_content
from sympy.polys.densetools import dmp_ground_content
from sympy.polys.densetools import dup_primitive
from sympy.polys.densetools import dmp_ground_primitive
from sympy.polys.densetools import dup_extract
from sympy.polys.densetools import dmp_ground_extract
from sympy.polys.densetools import dup_real_imag
from sympy.polys.densetools import dup_mirror
from sympy.polys.densetools import dup_scale
from sympy.polys.densetools import dup_shift
from sympy.polys.densetools import dup_transform
from sympy.polys.densetools import dup_compose
from sympy.polys.densetools import dmp_compose
from sympy.polys.densetools import dup_decompose
from sympy.polys.densetools import dmp_lift
from sympy.polys.densetools import dup_sign_variations
from sympy.polys.densetools import dup_clear_denoms
from sympy.polys.densetools import dmp_clear_denoms
from sympy.polys.densetools import dup_revert
from sympy.polys.euclidtools import dup_half_gcdex
from sympy.polys.euclidtools import dmp_half_gcdex
from sympy.polys.euclidtools import dup_gcdex
from sympy.polys.euclidtools import dmp_gcdex
from sympy.polys.euclidtools import dup_invert
from sympy.polys.euclidtools import dmp_invert
from sympy.polys.euclidtools import dup_euclidean_prs
from sympy.polys.euclidtools import dmp_euclidean_prs
from sympy.polys.euclidtools import dup_primitive_prs
from sympy.polys.euclidtools import dmp_primitive_prs
from sympy.polys.euclidtools import dup_inner_subresultants
from sympy.polys.euclidtools import dup_subresultants
from sympy.polys.euclidtools import dup_prs_resultant
from sympy.polys.euclidtools import dup_resultant
from sympy.polys.euclidtools import dmp_inner_subresultants
from sympy.polys.euclidtools import dmp_subresultants
from sympy.polys.euclidtools import dmp_prs_resultant
from sympy.polys.euclidtools import dmp_zz_modular_resultant
from sympy.polys.euclidtools import dmp_zz_collins_resultant
from sympy.polys.euclidtools import dmp_qq_collins_resultant
from sympy.polys.euclidtools import dmp_resultant
from sympy.polys.euclidtools import dup_discriminant
from sympy.polys.euclidtools import dmp_discriminant
from sympy.polys.euclidtools import dup_rr_prs_gcd
from sympy.polys.euclidtools import dup_ff_prs_gcd
from sympy.polys.euclidtools import dmp_rr_prs_gcd
from sympy.polys.euclidtools import dmp_ff_prs_gcd
from sympy.polys.euclidtools import dup_zz_heu_gcd
from sympy.polys.euclidtools import dmp_zz_heu_gcd
from sympy.polys.euclidtools import dup_qq_heu_gcd
from sympy.polys.euclidtools import dmp_qq_heu_gcd
from sympy.polys.euclidtools import dup_inner_gcd
from sympy.polys.euclidtools import dmp_inner_gcd
from sympy.polys.euclidtools import dup_gcd
from sympy.polys.euclidtools import dmp_gcd
from sympy.polys.euclidtools import dup_rr_lcm
from sympy.polys.euclidtools import dup_ff_lcm
from sympy.polys.euclidtools import dup_lcm
from sympy.polys.euclidtools import dmp_rr_lcm
from sympy.polys.euclidtools import dmp_ff_lcm
from sympy.polys.euclidtools import dmp_lcm
from sympy.polys.euclidtools import dmp_content
from sympy.polys.euclidtools import dmp_primitive
from sympy.polys.euclidtools import dup_cancel
from sympy.polys.euclidtools import dmp_cancel
from sympy.polys.factortools import dup_trial_division
from sympy.polys.factortools import dmp_trial_division
from sympy.polys.factortools import dup_zz_mignotte_bound
from sympy.polys.factortools import dmp_zz_mignotte_bound
from sympy.polys.factortools import dup_zz_hensel_step
from sympy.polys.factortools import dup_zz_hensel_lift
from sympy.polys.factortools import dup_zz_zassenhaus
from sympy.polys.factortools import dup_zz_irreducible_p
from sympy.polys.factortools import dup_cyclotomic_p
from sympy.polys.factortools import dup_zz_cyclotomic_poly
from sympy.polys.factortools import dup_zz_cyclotomic_factor
from sympy.polys.factortools import dup_zz_factor_sqf
from sympy.polys.factortools import dup_zz_factor
from sympy.polys.factortools import dmp_zz_wang_non_divisors
from sympy.polys.factortools import dmp_zz_wang_test_points
from sympy.polys.factortools import dmp_zz_wang_lead_coeffs
from sympy.polys.factortools import dup_zz_diophantine
from sympy.polys.factortools import dmp_zz_diophantine
from sympy.polys.factortools import dmp_zz_wang_hensel_lifting
from sympy.polys.factortools import dmp_zz_wang
from sympy.polys.factortools import dmp_zz_factor
from sympy.polys.factortools import dup_ext_factor
from sympy.polys.factortools import dmp_ext_factor
from sympy.polys.factortools import dup_gf_factor
from sympy.polys.factortools import dmp_gf_factor
from sympy.polys.factortools import dup_factor_list
from sympy.polys.factortools import dup_factor_list_include
from sympy.polys.factortools import dmp_factor_list
from sympy.polys.factortools import dmp_factor_list_include
from sympy.polys.factortools import dup_irreducible_p
from sympy.polys.factortools import dmp_irreducible_p
from sympy.polys.orthopolys import dup_jacobi
from sympy.polys.orthopolys import dup_gegenbauer
from sympy.polys.orthopolys import dup_chebyshevt
from sympy.polys.orthopolys import dup_chebyshevu
from sympy.polys.orthopolys import dup_hermite
from sympy.polys.orthopolys import dup_legendre
from sympy.polys.orthopolys import dup_laguerre
from sympy.polys.orthopolys import dup_spherical_bessel_fn
from sympy.polys.orthopolys import dup_spherical_bessel_fn_minus
from sympy.polys.rootisolation import dup_sturm
from sympy.polys.rootisolation import dup_root_upper_bound
from sympy.polys.rootisolation import dup_root_lower_bound
from sympy.polys.rootisolation import dup_step_refine_real_root
from sympy.polys.rootisolation import dup_inner_refine_real_root
from sympy.polys.rootisolation import dup_outer_refine_real_root
from sympy.polys.rootisolation import dup_refine_real_root
from sympy.polys.rootisolation import dup_inner_isolate_real_roots
from sympy.polys.rootisolation import dup_inner_isolate_positive_roots
from sympy.polys.rootisolation import dup_inner_isolate_negative_roots
from sympy.polys.rootisolation import dup_isolate_real_roots_sqf
from sympy.polys.rootisolation import dup_isolate_real_roots
from sympy.polys.rootisolation import dup_isolate_real_roots_list
from sympy.polys.rootisolation import dup_count_real_roots
from sympy.polys.rootisolation import dup_count_complex_roots
from sympy.polys.rootisolation import dup_isolate_complex_roots_sqf
from sympy.polys.rootisolation import dup_isolate_all_roots_sqf
from sympy.polys.rootisolation import dup_isolate_all_roots
from sympy.polys.specialpolys import dmp_fateman_poly_F_1
from sympy.polys.specialpolys import dmp_fateman_poly_F_2
from sympy.polys.specialpolys import dmp_fateman_poly_F_3
from sympy.polys.sqfreetools import dup_sqf_p
from sympy.polys.sqfreetools import dmp_sqf_p
from sympy.polys.sqfreetools import dup_sqf_norm
from sympy.polys.sqfreetools import dmp_sqf_norm
from sympy.polys.sqfreetools import dup_gf_sqf_part
from sympy.polys.sqfreetools import dmp_gf_sqf_part
from sympy.polys.sqfreetools import dup_sqf_part
from sympy.polys.sqfreetools import dmp_sqf_part
from sympy.polys.sqfreetools import dup_gf_sqf_list
from sympy.polys.sqfreetools import dmp_gf_sqf_list
from sympy.polys.sqfreetools import dup_sqf_list
from sympy.polys.sqfreetools import dup_sqf_list_include
from sympy.polys.sqfreetools import dmp_sqf_list
from sympy.polys.sqfreetools import dmp_sqf_list_include
from sympy.polys.sqfreetools import dup_gff_list
from sympy.polys.sqfreetools import dmp_gff_list

class IPolys(object):
    def wrap(self, element):
        from sympy.polys.rings import PolyElement
        if isinstance(element, PolyElement):
            if element.ring == self:
                return element
            else:
                raise NotImplementedError("domain conversions")
        else:
            return self.ground_new(element)

    def to_dense(self, element):
        return self.wrap(element).to_dense()

    def from_dense(self, element):
        return self.from_dict(dmp_to_dict(element, self.ngens-1, self.domain))

    def dup_add_term(self, f, c, i):
        return self.from_dense(dup_add_term(self.to_dense(f), c, i, self.domain))
    def dmp_add_term(self, f, c, i):
        return self.from_dense(dmp_add_term(self.to_dense(f), self.wrap(c).drop(0).to_dense(), i, self.ngens-1, self.domain))
    def dup_sub_term(self, f, c, i):
        return self.from_dense(dup_sub_term(self.to_dense(f), c, i, self.domain))
    def dmp_sub_term(self, f, c, i):
        return self.from_dense(dmp_sub_term(self.to_dense(f), self.wrap(c).drop(0).to_dense(), i, self.ngens-1, self.domain))
    def dup_mul_term(self, f, c, i):
        return self.from_dense(dup_mul_term(self.to_dense(f), c, i, self.domain))
    def dmp_mul_term(self, f, c, i):
        return self.from_dense(dmp_mul_term(self.to_dense(f), self.wrap(c).drop(0).to_dense(), i, self.ngens-1, self.domain))

    def dup_add_ground(self, f, c):
        return self.from_dense(dup_add_ground(self.to_dense(f), c, self.domain))
    def dmp_add_ground(self, f, c):
        return self.from_dense(dmp_add_ground(self.to_dense(f), c, self.ngens-1, self.domain))
    def dup_sub_ground(self, f, c):
        return self.from_dense(dup_sub_ground(self.to_dense(f), c, self.domain))
    def dmp_sub_ground(self, f, c):
        return self.from_dense(dmp_sub_ground(self.to_dense(f), c, self.ngens-1, self.domain))
    def dup_mul_ground(self, f, c):
        return self.from_dense(dup_mul_ground(self.to_dense(f), c, self.domain))
    def dmp_mul_ground(self, f, c):
        return self.from_dense(dmp_mul_ground(self.to_dense(f), c, self.ngens-1, self.domain))
    def dup_quo_ground(self, f, c):
        return self.from_dense(dup_quo_ground(self.to_dense(f), c, self.domain))
    def dmp_quo_ground(self, f, c):
        return self.from_dense(dmp_quo_ground(self.to_dense(f), c, self.ngens-1, self.domain))
    def dup_exquo_ground(self, f, c):
        return self.from_dense(dup_exquo_ground(self.to_dense(f), c, self.domain))
    def dmp_exquo_ground(self, f, c):
        return self.from_dense(dmp_exquo_ground(self.to_dense(f), c, self.ngens-1, self.domain))

    def dup_lshift(self, f, n):
        return self.from_dense(dup_lshift(self.to_dense(f), n, self.domain))
    def dup_rshift(self, f, n):
        return self.from_dense(dup_rshift(self.to_dense(f), n, self.domain))

    def dup_abs(self, f):
        return self.from_dense(dup_abs(self.to_dense(f), self.domain))
    def dmp_abs(self, f):
        return self.from_dense(dmp_abs(self.to_dense(f), self.ngens-1, self.domain))

    def dup_neg(self, f):
        return self.from_dense(dup_neg(self.to_dense(f), self.domain))
    def dmp_neg(self, f):
        return self.from_dense(dmp_neg(self.to_dense(f), self.ngens-1, self.domain))

    def dup_add(self, f, g):
        return self.from_dense(dup_add(self.to_dense(f), self.to_dense(g), self.domain))
    def dmp_add(self, f, g):
        return self.from_dense(dmp_add(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))

    def dup_sub(self, f, g):
        return self.from_dense(dup_sub(self.to_dense(f), self.to_dense(g), self.domain))
    def dmp_sub(self, f, g):
        return self.from_dense(dmp_sub(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))

    def dup_add_mul(self, f, g, h):
        return self.from_dense(dup_add_mul(self.to_dense(f), self.to_dense(g), self.to_dense(h), self.domain))
    def dmp_add_mul(self, f, g, h):
        return self.from_dense(dmp_add_mul(self.to_dense(f), self.to_dense(g), self.to_dense(h), self.ngens-1, self.domain))
    def dup_sub_mul(self, f, g, h):
        return self.from_dense(dup_sub_mul(self.to_dense(f), self.to_dense(g), self.to_dense(h), self.domain))
    def dmp_sub_mul(self, f, g, h):
        return self.from_dense(dmp_sub_mul(self.to_dense(f), self.to_dense(g), self.to_dense(h), self.ngens-1, self.domain))

    def dup_mul(self, f, g):
        return self.from_dense(dup_mul(self.to_dense(f), self.to_dense(g), self.domain))
    def dmp_mul(self, f, g):
        return self.from_dense(dmp_mul(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))

    def dup_sqr(self, f):
        return self.from_dense(dup_sqr(self.to_dense(f), self.domain))
    def dmp_sqr(self, f):
        return self.from_dense(dmp_sqr(self.to_dense(f), self.ngens-1, self.domain))
    def dup_pow(self, f, n):
        return self.from_dense(dup_pow(self.to_dense(f), n, self.domain))
    def dmp_pow(self, f, n):
        return self.from_dense(dmp_pow(self.to_dense(f), n, self.ngens-1, self.domain))

    def dup_pdiv(self, f, g):
        q, r = dup_pdiv(self.to_dense(f), self.to_dense(g), self.domain)
        return (self.from_dense(q), self.from_dense(r))
    def dup_prem(self, f, g):
        return self.from_dense(dup_prem(self.to_dense(f), self.to_dense(g), self.domain))
    def dup_pquo(self, f, g):
        return self.from_dense(dup_pquo(self.to_dense(f), self.to_dense(g), self.domain))
    def dup_pexquo(self, f, g):
        return self.from_dense(dup_pexquo(self.to_dense(f), self.to_dense(g), self.domain))

    def dmp_pdiv(self, f, g):
        q, r = dmp_pdiv(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return (self.from_dense(q), self.from_dense(r))
    def dmp_prem(self, f, g):
        return self.from_dense(dmp_prem(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))
    def dmp_pquo(self, f, g):
        return self.from_dense(dmp_pquo(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))
    def dmp_pexquo(self, f, g):
        return self.from_dense(dmp_pexquo(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))

    def dup_rr_div(self, f, g):
        q, r = dup_rr_div(self.to_dense(f), self.to_dense(g), self.domain)
        return (self.from_dense(q), self.from_dense(r))
    def dmp_rr_div(self, f, g):
        q, r = dmp_rr_div(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return (self.from_dense(q), self.from_dense(r))
    def dup_ff_div(self, f, g):
        q, r = dup_ff_div(self.to_dense(f), self.to_dense(g), self.domain)
        return (self.from_dense(q), self.from_dense(r))
    def dmp_ff_div(self, f, g):
        q, r = dmp_ff_div(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return (self.from_dense(q), self.from_dense(r))

    def dup_div(self, f, g):
        q, r = dup_div(self.to_dense(f), self.to_dense(g), self.domain)
        return (self.from_dense(q), self.from_dense(r))
    def dup_rem(self, f, g):
        return self.from_dense(dup_rem(self.to_dense(f), self.to_dense(g), self.domain))
    def dup_quo(self, f, g):
        return self.from_dense(dup_quo(self.to_dense(f), self.to_dense(g), self.domain))
    def dup_exquo(self, f, g):
        return self.from_dense(dup_exquo(self.to_dense(f), self.to_dense(g), self.domain))

    def dmp_div(self, f, g):
        q, r = dmp_div(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return (self.from_dense(q), self.from_dense(r))
    def dmp_rem(self, f, g):
        return self.from_dense(dmp_rem(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))
    def dmp_quo(self, f, g):
        return self.from_dense(dmp_quo(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))
    def dmp_exquo(self, f, g):
        return self.from_dense(dmp_exquo(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))

    def dup_max_norm(self, f):
        return dup_max_norm(self.to_dense(f), self.domain)
    def dmp_max_norm(self, f):
        return dmp_max_norm(self.to_dense(f), self.ngens-1, self.domain)

    def dup_l1_norm(self, f):
        return dup_l1_norm(self.to_dense(f), self.domain)
    def dmp_l1_norm(self, f):
        return dmp_l1_norm(self.to_dense(f), self.ngens-1, self.domain)

    def dup_expand(self, polys):
        return self.from_dense(dup_expand(map(self.to_dense, polys), self.domain))
    def dmp_expand(self, polys):
        return self.from_dense(dmp_expand(map(self.to_dense, polys), self.ngens-1, self.domain))

    def dup_LC(self, f):
        return dup_LC(self.to_dense(f), self.domain)
    def dmp_LC(self, f):
        LC = dmp_LC(self.to_dense(f), self.domain)
        if isinstance(LC, list):
            return self.drop(0).from_dense(LC)
        else:
            return LC
    def dup_TC(self, f):
        return dup_TC(self.to_dense(f), self.domain)
    def dmp_TC(self, f):
        TC = dmp_TC(self.to_dense(f), self.domain)
        if isinstance(TC, list):
            return self.drop(0).from_dense(TC)
        else:
            return TC

    def dmp_ground_LC(self, f):
        return dmp_ground_LC(self.to_dense(f), self.ngens-1, self.domain)
    def dmp_ground_TC(self, f):
        return dmp_ground_TC(self.to_dense(f), self.ngens-1, self.domain)

    def dup_degree(self, f):
        return dup_degree(self.to_dense(f))
    def dmp_degree(self, f):
        return dmp_degree(self.to_dense(f), self.ngens-1)
    def dmp_degree_in(self, f, j):
        return dmp_degree_in(self.to_dense(f), j, self.ngens-1)
    """
    def dmp_degree_list(f, u):
    def dup_strip(f):
    def dmp_strip(f, u):
    def dmp_validate(f, K=None):
    def dup_reverse(f):
    def dup_copy(f):
    def dmp_copy(f, u):
    def dup_to_tuple(f):
    def dmp_to_tuple(f, u):
    def dup_normal(f, K):
    def dmp_normal(f):
    def dup_convert(f, K0, K1):
    def dmp_convert(f, u, K0, K1):
    def dup_from_sympy(f, K):
    def dmp_from_sympy(f, u, K):
    def dup_nth(f, n, K):
    def dmp_nth(f, n, u, K):
    def dmp_ground_nth(f, N, u, K):
    def dmp_zero_p(f, u):
    def dmp_zero(u):
    def dmp_one_p(f, u, K):
    def dmp_one(u, K):
    def dmp_ground_p(f, c, u):
    def dmp_ground(c, u):
    def dmp_zeros(n, u, K):
    def dmp_grounds(c, n, u):
    def dmp_negative_p(f, u, K):
    def dmp_positive_p(f, u, K):
    def dup_from_dict(f, K):
    def dup_from_raw_dict(f, K):
    def dmp_from_dict(f, u, K):
    def dup_to_dict(f, K=None, zero=False):
    def dup_to_raw_dict(f, K=None, zero=False):
    def dmp_to_dict(f, u, K=None, zero=False):
    def dmp_swap(f, i, j, u, K):
    def dmp_permute(f, P, u, K):
    def dmp_nest(f, l, K):
    def dmp_raise(f, l, u, K):
    def dup_deflate(f, K):
    def dmp_deflate(f, u, K):
    def dup_multi_deflate(polys, K):
    def dmp_multi_deflate(polys, u, K):
    def dup_inflate(f, m, K):
    def dmp_inflate(f, M, u, K):
    def dmp_exclude(f, u, K):
    def dmp_include(f, J, u, K):
    def dmp_inject(f, u, K, front=False):
    def dmp_eject(f, u, K, front=False):
    def dup_terms_gcd(f, K):
    def dmp_terms_gcd(f, u, K):
    def dmp_list_terms(f, u, K, order=None):
    def dup_apply_pairs(f, g, h, args, K):
    def dmp_apply_pairs(f, g, h, args, u, K):
    def dup_slice(f, m, n, K):
    def dmp_slice(f, m, n, u, K):
    def dmp_slice_in(f, m, n, j, u, K):
    def dup_random(n, a, b, K):
    """

    def dup_integrate(self, f, m):
        return self.from_dense(dup_integrate(self.to_dense(f), m, self.domain))
    def dmp_integrate(self, f, m):
        return self.from_dense(dmp_integrate(self.to_dense(f), m, self.ngens-1, self.domain))


    def dup_diff(self, f, m):
        return self.from_dense(dup_diff(self.to_dense(f), m, self.domain))
    def dmp_diff(self, f, m):
        return self.from_dense(dmp_diff(self.to_dense(f), m, self.ngens-1, self.domain))

    def dmp_diff_in(self, f, m, j):
        return self.from_dense(dmp_diff_in(self.to_dense(f), m, j, self.ngens-1, self.domain))
    def dmp_integrate_in(self, f, m, j):
        return self.from_dense(dmp_integrate_in(self.to_dense(f), m, j, self.ngens-1, self.domain))

    def dup_eval(self, f, a):
        return dup_eval(self.to_dense(f), a, self.domain)
    def dmp_eval(self, f, a):
        result = dmp_eval(self.to_dense(f), a, self.ngens-1, self.domain)
        return self.drop(0).from_dense(result)

    def dmp_eval_in(self, f, a, j):
        result = dmp_eval_in(self.to_dense(f), a, j, self.ngens-1, self.domain)
        return self.drop(j).from_dense(result)
    def dmp_diff_eval_in(self, f, m, a, j):
        result = dmp_diff_eval_in(self.to_dense(f), m, a, j, self.ngens-1, self.domain)
        return self.drop(j).from_dense(result)

    def dmp_eval_tail(self, f, A):
        result = dmp_eval_tail(self.to_dense(f), A, self.ngens-1, self.domain)
        if isinstance(result, list):
            return self[:-len(A)].from_dense(result)
        else:
            return result

    def dup_trunc(self, f, p):
        return self.from_dense(dup_trunc(self.to_dense(f), p, self.domain))
    def dmp_trunc(self, f, g):
        return self.from_dense(dmp_trunc(self.to_dense(f), self.drop(0).to_dense(g), self.ngens-1, self.domain))
    def dmp_ground_trunc(self, f, p):
        return self.from_dense(dmp_ground_trunc(self.to_dense(f), p, self.ngens-1, self.domain))

    def dup_monic(self, f):
        return self.from_dense(dup_monic(self.to_dense(f), self.domain))
    def dmp_ground_monic(self, f):
        return self.from_dense(dmp_ground_monic(self.to_dense(f), self.ngens-1, self.domain))

    def dup_extract(self, f, g):
        c, F, G = dup_extract(self.to_dense(f), self.to_dense(g), self.domain)
        return (c, self.from_dense(F), self.from_dense(G))
    def dmp_ground_extract(self, f, g):
        c, F, G = dmp_ground_extract(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return (c, self.from_dense(F), self.from_dense(G))

    def dup_real_imag(self, f):
        p, q = dup_real_imag(self.wrap(f).drop(1).to_dense(), self.domain)
        return (self.from_dense(p), self.from_dense(q))

    def dup_mirror(self, f):
        return self.from_dense(dup_mirror(self.to_dense(f), self.domain))
    def dup_scale(self, f, a):
        return self.from_dense(dup_scale(self.to_dense(f), a, self.domain))
    def dup_shift(self, f, a):
        return self.from_dense(dup_shift(self.to_dense(f), a, self.domain))
    def dup_transform(self, f, p, q):
        return self.from_dense(dup_transform(self.to_dense(f), self.to_dense(p), self.to_dense(q), self.domain))

    def dup_compose(self, f, g):
        return self.from_dense(dup_compose(self.to_dense(f), self.to_dense(g), self.domain))
    def dmp_compose(self, f, g):
        return self.from_dense(dmp_compose(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))

    def dup_decompose(self, f):
        components = dup_decompose(self.to_dense(f), self.domain)
        return map(self.from_dense, components)

    def dmp_lift(self, f):
        result = dmp_lift(self.to_dense(f), self.ngens-1, self.domain)
        return self.to_ground().from_dense(result)

    def dup_sign_variations(self, f):
        return dup_sign_variations(self.to_dense(f), self.domain)

    def dup_clear_denoms(self, f, convert=False):
        c, F = dup_clear_denoms(self.to_dense(f), self.domain, convert=convert)
        if convert: ring = self.clone(domain=self.domain.get_ring())
        else: ring = self
        return (c, ring.from_dense(F))
    def dmp_clear_denoms(self, f, convert=False):
        c, F = dmp_clear_denoms(self.to_dense(f), self.ngens-1, self.domain, convert=convert)
        if convert: ring = self.clone(domain=self.domain.get_ring())
        else: ring = self
        return (c, ring.from_dense(F))

    def dup_revert(self, f, n):
        return self.from_dense(dup_revert(self.to_dense(f), n, self.domain))

    def dup_half_gcdex(self, f, g):
        s, h = dup_half_gcdex(self.to_dense(f), self.to_dense(g), self.domain)
        return (self.from_dense(s), self.from_dense(h))
    def dmp_half_gcdex(self, f, g):
        s, h = dmp_half_gcdex(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return (self.from_dense(s), self.from_dense(h))
    def dup_gcdex(self, f, g):
        s, t, h = dup_gcdex(self.to_dense(f), self.to_dense(g), self.domain)
        return (self.from_dense(s), self.from_dense(t), self.from_dense(h))
    def dmp_gcdex(self, f, g):
        s, t, h = dmp_gcdex(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return (self.from_dense(s), self.from_dense(t), self.from_dense(h))

    def dup_invert(self, f, g):
        return self.from_dense(dup_invert(self.to_dense(f), self.to_dense(g), self.domain))
    def dmp_invert(self, f, g):
        return self.from_dense(dmp_invert(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))

    def dup_euclidean_prs(self, f, g):
        prs = dup_euclidean_prs(self.to_dense(f), self.to_dense(g), self.domain)
        return map(self.from_dense, prs)
    def dmp_euclidean_prs(self, f, g):
        prs = dmp_euclidean_prs(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return map(self.from_dense, prs)
    def dup_primitive_prs(self, f, g):
        prs = dup_primitive_prs(self.to_dense(f), self.to_dense(g), self.domain)
        return map(self.from_dense, prs)
    def dmp_primitive_prs(self, f, g):
        prs = dmp_primitive_prs(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return map(self.from_dense, prs)

    def dup_inner_subresultants(self, f, g):
        prs, beta, delta = dup_inner_subresultants(self.to_dense(f), self.to_dense(g), self.domain)
        return (map(self.from_dense, prs), beta, delta)
    def dmp_inner_subresultants(self, f, g):
        prs, beta, delta = dmp_inner_subresultants(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return (map(self.from_dense, prs), beta, delta)

    def dup_subresultants(self, f, g):
        prs = dup_subresultants(self.to_dense(f), self.to_dense(g), self.domain)
        return map(self.from_dense, prs)
    def dmp_subresultants(self, f, g):
        prs = dmp_subresultants(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return map(self.from_dense, prs)

    def dup_prs_resultant(self, f, g):
        res, prs = dup_prs_resultant(self.to_dense(f), self.to_dense(g), self.domain)
        return (res, map(self.from_dense, prs))
    def dmp_prs_resultant(self, f, g):
        res, prs = dmp_prs_resultant(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return (self.drop(0).from_dense(res), map(self.from_dense, prs))

    def dmp_zz_modular_resultant(self, f, g, p):
        res = dmp_zz_modular_resultant(self.to_dense(f), self.to_dense(g), self.domain_new(p), self.ngens-1, self.domain)
        return self.drop(0).from_dense(res)
    def dmp_zz_collins_resultant(self, f, g):
        res = dmp_zz_collins_resultant(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.drop(0).from_dense(res)
    def dmp_qq_collins_resultant(self, f, g):
        res = dmp_qq_collins_resultant(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.drop(0).from_dense(res)

    def dup_resultant(self, f, g): #, includePRS=False):
        return dup_resultant(self.to_dense(f), self.to_dense(g), self.domain) #, includePRS=includePRS)
    def dmp_resultant(self, f, g): #, includePRS=False):
        res = dmp_resultant(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain) #, includePRS=includePRS)
        if isinstance(res, list):
            return self.drop(0).from_dense(res)
        else:
            return res

    def dup_discriminant(self, f):
        return dup_discriminant(self.to_dense(f), self.domain)
    def dmp_discriminant(self, f):
        disc = dmp_discriminant(self.to_dense(f), self.ngens-1, self.domain)
        if isinstance(disc, list):
            return self.drop(0).from_dense(disc)
        else:
            return disc

    def dup_rr_prs_gcd(self, f, g):
        H, F, G = dup_rr_prs_gcd(self.to_dense(f), self.to_dense(g), self.domain)
        return (self.from_dense(H), self.from_dense(F), self.from_dense(G))
    def dup_ff_prs_gcd(self, f, g):
        H, F, G = dup_ff_prs_gcd(self.to_dense(f), self.to_dense(g), self.domain)
        return (self.from_dense(H), self.from_dense(F), self.from_dense(G))
    def dmp_rr_prs_gcd(self, f, g):
        H, F, G = dmp_rr_prs_gcd(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return (self.from_dense(H), self.from_dense(F), self.from_dense(G))
    def dmp_ff_prs_gcd(self, f, g):
        H, F, G = dmp_ff_prs_gcd(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return (self.from_dense(H), self.from_dense(F), self.from_dense(G))
    def dup_zz_heu_gcd(self, f, g):
        H, F, G = dup_zz_heu_gcd(self.to_dense(f), self.to_dense(g), self.domain)
        return (self.from_dense(H), self.from_dense(F), self.from_dense(G))
    def dmp_zz_heu_gcd(self, f, g):
        H, F, G = dmp_zz_heu_gcd(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return (self.from_dense(H), self.from_dense(F), self.from_dense(G))
    def dup_qq_heu_gcd(self, f, g):
        H, F, G = dup_qq_heu_gcd(self.to_dense(f), self.to_dense(g), self.domain)
        return (self.from_dense(H), self.from_dense(F), self.from_dense(G))
    def dmp_qq_heu_gcd(self, f, g):
        H, F, G = dmp_qq_heu_gcd(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return (self.from_dense(H), self.from_dense(F), self.from_dense(G))
    def dup_inner_gcd(self, f, g):
        H, F, G = dup_inner_gcd(self.to_dense(f), self.to_dense(g), self.domain)
        return (self.from_dense(H), self.from_dense(F), self.from_dense(G))
    def dmp_inner_gcd(self, f, g):
        H, F, G = dmp_inner_gcd(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return (self.from_dense(H), self.from_dense(F), self.from_dense(G))
    def dup_gcd(self, f, g):
        H = dup_gcd(self.to_dense(f), self.to_dense(g), self.domain)
        return self.from_dense(H)
    def dmp_gcd(self, f, g):
        H = dmp_gcd(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.from_dense(H)
    def dup_rr_lcm(self, f, g):
        H = dup_rr_lcm(self.to_dense(f), self.to_dense(g), self.domain)
        return self.from_dense(H)
    def dup_ff_lcm(self, f, g):
        H = dup_ff_lcm(self.to_dense(f), self.to_dense(g), self.domain)
        return self.from_dense(H)
    def dup_lcm(self, f, g):
        H = dup_lcm(self.to_dense(f), self.to_dense(g), self.domain)
        return self.from_dense(H)
    def dmp_rr_lcm(self, f, g):
        H = dmp_rr_lcm(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.from_dense(H)
    def dmp_ff_lcm(self, f, g):
        H = dmp_ff_lcm(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.from_dense(H)
    def dmp_lcm(self, f, g):
        H = dmp_lcm(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.from_dense(H)

    def dup_content(self, f):
        cont = dup_content(self.to_dense(f), self.domain)
        return cont
    def dup_primitive(self, f):
        cont, prim = dup_primitive(self.to_dense(f), self.domain)
        return cont, self.from_dense(prim)

    def dmp_content(self, f):
        cont = dmp_content(self.to_dense(f), self.ngens-1, self.domain)
        if isinstance(cont, list):
            return self.drop(0).from_dense(cont)
        else:
            return cont
    def dmp_primitive(self, f):
        cont, prim = dmp_primitive(self.to_dense(f), self.ngens-1, self.domain)
        if isinstance(cont, list):
            return (self.drop(0).from_dense(cont), self.from_dense(prim))
        else:
            return (cont, self.from_dense(prim))

    def dmp_ground_content(self, f):
        cont = dmp_ground_content(self.to_dense(f), self.ngens-1, self.domain)
        return cont
    def dmp_ground_primitive(self, f):
        cont, prim = dmp_ground_primitive(self.to_dense(f), self.ngens-1, self.domain)
        return (cont, self.from_dense(prim))

    def dup_cancel(self, f, g, include=True):
        result = dup_cancel(self.to_dense(f), self.to_dense(g), self.domain, include=include)
        if not include:
            cf, cg, F, G = result
            return (cf, cg, self.from_dense(F), self.from_dense(G))
        else:
            F, G = result
            return (self.from_dense(F), self.from_dense(G))
    def dmp_cancel(self, f, g, include=True):
        result = dmp_cancel(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain, include=include)
        if not include:
            cf, cg, F, G = result
            return (cf, cg, self.from_dense(F), self.from_dense(G))
        else:
            F, G = result
            return (self.from_dense(F), self.from_dense(G))

    def dup_trial_division(self, f, factors):
        factors = dup_trial_division(self.to_dense(f), map(self.to_dense, factors), self.domain)
        return [ (self.from_dense(g), k) for g, k in factors ]
    def dmp_trial_division(self, f, factors):
        factors = dmp_trial_division(self.to_dense(f), map(self.to_dense, factors), self.ngens-1, self.domain)
        return [ (self.from_dense(g), k) for g, k in factors ]

    def dup_zz_mignotte_bound(self, f):
        return dup_zz_mignotte_bound(self.to_dense(f), self.domain)
    def dmp_zz_mignotte_bound(self, f):
        return dmp_zz_mignotte_bound(self.to_dense(f), self.ngens-1, self.domain)

    def dup_zz_hensel_step(self, m, f, g, h, s, t):
        D = self.to_dense
        G, H, S, T = dup_zz_hensel_step(m, D(f), D(g), D(h), D(s), D(t), self.domain)
        return (self.from_dense(G), self.from_dense(H), self.from_dense(S), self.from_dense(T))
    def dup_zz_hensel_lift(self, p, f, f_list, l):
        D = self.to_dense
        polys = dup_zz_hensel_lift(p, D(f), map(D, f_list), l, self.domain)
        return map(self.from_dense, polys)

    def dup_zz_zassenhaus(self, f):
        factors = dup_zz_zassenhaus(self.to_dense(f), self.domain)
        return [ (self.from_dense(g), k) for g, k in factors ]

    def dup_zz_irreducible_p(self, f):
        return dup_zz_irreducible_p(self.to_dense(f), self.domain)
    def dup_cyclotomic_p(self, f, irreducible=False):
        return dup_cyclotomic_p(self.to_dense(f), self.domain, irreducible=irreducible)
    def dup_zz_cyclotomic_poly(self, n):
        F = dup_zz_cyclotomic_poly(n, self.domain)
        return self.from_dense(F)
    def dup_zz_cyclotomic_factor(self, f):
        result = dup_zz_cyclotomic_factor(self.to_dense(f), self.domain)
        if result is None:
            return result
        else:
            return map(self.from_dense, result)

    # E: List[ZZ], cs: ZZ, ct: ZZ
    def dmp_zz_wang_non_divisors(self, E, cs, ct):
        return dmp_zz_wang_non_divisors(E, cs, ct, self.domain)

    # f: Poly, T: List[(Poly, int)], ct: ZZ, A: List[ZZ]
    #def dmp_zz_wang_test_points(f, T, ct, A):
    #   dmp_zz_wang_test_points(self.to_dense(f), T, ct, A, self.ngens-1, self.domain)

    # f: Poly, T: List[(Poly, int)], cs: ZZ, E: List[ZZ], H: List[Poly], A: List[ZZ]
    def dmp_zz_wang_lead_coeffs(self, f, T, cs, E, H, A):
        mv = self[1:]
        T = [ (mv.to_dense(t), k) for t, k in T ]
        uv = self[:1]
        H = map(uv.to_dense, H)
        f, HH, CC = dmp_zz_wang_lead_coeffs(self.to_dense(f), T, cs, E, H, A, self.ngens-1, self.domain)
        return self.from_dense(f), map(uv.from_dense, HH), map(mv.from_dense, CC)

    # f: List[Poly], m: int, p: ZZ
    def dup_zz_diophantine(self, F, m, p):
        result = dup_zz_diophantine(map(self.to_dense, F), m, p, self.domain)
        return map(self.from_dense, result)

    # f: List[Poly], c: List[Poly], A: List[ZZ], d: int, p: ZZ
    def dmp_zz_diophantine(self, F, c, A, d, p):
        result = dmp_zz_diophantine(map(self.to_dense, F), self.to_dense(c), A, d, p, self.ngens-1, self.domain)
        return map(self.from_dense, result)

    # f: Poly, H: List[Poly], LC: List[Poly], A: List[ZZ], p: ZZ
    def dmp_zz_wang_hensel_lifting(self, f, H, LC, A, p):
        uv = self[:1]
        mv = self[1:]
        H = map(uv.to_dense, H)
        LC = map(mv.to_dense, LC)
        result = dmp_zz_wang_hensel_lifting(self.to_dense(f), H, LC, A, p, self.ngens-1, self.domain)
        return map(self.from_dense, result)

    def dmp_zz_wang(self, f, mod=None, seed=None):
        factors = dmp_zz_wang(self.to_dense(f), self.ngens-1, self.domain, mod=mod, seed=seed)
        return [ self.from_dense(g) for g in factors ]

    def dup_zz_factor_sqf(self, f):
        coeff, factors = dup_zz_factor_sqf(self.to_dense(f), self.domain)
        return (coeff, [ self.from_dense(g) for g in factors ])

    def dup_zz_factor(self, f):
        coeff, factors = dup_zz_factor(self.to_dense(f), self.domain)
        return (coeff, [ (self.from_dense(g), k) for g, k in factors ])
    def dmp_zz_factor(self, f):
        coeff, factors = dmp_zz_factor(self.to_dense(f), self.ngens-1, self.domain)
        return (coeff, [ (self.from_dense(g), k) for g, k in factors ])

    def dup_ext_factor(self, f):
        coeff, factors = dup_ext_factor(self.to_dense(f), self.domain)
        return (coeff, [ (self.from_dense(g), k) for g, k in factors ])
    def dmp_ext_factor(self, f):
        coeff, factors = dmp_ext_factor(self.to_dense(f), self.ngens-1, self.domain)
        return (coeff, [ (self.from_dense(g), k) for g, k in factors ])

    def dup_gf_factor(self, f):
        coeff, factors = dup_gf_factor(self.to_dense(f), self.domain)
        return (coeff, [ (self.from_dense(g), k) for g, k in factors ])
    def dmp_gf_factor(self, f):
        coeff, factors = dmp_gf_factor(self.to_dense(f), self.ngens-1, self.domain)
        return (coeff, [ (self.from_dense(g), k) for g, k in factors ])

    def dup_factor_list(self, f):
        coeff, factors = dup_factor_list(self.to_dense(f), self.domain)
        return (coeff, [ (self.from_dense(g), k) for g, k in factors ])
    def dup_factor_list_include(self, f):
        factors = dup_factor_list_include(self.to_dense(f), self.domain)
        return [ (self.from_dense(g), k) for g, k in factors ]

    def dmp_factor_list(self, f):
        coeff, factors = dmp_factor_list(self.to_dense(f), self.ngens-1, self.domain)
        return (coeff, [ (self.from_dense(g), k) for g, k in factors ])
    def dmp_factor_list_include(self, f):
        factors = dmp_factor_list_include(self.to_dense(f), self.ngens-1, self.domain)
        return [ (self.from_dense(g), k) for g, k in factors ]

    def dup_irreducible_p(self, f):
        return dup_irreducible_p(self.to_dense(f), self.domain)
    def dmp_irreducible_p(self, f):
        return dmp_irreducible_p(self.to_dense(f), self.ngens-1, self.domain)

    def dup_sturm(self, f):
        seq = dup_sturm(self.to_dense(f), self.domain)
        return map(self.from_dense, seq)

    def dup_sqf_p(self, f):
        return dup_sqf_p(self.to_dense(f), self.domain)
    def dmp_sqf_p(self, f):
        return dmp_sqf_p(self.to_dense(f), self.ngens-1, self.domain)

    def dup_sqf_norm(self, f):
        s, F, R = dup_sqf_norm(self.to_dense(f), self.domain)
        return (s, self.from_dense(F), self.to_ground().from_dense(R))
    def dmp_sqf_norm(self, f):
        s, F, R = dmp_sqf_norm(self.to_dense(f), self.ngens-1, self.domain)
        return (s, self.from_dense(F), self.to_ground().from_dense(R))

    def dup_gf_sqf_part(self, f):
        return self.from_dense(dup_gf_sqf_part(self.to_dense(f), self.domain))
    def dmp_gf_sqf_part(self, f):
        return self.from_dense(dmp_gf_sqf_part(self.to_dense(f), self.domain))
    def dup_sqf_part(self, f):
        return self.from_dense(dup_sqf_part(self.to_dense(f), self.domain))
    def dmp_sqf_part(self, f):
        return self.from_dense(dmp_sqf_part(self.to_dense(f), self.ngens-1, self.domain))

    def dup_gf_sqf_list(self, f, all=False):
        coeff, factors = dup_gf_sqf_list(self.to_dense(f), self.domain, all=all)
        return (coeff, [ (self.from_dense(g), k) for g, k in factors ])
    def dmp_gf_sqf_list(self, f, all=False):
        coeff, factors = dmp_gf_sqf_list(self.to_dense(f), self.ngens-1, self.domain, all=all)
        return (coeff, [ (self.from_dense(g), k) for g, k in factors ])

    def dup_sqf_list(self, f, all=False):
        coeff, factors = dup_sqf_list(self.to_dense(f), self.domain, all=all)
        return (coeff, [ (self.from_dense(g), k) for g, k in factors ])
    def dup_sqf_list_include(self, f, all=False):
        factors = dup_sqf_list_include(self.to_dense(f), self.domain, all=all)
        return [ (self.from_dense(g), k) for g, k in factors ]
    def dmp_sqf_list(self, f, all=False):
        coeff, factors = dmp_sqf_list(self.to_dense(f), self.ngens-1, self.domain, all=all)
        return (coeff, [ (self.from_dense(g), k) for g, k in factors ])
    def dmp_sqf_list_include(self, f, all=False):
        factors = dmp_sqf_list_include(self.to_dense(f), self.ngens-1, self.domain, all=all)
        return [ (self.from_dense(g), k) for g, k in factors ]

    def dup_gff_list(self, f):
        factors = dup_gff_list(self.to_dense(f), self.domain)
        return [ (self.from_dense(g), k) for g, k in factors ]
    def dmp_gff_list(self, f):
        factors = dmp_gff_list(self.to_dense(f), self.ngens-1, self.domain)
        return [ (self.from_dense(g), k) for g, k in factors ]

    def dup_root_upper_bound(f, K):
        return dup_root_upper_bound(self.to_dense(f), self.domain)
    def dup_root_lower_bound(f, K):
        return dup_root_lower_bound(self.to_dense(f), self.domain)

    def dup_step_refine_real_root(self, f, M, fast=False):
        return dup_step_refine_real_root(self.to_dense(f), M, self.domain, fast=fast)
    def dup_inner_refine_real_root(self, f, M, eps=None, steps=None, disjoint=None, fast=False, mobius=False):
        return dup_inner_refine_real_root(self.to_dense(f), M, self.domain, eps=eps, steps=steps, disjoint=disjoint, fast=fast, mobius=mobius)
    def dup_outer_refine_real_root(self, f, s, t, eps=None, steps=None, disjoint=None, fast=False):
        return dup_outer_refine_real_root(self.to_dense(f), s, t, self.domain, eps=eps, steps=steps, disjoint=disjoint, fast=fast)
    def dup_refine_real_root(self, f, s, t, eps=None, steps=None, disjoint=None, fast=False):
        return dup_refine_real_root(self.to_dense(f), s, t, self.domain, eps=eps, steps=steps, disjoint=disjoint, fast=fast)
    def dup_inner_isolate_real_roots(self, f, eps=None, fast=False):
        return dup_inner_isolate_real_roots(self.to_dense(f), self.domain, eps=eps, fast=fast)
    def dup_inner_isolate_positive_roots(self, f, eps=None, inf=None, sup=None, fast=False, mobius=False):
        return dup_inner_isolate_positive_roots(self.to_dense(f), self.domain, eps=eps, inf=inf, sup=sup, fast=fast, mobius=mobius)
    def dup_inner_isolate_negative_roots(self, f, inf=None, sup=None, eps=None, fast=False, mobius=False):
        return dup_inner_isolate_negative_roots(self.to_dense(f), self.domain, inf=inf, sup=sup, eps=eps, fast=fast, mobius=mobius)
    def dup_isolate_real_roots_sqf(self, f, eps=None, inf=None, sup=None, fast=False, blackbox=False):
        return dup_isolate_real_roots_sqf(self.to_dense(f), self.domain, eps=eps, inf=inf, sup=sup, fast=fast, blackbox=blackbox)
    def dup_isolate_real_roots(self, f, eps=None, inf=None, sup=None, basis=False, fast=False):
        return dup_isolate_real_roots(self.to_dense(f), self.domain, eps=eps, inf=inf, sup=sup, basis=basis, fast=fast)
    def dup_isolate_real_roots_list(self, polys, eps=None, inf=None, sup=None, strict=False, basis=False, fast=False):
        return dup_isolate_real_roots_list(map(self.to_dense, polys), self.domain, eps=eps, inf=inf, sup=sup, strict=strict, basis=basis, fast=fast)
    def dup_count_real_roots(self, f, inf=None, sup=None):
        return dup_count_real_roots(self.to_dense(f), self.domain, inf=inf, sup=sup)
    def dup_count_complex_roots(self, f, inf=None, sup=None, exclude=None):
        return dup_count_complex_roots(self.to_dense(f), self.domain, inf=inf, sup=sup, exclude=exclude)
    def dup_isolate_complex_roots_sqf(self, f, eps=None, inf=None, sup=None, blackbox=False):
        return dup_isolate_complex_roots_sqf(self.to_dense(f), self.domain, eps=eps, inf=inf, sup=sup, blackbox=blackbox)
    def dup_isolate_all_roots_sqf(self, f, eps=None, inf=None, sup=None, fast=False, blackbox=False):
        return dup_isolate_all_roots_sqf(self.to_dense(f), self.domain, eps=eps, inf=inf, sup=sup, fast=fast, blackbox=blackbox)
    def dup_isolate_all_roots(self, f, eps=None, inf=None, sup=None, fast=False):
        return dup_isolate_all_roots(self.to_dense(f), self.domain, eps=eps, inf=inf, sup=sup, fast=fast)

    def fateman_poly_F_1(self):
        return tuple(map(self.from_dense, dmp_fateman_poly_F_1(self.ngens-1, self.domain)))
    def fateman_poly_F_2(self):
        return tuple(map(self.from_dense, dmp_fateman_poly_F_2(self.ngens-1, self.domain)))
    def fateman_poly_F_3(self):
        return tuple(map(self.from_dense, dmp_fateman_poly_F_3(self.ngens-1, self.domain)))
