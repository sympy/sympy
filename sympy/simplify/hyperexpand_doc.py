""" This module cooks up a docstring when imported. Its only purpose is to
    be displayed in the sphinx documentation. """
from __future__ import annotations

def _generate_doc():
    from sympy import Eq, hyper
    from sympy.printing.latex import latex
    from sympy.simplify.hyperexpand import FormulaCollection
    from sympy.core.function import AppliedUndef, Function

    FUNCTION_LINKS = {
        # Special functions - gamma family
        'gamma': ':func:`~sympy.functions.special.gamma_functions.gamma`',
        'lowergamma': ':func:`~sympy.functions.special.gamma_functions.lowergamma`',
        'uppergamma': ':func:`~sympy.functions.special.gamma_functions.uppergamma`',

        # Bessel functions
        'besselj': ':func:`~sympy.functions.special.bessel.besselj`',
        'besseli': ':func:`~sympy.functions.special.bessel.besseli`',

        # Error functions
        'erf': ':func:`~sympy.functions.special.error_functions.erf`',
        'expint': ':func:`~sympy.functions.special.error_functions.expint`',
        'Ei': ':func:`~sympy.functions.special.error_functions.Ei`',
        'Ci': ':func:`~sympy.functions.special.error_functions.Ci`',
        'Si': ':func:`~sympy.functions.special.error_functions.Si`',
        'Shi': ':func:`~sympy.functions.special.error_functions.Shi`',
        'Chi': ':func:`~sympy.functions.special.error_functions.Chi`',
        'fresnels': ':func:`~sympy.functions.special.error_functions.fresnels`',
        'fresnelc': ':func:`~sympy.functions.special.error_functions.fresnelc`',

        # Elliptic integrals
        'elliptic_k': ':func:`~sympy.functions.special.elliptic_integrals.elliptic_k`',
        'elliptic_e': ':func:`~sympy.functions.special.elliptic_integrals.elliptic_e`',

        # Elementary functions - exponential
        'exp': ':func:`~sympy.functions.elementary.exponential.exp`',
        'log': ':func:`~sympy.functions.elementary.exponential.log`',

        # miscellaneous
        'sqrt': ':func:`~sympy.functions.elementary.miscellaneous.sqrt`',
        'root': ':func:`~sympy.functions.elementary.miscellaneous.root`',

        # trigonometric
        'cos': ':func:`~sympy.functions.elementary.trigonometric.cos`',
        'sin': ':func:`~sympy.functions.elementary.trigonometric.sin`',

        # hyperbolic
        'cosh': ':func:`~sympy.functions.elementary.hyperbolic.cosh`',
        'sinh': ':func:`~sympy.functions.elementary.hyperbolic.sinh`',

        # Combinatorial functions
        'RisingFactorial': ':func:`~sympy.functions.combinatorial.factorials.RisingFactorial`',
        'factorial': ':func:`~sympy.functions.combinatorial.factorials.factorial`',

        # Other special functions
        'lerchphi': ':func:`~sympy.functions.special.zeta_functions.lerchphi`',

    }

    def get_special_functions_used(expr):
        """Extract special function names from an expression."""
        func_names = set()
        for atom in expr.atoms(Function):
            if not isinstance(atom, AppliedUndef):
                func_name = atom.func.__name__
                if func_name in FUNCTION_LINKS:
                    func_names.add(func_name)
        return sorted(func_names)

    c = FormulaCollection()

    lines = []

    for f in c.formulae:
        obj = Eq(hyper(f.func.ap, f.func.bq, f.z),
                f.closed_form.rewrite('nonrepsmall'))
        lines.append(".. math::")
        lines.append("")
        lines.append("    %s" % latex(obj))
        lines.append("")
        # Add links to special functions used
        special_funcs = get_special_functions_used(f.closed_form)
        if special_funcs:
            lines.append("")
            n = len(special_funcs)
            if n == 1:
                lines.append("This formula involves %s." % FUNCTION_LINKS[special_funcs[0]])
            elif n == 2:
                lines.append("This formula involves %s and %s." % (
                    FUNCTION_LINKS[special_funcs[0]],
                    FUNCTION_LINKS[special_funcs[1]]
                ))
            else:
                func_links = ', '.join(FUNCTION_LINKS[fn] for fn in special_funcs[:-1])
                lines.append("This formula involves %s, and %s."% (
                    func_links,
                    FUNCTION_LINKS[special_funcs[-1]]
                ))
            lines.append("")
    return '\n'.join(lines)

__doc__ = _generate_doc()
