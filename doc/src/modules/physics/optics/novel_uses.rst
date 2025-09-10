=========================
Novel Uses of Optics Features
=========================

This document details novel uses for the newly implemented features in the `sympy.physics.optics` module. These features include `geometric_conj_af`, `geometric_conj_bf`, `gaussian_conj`, and `conjugate_gauss_beams`.

Geometric Conjugation Functions
-------------------------------

The `geometric_conj_af` and `geometric_conj_bf` functions are used to find the conjugate distances in geometrical optics. These functions are useful in designing optical systems where precise control over object and image distances is required.

Example:
    ```python
    from sympy.physics.optics.gaussopt import geometric_conj_af, geometric_conj_bf
    from sympy import symbols

    a, f = symbols('a f')
    b = geometric_conj_af(a, f)
    print(b)  # Output: a*f/(a - f)

    b, f = symbols('b f')
    a = geometric_conj_bf(b, f)
    print(a)  # Output: b*f/(b - f)
    ```

Gaussian Beam Conjugation
-------------------------

The `gaussian_conj` function is used to find the conjugate parameters for Gaussian beams. This function is useful in designing laser systems where precise control over beam parameters is required.

Example:
    ```python
    from sympy.physics.optics.gaussopt import gaussian_conj
    from sympy import symbols

    s_in, z_r_in, f = symbols('s_in z_r_in f')
    s_out, z_r_out, m = gaussian_conj(s_in, z_r_in, f)
    print(s_out, z_r_out, m)
    ```

Conjugate Gaussian Beams
------------------------

The `conjugate_gauss_beams` function is used to find the optical setup that conjugates the object and image waists of Gaussian beams. This function is useful in designing optical systems where precise control over beam waists is required.

Example:
    ```python
    from sympy.physics.optics.gaussopt import conjugate_gauss_beams
    from sympy import symbols

    wavelen, waist_in, waist_out, f = symbols('wavelen waist_in waist_out f')
    s_in, s_out, f = conjugate_gauss_beams(wavelen, waist_in, waist_out, f=f)
    print(s_in, s_out, f)
    ```

References
----------

For more information on these functions, refer to the `sympy.physics.optics` module documentation.
