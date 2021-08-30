=============================================
Solved Textbook Problems using Control Module
=============================================

Given below, are some comprehensive textbook examples to demonstrate the possible use cases
of the Control Module.

Example 1
---------

        $g(t) = \begin{bmatrix}
        (1-t)e^{-t} & e^{-2t} \\
        -e^{-t}+5e^{-2t} & \left(-3\sqrt{3}\sin\left(\frac{\sqrt{3}t}{2}\right)+\cos\left(\frac{\sqrt{3}t}{2}\right)\right)e^{-\frac{t}{2}}
        \end{bmatrix}$

With respect to the signal matrix in time-domain, also known as the *impulse response matrix* (**g(t)**) given above, find -

* The system matrix (Transfer Function Matrix) in the Laplace domain (**g(t)** → **G(s)**).
* The number of input and output signals in the system.
* Whether the system matrix is **proper**. If so, is it **strictly proper**? Also, check whether it is **bi-proper** in case it is proper.
* Whether the given MIMO system is **IO-Coupled**. *(When each input affects many outputs, the system is called IO-Coupled)*
* **Poles** and **Zeros** of the system elements (individual TFs in TFM) in the Laplace domain. *(Note: The actual poles and zeros of a MIMO system are NOT the poles and zeros of the individual elements of the transfer function matrix)*
* Plot the **unit step response** of the Transfer Function corresponding to the **1st input** to **1st output** of the **G(s)** matrix.
* The equivalent system matrix having **P(s)** and **C(s)** arranged in a series configuration *(Values of P(s) and C(s) are provided below)*. Compute the equivalent system matrix when the order of blocks is reversed *(i.e. C(s) then P(s))*.

        $P(s) = \begin{bmatrix}
        \frac{1}{s} & \frac{2}{s+2} \\
        0 & 3
        \end{bmatrix}$

        $C(s) = \begin{bmatrix}
        1 & 1 \\
        2 & 2
        \end{bmatrix}$

* Also, find the **equivalent closed-loop system** *(or the ratio v/u from the block diagram given below)* for the system (negative-feedback loop) having **C(s)** as the **controller** and **P(s)** as **plant** *(Refer to the block diagram given below)*.

        .. image:: https://user-images.githubusercontent.com/53227127/120820301-0b368c80-c572-11eb-84c4-e372190cf0bd.png

Solution-

    >>> from sympy import *
    >>> from sympy.physics.control import lti
    >>> from sympy.integrals.transforms import laplace_transform, inverse_laplace_transform
    >>> from sympy.abc import s, t

    Subpart 1

    >>> g =  Matrix([[exp(-t)*(1 - t), exp(-2*t)], [5*exp((-2*t))-exp((-t)), (cos((sqrt(3)*t)/2) - 3*sqrt(3)*sin((sqrt(3)*t)/2))*exp(-t/2)]])  # doctest: +SKIP
    >>> G = g.applyfunc(lambda a: laplace_transform(a, t, s)[0])  # doctest: +SKIP
    >>> pprint(G)  # doctest: +SKIP
    ⎡      s              1     ⎤
    ⎢ ────────────      ─────   ⎥
    ⎢  2                s + 2   ⎥
    ⎢ s  + 2⋅s + 1              ⎥
    ⎢                           ⎥
    ⎢    4⋅s + 3        s - 4   ⎥
    ⎢───────────────  ──────────⎥
    ⎢(s + 1)⋅(s + 2)   2        ⎥
    ⎣                 s  + s + 1⎦

    Subpart 2

    >>> G = to_transfer_function_matrix(G, s)  # doctest: +SKIP
    >>> type(G)  # doctest: +SKIP
    sympy.physics.control.lti.TransferFunctionMatrix
    >>> type(G[0])  # doctest: +SKIP
    sympy.physics.control.lti.TransferFunction
    >>> print(f'Inputs = {G.num_inputs}, Outputs = {G.num_outputs}')  # doctest: +SKIP
    Inputs = 2, Outputs = 2

    Subpart 3

    >>> G.is_proper  # doctest: +SKIP
    True
    >>> G.is_strictly_proper  # doctest: +SKIP
    True
    >>> G.is_biproper  # doctest: +SKIP
    False

    Subpart 4

    >>> G.is_io_coupled  # doctest: +SKIP
    True

    Subpart 5

    >>> G.poles()  # doctest: +SKIP
    [[[-1, -1], [-2]],[[-1, -2], [-1 + I*sqrt(3)/2, -1 - I*sqrt(3)/2]]]
    >>> G.zeros()  # doctest: +SKIP
    [[[0], []],[-3/4, 4]]

    Subpart 6

    >>> tf1 = G[0, 0]  # doctest: +SKIP
    >>> pprint(tf1)  # doctest: +SKIP
        s      
    ────────────
    2          
    s  + 2⋅s + 1
    >>> unit_step_response(tf1)  # doctest: +SKIP


    .. image:: https://user-images.githubusercontent.com/53227127/120837901-e5ff4980-c584-11eb-8b97-211680cfceb9.png


    Subpart 7

    >>> exp1 = 1/s  # doctest: +SKIP
    >>> exp2 = 2/(2+s)  # doctest: +SKIP
    >>> exp2 = 0  # doctest: +SKIP
    >>> exp4 = 3  # doctest: +SKIP
    >>> P = TransferFunctionMatrix([[exp1, exp2],[exp3, exp4]], var=s)  # doctest: +SKIP
    >>> C = TransferFunctionMatrix([[1, 1],[2, 2]], var=s)  # doctest: +SKIP
    # Series equivalent, considering (Input)→[P]→[C]→(Output). Note that order of matrix multiplication is opposite to the order in which the elements are arranged.
    >>> pprint(C*P)  # doctest: +SKIP
    ⎡1  3⋅s + 8 ⎤
    ⎢─  ─────── ⎥
    ⎢s   s + 2  ⎥
    ⎢           ⎥
    ⎢2  6⋅s + 16⎥
    ⎢─  ────────⎥
    ⎣s   s + 2  ⎦
    # Series equivalent, considering (Input)→[C]→[P]→(Output).
    >>> pprint(P*C)  # doctest: +SKIP
    ⎡5⋅s + 2   5⋅s + 2 ⎤
    ⎢────────  ────────⎥
    ⎢ 2         2      ⎥
    ⎢s  + 2⋅s  s  + 2⋅s⎥
    ⎢                  ⎥
    ⎣   6         6    ⎦

    Subpart 8

    >>> tfm_feedback = Feedback(P, C, type='neg')  # doctest: +SKIP
    >>> pprint(tfm_feedback.doit())  # ((I + P*C)**-1)*P  # doctest: +SKIP 
    ⎡       2                                      ⎛   2       ⎞      ⎤
    ⎢    7⋅s  + 14⋅s        3⋅(-5⋅s - 2)         2⋅⎝7⋅s  + 14⋅s⎠      ⎥
    ⎢───────────────────  ─────────────── + ───────────────────────── ⎥
    ⎢  ⎛   2           ⎞     2                      ⎛   2           ⎞ ⎥
    ⎢s⋅⎝7⋅s  + 19⋅s + 2⎠  7⋅s  + 19⋅s + 2   (s + 2)⋅⎝7⋅s  + 19⋅s + 2⎠ ⎥
    ⎢                                                                 ⎥
    ⎢        2              ⎛ 2          ⎞         ⎛     2       ⎞    ⎥
    ⎢   - 6⋅s  - 12⋅s     3⋅⎝s  + 7⋅s + 2⎠       2⋅⎝- 6⋅s  - 12⋅s⎠    ⎥
    ⎢───────────────────  ──────────────── + ─────────────────────────⎥
    ⎢  ⎛   2           ⎞     2                       ⎛   2           ⎞⎥
    ⎣s⋅⎝7⋅s  + 19⋅s + 2⎠  7⋅s  + 19⋅s + 2    (s + 2)⋅⎝7⋅s  + 19⋅s + 2⎠⎦
