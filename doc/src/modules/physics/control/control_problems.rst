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
* **Poles** and **Zeros** of the system elements (individual TFs in TFM) in the Laplace domain *(Note: The actual poles and zeros of a MIMO system are NOT the poles and zeros of the individual elements of the transfer function matrix)*. Also, visualise the poles and zeros of the individual transfer function corresponding to the **1st input** and **1st output** of the **G(s)** matrix.
* Plot the **unit step response** of the individual Transfer Function corresponding to the **1st input** and **1st output** of the **G(s)** matrix.
* Analyse the Bode magnitude and phase plot of the Transfer Function corresponding to **1st input** and **2nd output** of the **G(s)** matrix.
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

    >>> # Imports
    >>> from sympy import Matrix, laplace_transform, inverse_laplace_transform, exp, cos, sqrt, sin, pprint
    >>> from sympy.abc import s, t
    >>> from sympy.physics.control import *

    Subpart 1

    >>> g =  Matrix([[exp(-t)*(1 - t), exp(-2*t)], [5*exp((-2*t))-exp((-t)), (cos((sqrt(3)*t)/2) - 3*sqrt(3)*sin((sqrt(3)*t)/2))*exp(-t/2)]])
    >>> G = g.applyfunc(lambda a: laplace_transform(a, t, s)[0])
    >>> pprint(G)
    [      s              1     ]
    [ ------------      -----   ]
    [  2                s + 2   ]
    [ s  + 2*s + 1              ]
    [                           ]
    [    4*s + 3        s - 4   ]
    [---------------  ----------]
    [(s + 1)*(s + 2)   2        ]
    [                 s  + s + 1]

    Subpart 2

    >>> G = TransferFunctionMatrix.from_Matrix(G, s)
    >>> type(G)
	<class 'sympy.physics.control.lti.TransferFunctionMatrix'>
    >>> type(G[0])
	<class 'sympy.physics.control.lti.TransferFunction'>
    >>> print(f'Inputs = {G.num_inputs}, Outputs = {G.num_outputs}')
    Inputs = 2, Outputs = 2

    Subpart 3

	>>> G.elem_poles()
	[[[-1, -1], [-2]], [[-2, -1], [-1/2 - sqrt(3)*I/2, -1/2 + sqrt(3)*I/2]]]
	>>> G.elem_zeros()
	[[[0], []], [[-3/4], [4]]]
	>>> pole_zero_plot(G[0, 0])   # doctest: +SKIP

    .. image:: https://user-images.githubusercontent.com/53227127/133929754-c3253709-d766-44d6-a472-4623b89beb35.png

    Subpart 4

    >>> tf1 = G[0, 0]
    >>> pprint(tf1)
         s      
    ------------
     2          
    s  + 2*s + 1
    >>> step_response_plot(tf1)  # doctest: +SKIP

    .. image:: https://user-images.githubusercontent.com/53227127/133929755-694da2b1-4097-4e28-9a61-00bc055c5e7f.png

    Subpart 5

    >>> tf2 = G[0, 1]
    >>> bode_magnitude_plot(tf2)  # doctest: +SKIP

    .. image:: https://user-images.githubusercontent.com/53227127/133929749-cb150528-ed5a-41dc-89dc-6b45eff22f4f.png

    >>> bode_phase_plot(tf2)  # doctest: +SKIP

    .. image:: https://user-images.githubusercontent.com/53227127/133929753-def77202-2d05-47eb-8e09-c1f049d60186.png

    Subpart 6

	>>> P_mat = Matrix([[1/s, 2/(2+s)], [0, 3]])
	>>> C_mat = Matrix([[1, 1], [2, 2]])

    >>> P = TransferFunctionMatrix.from_Matrix(P_mat, var=s)
    >>> C = TransferFunctionMatrix.from_Matrix(C_mat, var=s)
    >>> # Series equivalent, considering (Input)→[P]→[C]→(Output). Note that order of matrix multiplication is opposite to the order in which the elements are arranged.
    >>> pprint(C*P)
    [1  1]    [1    2  ]   
    [-  -]    [-  -----]   
    [1  1]    [s  s + 2]   
    [    ]   *[        ]   
    [2  2]    [0    3  ]   
    [-  -]    [-    -  ]   
    [1  1]{t} [1    1  ]{t}
    >>> # Series equivalent, considering (Input)→[C]→[P]→(Output).
    >>> pprint(P*C)
    [1    2  ]    [1  1]   
    [-  -----]    [-  -]   
    [s  s + 2]    [1  1]   
    [        ]   *[    ]   
    [0    3  ]    [2  2]   
    [-    -  ]    [-  -]   
    [1    1  ]{t} [1  1]{t}
    
	>>> pprint((P*C).doit())  # doctest: +SKIP
	[ 5*s + 2    5*s + 2 ]   
	[---------  ---------]   
	[s*(s + 2)  s*(s + 2)]   
	[                    ]   
	[    6          6    ]   
	[    -          -    ]   
	[    1          1    ]{t}

    Subpart 7

    >>> tfm_feedback = MIMOFeedback(P, C, sign=-1)
    >>> pprint(tfm_feedback.doit())  # ((I + P*C)**-1)*P
    [    7*s + 14          -s - 6    ]   
    [---------------  ---------------]   
    [   2                2           ]   
    [7*s  + 19*s + 2  7*s  + 19*s + 2]   
    [                                ]   
    [                     2          ]   
    [   -6*s - 12      3*s  + 9*s + 6]   
    [---------------  ---------------]   
    [   2                2           ]   
    [7*s  + 19*s + 2  7*s  + 19*s + 2]{t}
