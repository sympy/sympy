.. raw:: html

    <script type="text/javascript" >
        MathJax.Hub.Config({
            TeX: { equationNumbers: { autoNumber: "AMS" } }
        });
    </script>

.. role:: red
   :class: color:red

*****************
Geometric Algebra
*****************

:Author: Alan Bromborsky

.. |release| replace:: 0.10

.. % Complete documentation on the extended LaTeX markup used for Python
.. % documentation is available in ``Documenting Python'', which is part
.. % of the standard documentation for Python.  It may be found online
.. % at:
.. %
.. % http://www.python.org/doc/current/doc/doc.html
.. % \lstset{language=Python}
.. % \input{macros}
.. % This is a template for short or medium-size Python-related documents,
.. % mostly notably the series of HOWTOs, but it can be used for any
.. % document you like.
.. % The title should be descriptive enough for people to be able to find
.. % the relevant document.

.. % Increment the release number whenever significant changes are made.
.. % The author and/or editor can define 'significant' however they like.

.. % At minimum, give your name and an email address.  You can include a
.. % snail-mail address if you like.

.. % This makes the Abstract go on a separate page in the HTML version;
.. % if a copyright notice is used, it should go immediately after this.
.. %
.. % \ifhtml
.. % \chapter*{Front Matter\label{front}}
.. % \fi
.. % Copyright statement should go here, if needed.
.. % ...
.. % The abstract should be a paragraph or two long, and describe the

.. % scope of the document.

.. topic:: Abstract

   This document describes the implementation, installation and use of a
   geometric algebra module written in
   python that utilizes the sympy symbolic algebra library.  The python
   module ga has been developed for coordinate free calculations using
   the operations (geometric, outer, and inner products etc.) of geometric algebra.
   The operations can be defined using a completely arbitrary metric defined
   by the inner products of a set of arbitrary vectors or the metric can be
   restricted to enforce orthogonality and signature constraints on the set of
   vectors.  Additionally, a metric that is a function of a coordinate set can
   be defined so that a geometric algebra over a manifold can be implemented.
   Geometric algebras over submanifolds of the base manifold are also supported as
   well as linear multivector differential operators and linear transformations.
   In addition the module includes the geometric, outer (curl) and inner
   (div) derivatives. Tensors are included in the module as multilinear
   functions of vectors with contraction and covariant differentiation
   defined without the need of component indices.  For latex output a
   latex distribution must be installed.  A more detail description of the
   module and the mathematics behind it is at :download:`GA <./pdf/GA.pdf>`.

.. module:: sympy.galgebra.ga


Using LaTeX, Enhanced Printing, and Numerics
============================================

For developing code using *sympy* and the *GA* module I strongly suggested that you
go to <http://www.geany.org/Download/Releases> and install the version of the *geany*
editor appropriate for your system.  This editor allows you to execute python scripts
from within the editor, incorporate pep8 checking of python code, and use the ansicon
windows package to allow enhanced terminal printing on windows systems.

Additionally, *geany* also recongizes rst files and can easily be configured to compile
the rst files to html if you are documenting your code using Python-Sphinx.

In order to use the LaTeX printing options in the *GA* module install a LaTeX
distribution as follows:

    #. To install texlive in ubuntu goto console and enter *sudo apt-get install texlive-full*.

    #. To install miktex in windows goto <http://miktex.org/> and follow instructions for full installation.

    #. To install mactex in OSX goto <http://tug.org/mactex/> and follow instructions for full installation.

To use the numerical evaluation options install python-numpy:

    #. Install python-nympy if you want to calculate numerical matrix functons (determinant, inverse, eigenvalues, etc.).
       For windows go to <http://sourceforge.net/projects/numpy/files/NumPy/1.6.2/> and install the distribution of numpy
       appropriate for your system.

    #. For OSX go to <http://sourceforge.net/projects/numpy/files/NumPy/1.6.1/>.

If you wish to use "enhance_print" on windows install ansicon as follows:

        #. Go to <https://github.com/adoxa/ansicon/downloads> and download "ansicon"
        #. In the Edit -> Preferences -> Tools menu of "geany" enter into the Terminal input the full path of "ansicon.exe"

To use extra options *crop*, *png*, and *prog* in LaTeX printing goto section on LaTeX printing for instructions and limitations.


Module Components
=================

Instantiating a Geometric Algebra
---------------------------------

A geometric algebra is instantiated with

   *sympy.galgebra.Ga(basis, g=None, coords=None, norm=False, debug=False, X=None)*

   The *basis* and *g* parameters were described in `GA[page 10] <../../_downloads/GA.pdf#page=11>`_.

   If *debug=True* the data structure required to initialize the *Ga* class
   are printed out.

   *coords* is a tuple of *sympy* symbols equal in length to
   the number of basis vectors.  These symbols are used as the arguments of a
   multivector field as a function of position and for calculating the derivatives
   of a multivector field. Additionally, *Ga()* calculates the pseudo scalar,
   :math:`I` and pseudo scalar inverse :math:`I^{-1}` and makes them available to the programmer as *MV.I* and *MV.Iinv*.
   For the case of instantiating a 3-d geometric algebra in spherical coordinates we could use

       .. code-block:: python

          (r, th, phi) = coords = symbols('r,theta,phi', real=True)
          basis = 'e_r e_theta e_phi'
          g = [1, r**2, r**2*sin(th)**2]
          sp3d = Ga(basis,g=g,coords=coords,norm=True)

   The input :math:`X` allows the metric to be input as a vector manifold. :math:`X`
   is a list of functions of *coords* dimension, :math:`m`, equal to or greater than
   the number of coordinates. If *g=None* it is assumed that *X* is a vector in an
   :math:`m`-dimensional orthonormal Euclidean vector space.  If it is wished
   the embedding vector space to be non-Euclidean that condition is specified with
   *g*.  For example if we wish the embedding space to be a 5-dimensional Minkowski
   space then *g=[-1,1,1,1,1]*.  Then the *Ga* class uses *X* to calculate the
   manifold basis vectors as a function of the coordinates and from them the metric
   tensor.

   If *norm=True* the basis vectors of the manifold are normalized so that the
   absolute values of the squares of the basis vectors are one.  It is suggested
   that one only use this option for diagonal metric tensors, and even there do so
   with caution, due to the possible
   problems with taking the square root of a general *sympy* expression (one that has an
   unknown sign).

   In addition to the basis vectors, if coordinates are defined for the geometric algebra, the
   left and right geometric derivative operators are calculated and accessed with the *Ga*
   member function *grads()*.

   *Ga.grads()*

       *Ga.grads()* returns a tuple with the left and right geometric derivative operators. A
       typical usage would be

    .. code-block:: python

        (grad,rgrad) = sp3d.grads()

   for the spherical 3-d geometric algebra. The left derivative *grad* :math:`= \nabla` and the
   right derivative *rgrad* :math:`= \bar{\nabla}` are explained in
   `GA[page 23] <../../_downloads/GA.pdf#page=24>`_. Again
   the names *grad* and *rgrad* are whatever the user chooses them to be.

   an alternative instantiation method is

   *Ga.build(basis, g=None, coords=None, X=None, norm=False, debug=False)*

     The input parameters for *Ga.build()* are the same as for *Ga()*.  The difference is
     that in addition to returning the geometric algebra *Ga.build()* returns the basis vectors
     at the same time. Using *Ga.build()* in the previous example gives

     .. code-block:: python

       (r, th, phi) = coords = symbols('r,theta,phi', real=True)
       basis = 'e_r e_theta e_phi'
       g = [1, r**2, r**2*sin(th)**2]
       (sp3d,er,eth,ephi) = Ga.build(basis,g=g,coord=coords,norm=True)

   To access the pseudo scalar of the geometric algebra use the member function *I()*.

   *Ga.I()*

       *Ga.I()* returns the normalized pseudo scalar :math:`\left (\left | {I^{2}}\right |=1\right )` for the
       geometric algebra. For example :math:`I` = *o3d.I()* for the *o3d* geometric
       algebra.

   In general we have defined member fuctions of the *Ga* class that will instantiate objects
   of other classes since the objects of the other classes are all associated with a particular
   geometric algebra object.  Thus we have

    .. image:: tables/class_objs.png
        :width: 400px
        :align: center

   for the instantiation of various objects from the *Ga* class.  This means that in order to
   instantiate any of these objects we need only to import *Ga* into our program.


Instantiating a Multivector
---------------------------

    Since we need to associate each multivector with the geometric algebra that contains it
    we use a member function of *Ga* to instantiate every multivector (There is a
    multivector class, *Mv*, but in order the insure that every multivector is associated
    with the correct geometric algebra we always use the member function *Ga.mv* to instantiate
    the multivector.)  The multivector is instantiated with:

    *Ga.mv(name, mode, f=False)*

        As an example of both instantiating a geometric algebra and multivectors consider the
        following code fragment for a 3-d Euclidean geometric algebra.

        .. code-block:: python

            from sympy import symbols
            from ga import Ga
            (x, y, z) = coords = symbols('x,y,z',real=True)
            o3d = Ga('e_x e_y e_z', g=[1,1,1], coords=coords)
            (ex, ey, ez) = o3d.mv()
            V = o3d.mv('V','vector',f=True)

        First consider the multivector instantiation *V = o3d.mv('V','vector',f=True)*.  Here
        a 3-dimensional multivector field that is a function of *x*, *y*, and *z* (*f=True*) is
        being instantiated.  If latex output were used (to be discussed later) the multivector
        *V* would be displayed as

        .. math::
          :nowrap:

          \begin{equation}
            A^{x}\boldsymbol{e}_{x} + A^{y}\boldsymbol{e}_{y} + A^{z}\boldsymbol{e}_{z}
          \end{equation}

        Where the coefficients of the basis vectors are generalized *sympy* functions of the
        coordinates.  The superscripts (Denoted in text output by *A__x*, etc. so
        that for text output *A* would be printed as *A__x\*e_x+A__y\*e_y+A__z\*e_z*) are formed
        from the coordinate symbols or if there are no coordinates from the subscripts of
        the basis vectors.  The types of name and modes available for multivector instantiation are

            .. image:: tables/instanciate_mv.png
                :width: 750px
                :align: center

        Line 5 of the previous listing illustrates the case of using the *mv* member function with
        no arguments. The code does not return a multivector, but rather a tuple or the basis vectors of the geometric algebra *o3d*.
        The elements of the tuple then can
        be used to construct multivectors, or multivector fields through the operations
        of addition, subtraction, multiplication (geometric, inner, and outer products and left and right contraction).
        As an example we could construct the vector function

        .. code-block:: python

            F = x**2*ex + z*ey + x*y*ez

        or the bivector function

        .. code-block:: python

            B = z*(ex^ey) + y*(ey^ez) + y*(ex^ez).

    If one wished to calculate the left and right geometric derivatives of *F* and *B* the required code would be

    .. code-block:: python

        (grad,rgrad) = o3d.grads()
        dF = grad*F
        dB = grad*B
        dFr = F*rgrad
        dBr = B*rgrad

    *dF*, *dB*, *dFr*, and *dBr* are all multivector functions. For the code where the order of the operations are
    reversed

    .. code-block:: python

        (grad,rgrad) = o3d.grads()
        dFop = F*grad
        dBop = B*grad
        dFrop = rgrad*F
        dBrop = rgrad*B

    *dFop*, *dBop*, *dFrop*, and *dBrop* are all multivector differential operators (again see
    `GA[page 25] <../../_downloads/GA.pdf#page=26>`_).


Basic Multivector Class Functions
---------------------------------

    *components(self)*

       Return list of multivectors correspoinding to each base blade of multivector.

    *convert_to_blades(self)*

       Convert multivector from the base representation to the blade representation.
       If multivector is already in blade representation nothing is done.


    *convert_from_blades(self)*

       Convert multivector from the blade representation to the base representation.
       If multivector is already in base representation nothing is done.


    *diff(self,var)*

       Calculate derivative of each multivector coefficient with resepect to
       variable *var* and form new multivector from coefficients.


    *dual(self)*

       Return dual of multivector which is multivector left multiplied by
       pseudoscalar *Mv.i* (Hestenes,p22).

    *even(self)*

       Return the even grade components of the multivector.


    *exp(self,hint='+')*

        Return exponential of a multivector :math:`A` if :math:`A^{2}` is a scalar (if :math:`A^{2}` is not a scalar an
        error message is generated).  If :math:`A` is the multivector then :math:`\boldsymbol{e}^{A}` is returned
        where the default *hint*, *+*, assumes :math:`A^{2} > 0` so that


        .. math::
          :nowrap:

            \begin{equation}
                    \boldsymbol{e}^{A} = \cosh\sqrt{A^{2}}+\sinh\sqrt{A^{2}}\left (\frac{A}{\sqrt{A^{2}}}\right )
            \end{equation}

        If the mode is not *+* then :math:`A^{2} < 0` is assumed so that


        .. math::
          :nowrap:

            \begin{equation}
                    \boldsymbol{e}^{A} = \cos \sqrt{-A^{2}}+\sin\sqrt{-A^{2}}\left (\frac{A}{\sqrt{-A^{2}}}\right ).
            \end{equation}

        The hint is required for symbolic multivectors :math:`A` since in general *sympy* cannot determine if
        :math:`A^{2}` is positive or negative.  If :math:`A` is purely numeric the hint is ignored.


    *expand(self)*

       Return multivector in which each coefficient has been expanded using
       sympy *expand()* function.


    *factor(self)*

       Apply the sympy *factor* function to each coefficient of the multivector.


    *Fmt(self, fmt=1,title=None)*

        Function to print multivectors in different formats where

            .. image :: tables/fmt_opts.png
                :width: 400px
                :align: center

        *title* appends a title string to the beginning of the output.  An equal sign in
        the title string is not required, but is added as a default.


    *func(self,fct)*

       Apply the *sympy* scalar function *fct* to each coefficient of the multivector.


    *grade(self,igrade=0)*

        Return a multivector that consists of the part of the multivector of
        grade equal to *igrade*.  If the multivector has no *igrade* part
        return a zero multivector.


    *inv(self)*

       Return the inverse of the multivector :math:`M` (*M.inv()*) if :math:`MM^{\dagger}` is a nonzero
       scalar.  If :math:`MM^{\dagger}`
       is not a scalar the program exits with an error message.


    *norm(self)*

       Return the norm of the multivector :math:`M` (*M.norm()*) defined by :math:`\sqrt{MM^{\dagger}}` if
       :math:`MM^{\dagger}` is a
       scalar (a sympy scalar
       is returned).  If :math:`MM^{\dagger}` is not a scalar the program exits with an error message.


    *norm2(self)*

       Return the square of the norm of the multivector :math:`M` (*M.norm2()*) defined by :math:`MM^{\dagger}`
       if :math:`MM^{\dagger}`
       is a scalar (a sympy scalar
       is returned).  If :math:`MM^{\dagger}` is not a scalar the program exits with an error message.


    *proj(self,bases_lst)*

       Return the projection of the multivector :math:`M` (*M.proj(bases_lst)*) onto the subspace defined by the list of bases
       (*bases_lst*).


    *scalar(self)*

        Return the coefficient (sympy scalar) of the scalar part of a
        multivector.


    *simplify(self,mode=simplify)*

       *mode* is a sympy simplification function of a list/tuple of sympy
       simplification functions that are applied in sequence (if more than
       one function) each coefficient of the multivector.  For example if
       we wished to applied *trigsimp* and *ratsimp* sympy functions to the
       mulitvector *F* the code would be


       .. code-block:: python

          Fsimp = F.simplify(mode=[trigsimp,ratsimp]).


       Actually *simplify* could be used to apply any scalar sympy function to
       the coefficients of the multivector.


    *subs(self,x)*

       Return multivector where sympy *subs* function has been applied to each
       coefficient of multivector for argument dictionary/list *x*.


    *rev(self)*

       Return the reverse of the multivector.


    *set_coef(self,grade,base,value)*

       Set the multivector coefficient of index *(grade,base)* to *value*.


    *trigsimp(self,\*\*kwargs)*

       Apply the sympy trignometric simplification function *trigsimp* to
       each coefficient of the multivector. *\*\*kwargs* are the arguments of
       *trigsimp*.  See sympy documentation on *trigsimp* for more information.

Basic Multivector Functions
---------------------------

    *Com(A,B)*

       Calulate commutator of multivectors *A* and *B*.  Returns *(AB-BA)/2*.

    *GAeval(s,pstr=False)*

       Returns multivector expression for string *s* with operator precedence for
       string *s* defined by inputs to function *def_prec()*.  if *pstr=True*
       *s* and *s* with parenthesis added to enforce operator precedence are printed.

    *Nga(x,prec=5)*

       If *x* is a multivector with coefficients that contain floating point numbers, *Nga()*
       rounds all these numbers to a precision of *prec* and returns the rounded multivector.

    *ReciprocalFrame(basis,mode='norm')*

       If *basis* is a list/tuple of vectors, *ReciprocalFrame()* returns a tuple of reciprocal
       vectors.  If *mode=norm* the vectors are normalized.  If *mode* is anything other than
       *norm* the vectors are unnormalized and the normalization coefficient is added to the
       end of the tuple.  One must divide by the coefficient to normalize the vectors.

    *ScalarFunction(TheFunction)*

       If *TheFuction* is a real *sympy* fuction a scalar multivector function is returned.

    *cross(v1,v2)*

       If *v1* and *v2* are 3-dimensional euclidean vectors the vector cross product is
       returned, :math:`v_{1}\times v_{2} = -I\left ( v_{1}\wedge v_{2} \right )`.

    *def_prec(gd,op_ord='<>|,^,\*')*

       This is used with the *GAeval()* function to evaluate a string representing a multivector
       expression with a revised operator precedence.  *def_prec()* redefines the operator
       precedence for multivectors. *def_prec()* must be called in the main program an the
       argument *gd* must be *globals()*.  The argument *op_ord* defines the order of operator
       precedence from high to low with groups of equal precedence separated by commas. the default
       precedence *op_ord='<>|,^,\*'* is that used by Hestenes.

    *dual(M)*

       Return the dual of the multivector *M*, math:`MI^{-1}`.

    *inv(B)*

       If for the multivector :math:`B`, :math:`BB^{\dagger}` is a nonzero scalar, return
       :math:`B^{-1} = B^{\dagger}/(BB^{\dagger})`.

    *proj(B,A)*

       Project blade *A* on blade *B* returning :math:`\left ( A\lfloor B\right ) B^{-1}`.

    *refl(B,A)*

       Reflect blade *A* in blade *B*. If *r* is grade of *A* and *s* is grade of *B*
       returns :math:`(-1)^{s(r+1)}BAB^{-1}`.

    *rot(itheta,A)*

       Rotate blade *A* by 2-blade *itheta*.  Is is assumed that *itheta\*itheta > 0* so that
       the rotation is Euclidian and not hyperbolic so that the angle of
       rotation is *theta = itheta.norm()*.  Then in 3-dimensional Euclidian space. *theta* is the angle of rotation (scalar in radians) and
       *n* is the vector axis of rotation.  Returned is the rotor *cos(theta)+sin(theta)*N* where *N* is
       the normalized dual of *n*.


Linear Transformations
----------------------

    The mathematical background for linear transformations is in `GA[page 26] <../../_downloads/GA.pdf#page=27>`_.
    Linear transformations on the tangent space of
    the manifold are instantiated with the *Ga* member function *lt* (the actual class being instantiated is *Lt*) as shown in
    lines 12, 20, 26, and 44 of the code (Ltrans.py)

    .. literalinclude:: pysource/Ltrans.py
       :language: python


    In all of the examples in the above code the default instantiation is used which produces a general (all the
    coefficients of the linear transformation are symbolic constants) linear transformation. Note that to
    instantiate linear transformations
    coordinates, :math:`\left \{\boldsymbol{e}_{i} \right \}`, must be defined when the geometric algebra associated with the
    linear transformation is instantiated.
    This is due to the naming conventions of the general linear transformation (coordinate names are used)
    and for the calculation
    of the trace of the linear transformation which requires taking a divergence.
    To instantiate a specific linear transformation
    the usage of *lt()* is  *Ga.lt(M,f=False)*.

    *M* is an expression that can define the coefficients of the linear transformation in various ways
    defined as follows.

    .. image:: tables/lt_ops.png
        :width: 800px
        :align: center

    *f* is *True* or *False*. If *True* the symbolic coefficients of the general linear transformation are
    instantiated as functions of the coordinates.

    The different methods of instantiation are demonstrated in the code (*LtransInst.py*)

    .. literalinclude:: pysource/LtransInst.py
        :language: python

    with output :download:`Linear Transformation Instantiation <./pysource/LtransInst.pdf>`.

    The member functions of the *Lt* class are

    *Lt(A)*

        Returns the image of the multivector :math:`A` under the linear transformation :math:`L` where
        :math:`L(A)` is defined by the
        linearity of :math:`L`, the vector values :math:`L\left ( e_{i}\right )`, and the definition
        :math:`L \left ( e_{i_{1}}\wedge\dots\wedge e_{i_{r}}\right ) = L\left ( e_{i_{1}} \right )
        \wedge\dots\wedge L\left ( e_{i_{r}}\right )`.

    *Lt.det()*

        Returns the determinant (a scalar) of the linear transformation, :math:`L`, defined by
        :math:`\det (L)I = L(I)`.

    *Lt.adj()*

        Returns the adjoint (a linear transformation) of the linear transformation, :math:`L`, defined by
        :math:`a\cdot L(b) = b\cdot \bar{L}(a)` where
        :math:`a` and :math:`b` are any two vectors in the tangent space and :math:`\bar{L}` is the adjoint of :math:`L`.

    *Lt.tr()*

        Returns the trace (a scalar) of the linear transformation, :math:`L`, defined by
        :math:`\operatorname{tr}\left ( L\right )=\nabla_{a}\cdot L\left ( a \right )` where :math:`a` is a vector in the tangent space.

    *Lt.matrix()*

        Returns the matrix representation (sympy Matrix) of the linear transformation, :math:`L`, defined by
        :math:`L\left ( e_{i}\right ) = L_{ij} e_{j}` where :math:`L_{ij}` is the matrix representation.

    The *Ltrans.py* demonstrate the use of the various *Lt* member functions and operators. The operators that can
    be used with
    linear transformations are *+*, *-*, and *\**. If :math:`A` and :math:`B` are linear transformations,
    :math:`V` a multivector,
    and :math:`\alpha` a
    scalar then :math:`(A\pm B)(V) = A(V)\pm B(V)`, :math:`(AB)(V) = A(B(V))`, and
    :math:`(\alpha A)(V) = \alpha A(V)`.

    The ouput of *Ltrans.py* is :download:`Linear Transformations <./pysource/Ltrans.pdf>`.

Differential Operators
----------------------

    For the mathematical treatment of linear mulivector differential operators `GA[page 23] <../../_downloads/GA.pdf#page=24>`_.
    There is a differential
    operator class *Dop*. However, one never needs to use it directly.  The operators are constructed from linear
    combinations of multivector products of the differential operators *Ga.grad* and *Ga.rgrad* as shown in the following code for
    both orthogonal rectangular and spherical 3-d coordinate systems.

    .. literalinclude:: pysource/dop.py
        :language: python

    The output of this code is :download:`Differential Operators <./pysource/dop.pdf>`.

    .. warning::

        Since *Ga.grad* operates on multivectors and differential operators on it's right and *Ga.rgrad*
        operates on multivectors and differential operators on it's left one can never multiply (any
        multivector multiplication operator) a left and a right differential operator together since the
        result is not well defined and will generate an error message.  A right operator times a right
        operator gives a right operator and a left operator times a left operator gives a left operator.

    The printed output of a *Dop* (text or latex) is the sum of the products of scalar differential operators
    times the base-blades of the geometric algebra. If the *Dop* operates to the right the base-blades are
    printed to the left of the scalar differential operator, if the *Dop* operates to the left the base-blades are
    printed to the right of the scalar differential operator.  For examples see the first two lines in
    :download:`Differential Operators <./pysource/dop.pdf>`.

    The various derivatives of a multivector function is accomplished by
    multiplying the gradient differential operator with the function.  The gradient
    differential operation is returned by the *Ga.mv()* function if coordinates
    are defined.  For example if we have for a 3-D vector space

    .. code-block:: python

        X = (x,y,z) = symbols('x y z')
        o3d = Ga('e*x|y|z',metric='[1,1,1]',coords=X)
        (ex,ey,ez) = o3d.mv()
        (grad,rgrad) = o3d.grads()


    Then the gradient operator is *grad* (actually the user can give
    it any name he wants to).  Then the derivatives of the multivector
    function *F = o3d.mv('F','mv',f=True)* are given by multiplying by the
    left geometric derivative operator and the right geometric derivative operator
    *grad* :math:`= \nabla` and *rgrad* :math:`= \bar{\nabla}`.  Another option
    is to use the gradient operator members of the geometric algebra directly where we have
    :math:`\nabla =` *o3d.grad* and :math:`\bar{\nabla} =` *o3d.rgrad*.

    .. math::
      :nowrap:

      \begin{align*}
            \nabla F &=  \mbox{grad $*$ F} \\
            F \bar{\nabla} &=  \mbox{F $*$ rgrad} \\
            \nabla \wedge F &=  \mbox{grad ^ F} \\
            F \wedge \bar{\nabla} &=  \mbox{F ^ rgrad} \\
            \nabla \cdot F &=  \mbox{grad $|$ F} \\
            F \cdot \bar{\nabla} &=  \mbox{F $|$ rgrad} \\
            \nabla \lfloor F &=  \mbox{grad $<$ F} \\
            F \lfloor \bar{\nabla} &=  \mbox{F $<$ rgrad} \\
            \nabla \rfloor F &=  \mbox{grad $>$ F} \\
            F \rfloor \bar{\nabla} &= \mbox{F $>$ rgrad}
      \end{align*}


    The preceding code block gives examples of all possible multivector
    derivatives of the multivector function *F* where the operation returns
    a multivector function. The complementary operations

        .. math::
          :nowrap:

          \begin{align*}
                F \nabla &=  \mbox{F $*$ grad} \\
                \bar{\nabla} F &=  \mbox{rgrad $*$ F} \\
                F \wedge \nabla &=  \mbox{F ^ grad} \\
                \bar{\nabla} \wedge F &=  \mbox{rgrad ^ F} \\
                F \cdot \nabla &=  \mbox{F $|$ grad} \\
                \bar{\nabla}\cdot F &=  \mbox{rgrad $|$ F} \\
                F \lfloor \nabla &=  \mbox{F $<$ grad} \\
                \bar{\nabla} \lfloor F &=  \mbox{rgrad $<$ F} \\
                F \rfloor \nabla &=  \mbox{F $>$ grad} \\
                \bar{\nabla} \rfloor F &= \mbox{rgrad $>$ F}
          \end{align*}

    all return multivector linear differential operators.

    One additional class function of utility is

    *components(self)*

       Return list of differential operators corresponding to each base blade in *self*
       differential operator.  For example ``o3d.grad.components()`` would print as
       :math:`\left [ \boldsymbol{e}_{x}\frac{\partial}{\partial x}, \boldsymbol{e}_{y}\frac{\partial}{\partial y},\boldsymbol{e}_{z}\frac{\partial}{\partial z}\right ]`,
       while ``o3d.rgrad.components()`` would print as
       :math:`\left [ \frac{\partial}{\partial x}\boldsymbol{e}_{x}, \frac{\partial}{\partial y}\boldsymbol{e}_{y},\frac{\partial}{\partial z}\boldsymbol{e}_{z}\right ]`.


Submanifolds
------------

    In general the geometric algebra that the user defines exists on the tangent space of
    a manifold.  The submanifold class, *Sm*, is derived from
    the *Ga* class and allows one
    to define a submanifold of a manifold by defining a coordinate mapping between the submanifold
    coordinates and the manifold coordinates.  What is returned is the geometric
    algebra of the tangent space of the submanifold. The submanifold for a geometric algebra is
    instantiated with

    *Ga.sm(map,coords,root='e',norm=False)*

        To define the submanifold we must define a coordinate map from the coordinates of the submanifold to
        each of the coordinates of the base manifold.  Thus the arguments *map* and *coords* are
        respectively lists of functions and symbols.  The list of symbols, *coords*, are the coordinates of the
        submanifold and are of length equal to the dimension of the submanifold.  The list of functions, *map*,
        define the mapping from the coordinate space of the submanifold to the coordinate space of the
        base manifold.  The length of *map* is equal to the dimension of the base manifold and each function in
        *map* is a function of the coordinates of the submanifold. As a concrete example consider the
        following code.

        .. literalinclude:: pysource/manifold_src.py
            :language: python

        The base manifold, *sp3d*, is a 3-d Euclidian space using standard spherical coordinates. The submanifold
        *sph2d* of *sp3d* is a spherical surface of radius :math:`1`.  To take the sumanifold operation one step further
        the submanifold *cir1d* of *sph2d* is a circle in *sph2d* where the latitude of the circle is :math:`\pi/8`.

        In each case, for demonstration purposes, a scalar and vector function on each manifold is defined (*f* and *F*
        for the 2-d manifold and *h* and *H* for the 1-d manifold) and the geometric derivative of each function is taken.  The
        manifold mapping and the metric tensor for *cir1d* of *sph2d* are also shown. Note that if the submanifold
        basis vectors are not normalized the program output is :download:`Submanifold <./pysource/manifold_src.pdf>`.


Instantiating a Multi-linear Functions (Tensors)
------------------------------------------------

The mathematical background for multi-linear functions is in `GA[page 27] <../../_downloads/GA.pdf#page=28>`_.  Note
that the multilinear transformations represent *concrete* tensors <http://en.wikipedia.org/wiki/Multilinear_map> as
opposed to *abstract* tensors <http://en.wikipedia.org/wiki/Abstract_index_notation>.
To instantiate a multi-linear function use

*Mlt(self, f, Ga, nargs=None, fct=False)* where the arguments are

        .. image:: tables/tensor_inst.png
            :width: 600px
            :align: center

Examples of different types of instantiation are

    .. code-block:: python

        #TensorDef.py

        import sys
        from sympy import symbols,sin,cos
        from printer import Format,xpdf,Get_Program,Print_Function
        from ga import Ga
        from lt import Mlt

        coords = symbols('t x y z',real=True)
        (st4d,g0,g1,g2,g3) = Ga.build('gamma*t|x|y|z',g=[1,-1,-1,-1],coords=coords)

        A = st4d.mv('T','bivector')

        def TA(a1,a2):
            global A
            return A | (a1 ^ a2)

        T = Mlt(TA,st4d) # Define multi-linear function

Basic Multilinear Function Class Functions
------------------------------------------

If we can instantiate multilinear functions we can use all the multilinear function class functions as described as follows.
See `GA[page 27] <../../_downloads/GA.pdf#page=28>`_ for the mathematical description of each operation.

    *self(kargs)*

       Calling function to evaluates multilinear function for \T{kargs} list of vector arguments and returns a value.  Note that a
       sympy scalar is returned, \emph{not} a multilinear function.


    *self.contract(slot1,slot2)*

        Returns contraction of tensor between *slot1* and *slot2* where *slot1* is the index of the first vector argument and
        *slot2* is the index of the second vector argument of the tensor. For example if we have a rank two tensor, *T(a1,a2)*,
        then *T.contract(1,2)* is the contraction of *T*.  For this case since there are only two slots there can only be one
        contraction.


    *self.pdiff(slot)*

        Returns gradient of tensor, *T*, with respect to slot vector.  For example if the tensor is :math:`T \left ( a_{1},a_{2} \right )` then
        *T.pdiff(2)* is :math:`\nabla_{a_{2}}T`.  Since *T* is a scalar function, *T.pdiff(2)* is a vector function.


    *self.cderiv()*

        Returns covariant derivative of tensor field. If *T* is a tensor of rank :math:`k` then *T.cderiv()* is a tensor of rank :math:`k+1`.
        The operation performed is defined in `GA[page 30] <../../_downloads/GA.pdf#page=31>`_.

Standard Printing
-----------------

    Printing of multivectors is handled by the module *printer* which contains
    a string printer class derived from the sympy string printer class and a latex
    printer class derived from the sympy latex printer class.  Additionally, there
    is an *Eprint* class that enhances the console output of sympy to make
    the printed output multivectors, functions, and derivatives more readable.
    *Eprint* requires an ansi console such as is supplied in linux/osx or the
    program *ansicon* (github.com/adoxa/ansicon) for windows which replaces *cmd.exe*.

    For a windows user the simplest way to implement ansicon is to use the *geany*
    editor and in the Edit :math:`\rightarrow` Preferences :math:`\rightarrow` Tools menu replace *cmd.exe* with
    *ansicon.exe* (be sure to supply the path to *ansicon*).

    If *Eprint* is called in a program (linux) when multivectors are printed
    the basis blades or bases are printed in bold text, functions are printed in red,
    and derivative operators in green.

    For formatting the multivector output there is the member function

    *Fmt(self,fmt=1,title=None)*

        *Fmt* is used to control how the multivector is printed with the argument
        *fmt*.  If *fmt=1* the entire multivector is printed on one line.  If
        *fmt=2* each grade of the multivector is printed on one line.  If *fmt=3*
        each component (base) of the multivector is printed on one line.  If a
        *title* is given then *title=multivector* is printed.  If the usual print
        command is used the entire multivector is printed on one line.

Latex Printing
--------------

    For latex printing one uses one functions from the *ga* module and one
    function from the *printer* module.  The
    functions are

    *Format(Fmode=True,Dmode=True,ipy=False)*

       This function from the *ga* module turns on latex printing with the
       following options

        .. image:: tables/latex_prnt.png
            :width: 700px
            :align: center

    *xpdf(filename=None,paper=(14,11),crop=False,png=False,prog=False,debug=False)*

       This function from the *printer* module post-processes the output captured from
       print statements, writes the resulting latex strings to the file *filename},
       processes the file with pdflatex, and displays the resulting pdf file.   All latex files except
       the pdf file are deleted. If *debug = True* the file *filename* is printed to
       standard output for debugging purposes and *filename* (the tex file) is saved.  If *filename* is not entered the default
       filename is the root name of the python program being executed with *.tex* appended.
       The *paper* option defines the size of the paper sheet for latex. The format for the *paper* input is

        .. image:: tables/latex_paper.png
            :width: 600px
            :align: center

       The default of *paper=(14,11)* was chosen so that long multivector expressions would not be truncated on the display.

       If the *crop* input is *True* the linux *pdfcrop* program is used to crop the pdf output (if output is one page).  This only works
       for linux installations (where *pdfcrop* is installed).

       If the *png* input is *True* the imagemagick convert program is used to output a portable network graphics format file.  This only
       work in linux installations where imagemagick is installed and the the script *Pdf2Png* (copy in galgebra/utilities directory) is
       in the executable path.

       If the *prog* input is *True* the python program that produces the \LaTeX output is printed to the pdf output first using the
       lstlisting package in the \LaTeX distribution.

       If the *debug* input is *True* the tex file generated is printed to the console.

       The *xpdf* function requires that latex and a pdf viewer be installed on
       the computer.

       *xpdf* is not required when printing latex in IPython notebook.


    As an example of using the latex printing options when the following code is
    executed.

    .. literalinclude:: pysource/print_example1.py
        :language: python

    The output of the code is :download:`LaTeX output example <./pysource/print_example1.pdf>` is displayed

    For the cases of derivatives the code is

    .. literalinclude:: pysource/print_example2.py
        :language: python

    and the latex displayed output is (:math:`f` is a scalar function) :download:`LaTeX output for derivatives <./pysource/print_example2.pdf>`.

    This example also demonstrates several other features of the latex printer.  In the
    case that strings are input into the latex printer such as *grad\*\\bm{A}*,
    *grad^\\bm{A}*, or *grad\*\\bm{A}'!.  The text symbols *grad*, *^*, *|*, and
    *\** are mapped by the *xpdf()* post-processor as follows if the string contains
    an *=*.

    .. image:: tables/latex_replc.png
        :width: 375px
        :align: center

    If the first character in the string to be printed is a *\%* none of the above substitutions
    are made before the latex processor is applied.  In general for the latex
    printer strings are assumed to be in a math environment (equation or
    align) unless the first character in the string is a *\#*. [#f3]_

      Except where noted the conventions for latex printing follow those of the
      latex printing module of sympy. This includes translating sympy variables
      with Greek name (such as *alpha*) to the equivalent Greek symbol
      (math:`\alpha`) for the purpose of latex printing.  Also a single
      underscore in the variable name (such as *X_j*) indicates a subscript
      (:math:`X_{j}`), and a double underscore (such as *X__k*) a
      superscript (:math:`X^{k}`).  The only other change with regard to the
      sympy latex printer is that matrices are printed full size (equation
      displaystyle).

For formatting the tensor (*Mlt*) output there is the member function

*Fmt(self,cnt=1,title=None)*

    *Fmt* is used to control how the tensor is printed with the argument
    *cnt*.  If *cnt=1* the each tensor component is printed on one line.  If
    *fmt=n* :math:`n` tensor components are printed on one line.  If a
    *title* is given then *title=tensor* is printed.  If the usual print
    command is used one tensor component is printed on one line. If *cnt* is
    greater or equal to the number of tensor components then the entire tensor
    is printer on one line.

Examples
========


Algebra
-------

BAC-CAB Formulas
^^^^^^^^^^^^^^^^

This example demonstrates the most general metric tensor

.. math::
  :nowrap:

  \begin{equation}
  g_{ij} = \left [ \begin{array}{cccc} \left ( a\cdot a\right )  & \left ( a\cdot b\right )  & \left ( a\cdot c\right )  & \left ( a\cdot d\right )  \\
  \left ( a\cdot b\right )  & \left ( b\cdot b\right )  & \left ( b\cdot c\right )  & \left ( b\cdot d\right )  \\
  \left ( a\cdot c\right )  & \left ( b\cdot c\right )  & \left ( c\cdot c\right )  & \left ( c\cdot d\right )  \\
  \left ( a\cdot d\right )  & \left ( b\cdot d\right )  & \left ( c\cdot d\right )  & \left ( d\cdot d\right )
  \end{array}\right ]
  \end{equation}

and how the *galgebra* module can be used to verify and expand geometric algebra identities consisting of relations between
the abstract vectors :math:`a`, :math:`b`, :math:`c`, and :math:`d`.

    .. literalinclude:: pysource/baccab.py
        :language: python

The preceeding code block also demonstrates the mapping of *\ **, *^*, and *|* to appropriate latex
symbols.

.. note::

  The :math:`\times` symbol is the commutator product of two multivectors, :math:`A\times B = (AB-BA)/2`.

The output of the code is :download:`BACCAB <./pysource/baccab.pdf>`

Reciprocal Frame
^^^^^^^^^^^^^^^^

The reciprocal frame of vectors with respect to the basis vectors is required
for the evaluation of the geometric dervative.  The following example demonstrates
that for the case of an arbitrary 3-dimensional Euclidian basis the reciprocal
basis vectors are correctly calculated.

    .. literalinclude:: pysource/recp_frame.py
        :language: python

The preceeding code also demonstrated the use of the *\%* directive in
printing a string so that *^* is treated literally and not translated
to *\\wedge*. Note that "%E^{2} =" is printed as ":math:`E^{2} =`"
and not as ":math:`E\wedge {2} =`".  The output is in :download:`Reciprocal Frame <./pysource/recp_frame.pdf>`.

The formulas derived for :math:`E1`, :math:`E2`, :math:`E3`, and :math:`E^{2}` could
also be applied to the numerical calculations of crystal properties.

Lorentz-Transformation
^^^^^^^^^^^^^^^^^^^^^^

A simple physics demonstration of geometric algebra is the derivation of
the Lorentz-Transformation.  In this demonstration a 2-dimensional
Minkowski space is defined and the Lorentz-Transformation is generated
from a rotation of a vector in the Minkowski space using the rotor
:math:`R`.

    .. literalinclude:: pysource/lorentz_trans.py
        :language: python

The preceding code also demonstrates how to use the sympy *subs* functions
to perform the hyperbolic half angle transformation.  The code also shows
the use of both the *#* and *\%* directives in the text string
``r"#%t\bm{\gamma_{t}}+x\bm{\gamma_{x}} = t'\bm{\gamma'_{t}}+x'\bm{\gamma'_{x}} = R\left ( t'\bm{\gamma_{t}}+x'\bm{\gamma_{x}}\right ) R^{\dagger}"``.
Both the *#* and *\%* are needed in this text string for two reasons.  First, the text string contains an *=* sign.  The latex preprocessor
uses this a key to combine the text string with a sympy expression to be printed after the text string.  The *#* is required to inform
the preprocessor that there is no sympy expression to follow.  Second, the *\%* is requires to inform the preprocessor that the text
string is to be displayed in latex math mode and not in text mode (if *#* is present the default latex mode is text mode unless
overridden by the *\%* directive).  The output of the code is in :download:`Lorentz Transformation <./pysource/lorentz_trans.pdf>`.

Linear Transformations
^^^^^^^^^^^^^^^^^^^^^^

Examples of linear transformations are produced by the following code -

    .. literalinclude:: pysource/Ltrans.py
        :language: python

Note the trace of the general basis 2-d case and the adjoint for the Minkowski 4-d case.  Also note
the results of the adjoint test for both cases.  The output is in :download:`Linear Transformation <./pysource/Ltrans.pdf>`.

Calculus
--------

Derivatives in Spherical Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following code shows how to use *galgebra* to perform calculations in spherical coordinates.
The gradient of a scalar function, :math:`f`, the divergence and curl
of a vector function, :math:`A`, and the exterior derivative (curl) of
a bivector function, :math:`B` are calculated.  Note that to get the
standard curl of a 3-dimension function the result is multiplied by
:math:`-I` the negative of the pseudoscalar.

.. note::

    In geometric calculus the operator :math:`\nabla^{2}` is well defined
    on its own as the geometic derivative of the geometric derivative.
    However, if needed we have for the vector function :math:`A` the relations
    (since :math:`\nabla\cdot A` is a scalar it's curl is equal to it's
    geometric derivative and it's divergence is zero) -

    .. math::
        :nowrap:

        \begin{align*}
        \nabla A =& \nabla\wedge A + \nabla\cdot A \\
        \nabla^{2} A =& \nabla\left ( {{\nabla\wedge A}} \right ) + \nabla\left ( {{\nabla\cdot A}} \right ) \\
        \nabla^{2} A =& \nabla\wedge\left ( {{\nabla\wedge A}} \right ) + \nabla\cdot\left ( {{\nabla\wedge A}} \right )
        +\nabla\wedge\left ( {{\nabla\cdot A}} \right ) + \nabla\cdot\left ( {{\nabla\cdot A}} \right ) \\
        \nabla^{2} A =& \nabla\wedge\left ( {{\nabla\wedge A}} \right ) + \left ( {{\nabla\cdot\nabla}} \right ) A
        - \nabla\left ( {{\nabla\cdot A}} \right ) + \nabla\left ( {{\nabla\cdot A}} \right ) \\
        \nabla^{2} A =& \nabla\wedge\nabla\wedge A + \left ( {{\nabla\cdot\nabla}} \right )A
        \end{align*}



    In the derivation we have used that :math:`\nabla\cdot\left ( {{\nabla\wedge A}} \right ) = \left ( {{\nabla\cdot\nabla}} \right )A - \nabla\left ( {{\nabla\cdot A}} \right )`
    which is implicit in the second *BAC-CAB* formula.
    No parenthesis is needed for the geometric curl of the curl (exterior derivative of exterior derivative)
    since the :math:`\wedge` operation is associative unlike the vector curl operator and :math:`\nabla\cdot\nabla` is the usual Laplacian
    operator.

.. literalinclude:: pysource/spherical.py
    :language: python

Results of the code are in :download:`Spherical Coordinates <./pysource/spherical.pdf>`.

Maxwell's Equations
^^^^^^^^^^^^^^^^^^^

The geometric algebra formulation of Maxwell's equations is deomonstrated
with the formalism developed in "Geometric Algebra for Physicists" [Doran]_.
In this formalism the signature of the metric is :math:`(1,-1,-1,-1)` and the
basis vectors are :math:`\gamma_{t}`, :math:`\gamma_{x}`, :math:`\gamma_{y}`,
and :math:`\gamma_{z}`.  The if :math:`\boldsymbol{E}` and :math:`\boldsymbol{B}` are the
normal electric and magnetic field vectors the electric and magnetic
bivectors are given by :math:`E = \boldsymbol{E}\gamma_{t}` and :math:`B = \boldsymbol{B}\gamma_{t}`.
The electromagnetic bivector is then :math:`F = E+IB` where
:math:`I = \gamma_{t}\gamma_{x}\gamma_{y}\gamma_{z}` is the pesudo-scalar
for the Minkowski space.  Note that the electromagnetic bivector is isomorphic
to the electromagnetic tensor.  Then if :math:`J` is the 4-current all of
Maxwell's equations are given by :math:`\boldsymbol{\nabla}F = J`.  For more details
see [Doran]_ chapter 7.

.. literalinclude:: pysource/maxwell.py
    :language: python

Code output is in :download:`Maxwell's Equations <./pysource/maxwell.pdf>`.

Dirac Equation
^^^^^^^^^^^^^^

In [Doran]_ equation 8.89 (page 283) is the geometric algebra formulation of the Dirac equation.  In this equation
:math:`\psi` is an 8-component real spinor which is to say that it is a multivector with sacalar, bivector, and
pseudo-vector components in the space-time geometric algebra (it consists only of even grade components).

.. literalinclude:: pysource/dirac.py
    :language: python

The equations displayed are the partial differential equations for each component of the Dirac equation
in rectangular coordinates we the driver for the equations is the 4-potential :math:`A`.  One utility
of these equations is to setup a numerical solver for the Dirac equation.  The code output is in :download:`Dirac Equation <./pysource/dirac.pdf>`.

Manifolds and Tensors
^^^^^^^^^^^^^^^^^^^^^

The example code for tensors on a manifold is

    .. literalinclude:: pysource/manifold.py
        :language: python

This code uses a unit sphere as the example manifold and shows the
geometric derivative and directional derivative operators.  Additionally,
rank 1 and rank 2 general tensors are defined.  The rank 2 tensor is
contracted and evaluated to show multilinear properties.  The covariant
derivatives of both tensors are caculated.  The results are in :download:`Manifold <./pysource/manifold.pdf>`

.. rubric:: Citations

.. [Doran]  `<http://www.mrao.cam.ac.uk/~cjld1/pages/book.htm>`_
    ``Geometric Algebra for Physicists`` by C. Doran and A. Lasenby, Cambridge
    University Press, 2003.

.. [Hestenes]  `<http://geocalc.clas.asu.edu/html/CA_to_GC.html>`_
    ``Clifford Algebra to Geometric Calculus`` by D.Hestenes and G. Sobczyk, Kluwer
    Academic Publishers, 1984.

.. [Macdonald] '<http://faculty.luther.edu/~macdonal>'_
   ``Linear and Geometric Algebra`` by Alan Macdonald, `<http://www.amazon.com/Alan-Macdonald/e/B004MB2QJQ>`_


.. rubric:: Footnotes

.. [#f0] In this case :math:`D_{B}^{j_{1}\dots j_{n}} = F` and :math:`\partial_{j_{1}\dots j_{n}} = 1`.

.. [#f1] Since :math:`\underline{T}` is linear we do not require :math:`I^{2} = \pm 1`.

.. [#f2] In this case :math:`y` is a vector in the tangent space and not a coordinate vector so that the
         basis vectors are {\em not} a function of :math:`y`.

.. [#f3] Preprocessing do not occur for the Ipython notebook and the string post processing commands
         *\%* and *\#* are not used in this case



