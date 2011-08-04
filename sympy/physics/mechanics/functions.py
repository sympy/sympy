__all__ = ['cross',
           'dot',
           'express',
           'outer',
           'inertia',
           'mechanics_printing',
           'mpprint',
           'mlatex',
           'kinematic_equations',
           'code_output']

from sympy.physics.mechanics.essential import (Vector, Dyad, ReferenceFrame,
                                               MechanicsStrPrinter,
                                               MechanicsPrettyPrinter,
                                               MechanicsLatexPrinter,
                                               dynamicsymbols)
from sympy.physics.mechanics.kane import Kane
from sympy import (sympify, symbols, numbered_symbols, cse, diff, sin, cos)

def cross(vec1, vec2):
    """Cross product convenience wrapper for Vector.cross(): \n"""
    if not isinstance(vec1, (Vector, Dyad)):
        raise TypeError('Cross product is between two vectors')
    return vec1 ^ vec2
cross.__doc__ += Vector.cross.__doc__

def dot(vec1, vec2):
    """Dot product convenience wrapper for Vector.dot(): \n"""
    if not isinstance(vec1, (Vector, Dyad)):
        raise TypeError('Dot product is between two vectors')
    return vec1 & vec2
dot.__doc__ += Vector.dot.__doc__

def express(vec, frame, frame2=None):
    """Express convenience wrapper for Vector.express(): \n"""
    if not isinstance(vec, (Vector, Dyad)):
        raise TypeError('Can only express Vectors')
    if isinstance(vec, Vector):
        return vec.express(frame)
    else:
        return vec.express(frame, frame2)

express.__doc__ += Vector.express.__doc__

def outer(vec1, vec2):
    """Outer prodcut convenience wrapper for Vector.outer():\n"""
    if not isinstance(vec1, Vector):
        raise TypeError('Outer product is between two Vectors')
    return vec1 | vec2
outer.__doc__ += Vector.express.__doc__

def inertia(frame, ixx, iyy, izz, ixy=0, iyz=0, izx=0):
    """Simple way to create inertia Dyad object.

    If you don't know what a Dyad is, just treat this like the inertia tensor.
    Then, do the easy thing and define it in a body-fixed frame.

    Parameters
    ==========
    frame : ReferenceFrame
        The frame the inertia is defined in
    ixx : Sympifyable
        the xx element in the inertia dyad
    iyy : Sympifyable
        the yy element in the inertia dyad
    izz : Sympifyable
        the zz element in the inertia dyad
    ixy : Sympifyable
        the xy element in the inertia dyad
    iyz : Sympifyable
        the yz element in the inertia dyad
    izx : Sympifyable
        the zx element in the inertia dyad

    Examples
    ========

    >>> from sympy.physics.mechanics import ReferenceFrame, inertia
    >>> N = ReferenceFrame('N')
    >>> inertia(N, 1, 2, 3)
    (N.x|N.x) + 2*(N.y|N.y) + 3*(N.z|N.z)

    """

    if not isinstance(frame, ReferenceFrame):
        raise TypeError('Need to define the inertia in a frame')
    ol  = sympify(ixx) * (frame.x | frame.x)
    ol += sympify(ixy) * (frame.x | frame.y)
    ol += sympify(izx) * (frame.x | frame.z)
    ol += sympify(ixy) * (frame.y | frame.x)
    ol += sympify(iyy) * (frame.y | frame.y)
    ol += sympify(iyz) * (frame.y | frame.z)
    ol += sympify(izx) * (frame.z | frame.x)
    ol += sympify(iyz) * (frame.z | frame.y)
    ol += sympify(izz) * (frame.z | frame.z)
    return ol

def mechanics_printing():
    """Sets up interactive printing for mechanics' derivatives.

    The main benefit of this is for printing of time derivatives;
    instead of displaying as Derivative(f(t),t), it will display f'
    This is only actually needed for when derivatives are present and are not
    in a physics.mechanics object.

    Examples
    ========

    >>> from sympy import Function, Symbol, diff
    >>> from sympy.physics.mechanics import mechanics_printing
    >>> f = Function('f')
    >>> t = Symbol('t')
    >>> x = Symbol('x')
    >>> diff(f(t), t)
    Derivative(f(t), t)
    >>> mechanics_printing()
    displayhook set
    >>> diff(f(t), t)
    f'
    >>> diff(f(x), x)
    Derivative(f(x), x)

    """

    import sys
    sys.displayhook = mprint
    return 'displayhook set'

def mprint(expr, **settings):
    """Function for printing of expressions generated in mechanics. """
    pr = MechanicsStrPrinter(settings)
    print pr.doprint(expr)

def mpprint(expr, **settings):
    """Function for pretty printing of expressions generated in mechanics.

    Mainly used for expressions not inside a vector; the output of running
    scripts and generating equations of motion.

    """

    mp = MechanicsPrettyPrinter(settings)
    print mp.doprint(expr)

def mlatex(expr, **settings):
    """Function for printing latex representation of mechanics objects.

    For latex representation of Vectors, Dyads, and dynamicsymbols.

    Examples
    ========

    >>> from sympy.physics.mechanics import mlatex, ReferenceFrame
    >>> N = ReferenceFrame('N')
    >>> mlatex(N.x + N.y)
    '\\\\mathbf{\\\\hat{n_{x}}} + \\\\mathbf{\\\\hat{n_{y}}}'

    """

    return MechanicsLatexPrinter(settings).doprint(expr)

def kinematic_equations(speeds, coords, rot_type, rot_order=''):
    """Gives equations relating the qdot's to u's for a rotation type.

    Supply rotation type and order as in orient. Speeds are assumed to be
    body-fixed; if we are defining the orientation of B in A using by rot_type,
    the angular velocity of B in A is assumed to be in the form: speed[0]*B.x +
    speed[1]*B.y + speed[2]*B.z

    Parameters
    ==========
    speeds : list of length 3
        The body fixed angular velocity measure numbers.
    coords : list of length 3 or 4
        The coordinates used to define the orientation of the two frames.
    rot_type : str
        The type of rotation used to create the equations. Body, Space, or
        Quaternion only
    rot_order : str
        If applicable, the order of a series of rotations.

    Examples
    ========

    >>> from sympy.physics.mechanics import dynamicsymbols, kinematic_equations
    >>> from sympy.physics.mechanics import mechanics_printing
    >>> mechanics_printing()
    displayhook set
    >>> u1, u2, u3 = dynamicsymbols('u1 u2 u3')
    >>> q1, q2, q3 = dynamicsymbols('q1 q2 q3')
    >>> kinematic_equations([u1,u2,u3],[q1,q2,q3],'body','313')
    [-(u1*sin(q3) + u2*cos(q3))/sin(q2) + q1', -u1*cos(q3) + u2*sin(q3) + q2', (u1*sin(q3) + u2*cos(q3))*cos(q2)/sin(q2) - u3 + q3']

    """

    # Code below is checking and sanitizing input
    approved_orders = ('123', '231', '312', '132', '213', '321', '121', '131',
                       '212', '232', '313', '323', '1', '2', '3', '')
    rot_order = str(rot_order).upper() # Now we need to make sure XYZ = 123
    rot_type  = rot_type.upper()
    rot_order = [i.replace('X', '1') for i in rot_order]
    rot_order = [i.replace('Y', '2') for i in rot_order]
    rot_order = [i.replace('Z', '3') for i in rot_order]
    rot_order = ''.join(rot_order)

    if not isinstance(speeds,(list, tuple)):
        raise TypeError('Need to supply speeds in a list')
    if len(speeds) != 3:
        raise TypeError('Need to supply 3 body-fixed speeds')
    if not isinstance(coords,(list, tuple)):
        raise TypeError('Need to supply coordinates in a list')
    if rot_type.lower() in ['body', 'space']:
        if rot_order not in approved_orders:
            raise ValueError('Not an acceptable rotation order')
        if len(coords) != 3:
            raise ValueError('Need 3 coordinates for body or space')
        # Actual hard-coded kinematic differential equations
        q1, q2, q3 = coords
        q1d, q2d, q3d = [diff(i, dynamicsymbols._t) for i in coords]
        w1, w2, w3 = speeds
        s1, s2, s3 = [sin(q1), sin(q2), sin(q3)]
        c1, c2, c3 = [cos(q1), cos(q2), cos(q3)]
        if rot_type.lower() == 'body':
            if rot_order == '123':
                return [q1d - (w1 * c1 - w2 * s3) / c2, q2d - w1 * s3 - w2 * c3,
                        q3d - (-w1 * c3 + w2 * s3) * s2 / c2 - w3]
            if rot_order == '231':
                return [q1d - (w2 * c3 - w3 * s3) / c2, q2d - w2 * s3 - w3 *
                        c3, q3d - w1 - (- w2 * c3 + w3 * s3) * s2 / c2]
            if rot_order == '312':
                return [q1d - (-w1 * s3 + w3 * c3) / c3, q2d - w1 * c3 - w3 *
                        s3, q3d - (w1 * s3 - w3 * c3) * s2 / c2 - w2]
            if rot_order == '132':
                return [q1d - (w1 * c3 + w3 * s3) / c2, q2d + w1 * s3 - w3 *
                        c3, q3d - (w1 * c3 + w3 * s3) * s2 / c2 - w2]
            if rot_order == '213':
                return [q1d - (q1 * s3 + w2 * c3) / c2, q2d - w1 * c3 + w2 *
                        s3, q3d - (w1 * s3 + w2 * c3) * s2 / c2 - w3]
            if rot_order == '321':
                return [q1d - (w2 * s3 + w3 * c3) / c2, q2d - w2 * c3 + w3 *
                        s3, q3d - w1 - (w2 * s3 + w3 * c3) * s2 / c2]
            if rot_order == '121':
                return [q1d - (w2 * s3 + w3 * c3) / s2, q2d - w2 * c3 + w3 *
                        s3, q3d - w1 + (w2 * s3 + w3 * c3) * c2 / s2]
            if rot_order == '131':
                return [q1d - (-w2 * c3 + w3 * s3) / s2, q2d - w2 * s3 - w3 *
                        c3, q3d - w1 - (w2 * c3 - w3 * s3) * c2 / s2]
            if rot_order == '212':
                return [q1d - (w1 * s3 - w3 * c3) / s2, q2d - w1 * c3 - w3 *
                        s3, q3d - (-w1 * s3 + w3 * c3) * c2 / s2 - w2]
            if rot_order == '232':
                return [q1d - (w1 * c3 + w3 * s3) / s2, q2d + w1 * s3 - w3 *
                        c3, q3d + (w1 * c3 + w3 * s3) * c2 / s2 - w2]
            if rot_order == '313':
                return [q1d - (w1 * s3 + w2 * c3) / s2, q2d - w1 * c3 + w2 *
                        s3, q3d + (w1 * s3 + w2 * c3) * c2 / s2 - w3]
            if rot_order == '323':
                return [q1d - (-w1 * c3 + w2 * s3) / s2, q2d - w1 * s3 - w2 *
                        c3, q3d - (w1 * c3 - w2 * s3) * c2 / s2 - w3]
        if rot_type.lower() == 'space':
            if rot_order == '123':
                return [q1d - w1 - (w2 * s1 + w3 * c1) * s2 / c2, q2d - w2 *
                        c1 + w3 * s1, q3d - (w2 * s1 + w3 * c1) / c2]
            if rot_order == '231':
                return [q1d - (w1 * c1 + w3 * s1) * s2 / c2 - w2, q2d + w1 *
                        s1 - w3 * c1, q3d - (w1 * c1 + w3 * s1) / c2]
            if rot_order == '312':
                return [q1d - (w1 * s1 + w2 * c1) * s2 / c2 - w3, q2d - w1 *
                        c1 + w2 * s1, q3d - (w1 * s1 + w2 * c1) / c2]
            if rot_order == '132':
                return [q1d - w1 - (-w2 * c1 + w3 * s1) * s2 / c2, q2d - w2 *
                        s1 + w3 * c1, q3d - (w2 * c1 - w3 * s1) / c2]
            if rot_order == '213':
                return [q1d - (w1 * s1 - w3 * c1) * s2 / c2 - w2, q2d - w1 *
                        c1 + w3 * s1, q3d - (-w1 * s1 + w3 * c1) / c2]
            if rot_order == '321':
                return [q1d - (-w1 * c1 + w2 * s1) * s2 / c2 - w3, q2d - w1 *
                        s1 - w2 * c1, q3d - (w1 * c1 - w2 * s1) / c2]
            if rot_order == '121':
                return [q1d - w1 + (w2 * s1 + w3 * c2) * c2 / s2, q2d - w2 *
                        c1 + w3 * s1, q3d - (w2 * s1 + w3 * c1) / s2]
            if rot_order == '131':
                return [q1d - w1 - (w2 * c1 - w3 * s1) * c2 / s2, q2d - w2 *
                        s1 - w3 * c1, q3d - (-w2 * c1 + w3 * s1) / s2]
            if rot_order == '212':
                return [q1d - (-w1 * s1 + w3 * c1) * c2 / s2 - w2, q2d - w1 *
                        c1 + w3 * s1, q3d - (w1 * s1 - w3 * c1) / s2]
            if rot_order == '232':
                return [q1d + (w1 * c1 + w3 * s1) * c2 / s2 - w2, q2d + w1 *
                        s1 - w3 * c1, q3d - (w1 * c1 + w3 * s1) / s2]
            if rot_order == '313':
                return [q1d + (w1 * s1 + w2 * c1) * c2 / s2 - w3, q2d - w1 *
                        c1 + w2 * s1, q3d - (w1 * s1 + w2 * c1) / s2]
            if rot_order == '323':
                return [q1d - (w1 * c1 - w2 * s1) * c2 / s2 - w3, q2d - w1 *
                        s1 - w2 * c1, q3d - (-w1 * c1 + w2 * s1) / s2]
    elif rot_type.lower() == 'quaternion':
        if rot_order != '':
            raise ValueError('Cannot have rotation order for quaternion')
        if len(coords) != 4:
            raise ValueError('Need 4 coordinates for quaternion')
        # Actual hard-coded kinematic differential equations
        e0, e1, e2, e3 = coords
        w = Matrix(speeds + [0])
        E = Matrix([[e0, -e3, e2, e1], [e3, e0, -e1, e2], [-e2, e1, e0, e3],
            [-e1, -e2, -e3, e0]])
        edots = Matrix([diff(i, dynamicsymbols._t) for i in coords])
        return (edots - 0.5 * w * E.t)[:]
    else:
        raise ValueError('Not an approved rotation type for this function')

def code_output(typ, KM, fname, params=[]):
    """Function for generating files with code to integrate.

    Parameters
    ==========
    typ : str
        Type of output file
    KM : Kane
        Kane object to get equations from
    fname : str
        The filename to write to code to, no extension please
    params : tuple
        A list of symbols needed by the rhs function

    """

    def _py_code_output(KM, fname, params):
        """Function for generating code for numerical integration.

        This is for SciPy's odeint function and uses NumPy.

        """

        mstr = MechanicsStrPrinter().doprint
        f = open("./%s.py" % fname, 'w')
        f.write('from scipy.integrate import odeint\n')
        f.write('from numpy import *\n')
        f.write('acos = arccos\n')
        f.write('asin = arcsin\n')
        f.write('atan = arctan\n\n')
        f.write('def rhs(y, t, params):\n')
        qsus = KM._q + KM._u
        kdd = KM.kindiffdict()
        MM = KM.mass_matrix_full.subs(kdd)
        fo = KM.forcing_full.subs(kdd)
        others = list(set(KM._find_dynamicsymbols(MM, qsus) +
                          KM._find_dynamicsymbols(fo, qsus)))
        if params == []:
            syms = list(set(KM._find_othersymbols(MM) +
                            KM._find_othersymbols(fo)))
        else:
            syms = params
        others2 = {}
        for i in others:
            stri = mstr(i)
            others2.update({i: symbols(stri.replace('\'','p'))})

        l1 = '    ' + mstr(qsus)[1:-1] + ' = y\n'
        f.write(l1)
        l2 = '    ' + mstr(syms)[1:-1] + ' = params\n'
        f.write(l2)

        for i in others2.keys():
            outl = '    ' + str(i) + ' = ?\n'
            f.write(outl)
        fo = fo.subs(others2)
        MM = MM.subs(others2)
        f.write('    M = zeros(%s)\n' % str(MM.shape))
        f.write('    fo = zeros(%s)\n' % str(fo.shape))
        ma, mb = MM.shape
        fa, fb = fo.shape
        for i in range(ma):
            for j in range(mb):
                f.write('    M[%s][%s] = %s\n' %
                        (i, j, mstr(MM.extract([i],[j]))[1:-1]))
        for i in range(fa):
            f.write('    fo[%s] = %s\n' % (i, mstr(fo[i])))
        f.write('    return linalg.solve(M, fo).T.flatten().tolist()\n\n')
        for i in syms:
            f.write('%s = ?\n' % i)
        f.write('params = (%s)\n' % (mstr(syms)[1:-1]))

    def _py_cse_code_output(KM, fname, params):
        """Function for generating code for numerical integration.

        This is for SciPy's odeint function and uses NumPy. Uses SymPy's cse() to
        generate shorter expressions, but takes longer.

        """

        mstr = MechanicsStrPrinter().doprint
        f = open("./%s.py" % fname, 'w')
        f.write('from scipy.integrate import odeint\n')
        f.write('from numpy import *\n')
        f.write('acos = arccos\n')
        f.write('asin = arcsin\n')
        f.write('atan = arctan\n\n')
        f.write('def rhs(y, t, params):\n')
        qsus = KM._q + KM._u
        kdd = KM.kindiffdict()
        MM = KM.mass_matrix_full.subs(kdd)
        fo = KM.forcing_full.subs(kdd)
        others = list(set(KM._find_dynamicsymbols(MM, qsus) +
                          KM._find_dynamicsymbols(fo, qsus)))
        if params == []:
            syms = list(set(KM._find_othersymbols(MM) +
                            KM._find_othersymbols(fo)))
        else:
            syms = params
        others2 = {}
        for i in others:
            stri = mstr(i)
            others2.update({i: symbols(stri.replace('\'','p'))})

        l1 = '    ' + mstr(qsus)[1:-1] + ' = y\n'
        f.write(l1)
        l2 = '    ' + mstr(syms)[1:-1] + ' = params\n'
        f.write(l2)

        for i in others2.keys():
            outl = '    ' + str(i) + ' = ?\n'
            f.write(outl)
        fo = fo.subs(others2)
        MM = MM.subs(others2)
        f.write('    M = zeros(%s)\n' % str(MM.shape))
        f.write('    fo = zeros(%s)\n' % str(fo.shape))
        ma, mb = MM.shape
        fa, fb = fo.shape
        # cse stuff
        z = numbered_symbols('z')
        exprs = fo[:] + MM[:]
        cseout = cse(exprs, z)
        mats = cseout[1]
        zs = cseout[0]
        fo = mats[:fa]
        MM = mats[fa:]
        for i in zs:
            f.write('    %s = %s\n' % (mstr(i[0]), mstr(i[1])))
        # end cse stuff
        for i in range(ma):
            for j in range(mb):
                # line below also formatted differently for cse
                f.write('    M[%s][%s] = %s\n' %
                        (i, j, mstr(MM[i * ma + j])))
        for i in range(fa):
            f.write('    fo[%s] = %s\n' % (i, mstr(fo[i])))
        f.write('    return linalg.solve(M, fo).T.flatten().tolist()\n\n')
        for i in syms:
            f.write('%s = ?\n' % i)
        f.write('params = (%s)\n' % (mstr(syms)[1:-1]))

    def _matlab_code_output(KM, fname, params):
        """Function for generating code for numerical integration.

        This is for MATLAB's ode45 function.

        """

        mstr = MechanicsStrPrinter().doprint
        qsus = KM._q + KM._u
        kdd = KM.kindiffdict()
        MM = KM.mass_matrix_full.subs(kdd)
        fo = KM.forcing_full.subs(kdd)
        others = list(set(KM._find_dynamicsymbols(MM, qsus) +
                          KM._find_dynamicsymbols(fo, qsus)))
        if params == []:
            syms = list(set(KM._find_othersymbols(MM) +
                            KM._find_othersymbols(fo)))
        else:
            syms = params
        others2 = {}
        for i in others:
            stri = mstr(i)
            others2.update({i: symbols(stri.replace('\'','p'))})

        syms_str = mstr(syms)[1:-1]

        f = open("./%s.m" % fname, 'w')
        f.write('function ydot = %s(t, y, pars):\n' % fname)
        for i, v in enumerate(qsus):
            f.write('%s = y(%s);\n' % (mstr(v), i + 1))
        for i in others2.keys():
            f.write(mstr(i) + ' = ?;\n')
        for i, v in enumerate(syms):
            f.write('%s = pars(%s)' % (mstr(v), i + 1))
        fo = fo.subs(others2)
        MM = MM.subs(others2)
        f.write('M = zeros%s;\n' % str(MM.shape))
        f.write('fo = zeros%s;\n' % str(fo.shape))
        ma, mb = MM.shape
        fa, fb = fo.shape
        for i in range(ma):
            for j in range(mb):
                f.write('M(%s, %s) = %s;\n' %
                    (i + 1, j + 1, (mstr(MM.extract([i],
                        [j]))[1:-1]).replace('**', '^')))
        for i in range(fa):
            f.write('fo(%s) = %s;\n' % (i + 1, mstr(fo[i]).replace('**', '^')))
        f.write('ydot = M\\fo;\nend\n')
        for i in syms:
            f.write('%s = ?\n' % i)
        f.write('%' + 'Order of symbols for rhs function: %s\n' % (syms_str))
        f.write('% ' + 'might be helpful: pars = [%s];' % syms_str)

    def _matlab_cse_code_output(KM, fname, params):
        """Function for generating code for numerical integration.

        This is for MATLAB's ode45 function.

        """

        mstr = MechanicsStrPrinter().doprint
        qsus = KM._q + KM._u
        kdd = KM.kindiffdict()
        MM = KM.mass_matrix_full.subs(kdd)
        fo = KM.forcing_full.subs(kdd)
        others = list(set(KM._find_dynamicsymbols(MM, qsus) +
                          KM._find_dynamicsymbols(fo, qsus)))
        if params == []:
            syms = list(set(KM._find_othersymbols(MM) +
                            KM._find_othersymbols(fo)))
        else:
            syms = params
        others2 = {}
        for i in others:
            stri = mstr(i)
            others2.update({i: symbols(stri.replace('\'','p'))})

        syms_str = mstr(syms)[1:-1]

        f = open("./%s.m" % fname, 'w')
        f.write('function ydot = %s(t, y, pars):\n' % fname)
        for i, v in enumerate(qsus):
            f.write('%s = y(%s);\n' % (mstr(v), i + 1))
        for i in others2.keys():
            f.write(mstr(i) + ' = ?;\n')
        for i, v in enumerate(syms):
            f.write('%s = pars(%s)' % (mstr(v), i + 1))
        fo = fo.subs(others2)
        MM = MM.subs(others2)
        f.write('M = zeros%s;\n' % str(MM.shape))
        f.write('fo = zeros%s;\n' % str(fo.shape))
        ma, mb = MM.shape
        fa, fb = fo.shape

        # cse stuff

        z = numbered_symbols('z')
        exprs = fo[:] + MM[:]
        cseout = cse(exprs, z)
        mats = cseout[1]
        zs = cseout[0]
        fo = mats[:fa]
        MM = mats[fa:]
        for i in zs:
            f.write('%s = %s\n' % (mstr(i[0]), mstr(i[1]).replace('**', '^')))
        # end cse stuff

        for i in range(ma):
            for j in range(mb):
                f.write('M[%s][%s] = %s\n' %
                        (i, j, mstr(MM[i * ma + j]).replace('**', '^')))
        for i in range(fa):
            f.write('fo(%s) = %s;\n' % (i + 1, mstr(fo[i]).replace('**', '^')))
        f.write('ydot = M\\fo;\nend\n')
        for i in syms:
            f.write('%s = ?\n' % i)
        f.write('% ' + 'Order of symbols for rhs function: %s\n' % (syms_str))
        f.write('% ' + 'might be helpful: pars = [%s];' % syms_str)

    if typ == 'py':
        _py_code_output(KM, fname, params)
    elif typ == 'py_cse':
        _py_cse_code_output(KM, fname, params)
    elif typ == 'matlab':
        _matlab_code_output(KM, fname, params)
    elif typ == 'matlab_cse':
        _matlab_cse_code_output(KM, fname, params)
    else:
        raise ValueError('Unsupported code output type')

if __name__ == "__main__":
    import doctest
    doctest.testmod()
