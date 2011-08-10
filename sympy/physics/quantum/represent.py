"""Logic for representing operators in state in various bases.

TODO:
* Get represent working with continuous hilbert spaces.
* Document default basis functionality.
* Fix representations of Pow for continuous bases (see note in code)
"""


from sympy import Add, Expr, I, integrate, Mul, oo, Pow, Symbol
from sympy.functions import conjugate, DiracDelta
from sympy.utilities import default_sort_key
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.commutator import Commutator
from sympy.physics.quantum.anticommutator import AntiCommutator
from sympy.physics.quantum.innerproduct import InnerProduct
from sympy.physics.quantum.qexpr import QExpr
from sympy.physics.quantum.tensorproduct import TensorProduct
from sympy.physics.quantum.matrixutils import flatten_scalar
from sympy.physics.quantum.state import KetBase, BraBase, StateBase, Wavefunction
from sympy.physics.quantum.operator import (
    Operator, HermitianOperator, OuterProduct, DifferentialOperator
)
from sympy.physics.quantum.qapply import qapply
from sympy.physics.quantum.operatorset import operators_to_state, state_to_operators

__all__ = [
    'represent',
    'innerproduct_helper',
    'expectation_helper',
    'integrate_result',
    'get_basis_state',
    'enumerate_states'
]

#-----------------------------------------------------------------------------
# Represent
#-----------------------------------------------------------------------------

def _sympy_to_scalar(e):
    """Convert from a sympy scalar to a Python scalar."""
    if isinstance(e, Expr):
        if e.is_Integer:
            return int(e)
        elif e.is_Float:
            return float(e)
        elif e.is_Rational:
            return float(e)
        elif e.is_Number or e.is_NumberSymbol or e == I:
            return complex(e)
    raise TypeError('Expected number, got: %r' % e)


def represent(expr, **options):
    """Represent the quantum expression in the given basis.

    In quantum mechanics abstract states and operators can be represented in
    various basis sets. Under this operation the follow transforms happen:

    * Ket -> column vector or function
    * Bra -> row vector of function
    * Operator -> matrix or differential operator

    This function is the top-level interface for this action.

    This function walks the sympy expression tree looking for ``QExpr``
    instances that have a ``_represent`` method. This method is then called
    and the object is replaced by the representation returned by this method.
    By default, the ``_represent`` method will dispatch to other methods
    that handle the representation logic for a particular basis set. The
    naming convention for these methods is the following::

        def _represent_FooBasis(self, e, basis, **options)

    This function will have the logic for representing instances of its class
    in the basis set having a class named ``FooBasis``.

    The function also keeps track of an 'index' option, which
    represents the index of dummy kets that are inserted into the
    expression. This index defaults to 1 at the start of the
    representation of the expression. If the user specifies an index,
    represent will start counting indices at the passed index.

    Parameters
    ==========

    expr  : Expr
        The expression to represent.
    basis : Operator, basis set
        An object that contains the information about the basis set. If an
        operator is used, the basis is assumed to be the orthonormal
        eigenvectors of that operator. In general though, the basis argument
        can be any object that contains the basis set information.
    options : dict
        Key/value pairs of options that are passed to the underlying method
        that does finds the representation. These options can be used to
        control how the representation is done. For example, this is where
        the size of the basis set would be set.

    Returns
    =======

    e : Expr
        The sympy expression of the represented quantum expression.

    Examples
    ========

    Here we subclass ``Operator`` and ``Ket`` to create the z-spin operator
    and its spin 1/2 up eigenstate. By definining the ``_represent_SzOp``
    method, the ket can be represented in the z-spin basis.

        >>> from sympy.physics.quantum import Operator, represent, Ket
        >>> from sympy import Matrix

        >>> class SzUpKet(Ket):
        ...     def _represent_SzUpKet(self, basis, **options):
        ...         return Matrix([1,0])
        ...
        >>> class SzOp(Operator):
        ...     pass
        ...
        >>> sz = SzOp('Sz')
        >>> up = SzUpKet('up')
        >>> represent(up, basis=up)
        [1]
        [0]

    Here we see an example of representations in a continuous
    basis. We see that the result of representing various combinations
    of cartesian position operators and kets give us continuous
    expressions involving DiracDelta functions.

        >>> from sympy.physics.quantum.cartesian import XOp, XKet, XBra
        >>> X = XOp()
        >>> x = XKet()
        >>> y = XBra('y')
        >>> represent(X*x)
        DiracDelta(x - x_2)*Wavefunction(x, x)
        >>> represent(X*x*y)
        DiracDelta(x - x_3)*DiracDelta(-x_1 + y)*Wavefunction(x, x)
    """

    rep_expr, basis = _represent_helper(expr, **options)

    return rep_expr

def _represent_helper(expr, **options):
    """A recursively called helper function for represent. Handles all of the
    logic described in the represent docstring, but it returns the basis used
    for representing in addition to the represented expression"""

    format = options.get('format', 'sympy')

    if not 'basis' in options or not isinstance(options['basis'], StateBase):
        basis = get_basis_state(expr, **options)
        options['basis'] = basis

    basis = options['basis']

    if isinstance(expr, QExpr) and not isinstance(expr, OuterProduct):
        return (expr._represent(**options), basis)
    elif isinstance(expr, Add):
        result = represent(expr.args[0], **options)
        for args in expr.args[1:]:
            # scipy.sparse doesn't support += so we use plain = here.
            result = result + represent(args, **options)
        return (result, basis)
    elif isinstance(expr, Pow):
        #NOTE: This will not currently work with continuous operators
        #For example, represent(X**2) != represent(X)**2
        base, exp = expr.as_base_exp()
        if format == 'numpy' or format == 'scipy.sparse':
            exp = _sympy_to_scalar(exp)
        return (represent(base, **options)**exp, basis)
    elif isinstance(expr, TensorProduct):
        new_args = [represent(arg, **options) for arg in expr.args]
        return (TensorProduct(*new_args), basis)
    elif isinstance(expr, Dagger):
        return (Dagger(represent(expr.args[0], **options)), basis)
    elif isinstance(expr, Commutator):
        A = expr.args[0]
        B = expr.args[1]
        return (represent(A*B, **options) - represent(B*A, **options), basis)
    elif isinstance(expr, AntiCommutator):
        A = expr.args[0]
        B = expr.args[1]
        return (represent(A*B, **options) + represent(B*A, **options), basis)
    elif isinstance(expr, InnerProduct):
        return (represent(Mul(expr.bra,expr.ket), **options), basis)
    elif not (isinstance(expr, Mul) or isinstance(expr, OuterProduct)):
        # For numpy and scipy.sparse, we can only handle numerical prefactors.
        if format == 'numpy' or format == 'scipy.sparse':
            return (_sympy_to_scalar(expr), basis)
        return (expr, basis)

    if not (isinstance(expr, Mul) or isinstance(expr, OuterProduct)):
        raise TypeError('Mul expected, got: %r' % expr)

    if "index" in options:
        options["index"] += 1
    else:
        options["index"] = 1

    if not "unities" in options:
        options["unities"] = []

    if not "delta_unities" in options:
        options["delta_unities"] = []

    result = represent(expr.args[-1], **options)
    last_arg = expr.args[-1]

    for arg in reversed(expr.args[:-1]):
        unity = None
        if isinstance(last_arg, Operator):
            options["index"] += 1
            unity = options["index"]
        elif isinstance(last_arg, BraBase) and isinstance(arg, KetBase):
            options["index"] += 1
        elif isinstance(last_arg, KetBase) and isinstance(arg, Operator):
            unity = options["index"]
        elif isinstance(last_arg, KetBase) and isinstance(arg, BraBase):
            unity = options["index"]

        tmp_result = represent(arg, **options)
        if unity is not None:
            if isinstance(tmp_result, Expr) and tmp_result.has(DiracDelta):
                options["delta_unities"].append(unity)
            else:
                options["unities"].append(unity)

        result = tmp_result*result
        last_arg = arg

    # All three matrix formats create 1 by 1 matrices when inner products of
    # vectors are taken. In these cases, we simply return a scalar.
    result = flatten_scalar(result)

    should_apply = options.pop('qapply', True)
    should_integrate = options.pop('integrate', True)
    should_wrap = options.pop('wrap_wf', True)

    if should_integrate:
        result = integrate_result(expr, result, delta=True, **options)

    if should_apply:
        result = do_qapply(result)

    unwrapped_res, unwrapped_vars = _unwrap_wf(result)
    vars_beforeint = set(unwrapped_vars)

    if should_integrate:
        result = integrate_result(expr, unwrapped_res, **options)

    if not should_wrap:
        return (result, basis)

    if len(unwrapped_vars) != 0:
        vars_afterint = result.free_symbols
        intersect = list(vars_beforeint.intersection(vars_afterint))
        intersect.sort(key=default_sort_key)

        if len(intersect) == 0 and len(vars_afterint) != 0:
            #If none of the original variables are left in the integrated
            #expression, then make the Wf variable the first coordinate that we
            #inserted
            vars_afterint = list(vars_afterint)
            vars_afterint.sort(key=default_sort_key)
            intersect.append(vars_afterint[0])

        unwrapped_vars = intersect

        if len(unwrapped_vars) != 0:
            result = _rewrap_wf(result, unwrapped_vars, **options)

    return (result, basis)

def innerproduct_helper(expr, **options):
    """
    Returns an innerproduct like representation (e.g. <x'|x>) for the given
    state. This is meant to be a helper function for the internal _represent
    methods when they want to form a standard <x'|x> type reprsentation. As this
    is a helper function, the basis must be specified in the options or an error
    will be raised. We also expect that any basis passed to this will be a state
    instance, as this is the convention specified for representations.

    Attempts to calculate inner product with a bra from the specified
    basis. Should only be passed an instance of KetBase or BraBase, and the
    basis specified in option must be a State instance.

    If a ``wrap_wavefunction`` option is passed and is True, the function will
    wrap the result in a Wavefunction object before returning. This option
    should be used if the basis is continuous.

    Parameters
    ==========

    expr : KetBase or BraBase
        The expression to be represented

    Examples
    ========

    >>> from sympy.physics.quantum.represent import innerproduct_helper
    >>> from sympy.physics.quantum.cartesian import XOp, XKet, PxOp, PxKet
    >>> innerproduct_helper(XKet(), basis=XKet())
    DiracDelta(x - x_1)
    >>> rep_innerproduct(XKet(), basis=PxOp())
    sqrt(2)*exp(-I*px_1*x/hbar)/(2*sqrt(hbar)*sqrt(pi))
    >>> rep_innerproduct(PxKet(), basis=XOp())
    sqrt(2)*exp(I*px*x_1/hbar)/(2*sqrt(hbar)*sqrt(pi))

    """

    if not isinstance(expr, (KetBase, BraBase)):
        raise TypeError("expr passed is not a Bra or Ket")

    basis = options.pop("basis", None)

    if basis is None:
        raise NotImplementedError("Basis not specified!")

    if not isinstance(basis, StateBase):
        raise NotImplementedError("Specified basis is not a State!")

    if not "index" in options:
        options["index"] = 1

    basis_kets = enumerate_states(basis, options["index"], 2)

    if isinstance(expr, BraBase):
        bra = expr
        ket =  (basis_kets[1] if basis_kets[0].dual == expr else basis_kets[0])
    else:
        bra = (basis_kets[1].dual if basis_kets[0] == expr else basis_kets[0].dual)
        ket = expr

    prod = InnerProduct(bra, ket)
    result = prod.doit()

    format = options.get('format', 'sympy')
    ret = expr._format_represent(result, format)

    wrap = options.pop('wrap_wavefunction', False)
    delta = options.pop('keep_delta', False)
    #TODO: Insert proper limits
    if wrap:
        # The label of the original expr will be used as the free coordinate
        return Wavefunction(ret, *expr.label, keep_delta=delta, **options)
    else:
        return ret

def expectation_helper(expr, **options):
    """
    Returns an ``<x'|A|x>`` type representation for the given operator.

    This is meant to be a helper function for the internal _represent methods
    when they want to form a standard <x'|A|x> type reprsentation. As it is a
    helper function, the basis must be specified in the options.  We also expect
    that any basis passed to this will be a state instance, as this is the
    convention specified for representations.

    If a ``wrap_wavefunction`` option is passed and is True, the function will
    wrap the result in a Wavefunction object before returning. This option
    should be used if the basis is continuous.

    Parameters
    ==========

    expr : Operator
        Operator to be represented in the specified basis

    Examples
    ========

    >>> from sympy.physics.quantum.cartesian import XOp, XKet, PxOp, PxKet
    >>> from sympy.physics.quantum.represent import expectation_helper
    >>> expectation_helper(XOp(), basis=XKet())
    x_1*DiracDelta(x_1 - x_2)
    >>> expectation_helper(XOp(), basis=PxKet())
    <px_2|*X*|px_1>
    >>> expectation_helper(XOp(), basis=PxKet())
    <px_2|*X*|px_1>

    """

    if not "index" in options:
        options["index"] = 1

    if not isinstance(expr, Operator):
        raise TypeError("The passed expression is not an operator")

    basis_state = options.pop('basis', None)

    if basis_state is None or not isinstance(basis_state, StateBase):
        raise NotImplementedError("Could not get basis kets for this operator")

    basis_kets = enumerate_states(basis_state, options["index"], 2)

    bra = basis_kets[1].dual
    ket = basis_kets[0]

    ret = qapply(bra*expr*ket)

    wrap = options.pop('wrap_wavefunction', False)
    delta = options.pop('keep_delta', False)
    #TODO: Insert proper limits
    if wrap:
        # We consider the ket coordinate to be the free variable, and
        # the bra coordinate a dummy constant
        return Wavefunction(ret, *ket.label, keep_delta=delta, **options)
    else:
        return ret

def integrate_result(orig_expr, result, **options):
    """
    Returns the result of integrating over any unities ``(|x><x|)`` in
    the given expression. Intended for integrating over the result of
    representations in continuous bases.

    This function integrates over any unities that may have been
    inserted into the quantum expression and returns the result.
    It uses the interval of the Hilbert space of the basis state
    passed to it in order to figure out the limits of integration.
    The unities option must be
    specified for this to work.

    Note: This is mostly used internally by represent(). Examples are
    given merely to show the use cases.

    Parameters
    ==========

    orig_expr : quantum expression
        The original expression which was to be represented

    result: Expr
        The resulting representation that we wish to integrate over

    Examples
    ========

    >>> from sympy import symbols, DiracDelta
    >>> from sympy.physics.quantum.represent import integrate_result
    >>> from sympy.physics.quantum.cartesian import XOp, XKet
    >>> x_ket = XKet()
    >>> X_op = XOp()
    >>> x, x_1, x_2 = symbols('x, x_1, x_2')
    >>> integrate_result(X_op*x_ket, x*DiracDelta(x-x_1)*DiracDelta(x_1-x_2))
    x*DiracDelta(x - x_1)*DiracDelta(x_1 - x_2)
    >>> integrate_result(X_op*x_ket, x*DiracDelta(x-x_1)*DiracDelta(x_1-x_2), basis=XKet, unities=[1])
    x*DiracDelta(x - x_2)

    """
    if not isinstance(result, Expr):
        return result

    if not "basis" in options:
        arg = orig_expr.args[-1]
        options["basis"] = get_basis_state(arg, **options)
    elif not isinstance(options["basis"], StateBase):
        options["basis"] = get_basis_state(orig_expr, **options)

    basis = options.pop("basis", None)

    if basis is None:
        return result

    delta = options.pop("delta", False)
    if delta:
        unities = options.pop("delta_unities", [])
    else:
        unities = options.pop("unities", [])

    if len(unities) == 0:
        return result

    kets = enumerate_states(basis, unities)
    coords = [k.label[0] for k in kets]

    for coord in coords:
        if coord in result.free_symbols:
            #TODO: Add support for sets of operators
            basis_op = state_to_operators(basis)
            start = basis_op.hilbert_space.interval.start
            end = basis_op.hilbert_space.interval.end
            result = integrate(result, (coord, start, end))

    return result

def get_basis_state(expr, **options):
    """
    Returns a basis state instance corresponding to the basis
    specified in options=s. If no basis is specified, the function
    will try to call the internal _get_default_basis function.

    There are three behaviors:

    1) The basis specified in options is already an instance of
    StateBase. If this is the case, it is simply returned. If the
    class is specified but not an instance, a default instance is returned.

    2) The basis specified is an operator or set of operators. If this
    is the case, the operator_to_state mapping method is used.

    3) No basis is specified. If expr is a state, then a default
    instance of its class is returned.
    If expr is an operator, then it is mapped to the corresponding state.
    If it is neither, then we cannot obtain the basis state.

    If the basis cannot be mapped, then it is not changed.

    This will be called from within represent, and represent will
    only pass QExpr's.

    TODO (?): Support for Muls and other types of expressions?

    Parameters
    ==========

    expr : Operator or StateBase
        Expression whose basis is sought

    Examples
    ========

    >>> from sympy.physics.quantum.represent import get_basis_state
    >>> from sympy.physics.quantum.cartesian import XOp, XKet, PxOp, PxKet
    >>> x = XKet()
    >>> X = XOp()
    >>> get_basis_state(x)
    |x>
    >>> get_basis_state(X)
    |x>
    >>> get_basis_state(x, basis=PxOp())
    |px>
    >>> get_basis_state(x, basis=PxKet)
    |px>

    """

    basis = options.pop("basis", None)

    if basis is None:
        return _make_default(_find_default_basis(expr, **options))
    elif (isinstance(basis, Operator) or isinstance(basis, set) or \
          (not isinstance(basis, StateBase) and issubclass(basis, Operator))):
        state = operators_to_state(basis)
        if state is None:
            return None
        elif isinstance(state, StateBase):
            return state
        else:
            return _make_default(state)
    elif isinstance(basis, StateBase):
        return basis
    elif issubclass(basis, StateBase):
        return _make_default(basis)
    else:
        return None


def _find_default_basis(expr, **options):
    if isinstance(expr, (Mul, Add, Pow, TensorProduct, Commutator, \
                         AntiCommutator, InnerProduct, OuterProduct)):
        for arg in expr.args:
            basis = _find_default_basis(arg, **options)
            if basis is not None:
                return basis

        return None
    elif isinstance(expr, QExpr):
        return expr._get_default_basis(**options)
    else:
        return None

def _make_default(expr):
    try:
        expr = expr()
    except Exception:
        return expr

    return expr

def enumerate_states(*args, **options):
    """
    Returns instances of the given state with dummy indices appended

    Operates in two different modes:

    1) Two arguments are passed to it. The first is the base state
    which is to be indexed, and the second argument is a list of
    indices to append.

    2) Three arguments are passed. The first is again the base state
    to be indexed. The second is the start index for counting.
    The final argument is the number of kets you wish to receive.

    Tries to call state._enumerate_state. If this fails, returns an empty list

    Parameters
    ==========

    args : list
        See list of operation modes above for explanation

    Examples
    ========

    >>> from sympy.physics.quantum.cartesian import XBra, XKet
    >>> from sympy.physics.quantum.represent import enumerate_states
    >>> test = XKet('foo')
    >>> enumerate_states(test, 1, 3)
    [|foo_1>, |foo_2>, |foo_3>]
    >>> test2 = XBra('bar')
    >>> enumerate_states(test2, [4, 5, 10])
    [<bar_4|, <bar_5|, <bar_10|]

    """

    state = args[0]

    if not (len(args) == 2 or len(args) == 3):
        raise NotImplementedError("Wrong number of arguments!")

    if not isinstance(state, StateBase):
        raise TypeError("First argument is not a state!")

    if len(args) == 3:
        num_states = args[2]
        options['start_index'] = args[1]
    else:
        num_states = len(args[1])
        options['index_list'] = args[1]

    try:
        ret = state._enumerate_state(num_states, **options)
    except NotImplementedError:
        ret = []

    return ret

def do_qapply(expr):
    if not isinstance(expr, Mul):
        return expr

    if expr.has(Wavefunction) and expr.has(DifferentialOperator):
        return qapply(expr)
    else:
        return expr

def _unwrap_wf(expr):
    if not isinstance(expr, Expr):
        return (expr, [])

    new_args = [None for arg in expr.args]
    unwrapped_vars = []
    has_wf = False
    for i,arg in enumerate(expr.args):
        if isinstance(arg, Wavefunction):
            new_args[i] = arg.expr
            for v in arg.variables:
                unwrapped_vars.append(v)
            has_wf = True
        else:
            new_args[i] = arg

    if not has_wf:
        return (expr, [])
    else:
        return (expr.__class__(*new_args), unwrapped_vars)

def _rewrap_wf(expr, unwrapped_vars, **options):
    if not isinstance(expr, Expr):
        return expr

    new_args = []
    tmp_args = []

    delta = options.pop("keep_delta", False)

    if not expr.has(DifferentialOperator):
        return Wavefunction(expr, *unwrapped_vars, keep_delta=delta)

    for arg in expr.args:
        if not isinstance(arg, DifferentialOperator):
            tmp_args.append(arg)
        else:
            expr_part = expr.__class__(*tmp_args)
            wf = Wavefunction(expr_part, *unwrapped_vars, keep_delta=delta)
            new_args.append(wf)
            new_args.append(arg)
            tmp_args = []

    if len(tmp_args) != 0:
        wf = Wavefunction(expr.__class__(*tmp_args), *unwrapped_vars, keep_delta=delta)
        new_args.append(wf)

    return expr.__class__(*new_args)
