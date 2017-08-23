from __future__ import (absolute_import, division, print_function)

from sympy import And, Gt, Lt, Abs, Dummy, oo, Tuple, Symbol, Function, Pow
from sympy.codegen.ast import (
    Assignment, AddAugmentedAssignment, CodeBlock, Declaration, FunctionDefinition,
    Print, Return, Scope, Statement, While, Variable, Pointer, real
)


def _decl_real(arg, value=None):
    return Declaration(Variable(arg, type=real), value)


def newtons_method(expr, wrt, atol=1e-12, delta=None, debug=False,
                   itermax=None, counter=None):
    """ Generates an AST for Newton-Raphson method (a root-finding algorithm).

    Returns an abstract syntax tree (AST) based on ``sympy.codegen.ast`` for Netwon's
    method of root-finding based on the following formula:

    .. math::

        x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}

    Parameters
    ----------
    expr : expression
    wrt : Symbol
        With respect to, i.e. what is the variable

    Examples
    --------
    >>> from sympy import symbols, cos
    >>> from sympy.codegen.ast import Assignment
    >>> from sympy.codegen.algorithms import newtons_method
    >>> x, dx, atol = symbols('x dx atol')
    >>> expr = cos(x) - x**3
    >>> algo = newtons_method(expr, x, atol, dx)
    >>> algo.has(Assignment(dx, -expr/expr.diff(x)))
    True

    References
    ==========
    https://en.wikipedia.org/wiki/Newton%27s_method

    """

    if delta is None:
        delta = Dummy()
        Wrapper = Scope
        name_d = 'delta'
    else:
        Wrapper = lambda x: x
        name_d = delta.name

    delta_expr = -expr/expr.diff(wrt)
    body = [Assignment(delta, delta_expr), AddAugmentedAssignment(wrt, delta)]
    if debug:
        prnt = Statement(Print([wrt, delta], r"{0}=%12.5g {1}=%12.5g\n".format(wrt.name, name_d)))
        body = [body[0], prnt] + body[1:]
    req = Gt(Abs(delta), atol)
    declars = [Statement(_decl_real(delta, oo))]
    if itermax is not None:
        counter = counter or Dummy(integer=True)
        v_counter = Variable.deduced(counter, 0)
        declars.append(Declaration(v_counter))
        body.append(AddAugmentedAssignment(counter, 1))
        req = And(req, Lt(counter, itermax))
    whl = While(req, CodeBlock(*body))
    blck = declars + [whl]
    return Wrapper(CodeBlock(*blck))


def _symbol_of(arg):
    if isinstance(arg, Declaration):
        arg = arg.variable
    if isinstance(arg, (Variable, Pointer)):
        arg = arg.symbol
    return arg


def newtons_method_function(expr, wrt, params=None, func_name="newton", attrs=Tuple(), **kwargs):
    """ Generates an AST for a function implementing the Newton-Raphson method.

    See also
    ========
    - sympy.codegen.ast.newtons_method

    """
    if params is None:
        params = (params,)
    pointer_subs = {p.symbol: Symbol('(*%s)' % p.symbol.name)
                    for p in params if isinstance(p, Pointer)}
    delta = kwargs.pop('delta', None)
    if delta is None:
        delta = Symbol('d_' + wrt.name)
        if expr.has(delta):
            delta = None  # will use Dummy
    algo = newtons_method(expr, wrt, delta=delta, **kwargs).xreplace(pointer_subs)
    if isinstance(algo, Scope):
        algo = algo.body
    not_in_params = expr.free_symbols.difference(set(_symbol_of(p) for p in params))
    if not_in_params:
        raise ValueError("Missing symbols in params: %s" % ', '.join(map(str, not_in_params)))
    declars = tuple(map(Variable.deduced, params))
    body = CodeBlock(algo, Statement(Return(wrt)))
    return FunctionDefinition(real, func_name, declars, body, attrs=attrs)
